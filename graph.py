from Bio import SeqIO, pairwise2
from Bio.Blast import NCBIWWW, NCBIXML
import requests
import io
import difflib
from IPython.display import HTML, display
import streamlit as st
from io import StringIO

file_path = "/content/user_sequence.fasta"  # Your FASTA file path

def read_user_sequence(path):
    record = SeqIO.read(path, "fasta")
    seq = str(record.seq).strip().replace("\n", "").replace("\r", "")
    print(f"✅ Uploaded sequence length: {len(seq)}")
    return seq, record.id

def detect_sequence_type(seq):
    seq_u = seq.upper()
    prot = set("ACDEFGHIKLMNPQRSTVWY")
    nuc = set("ACGTU")
    s = set(seq_u)
    if s.issubset(nuc) or len(s & nuc) / len(s) > 0.95:
        print("🧪 Detected type: Nucleotide")
        return "nucleotide"
    elif s.issubset(prot):
        print("🧪 Detected type: Protein")
        return "protein"
    else:
        raise ValueError("Unknown sequence type")

def run_blast(seq, seq_type, db):
    prog = "blastp" if seq_type == "protein" else "blastn"
    print(f"Running {prog.upper()} against {db}")
    handle = NCBIWWW.qblast(prog, db, seq, format_type="XML", query_believe_defline=False)
    record = NCBIXML.read(handle)
    top = record.alignments[0]
    return top.accession, top.hit_def

def fetch_ncbi(accession, seq_type):
    db = "protein" if seq_type == "protein" else "nucleotide"
    print(f"Fetching NCBI {db} accession: {accession}")
    url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
    params = {"db": db, "id": accession, "rettype": "fasta", "retmode": "text"}
    resp = requests.get(url, params=params)
    resp.raise_for_status()
    print("Fetched raw data:", resp.text.splitlines()[0])
    fio = io.StringIO(resp.text)
    rec = SeqIO.read(fio, "fasta")
    return str(rec.seq), rec.id

def find_mutations(ref, seq):
    aln = pairwise2.align.globalms(ref, seq, 2, -1, -2, -0.5)[0]
    return aln.seqA, aln.seqB

def calculate_stats(aln_seqA, aln_seqB):
    length = len(aln_seqA)
    identity = sum(a == b and a != '-' for a, b in zip(aln_seqA, aln_seqB))
    gaps = sum(a == '-' or b == '-' for a, b in zip(aln_seqA, aln_seqB))
    similarity = identity  # Simplified similarity
    return length, identity, similarity, gaps

# Change print_emboss_style_html to return HTML string
def emboss_style_html(name1, name2, matrix, gap_open, gap_extend, aln_seqA, aln_seqB, label1="Ref", label2="User Seq"):
    length, identity, similarity, gaps = calculate_stats(aln_seqA, aln_seqB)

    html_output = f"""
    <pre>
# 1: {label1}
# 2: {label2}
# Matrix: {matrix}
# Gap_penalty: {gap_open}
# Extend_penalty: {gap_extend}
#
# Length: {length}
# Identity:     {identity}/{length} ({identity/length*100:.1f}%)
# Similarity:   {similarity}/{length} ({similarity/length*100:.1f}%)
# Gaps:        {gaps}/{length} ({gaps/length*100:.1f}%)
#
"""

    line_width = 60
    ref_pos = user_pos = 1
    label_width = max(len(label1), len(label2)) + 2

    for i in range(0, length, line_width):
        ref_block = aln_seqA[i:i+line_width]
        user_block = aln_seqB[i:i+line_width]
        match_line = ""

        ref_colored = ""
        user_colored = ""

        for a, b in zip(ref_block, user_block):
            if a == b and a != '-':
                match_line += "|"
                ref_colored += a
                user_colored += b
            else:
                match_line += " "
                if a != '-' and b != '-':
                    ref_colored += f"<span style='background-color:#f88'>{a}</span>"
                    user_colored += f"<span style='background-color:#8f8'>{b}</span>"
                else:
                    ref_colored += a
                    user_colored += b

        def count_residues(seq):
            return sum(1 for c in seq if c != '-')

        ref_end = ref_pos + count_residues(ref_block) - 1
        user_end = user_pos + count_residues(user_block) - 1

        html_output += (
            f"{label2:<{label_width}}{user_pos:>4} {user_colored} {user_end}\n"
            f"{'':<{label_width}}     {match_line}\n"
            f"{label1:<{label_width}}{ref_pos:>4} {ref_colored} {ref_end}\n\n"
        )

        ref_pos = ref_end + 1
        user_pos = user_end + 1

    html_output += "</pre>"
    return html_output

def list_mutations(ref_seq, user_seq):
    matcher = difflib.SequenceMatcher(None, ref_seq, user_seq)
    mutations = []
    ref_pos = 0
    for tag, i1, i2, j1, j2 in matcher.get_opcodes():
        if tag == "equal":
            ref_pos += i2 - i1
        else:
            for off in range(max(i2 - i1, j2 - j1)):
                r = ref_seq[i1 + off] if i1 + off < i2 else '-'
                u = user_seq[j1 + off] if j1 + off < j2 else '-'
                if r != u:
                    mutations.append((ref_pos + 1, r, u))
                if r != '-': ref_pos += 1
    if mutations:
        print("📌 Mutations found (Position: Reference → User):")
        for pos, r, u in mutations:
            print(f"{pos}: {r} → {u}")
    else:
        print("✅ No mutations detected.")

def run_pipeline(path):
    user_seq, user_id = read_user_sequence(path)
    seq_type = detect_sequence_type(user_seq)
    accession, hit_title = run_blast(user_seq, seq_type, "nr" if seq_type == "protein" else "nt")
    ref_seq, ref_id = fetch_ncbi(accession, seq_type)
    aln_seqA, aln_seqB = find_mutations(ref_seq, user_seq)

    print_emboss_style_html(
        name1=ref_id,
        name2=user_id,
        matrix="EBLOSUM62",
        gap_open=10.0,
        gap_extend=0.5,
        aln_seqA=aln_seqA,
        aln_seqB=aln_seqB,
        label1=f"Ref ({ref_id})",
        label2=f"User Seq ({user_id})"
    )

    list_mutations(aln_seqA, aln_seqB)

def run_pipeline_streamlit(path):
    st.header("🔬 Fastalysis Sequence Analysis Dashboard")
    user_seq, user_id = read_user_sequence(path)
    seq_type = detect_sequence_type(user_seq)
    accession, hit_title = run_blast(user_seq, seq_type, "nr" if seq_type == "protein" else "nt")
    ref_seq, ref_id = fetch_ncbi(accession, seq_type)
    aln_seqA, aln_seqB = find_mutations(ref_seq, user_seq)

    st.subheader("📄 Sequence Info")
    st.markdown(f"- **User Sequence ID:** `{user_id}`")
    st.markdown(f"- **Detected Type:** `{seq_type}`")
    st.markdown(f"- **Top BLAST Hit:** `{hit_title}` (`{accession}`)")
    st.markdown(f"- **Reference Sequence ID:** `{ref_id}`")

    st.subheader("🧬 Alignment Summary")
    length, identity, similarity, gaps = calculate_stats(aln_seqA, aln_seqB)
    st.table({
        "Length": [length],
        "Identity": [f"{identity}/{length} ({identity/length*100:.1f}%)"],
        "Similarity": [f"{similarity}/{length} ({similarity/length*100:.1f}%)"],
        "Gaps": [f"{gaps}/{length} ({gaps/length*100:.1f}%)"]
    })

    st.subheader("🔗 Alignment Visualization")
    html_output = emboss_style_html(
        name1=ref_id,
        name2=user_id,
        matrix="EBLOSUM62",
        gap_open=10.0,
        gap_extend=0.5,
        aln_seqA=aln_seqA,
        aln_seqB=aln_seqB,
        label1=f"Ref ({ref_id})",
        label2=f"User Seq ({user_id})"
    )
    # Use st.components.v1.html for better HTML rendering
    st.components.v1.html(
        f"""
        <div style="font-family:monospace; font-size:16px;">
        {html_output}
        </div>
        """,
        height=600,
        scrolling=True
    )

    st.subheader("🧬 Mutation List")
    matcher = difflib.SequenceMatcher(None, aln_seqA, aln_seqB)
    mutations = []
    ref_pos = 0
    for tag, i1, i2, j1, j2 in matcher.get_opcodes():
        if tag == "equal":
            ref_pos += i2 - i1
        else:
            for off in range(max(i2 - i1, j2 - j1)):
                r = aln_seqA[i1 + off] if i1 + off < i2 else '-'
                u = aln_seqB[j1 + off] if j1 + off < j2 else '-'
                if r != u:
                    mutations.append((ref_pos + 1, r, u))
                if r != '-': ref_pos += 1
    if mutations:
        st.write("📌 Mutations found (Position: Reference → User):")
        st.table([{"Position": pos, "Reference": r, "User": u} for pos, r, u in mutations])
    else:
        st.success("✅ No mutations detected.")

if __name__ == "__main__":
    st.set_page_config(page_title="Fastalysis Dashboard", layout="wide")
    st.sidebar.header("Upload FASTA File")
    uploaded_file = st.sidebar.file_uploader("Choose a FASTA file", type=["fasta"])
    if uploaded_file:
        # Save uploaded file to a temp location
        temp_path = "user_sequence.fasta"
        with open(temp_path, "wb") as f:
            f.write(uploaded_file.getbuffer())
        run_pipeline_streamlit(temp_path)
    else:
        st.info("Please upload a FASTA file to start analysis.")
