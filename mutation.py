from Bio import SeqIO, pairwise2
from Bio.Blast import NCBIWWW, NCBIXML
import requests
import io
import difflib
from IPython.display import HTML, display
import streamlit as st
from io import StringIO
from Bio import Entrez
from tabulate import tabulate
import plotly.graph_objects as go

Entrez.email = "your_email@example.com"  # Replace with your real email

def fetch_sequence_by_accession(accession, seq_type):
    db = "nuccore" if seq_type == "nucleotide" else "protein"
    st.write(f"⬇️  Fetching user sequence for {accession} from NCBI ({db})...")
    handle = Entrez.efetch(db=db, id=accession, rettype="fasta", retmode="text")
    record = SeqIO.read(handle, "fasta")
    return str(record.seq), record.id

def detect_sequence_type(seq):
    seq_u = seq.upper()
    prot = set("ACDEFGHIKLMNPQRSTVWY")
    nuc = set("ACGTU")
    s = set(seq_u)
    if s.issubset(nuc) or len(s & nuc) / len(s) > 0.95:
        return "nucleotide"  # Add missing return
    elif s.issubset(prot):
        return "protein"     # Add missing return
    else:
        return "unknown"     # Add missing return

def run_blast(seq, seq_type, db):
    prog = "blastp" if seq_type == "protein" else "blastn"
    print(f"Running {prog.upper()} against {db}")
    handle = NCBIWWW.qblast(prog, db, seq, format_type="XML", query_believe_defline=False)
    record = NCBIXML.read(handle)
    top = record.alignments[0]
    return top  # Add return statement

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
                elif a == '-' and b != '-':
                    ref_colored += "-"
                    user_colored += f"<span style='background-color:#8ff'>{b}</span>"
                elif a != '-' and b == '-':
                    ref_colored += f"<span style='background-color:#fbb'>{a}</span>"
                    user_colored += "-"

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

def emboss_style_html_doc(name1, name2, matrix, gap_open, gap_extend, aln_seqA, aln_seqB, label1="Ref", label2="User Seq"):
    length, identity, similarity, gaps = calculate_stats(aln_seqA, aln_seqB)

    html_output = f"""
    <div style="font-family:Consolas, 'Courier New', monospace; font-size:15px; white-space:pre;">
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

        ref_colored = ""
        user_colored = ""

        for a, b in zip(ref_block, user_block):
            if a == b and a != '-':
                ref_colored += f"<span>{a}</span>"
                user_colored += f"<span>{b}</span>"
            elif a != '-' and b != '-':
                # Substitution
                ref_colored += f"<span style='background-color:#f88'>{a}</span>"
                user_colored += f"<span style='background-color:#8f8'>{b}</span>"
            elif a == '-' and b != '-':
                # Insertion (in user)
                ref_colored += "<span style='background-color:#e0e0e0'>-</span>"
                user_colored += f"<span style='background-color:#8ff'>{b}</span>"
            elif a != '-' and b == '-':
                # Deletion (in user)
                ref_colored += f"<span style='background-color:#fbb'>{a}</span>"
                user_colored += "<span style='background-color:#e0e0e0'>-</span>"

        def count_residues(seq):
            return sum(1 for c in seq if c != '-')

        ref_end = ref_pos + count_residues(ref_block) - 1
        user_end = user_pos + count_residues(user_block) - 1

        html_output += (
            f"<span style='font-weight:bold'>{label2:<{label_width}}</span>{user_pos:>4} {user_colored} {user_end}<br>"
            f"<br>"  # Blank line after user seq
            f"<span style='font-weight:bold'>{label1:<{label_width}}</span>{ref_pos:>4} {ref_colored} {ref_end}<br>"
            f"<br>"  # Blank line after ref seq
        )

        ref_pos = ref_end + 1
        user_pos = user_end + 1

    html_output += "</div>"
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
    accession, hit_title = run_blast(user_szeq, seq_type, "nr" if seq_type == "protein" else "nt")
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

@st.cache_data(show_spinner="Running BLAST (cached)...")
def run_blast_against_refseq(seq, seq_type, hitlist_size=2):
    from Bio.Blast import NCBIWWW, NCBIXML
    program = "blastn" if seq_type == "nucleotide" else "blastp"
    db = "refseq_rna" if seq_type == "nucleotide" else "refseq_protein"
    result_handle = NCBIWWW.qblast(program, db, seq, hitlist_size=hitlist_size)
    blast_record = NCBIXML.read(result_handle)

    if not blast_record.alignments:
        return []

    hits = []
    for alignment in blast_record.alignments[:hitlist_size]:
        hsp = alignment.hsps[0]
        hit_info = {
            "accession": alignment.accession,
            "title": alignment.hit_def,
            "ref_aligned": hsp.sbjct,
            "user_aligned": hsp.query,
            "alignment_match_line": hsp.match,
            "ref_start": hsp.sbjct_start,
            "ref_end": hsp.sbjct_end,
            "score": hsp.score,
            "evalue": hsp.expect,
            "identity": hsp.identities,
            "align_len": hsp.align_length,
            "percent_identity": round((hsp.identities / hsp.align_length) * 100, 2),
            "query_coverage": round((hsp.align_length / len(seq)) * 100, 2)
        }
        hits.append(hit_info)
    return hits

def print_blast_hits_table(hits):
    st.subheader("📊 Top BLAST Hits Summary")
    table_data = []
    headers = ["Index", "Accession", "Identity (%)", "E-value", "Score", "Coverage (%)", "Title"]
    for i, hit in enumerate(hits):
        acc_link = f"https://www.ncbi.nlm.nih.gov/nucleotide/{hit['accession']}"
        acc_md = f"[{hit['accession']}]({acc_link})"
        row = [
            i + 1,
            acc_md,
            hit['percent_identity'],
            hit['evalue'],
            hit['score'],
            hit['query_coverage'],
            hit['title']
        ]
        table_data.append(row)
    table_md = tabulate(table_data, headers=headers, tablefmt="github")
    st.markdown(table_md, unsafe_allow_html=True)

def print_mutations_table(mutations):
    if not mutations:
        st.write(" - None")
        return

    st.subheader("🧬 Detected Mutations (Ref → User)")
    table_data = []
    headers = ["Position", "Ref", "User", "Mutation"]
    for pos, ref, user in mutations:
        mutation_str = f"c.{pos}{ref}>{user}"
        row = [pos, ref, user, mutation_str]
        table_data.append(row)
    table_md = tabulate(table_data, headers=headers, tablefmt="github")
    st.markdown(table_md, unsafe_allow_html=True)

def print_emboss_style_alignment(user_seq, ref_seq, match_line=None, block_size=60, use_pairwise=False):
    if use_pairwise:
        from Bio import pairwise2
        # Use Biopython's pairwise2 for local alignment
        alignments = pairwise2.align.localms(user_seq, ref_seq, 2, -1, -0.5, -0.1)
        best_alignment = alignments[0]
        aligned_user = best_alignment.seqA
        aligned_ref = best_alignment.seqB

        # Generate match line and track mutations
        match_line = []
        mutations = []
        for i, (u, r) in enumerate(zip(aligned_user, aligned_ref)):
            if u == r:
                match_line.append('|')
            else:
                match_line.append(' ')
                if u != '-' and r != '-':
                    mutations.append((i + 1, r, u))  # 1-based position
        match_line = ''.join(match_line)
    else:
        # Use BLAST-provided alignment
        aligned_user = user_seq
        aligned_ref = ref_seq
        if match_line is None:
            match_line = ''.join('|' if u == r else ' ' for u, r in zip(user_seq, ref_seq))
        mutations = []
        for i, (u, r) in enumerate(zip(aligned_user, aligned_ref)):
            if u != r and u != '-' and r != '-':
                mutations.append((i + 1, r, u))  # 1-based position

    # Print alignment
    st.subheader("🧬 Alignment")
    alignment_text = ""
    length = len(aligned_user)
    for start in range(0, length, block_size):
        end = min(start + block_size, length)
        alignment_text += f"User  {start+1:>5}  {aligned_user[start:end]}  {end}\n"
        alignment_text += f"      {'':5}  {match_line[start:end]}\n"
        alignment_text += f"Ref   {start+1:>5}  {aligned_ref[start:end]}  {end}\n\n"
    
    st.code(alignment_text, language=None)
    
    # Print compact mutation table
    print_mutations_table(mutations)

def analyze_hit(user_seq, hit, index, use_pairwise=False):
    st.subheader(f"🧬 Analyzing Hit {index + 1}: {hit['accession']}")
    st.write(f"📍 Alignment Region: Ref {hit['ref_start']}–{hit['ref_end']}")
    st.write(f"📈 Identity: {hit['percent_identity']}% | Coverage: {hit['query_coverage']}% | E-value: {hit['evalue']}")

    # --- Alignment Summary ---
    length, identity, similarity, gaps = calculate_stats(hit['ref_aligned'], hit['user_aligned'])
    st.subheader("📄 Alignment Summary")
    st.table({
        "Length": [length],
        "Identity": [f"{identity}/{length} ({identity/length*100:.1f}%)"],
        "Similarity": [f"{similarity}/{length} ({similarity/length*100:.1f}%)"],
        "Gaps": [f"{gaps}/{length} ({gaps/length*100:.1f}%)"]
    })

    # --- Alignment Visualization (HTML only) ---
    st.subheader("🔗 Alignment Visualization")
    html_output = emboss_style_html(
        name1=f"Ref ({hit['accession']})",
        name2="User Seq",
        matrix="EBLOSUM62",
        gap_open=10.0,
        gap_extend=0.5,
        aln_seqA=hit['ref_aligned'],
        aln_seqB=hit['user_aligned'],
        label1=f"Ref ({hit['accession']})",
        label2="User Seq"
    )
    html_output_doc = emboss_style_html_doc(
        name1=f"Ref ({hit['accession']})",
        name2="User Seq",
        matrix="EBLOSUM62",
        gap_open=10.0,
        gap_extend=0.5,
        aln_seqA=hit['ref_aligned'],
        aln_seqB=hit['user_aligned'],
        label1=f"Ref ({hit['accession']})",
        label2="User Seq"
    )
    st.components.v1.html(
        f"""
        <div style="font-family:monospace; font-size:16px;">
        {html_output}
        </div>
        """,
        height=600,
        scrolling=True
    )

    # --- Mutations Table ---
    mutations = []
    for i, (u, r) in enumerate(zip(hit['user_aligned'], hit['ref_aligned'])):
        if u != r and u != '-' and r != '-':
            mutations.append((i + 1, r, u))

    # Show mutation table in Streamlit (with expander)
    with st.expander("🧬 Detected Mutations (Ref → User)"):
        print_mutations_table(mutations)

    # --- Mutation Table as HTML ---
    mutation_table_html = ""
    if mutations:
        mutation_table_html = "<h3>Detected Mutations (Ref → User)</h3><table border='1' cellpadding='4' cellspacing='0'><tr><th>Position</th><th>Ref</th><th>User</th><th>Mutation</th></tr>"
        for pos, ref, user in mutations:
            mutation_str = f"c.{pos}{ref}&gt;{user}"
            mutation_table_html += f"<tr><td>{pos}</td><td>{ref}</td><td>{user}</td><td>{mutation_str}</td></tr>"
        mutation_table_html += "</table>"

    # --- Download Option (HTML .doc, fully formatted) ---
    html_doc = f"""
    <html>
    <head>
    <meta charset="UTF-8">
    <title>Alignment Summary</title>
    <style>
    body {{ font-family: Arial, sans-serif; }}
    pre {{ font-family: 'Consolas', 'Courier New', monospace; font-size: 15px; }}
    table {{ border-collapse: collapse; font-size: 15px; }}
    th, td {{ border: 1px solid #888; padding: 4px 8px; }}
    th {{ background: #e0e0e0; }}
    </style>
    </head>
    <body>
    <h2>Alignment Summary</h2>
    <b>Ref:</b> {hit['accession']}<br>
    <b>User Seq</b><br>
    <b>Identity:</b> {hit['percent_identity']}%<br>
    <b>Coverage:</b> {hit['query_coverage']}%<br>
    <b>E-value:</b> {hit['evalue']}<br>
    <h3>Alignment Visualization</h3>
    <div style="font-family:monospace; font-size:16px;">
    {html_output_doc}
    </div>
    {mutation_table_html}
    </body>
    </html>
    """
    st.download_button(
        label="⬇️ Download Alignment (.doc)",
        data=html_doc,
        file_name=f"alignment_{hit['accession']}.doc",
        mime="application/msword"
    )

    # --- Mutation Track Diagram ---
    st.subheader("🧬 Mutation Track Diagram")
    plot_mutation_track(
        hit['ref_aligned'],
        hit['user_aligned'],
        ref_label=f"Ref ({hit['accession']})",
        user_label="User Seq"
    )
    st.subheader("Interactive Sequence Differences")
    fig = plot_mutation_scatter(hit['ref_aligned'], hit['user_aligned'], ref_label=f"Ref ({hit['accession']})", user_label="User Seq")
    st.plotly_chart(fig, use_container_width=True)

def plot_mutation_scatter(ref_seq, user_seq, ref_label="Ref", user_label="User"):
    positions = []
    ref_bases = []
    user_bases = []
    hover_texts = []
    y_vals = []

    for i, (r, u) in enumerate(zip(ref_seq, user_seq)):
        if r != u:
            positions.append(i+1)
            ref_bases.append(r)
            user_bases.append(u)
            hover_texts.append(f"{ref_label}: {r} / {user_label}: {u}")
            y_vals.append("Overlap")

    fig = go.Figure()

    fig.add_trace(go.Scatter(
        x=positions,
        y=y_vals,
        mode='markers',
        marker=dict(color='purple', size=12, symbol='diamond'),
        name="Mutations",
        text=hover_texts,
        hovertemplate="Position: %{x}<br>%{text}<extra></extra>"
    ))

    fig.update_layout(
        title="Interactive Sequence Differences",
        xaxis_title="Position",
        yaxis=dict(
            tickvals=["Ref", "Overlap", "User"],
            ticktext=[ref_label, "Overlap", user_label],
            categoryorder='array',
            categoryarray=["Ref", "Overlap", "User"]
        ),
        legend=dict(x=1, y=1),
        height=400,
        margin=dict(l=60, r=60, t=60, b=60)
    )

    return fig

def plot_mutation_track(ref_seq, user_seq, ref_label="REF", user_label="ALT"):
    """Draw a genome-browser style mutation track diagram for two sequences."""
    positions = []
    ref_bases = []
    user_bases = []
    mutation_labels = []
    for i, (r, u) in enumerate(zip(ref_seq, user_seq)):
        if r != u:
            positions.append(i+1)
            ref_bases.append(r)
            user_bases.append(u)
            mutation_labels.append(f"{r}>{u} at {i+1}")

    # Main track for REF
    fig = go.Figure()
    fig.add_trace(go.Scatter(
        x=list(range(1, len(ref_seq)+1)),
        y=[1]*len(ref_seq),
        mode="lines",
        line=dict(color="black", width=3),
        name=ref_label,
        hoverinfo="skip"
    ))

    # Main track for ALT (user)
    fig.add_trace(go.Scatter(
        x=list(range(1, len(user_seq)+1)),
        y=[0.8]*len(user_seq),
        mode="lines",
        line=dict(color="red", width=3),
        name=user_label,
        hoverinfo="skip"
    ))

    # Mutation markers
    if positions:
        fig.add_trace(go.Scatter(
            x=positions,
            y=[1.1]*len(positions),
            mode="markers+text",
            marker=dict(color="orange", size=14, symbol="line-ns"),
            text=mutation_labels,
            textposition="top center",
            name="Mutation",
            hovertemplate="%{text}<extra></extra>"
        ))
        # Draw stems
        for x in positions:
            fig.add_shape(type="line", x0=x, x1=x, y0=0.82, y1=1.08, line=dict(color="orange", width=2))

    fig.update_layout(
        title="Mutation Track Diagram",
        xaxis_title="Sequence Position",
        yaxis=dict(
            tickvals=[0.8, 1],
            ticktext=[user_label, ref_label],
            showgrid=False,
            zeroline=False,
            range=[0.7, 1.2]
        ),
        showlegend=True,
        height=300,
        margin=dict(l=40, r=40, t=40, b=40),
        plot_bgcolor="white"
    )
    st.plotly_chart(fig, use_container_width=True)

def main():
    st.set_page_config(page_title="Mutation Analysis", layout="wide")
    st.title("🌟 Mutation Analysis Dashboard")
    
    # Sidebar inputs
    st.sidebar.header("Input Parameters")
    user_accession = st.sidebar.text_input("NCBI Accession", value="OR905937")
    seq_type = st.sidebar.selectbox("Sequence Type", ["nucleotide", "protein"])
    
    if user_accession:
        with st.spinner("Fetching sequence and running BLAST analysis..."):
            user_seq, _ = fetch_sequence_by_accession(user_accession, seq_type)
            blast_hits = run_blast_against_refseq(user_seq, seq_type, hitlist_size=2)
        
        if not blast_hits:
            return

        print_blast_hits_table(blast_hits)

        # 🔁 Automatically analyze top 1 hit (default)
        analyze_hit(user_seq, blast_hits[0], index=0)

        # 🔁 Button for second hit analysis
        if len(blast_hits) > 1:
            button_label = f"Analyze Second Top Hit ({blast_hits[1]['accession']})"
            button_html = f"""
                <style>
                .analyze-btn {{
                    background-color: #1976d2;
                    color: white;
                    padding: 0.5em 1.5em;
                    border: none;
                    border-radius: 6px;
                    font-size: 1.1em;
                    font-weight: bold;
                    cursor: pointer;
                    margin-top: 10px;
                    margin-bottom: 10px;
                }}
                .analyze-btn:hover {{
                    background-color: #0d47a1;
                }}
                </style>
                <form action="" method="post">
                    <button class="analyze-btn" type="submit">{button_label}</button>
                </form>
            """
            if st.markdown(button_html, unsafe_allow_html=True):
                analyze_hit(user_seq, blast_hits[1], index=1)

if __name__ == "__main__":
    main()
