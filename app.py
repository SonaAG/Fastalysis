from Bio import SeqIO, Entrez
from Bio.Blast import NCBIWWW, NCBIXML
import streamlit as st
import pandas as pd
from tabulate import tabulate
import re

Entrez.email = "your@email.com"  # Replace with your actual email

def load_user_sequence(path):
    record = SeqIO.read(path, "fasta")
    seq = str(record.seq)
    return seq, record.id

def detect_sequence_type(seq):
    protein_set = set("ACDEFGHIKLMNPQRSTVWY")
    nucleotide_set = set("ACGTU")
    seq_set = set(seq.upper())
    if seq_set.issubset(nucleotide_set) or len(seq_set & nucleotide_set) / len(seq_set) > 0.95:
        return "nucleotide"
    elif seq_set.issubset(protein_set):
        return "protein"
    else:
        raise ValueError("Unable to detect sequence type.")

def run_blast_top_hits(seq, seq_type, hitlist_size=5):
    program = "blastn" if seq_type == "nucleotide" else "blastp"
    db = "nt" if seq_type == "nucleotide" else "nr"
    result_handle = NCBIWWW.qblast(program, db, seq, hitlist_size=hitlist_size)
    blast_record = NCBIXML.read(result_handle)
    hits = []
    for alignment in blast_record.alignments:
        hsp = alignment.hsps[0]
        accession = alignment.accession
        title = alignment.hit_def
        score = hsp.score
        evalue = hsp.expect
        identities = hsp.identities
        align_len = hsp.align_length
        percent_identity = round((identities / align_len) * 100, 2)
        query_coverage = min(round((align_len / len(seq)) * 100, 2), 100.0)
        hits.append({
            "Accession": accession,
            "Title": title,
            "Identity (%)": percent_identity,
            "E-value": evalue,
            "Score": score,
            "Coverage (%)": query_coverage
        })
    return hits

def run_pipeline_streamlit(path):
    st.header("🔬 Fastalysis Sequence Analysis Dashboard")
    user_seq, user_id = load_user_sequence(path)
    seq_type = detect_sequence_type(user_seq)
    st.subheader("📄 Sequence Info")
    st.markdown(f"- **User Sequence ID:** `{user_id}`")
    st.markdown(f"- **Detected Type:** `{seq_type}`")
    st.markdown(f"- **User Sequence Length:** `{len(user_seq)}`")

    st.subheader("🔍 BLAST Top Hits (NCBI)")
    num_hits = st.sidebar.slider("Number of BLAST Hits to Display", min_value=1, max_value=20, value=5)
    with st.spinner("Running BLAST search and fetching top hits from NCBI..."):
        hits = run_blast_top_hits(user_seq, seq_type, hitlist_size=num_hits)
    if hits:
        st.success("BLAST search completed!")

        # Prepare data for tabulate
        table_data = []
        headers = ["Index", "Accession", "Identity (%)", "E-value", "Score", "Coverage (%)", "Title", "Download"]
        for idx, hit in enumerate(hits, 1):
            if seq_type == "nucleotide":
                acc_link = f"https://www.ncbi.nlm.nih.gov/nucleotide/{hit['Accession']}"
                fasta_link = f"https://www.ncbi.nlm.nih.gov/sviewer/viewer.fcgi?id={hit['Accession']}&db=nuccore&report=fasta"
            else:
                acc_link = f"https://www.ncbi.nlm.nih.gov/protein/{hit['Accession']}"
                fasta_link = f"https://www.ncbi.nlm.nih.gov/sviewer/viewer.fcgi?id={hit['Accession']}&db=protein&report=fasta"
            acc_md = f"[{hit['Accession']}]({acc_link})"
            download_md = f"[Fasta]({fasta_link})"
            row = [
                idx,
                acc_md,
                hit['Identity (%)'],
                hit['E-value'],
                hit['Score'],
                hit['Coverage (%)'],
                hit['Title'],
                download_md
            ]
            table_data.append(row)
        table_md = tabulate(table_data, headers=headers, tablefmt="github")
        st.markdown(table_md, unsafe_allow_html=True)

       
    else:
        st.error("No BLAST hits found.")

if __name__ == "__main__":
    st.set_page_config(page_title="Fastalysis Dashboard", layout="wide")
    st.sidebar.header("Upload FASTA File")
    uploaded_file = st.sidebar.file_uploader("Choose a FASTA file", type=["fasta"])
    if uploaded_file:
        temp_path = "user_sequence.fasta"
        with open(temp_path, "wb") as f:
            f.write(uploaded_file.getbuffer())
        run_pipeline_streamlit(temp_path)
    else:
        st.info("Please upload a FASTA file to start analysis.")
