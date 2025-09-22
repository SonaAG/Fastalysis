import pandas as pd
import streamlit as st
import pygwalker as pyg
from Bio import Entrez
from xml.etree import ElementTree as ET
from tabulate import tabulate
import requests
import json

Entrez.email = "your.email@example.com"  # Replace with your real email

def decode_iupac(code):
    """Decode IUPAC nucleotide codes"""
    iupac_codes = {
        "A": "A", "C": "C", "G": "G", "T": "T",
        "R": "A/G", "Y": "C/T", "S": "G/C", "W": "A/T",
        "K": "G/T", "M": "A/C", "B": "C/G/T", "D": "A/G/T",
        "H": "A/C/T", "V": "A/C/G", "N": "A/C/G/T"
    }
    return iupac_codes.get(code.upper(), code)

def extract_pubmed_ids_from_snp(snp_id):
    """Extract PubMed IDs from dbSNP entry"""
    try:
        handle = Entrez.efetch(db="snp", id=snp_id, rettype="xml")
        xml_data = handle.read()
        handle.close()
        
        root = ET.fromstring(xml_data)
        pubmed_ids = []
        
        # Look for citation elements
        for citation in root.iter():
            if 'pmid' in citation.tag.lower() or 'pubmed' in citation.tag.lower():
                if citation.text and citation.text.isdigit():
                    pubmed_ids.append(citation.text)
        
        return list(set(pubmed_ids))[:5]  # Return first 5 unique PMIDs
    except Exception as e:
        return []

def get_gene_symbol_from_refseq(accession):
    """Get gene symbol from RefSeq accession"""
    try:
        handle = Entrez.efetch(db="nucleotide", id=accession, rettype="gb", retmode="xml")
        records = Entrez.read(handle)
        handle.close()
        
        for record in records:
            if 'GBSeq_feature-table' in record:
                for feature in record['GBSeq_feature-table']:
                    if feature['GBFeature_key'] == 'gene':
                        for qualifier in feature.get('GBFeature_quals', []):
                            if qualifier['GBQualifier_name'] == 'gene':
                                return qualifier['GBQualifier_value']
        return None
    except Exception as e:
        return None

def fetch_dbsnp_variants_for_gene(gene_symbol, max_records=10):
    """
    Fetch dbSNP variant summaries for the given gene symbol.
    """
    try:
        handle = Entrez.esearch(db="snp", term=f"{gene_symbol}[Gene]", retmax=max_records)
        rec = Entrez.read(handle)
        handle.close()
        snp_ids = rec.get("IdList", [])
        variants = []
        for snp_id in snp_ids:
            sumh = Entrez.esummary(db="snp", id=snp_id)
            sums = Entrez.read(sumh)
            sumh.close()

            doc = sums["DocumentSummarySet"]["DocumentSummary"][0]
            rsid = f"rs{snp_id}"
            position = f"{doc.get('CHR', '?')}:{doc.get('CHRPOS', '?')}"
            allele_raw = doc.get("ALLELE", "-")
            change = decode_iupac(allele_raw)
            desc = doc.get("SNP_CLASS", "-")

            pubmed_ids = extract_pubmed_ids_from_snp(snp_id)

            variants.append({
                "Variant ID": rsid,
                "Position": position,
                "Change": change,
                "Description": desc,
                "PubMed IDs": pubmed_ids,
                "rsID": rsid,
                "Link": f"https://www.ncbi.nlm.nih.gov/snp/{rsid}"
            })
        return variants
    except Exception as e:
        print(f"Error fetching variants for gene {gene_symbol}: {e}")
        return []

def display_variants_table(variants):
    if not variants:
        st.warning("No variants found.")
        return

    table_data = []
    headers = ["Variant ID", "Position", "Change", "Description", "PubMed IDs", "rsID"]
    for variant in variants:
        pubmed_ids = variant.get('PubMed IDs', [])
        pubmed_str = "-" if not pubmed_ids else ", ".join(pubmed_ids)
        variant_id_md = f"[{variant['Variant ID']}]({variant['Link']})"
        row = [
            variant_id_md,
            variant['Position'],
            variant['Change'],
            variant['Description'],
            pubmed_str,
            variant['rsID']
        ]
        table_data.append(row)

    table_md = tabulate(table_data, headers=headers, tablefmt="github")
    st.markdown(table_md, unsafe_allow_html=True)

def main():
    st.set_page_config(page_title="Mutation Variant Table", layout="wide")
    st.title("🧬 Genetic Variant Table Generator")
    st.subheader("Generate detailed variant tables from RefSeq accessions or gene symbols")
    
    # Initialize session state
    if 'gene_symbol' not in st.session_state:
        st.session_state.gene_symbol = None
    if 'accession' not in st.session_state:
        st.session_state.accession = None
    
    # Sidebar controls
    st.sidebar.header("Search Parameters")
    search_method = st.sidebar.radio("Search Method", ["RefSeq Accession", "Gene Symbol"])
    
    if search_method == "RefSeq Accession":
        accession = st.sidebar.text_input("RefSeq Accession", value="NM_007294", 
                                         help="e.g., NM_007294 (BRCA1)")
        st.session_state.accession = accession
        
        if st.sidebar.button("🔍 Search by Accession"):
            with st.spinner("Extracting gene information..."):
                gene_symbol = get_gene_symbol_from_refseq(accession)
                if gene_symbol:
                    st.session_state.gene_symbol = gene_symbol
                    st.success(f"Found gene: **{gene_symbol}** for accession **{accession}**")
                else:
                    st.session_state.gene_symbol = None
                    st.error("Could not extract gene information from accession")
    else:
        gene_symbol = st.sidebar.text_input("Gene Symbol", value="BRCA1",
                                           help="e.g., BRCA1, TP53, EGFR")
        st.session_state.gene_symbol = gene_symbol
        st.session_state.accession = None
    
    max_variants = st.sidebar.slider("Maximum Variants", 10, 50, 20)
    
    # Show current gene symbol if available
    if st.session_state.gene_symbol:
        st.info(f"Current gene: **{st.session_state.gene_symbol}**")
    
    # Search button - now works with session state
    if st.sidebar.button("📊 Generate Variant Table"):
        if st.session_state.gene_symbol:
            with st.spinner(f"Fetching variants for {st.session_state.gene_symbol}..."):
                variants = fetch_dbsnp_variants_for_gene(st.session_state.gene_symbol, max_variants)
                
                if variants:
                    st.success(f"Found {len(variants)} variants for **{st.session_state.gene_symbol}**")
                    display_variants_table(variants)
                    
                    
                    # Download option
                    if st.button("💾 Download Variant Table"):
                        # Create downloadable text file
                        output = f"Variant Table for {st.session_state.gene_symbol}\n"
                        output += "=" * 80 + "\n\n"
                        
                        table_data = []
                        headers = ["Variant ID", "Position", "Change", "Description", "PubMed IDs", "rsID", "Link"]
                        
                        for variant in variants:
                            pubmed_str = "-" if not variant['pubmed_ids'] else ", ".join(variant['pubmed_ids'])
                            row = [
                                variant['variant_id'],
                                variant['position'],
                                variant['change'],
                                variant['description'],
                                pubmed_str,
                                variant['rsid'],
                                variant['link']
                            ]
                            table_data.append(row)
                        
                        output += tabulate(table_data, headers=headers, tablefmt="grid")
                        
                        st.download_button(
                            label="📋 Download Table",
                            data=output,
                            file_name=f"variant_table_{st.session_state.gene_symbol}.txt",
                            mime="text/plain"
                        )
                else:
                    st.warning(f"No variants found for gene: **{st.session_state.gene_symbol}**")
        else:
            st.error("Please search for a gene first using RefSeq accession or enter a gene symbol.")
    
    # Information panel
    with st.expander("ℹ️ About Variant Tables"):
        st.markdown("""
        **Variant Types:**
        - **SNV**: Single Nucleotide Variant (e.g., A/G)
        - **DEL**: Deletion variant
        - **DELINS**: Deletion-insertion variant
        
        **Position Format:**
        - Format: chromosome:chromosome:position
        - Example: 17:17:43126843 (chromosome 17, position 43126843)
        
        **Change Column:**
        - Shows reference/alternate alleles
        - Multiple alleles shown as A/C/G/T for complex variants
        
        **Links:**
        - Each RS ID links to dbSNP database
        - PubMed IDs link to research papers
        """)


if __name__ == "__main__":
    main()
