import streamlit as st
from Bio import SeqIO, Entrez
from tabulate import tabulate
from xml.etree import ElementTree as ET
import requests

# Set your email for Entrez
Entrez.email = "your.email@example.com"

def get_gene_symbol_from_accession(accession):
    """Extract gene name from nucleotide accession"""
    try:
        handle = Entrez.efetch(db="nucleotide", id=accession, rettype="gb", retmode="xml")
        records = Entrez.read(handle)
        handle.close()
        
        # Extract gene information from GenBank record
        for record in records:
            if 'GBSeq_feature-table' in record:
                for feature in record['GBSeq_feature-table']:
                    if feature['GBFeature_key'] == 'gene':
                        for qualifier in feature.get('GBFeature_quals', []):
                            if qualifier['GBQualifier_name'] == 'gene':
                                return qualifier['GBQualifier_value']
        return None
    except Exception as e:
        st.error(f"Error fetching gene info: {str(e)}")
        return None

def search_pubmed(term, max_results=10):
    """Search PubMed for relevant papers"""
    try:
        # Search for papers
        search_handle = Entrez.esearch(db="pubmed", term=term, retmax=max_results)
        search_results = Entrez.read(search_handle)
        search_handle.close()
        
        pmids = search_results["IdList"]
        if not pmids:
            return []
        
        # Fetch paper details
        fetch_handle = Entrez.efetch(db="pubmed", id=pmids, rettype="xml")
        papers = Entrez.read(fetch_handle)
        fetch_handle.close()
        
        results = []
        for paper in papers['PubmedArticle']:
            article = paper['MedlineCitation']['Article']
            
            # Extract title
            title = article.get('ArticleTitle', 'No title available')
            
            # Extract authors
            authors = []
            if 'AuthorList' in article:
                for author in article['AuthorList'][:3]:  # First 3 authors
                    if 'LastName' in author and 'ForeName' in author:
                        authors.append(f"{author['LastName']} {author['ForeName']}")
            author_str = ', '.join(authors) + (' et al.' if len(article.get('AuthorList', [])) > 3 else '')
            
            # Extract journal and year
            journal = article.get('Journal', {}).get('Title', 'Unknown Journal')
            pub_date = paper['MedlineCitation']['Article'].get('Journal', {}).get('JournalIssue', {}).get('PubDate', {})
            year = pub_date.get('Year', 'Unknown Year')
            
            # Extract abstract
            abstract = ""
            if 'Abstract' in article:
                abstract_texts = []
                for abstract_part in article['Abstract']['AbstractText']:
                    if isinstance(abstract_part, str):
                        abstract_texts.append(abstract_part)
                    else:
                        abstract_texts.append(str(abstract_part))
                abstract = ' '.join(abstract_texts)
            
            # Get PMID
            pmid = paper['MedlineCitation']['PMID']
            
            results.append({
                'pmid': pmid,
                'title': title,
                'authors': author_str,
                'journal': journal,
                'year': year,
                'abstract': abstract[:500] + '...' if len(abstract) > 500 else abstract
            })
        
        return results
    
    except Exception as e:
        st.error(f"Error searching PubMed: {str(e)}")
        return []

def display_papers_table(papers):
    """Display papers in a Markdown table using tabulate (like mu_table.py)"""
    if not papers:
        st.warning("No papers found.")
        return

    table_data = []
    headers = ["PMID", "Title", "Authors", "Journal (Year)", "Abstract"]
    for paper in papers:
        pmid_md = f"[{paper['pmid']}](https://pubmed.ncbi.nlm.nih.gov/{paper['pmid']})"
        title = paper['title']
        authors = paper['authors']
        journal_year = f"{paper['journal']} ({paper['year']})"
        abstract = paper['abstract']
        row = [
            pmid_md,
            title,
            authors,
            journal_year,
            abstract
        ]
        table_data.append(row)

    table_md = tabulate(table_data, headers=headers, tablefmt="github")
    st.markdown(table_md, unsafe_allow_html=True)

def main():
    st.set_page_config(page_title="PubMed Search by Accession", layout="wide")
    st.title("📚 PubMed Research Papers Search")
    st.subheader("Search research papers using NCBI accession numbers")
    
    # Sidebar inputs
    st.sidebar.header("Search Parameters")
    accession = st.sidebar.text_input("Enter NCBI Accession Number", 
                                     value="NM_007294", 
                                     help="e.g., NM_007294, NC_000017, OR905937")
    max_results = st.sidebar.slider("Maximum Results", min_value=5, max_value=50, value=20)
    
    # Additional search terms
    additional_terms = st.sidebar.text_input("Additional Search Terms (optional)", 
                                           help="e.g., mutation, variant, disease")
    
    if accession:
        with st.spinner("🔍 Searching for gene information and PubMed papers..."):
            # Get gene symbol from accession
            gene_symbol = get_gene_symbol_from_accession(accession)
            
            if gene_symbol:
                st.success(f"✅ Found gene: **{gene_symbol}** for accession **{accession}**")
                
                # Build search term
                search_term = gene_symbol
                if additional_terms:
                    search_term += f" AND {additional_terms}"
                
                st.info(f"Searching PubMed with term: **{search_term}**")
                
                # Search PubMed
                papers = search_pubmed(search_term, max_results)
                
                if papers:
                    st.subheader(f"📄 Found {len(papers)} Research Papers")
                    display_papers_table(papers)
                    
                    # Download option
                    if st.button("📥 Download Results as Text"):
                        text_content = f"PubMed Search Results for {gene_symbol} ({accession})\n"
                        text_content += "=" * 60 + "\n\n"
                        
                        for i, paper in enumerate(papers, 1):
                            text_content += f"{i}. {paper['title']}\n"
                            text_content += f"   Authors: {paper['authors']}\n"
                            text_content += f"   Journal: {paper['journal']} ({paper['year']})\n"
                            text_content += f"   PMID: {paper['pmid']}\n"
                            text_content += f"   Abstract: {paper['abstract']}\n\n"
                        
                        st.download_button(
                            label="💾 Download",
                            data=text_content,
                            file_name=f"pubmed_results_{gene_symbol}_{accession}.txt",
                            mime="text/plain"
                        )
                else:
                    st.warning(f"No papers found for gene: **{gene_symbol}**")
            else:
                st.warning(f"Could not extract gene information from accession: **{accession}**")
                
                # Try searching with accession directly
                st.info("Trying direct search with accession number...")
                papers = search_pubmed(accession, max_results)
                
                if papers:
                    st.subheader(f"📄 Found {len(papers)} Research Papers")
                    display_papers_table(papers)
                else:
                    st.error("No papers found with this accession number.")

if __name__ == "__main__":
    main()
