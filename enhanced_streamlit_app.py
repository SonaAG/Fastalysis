"""
Enhanced Streamlit Frontend for Genomics Research Assistant
Integrates with FastAPI backend and agent system
"""

import streamlit as st
import requests
import json
from typing import Dict, Optional
import plotly.express as px
import pandas as pd
from Bio import Entrez
from xml.etree import ElementTree as ET
from tabulate import tabulate

# Configuration
API_BASE_URL = "http://localhost:8000"
Entrez.email = "your.email@example.com"  # Replace with your real email

st.set_page_config(
    page_title="🧬 Genomics Research Assistant",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="expanded"
)

def init_session_state():
    """Initialize session state variables"""
    if "messages" not in st.session_state:
        st.session_state.messages = []
    if "analysis_results" not in st.session_state:
        st.session_state.analysis_results = {}
    if "current_sequence" not in st.session_state:
        st.session_state.current_sequence = ""

def call_api(endpoint: str, data: Dict) -> Dict:
    """Make API calls to FastAPI backend"""
    try:
        response = requests.post(f"{API_BASE_URL}/{endpoint}", json=data)
        if response.ok:
            return response.json()
        else:
            return {"status": "error", "message": f"API error: {response.status_code}"}
    except Exception as e:
        return {"status": "error", "message": str(e)}

# Variant lookup functions from mu_table.py
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
    """Fetch dbSNP variant summaries for the given gene symbol."""
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
        st.error(f"Error fetching variants for gene {gene_symbol}: {e}")
        return []

def display_variants_table(variants):
    """Display variants in a formatted table"""
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

# Literature search functions from pubmed.py
def search_pubmed(term, max_results=10):
    """Search PubMed for relevant papers"""
    try:
        handle = Entrez.esearch(db="pubmed", term=term, retmax=max_results)
        record = Entrez.read(handle)
        handle.close()
        return record.get("IdList", [])
    except Exception as e:
        st.error(f"Error searching PubMed: {e}")
        return []

def fetch_pubmed_details(pmids):
    """Fetch detailed information from PubMed using XML format for more reliable parsing"""
    if not pmids:
        return []
    
    try:
        # Use XML format for more reliable parsing
        handle = Entrez.efetch(db="pubmed", id=",".join(pmids), rettype="xml")
        records = Entrez.read(handle)
        handle.close()
        
        papers = []
        
        # Parse XML records
        for record in records['PubmedArticle']:
            try:
                article = record['MedlineCitation']['Article']
                pmid = str(record['MedlineCitation']['PMID'])
                
                # Get title
                title = article.get('ArticleTitle', 'No title available')
                
                # Get authors
                authors = []
                if 'AuthorList' in article and article['AuthorList']:
                    for author in article['AuthorList'][:3]:  # Limit to first 3 authors
                        if 'LastName' in author and 'ForeName' in author:
                            authors.append(f"{author['LastName']} {author['ForeName']}")
                        elif 'LastName' in author:
                            authors.append(author['LastName'])
                
                author_str = ', '.join(authors)
                if len(article.get('AuthorList', [])) > 3:
                    author_str += ' et al.'
                if not author_str:
                    author_str = 'Not available'
                
                # Get journal info
                journal = 'Not available'
                if 'Journal' in article:
                    journal_info = article['Journal']
                    journal_title = journal_info.get('Title', journal_info.get('ISOAbbreviation', ''))
                    year = 'Unknown'
                    if 'JournalIssue' in journal_info and 'PubDate' in journal_info['JournalIssue']:
                        pub_date = journal_info['JournalIssue']['PubDate']
                        year = pub_date.get('Year', 'Unknown')
                    journal = f"{journal_title} ({year})"
                
                # Get abstract
                abstract = 'No abstract available'
                if 'Abstract' in article and 'AbstractText' in article['Abstract']:
                    abstract_texts = article['Abstract']['AbstractText']
                    if isinstance(abstract_texts, list):
                        abstract = ' '.join([str(text) for text in abstract_texts])
                    else:
                        abstract = str(abstract_texts)
                
                papers.append({
                    'pmid': pmid,
                    'title': title,
                    'authors': author_str,
                    'journal': journal,
                    'abstract': abstract
                })
                
            except Exception as e:
                # Skip this record if there's an error parsing it
                continue
                
        return papers
        
    except Exception as e:
        st.error(f"Error fetching PubMed details: {e}")
        return []

def display_chat_interface():
    """Enhanced chat interface with context awareness"""
    st.header("🤖 Genomics Research Assistant")
    
    # Get available models
    try:
        models_response = requests.get(f"{API_BASE_URL}/models")
        if models_response.ok:
            available_models = models_response.json()["models"]
        else:
            available_models = ["kimi-k2", "llama-3.1-8b-instant"]
    except:
        available_models = ["kimi-k2", "llama-3.1-8b-instant"]
    
    # Model selection
    selected_model = st.selectbox(
        "Choose AI Model", 
        available_models,
        help="Select the AI model for conversations"
    )
    
    # Chat input with file context
    if st.session_state.current_sequence:
        seq_preview = st.session_state.current_sequence[:50] + "..." if len(st.session_state.current_sequence) > 50 else st.session_state.current_sequence
        st.info(f"📄 Current sequence available: {seq_preview} ({len(st.session_state.current_sequence)} chars)")
    
    user_input = st.chat_input("Ask me about genomics, genes, mutations, or request analysis...")
    
    # Quick action buttons
    if st.session_state.current_sequence:
        col1, col2, col3 = st.columns(3)
        with col1:
            if st.button("🤖 Analyze my sequence", key="chat_analyze"):
                user_input = f"Please analyze this sequence and identify what it is"
        with col2:
            if st.button("🔍 BLAST search", key="chat_blast"):
                user_input = f"Run a BLAST search on this sequence to find similar sequences"
        with col3:
            if st.button("📚 Find literature", key="chat_literature"):
                user_input = "Find relevant literature for this sequence and explain its significance"
    
    if user_input:
        # Add user message
        st.session_state.messages.append({"role": "user", "content": user_input})
        
        # Prepare context
        context = {}
        if st.session_state.current_sequence:
            context["sequence"] = st.session_state.current_sequence
        if st.session_state.analysis_results:
            context["previous_results"] = st.session_state.analysis_results
        
        # Call chat API with context
        chat_request = {
            "message": user_input,
            "model": selected_model,
            "context": context
        }
        
        with st.spinner("🧠 AI is thinking..."):
            response = call_api("chat", chat_request)
        
        if response["status"] == "success":
            chat_response = response["response"]
            
            # Handle different response types
            if isinstance(chat_response, dict):
                if "message" in chat_response:
                    st.session_state.messages.append({
                        "role": "assistant", 
                        "content": chat_response["message"]
                    })
                    
                    # If analysis was performed, also store the results
                    if chat_response.get("status") == "analysis_complete" and "analysis_data" in chat_response:
                        # Merge new analysis results
                        for analysis_type, data in chat_response["analysis_data"].items():
                            st.session_state.analysis_results[analysis_type] = {"status": "success", "data": data}
                        st.rerun()  # Refresh to show new results
                else:
                    st.session_state.messages.append({
                        "role": "assistant", 
                        "content": str(chat_response)
                    })
            else:
                st.session_state.messages.append({
                    "role": "assistant", 
                    "content": str(chat_response)
                })
        else:
            st.error(f"Chat error: {response.get('message', 'Unknown error')}")
    
    # Display conversation
    for message in st.session_state.messages:
        with st.chat_message(message["role"]):
            if isinstance(message["content"], str):
                st.markdown(message["content"])
            else:
                st.json(message["content"])

def parse_fasta_file(uploaded_file) -> str:
    """Parse uploaded FASTA file and extract sequence"""
    try:
        # Read file content
        file_content = uploaded_file.read().decode('utf-8')
        
        # Parse FASTA format
        lines = file_content.strip().split('\n')
        sequences = []
        current_sequence = ""
        headers = []
        
        for line in lines:
            line = line.strip()
            if line.startswith('>'):
                # Save previous sequence if exists
                if current_sequence:
                    sequences.append(current_sequence)
                    current_sequence = ""
                headers.append(line[1:])  # Remove '>'
            elif line:
                # Remove any whitespace and join sequence lines
                current_sequence += line.replace(' ', '').replace('\t', '')
        
        # Add last sequence
        if current_sequence:
            sequences.append(current_sequence)
        
        if not sequences:
            st.error("No valid sequences found in FASTA file")
            return ""
        
        # If multiple sequences, show selection
        if len(sequences) > 1:
            st.info(f"📄 Found {len(sequences)} sequences in file")
            selected_seq_idx = st.selectbox(
                "Choose sequence:",
                range(len(sequences)),
                format_func=lambda x: f"Sequence {x+1}: {headers[x][:50]}..." if x < len(headers) else f"Sequence {x+1}"
            )
            return sequences[selected_seq_idx]
        else:
            return sequences[0]
            
    except UnicodeDecodeError:
        st.error("File encoding error. Please ensure your FASTA file is in UTF-8 format.")
        return ""
    except Exception as e:
        st.error(f"Error parsing FASTA file: {str(e)}")
        return ""

def display_analysis_panel():
    """Analysis tools panel"""
    st.header("🔬 Analysis Tools")
    
    # Sequence input options
    st.subheader("Sequence Input")
    
    # File upload or text input tabs
    input_tab1, input_tab2 = st.tabs(["📁 Upload FASTA File", "📝 Text Input"])
    
    with input_tab1:
        uploaded_file = st.file_uploader(
            "Choose a FASTA file",
            type=['fasta', 'fa', 'fas', 'txt'],
            help="Upload a FASTA file containing your sequence"
        )
        
        if uploaded_file is not None:
            # Display file info
            st.info(f"📄 File: {uploaded_file.name} ({uploaded_file.size} bytes)")
            
            # Parse FASTA file
            sequence_from_file = parse_fasta_file(uploaded_file)
            
            if sequence_from_file:
                st.session_state.current_sequence = sequence_from_file
                
                # Show sequence preview
                with st.expander("🔍 Sequence Preview"):
                    preview_length = min(200, len(sequence_from_file))
                    st.code(f"{sequence_from_file[:preview_length]}{'...' if len(sequence_from_file) > preview_length else ''}")
                    st.caption(f"Total length: {len(sequence_from_file)} nucleotides/amino acids")
    
    with input_tab2:
        sequence_input = st.text_area(
            "Enter DNA/Protein Sequence",
            value=st.session_state.current_sequence,
            height=100,
            help="Paste your sequence here (FASTA format or raw sequence)"
        )
        
        if sequence_input != st.session_state.current_sequence:
            st.session_state.current_sequence = sequence_input
    
    # Show current sequence status
    if st.session_state.current_sequence:
        seq_length = len(st.session_state.current_sequence)
        st.success(f"✅ Sequence loaded ({seq_length} characters)")
        
        # Detect sequence type
        seq_upper = st.session_state.current_sequence.upper()
        if set(seq_upper).issubset(set("ATCGNU")):
            seq_type = "🧬 DNA/RNA"
        elif set(seq_upper).issubset(set("ACDEFGHIKLMNPQRSTVWY")):
            seq_type = "🧪 Protein"
        else:
            seq_type = "❓ Unknown"
        
        st.info(f"Detected type: {seq_type}")
    else:
        st.warning("⚠️ Please upload a FASTA file or enter a sequence")
    
    if st.session_state.current_sequence:
        # Analysis buttons
        col1, col2, col3, col4 = st.columns(4)
        
        with col1:
            if st.button("🔍 BLAST Search", use_container_width=True):
                with st.spinner("Running BLAST..."):
                    result = call_api("blast", {
                        "sequence": st.session_state.current_sequence,
                        "num_hits": 5
                    })
                    st.session_state.analysis_results["blast"] = result
        
        with col2:
            if st.button("🧬 Mutation Analysis", use_container_width=True):
                with st.spinner("Analyzing mutations..."):
                    result = call_api("mutation", {
                        "sequence": st.session_state.current_sequence,
                        "hit_index": 0
                    })
                    st.session_state.analysis_results["mutation"] = result
        
        with col3:
            st.markdown("**Genetic Variants Lookup**")
            gene_input = st.text_input("Gene Symbol/RefSeq", placeholder="e.g., BRCA1 or NM_007294")
            is_refseq = st.checkbox("Input is RefSeq accession", key="refseq_checkbox")
            
            if st.button("📊 Variants", use_container_width=True) and gene_input:
                with st.spinner("Looking up variants..."):
                    try:
                        gene_symbol = gene_input
                        # If RefSeq accession, convert to gene symbol first
                        if is_refseq:
                            with st.spinner("Converting RefSeq to gene symbol..."):
                                gene_symbol = get_gene_symbol_from_refseq(gene_input)
                                if gene_symbol:
                                    st.success(f"Found gene: **{gene_symbol}** for accession **{gene_input}**")
                                else:
                                    st.error("Could not extract gene information from accession")
                                    gene_symbol = None
                        
                        if gene_symbol:
                            variants = fetch_dbsnp_variants_for_gene(gene_symbol, max_records=20)
                            if variants:
                                result = {"status": "success", "data": {"variants": variants, "gene_symbol": gene_symbol}}
                            else:
                                result = {"status": "error", "message": f"No variants found for gene: {gene_symbol}"}
                        else:
                            result = {"status": "error", "message": "Could not determine gene symbol"}
                    except Exception as e:
                        result = {"status": "error", "message": str(e)}
                    st.session_state.analysis_results["variants"] = result
        
        with col4:
            if st.button("📚 Literature", use_container_width=True) and gene_input:
                with st.spinner("Searching literature..."):
                    try:
                        # Search PubMed for the gene input
                        search_term = gene_input
                        
                        pmids = search_pubmed(search_term, max_results=10)
                        
                        if pmids:
                            papers = fetch_pubmed_details(pmids[:5])  # Limit to 5 for display
                            
                            if papers:
                                result = {"status": "success", "data": {"papers": papers, "search_term": search_term}}
                            else:
                                result = {"status": "error", "message": "Found PMIDs but could not retrieve paper details"}
                        else:
                            result = {"status": "error", "message": f"No papers found for: {search_term}"}
                    except Exception as e:
                        result = {"status": "error", "message": f"Search error: {str(e)}"}
                    st.session_state.analysis_results["literature"] = result

def display_results():
    """Display analysis results"""
    if st.session_state.analysis_results:
        st.header("📈 Analysis Results")
        
        # Create tabs for different results
        result_types = list(st.session_state.analysis_results.keys())
        tabs = st.tabs([f"🔍 {rt.title()}" for rt in result_types])
        
        for i, (result_type, tab) in enumerate(zip(result_types, tabs)):
            with tab:
                result_data = st.session_state.analysis_results[result_type]
                
                if result_data.get("status") == "success":
                    data = result_data["data"]
                    
                    if result_type == "blast":
                        display_blast_results(data)
                    elif result_type == "mutation":
                        display_mutation_results(data)
                    elif result_type == "variants":
                        display_variant_results(data)
                    elif result_type == "literature":
                        display_literature_results(data)
                else:
                    st.error(f"Analysis failed: {result_data.get('message', 'Unknown error')}")

def display_blast_results(data: Dict):
    """Display BLAST results"""
    if "hits" in data and data["hits"]:
        st.subheader(f"🎯 BLAST Hits ({len(data['hits'])})")
        
        # Create DataFrame for better display
        hits_df = pd.DataFrame([
            {
                "Rank": i+1,
                "Accession": hit.get("accession", hit.get("Accession", "N/A")),
                "Title": (hit.get("title", hit.get("Title", ""))[:100] + "...") if len(hit.get("title", hit.get("Title", ""))) > 100 else hit.get("title", hit.get("Title", "")),
                "Identity %": hit.get("percent_identity", hit.get("Identity (%)", 0)),
                "E-value": hit.get("evalue", hit.get("E-value", 0)),
                "Score": hit.get("score", hit.get("Score", 0)),
                "Coverage %": hit.get("query_coverage", hit.get("Coverage (%)", 0))
            }
            for i, hit in enumerate(data["hits"])
        ])
        
        st.dataframe(hits_df, use_container_width=True)
        
        # Visualization - only if we have valid data
        valid_hits = hits_df[hits_df["Identity %"] > 0]
        if len(valid_hits) > 1:
            fig = px.bar(
                valid_hits, 
                x="Rank", 
                y="Identity %",
                title="BLAST Hit Identity Scores",
                color="Identity %",
                color_continuous_scale="viridis"
            )
            st.plotly_chart(fig, use_container_width=True)
        elif len(valid_hits) == 0:
            st.warning("⚠️ BLAST search completed but found no significant matches")
        else:
            st.info("Only one hit found - visualization needs multiple hits")
            
    else:
        st.error("❌ No BLAST hits found or BLAST search failed")
        st.info("Try with a different sequence or check if the sequence is valid")

def display_mutation_results(data: Dict):
    """Display mutation analysis results with rich visualizations and downloads"""
    if "mutations" in data:
        # Header with reference information
        st.subheader("🧬 Mutation Analysis Dashboard")
        
        # Reference sequence information
        if "ref_id" in data and "ref_title" in data:
            with st.expander("📍 Reference Sequence Information", expanded=True):
                st.info(f"**Accession:** {data['ref_id']}")
                st.info(f"**Description:** {data['ref_title']}")
                if "alignment" in data and "ref_start" in data["alignment"]:
                    st.info(f"**Alignment Region:** {data['alignment']['ref_start']} - {data['alignment']['ref_end']}")
        
        # Alignment statistics dashboard
        if "alignment" in data:
            st.subheader("📊 Alignment Quality Metrics")
            align_data = data["alignment"]
            
            col1, col2, col3, col4 = st.columns(4)
            with col1:
                st.metric("Length", f"{align_data.get('length', 'N/A')} bp", help="Alignment length")
            with col2:
                identity_val = align_data.get('identity', 0)
                total_len = align_data.get('length', 1)
                st.metric("Matches", f"{identity_val}/{total_len}", help="Identical positions")
            with col3:
                pct_id = align_data.get('percent_identity', 0)
                st.metric("Identity", f"{pct_id:.1f}%", help="Percentage identity")
            with col4:
                coverage = align_data.get('query_coverage', 0)
                st.metric("Coverage", f"{coverage:.1f}%", help="Query sequence coverage")
            
            # Quality indicator
            if pct_id >= 95:
                st.success("🎯 Excellent alignment quality")
            elif pct_id >= 85:
                st.info("✅ Good alignment quality") 
            else:
                st.warning("⚠️ Moderate alignment quality")
        
        # Add prominent sequence alignment visualization button
        if ("alignment" in data and "ref_aligned" in data["alignment"] 
            and "user_aligned" in data["alignment"]):
            
            st.markdown("---")
            st.subheader("🧬 Sequence Alignment Visualization")
            
            col1, col2 = st.columns(2)
            with col1:
                if st.button("🧬 View Detailed Alignment", key="main_alignment_btn", use_container_width=True):
                    st.session_state.show_alignment = True
            
            # Show alignment if button was clicked
            if st.session_state.get('show_alignment', False):
                with st.expander("🧬 Detailed Sequence Alignment Visualization", expanded=True):
                    alignment_html = generate_emboss_alignment_html(data)
                    st.components.v1.html(
                        f"""
                        <div style="font-family:monospace; font-size:15px; background: white; padding: 20px; border: 1px solid #ddd; border-radius: 5px;">
                        {alignment_html}
                        </div>
                        """,
                        height=600,
                        scrolling=True
                    )
        else:
            # Debug information - show what data is available
            st.markdown("---")
            st.error("🔍 **Debug Info:** Alignment visualization not available")
            st.write("**Available data keys:**", list(data.keys()))
            if "alignment" in data:
                st.write("**Alignment keys:**", list(data["alignment"].keys()))
            else:
                st.write("No 'alignment' key found in data")
            
            # Add basic alignment button even without detailed data
            if "alignment" in data:
                st.subheader("🧬 Basic Alignment View")
                if st.button("📊 Show Available Alignment Data", key="basic_alignment_btn"):
                    st.json(data["alignment"])
        
        # Mutation results
        if data["mutations"]:
            st.subheader(f"🔬 Detected Mutations ({len(data['mutations'])})")
            
            # Summary metrics
            total_mutations = len(data['mutations'])
            col1, col2 = st.columns(2)
            with col1:
                st.metric("Total Mutations", total_mutations)
            with col2:
                if "alignment" in data:
                    mutation_rate = (total_mutations / data['alignment'].get('length', 1)) * 100
                    st.metric("Mutation Rate", f"{mutation_rate:.2f}%")
            
            # Mutations table
            mutations_df = pd.DataFrame(data["mutations"])
            
            # Format the dataframe for better display
            if 'position' in mutations_df.columns:
                mutations_df = mutations_df.sort_values('position')
            
            st.dataframe(
                mutations_df,
                use_container_width=True,
                column_config={
                    "position": st.column_config.NumberColumn("Position", help="Position in alignment"),
                    "ref": st.column_config.TextColumn("Reference", help="Reference nucleotide/amino acid"),
                    "user": st.column_config.TextColumn("Query", help="Query sequence nucleotide/amino acid"),
                    "mutation": st.column_config.TextColumn("Mutation", help="Mutation notation")
                }
            )
            
            # Visualization tabs
            st.subheader("📈 Interactive Visualizations")
            viz_tab1, viz_tab2 = st.tabs(["🎯 Mutation Scatter", "📊 Mutation Track"])
            
            with viz_tab1:
                st.markdown("**Interactive Mutation Scatter Plot**")
                st.markdown("*Hover over points to see mutation details*")
                try:
                    import os
                    if "mutation_scatter_html" in data and os.path.exists(data["mutation_scatter_html"]):
                        with open(data["mutation_scatter_html"], 'r', encoding='utf-8') as f:
                            scatter_html = f.read()
                        st.components.v1.html(scatter_html, height=500, scrolling=True)
                    else:
                        st.warning("Scatter plot visualization not available")
                except Exception as e:
                    st.error(f"Error loading scatter plot: {e}")
            
            with viz_tab2:
                st.markdown("**Mutation Track Diagram**") 
                st.markdown("*Bar chart showing mutation positions and types*")
                try:
                    import os
                    if "mutation_track_html" in data and os.path.exists(data["mutation_track_html"]):
                        with open(data["mutation_track_html"], 'r', encoding='utf-8') as f:
                            track_html = f.read()
                        st.components.v1.html(track_html, height=400, scrolling=True)
                    else:
                        st.warning("Track diagram not available")
                except Exception as e:
                    st.error(f"Error loading track diagram: {e}")
            
            # Document Download Section (like mutation.py)
            st.subheader("📄 Download Alignment Document")
            
            # Generate HTML document with same format as mutation.py
            alignment_html_doc = generate_emboss_alignment_html_doc(data)
            
            # Create mutation table HTML (same format as mutation.py)
            mutation_table_html = ""
            if "mutations" in data and data["mutations"]:
                mutation_table_html = "<h3>Detected Mutations (Ref → User)</h3><table border='1' cellpadding='4' cellspacing='0'><tr><th>Position</th><th>Ref</th><th>User</th><th>Mutation</th></tr>"
                for mut in data["mutations"]:
                    pos = mut.get("position", "N/A")
                    ref = mut.get("ref", "N/A")
                    user = mut.get("user", "N/A")
                    mutation_str = f"c.{pos}{ref}&gt;{user}"
                    mutation_table_html += f"<tr><td>{pos}</td><td>{ref}</td><td>{user}</td><td>{mutation_str}</td></tr>"
                mutation_table_html += "</table>"
            else:
                mutation_table_html = "<h3>Detected Mutations</h3><p>No mutations detected - perfect sequence match!</p>"
            
            # Create HTML document (exact same format as mutation.py)
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
            <b>Ref:</b> {data.get('ref_id', 'Reference')}<br>
            <b>User Seq:</b> {data.get('user_id', 'User Sequence')}<br>
            <b>Identity:</b> {data.get('percent_identity', 'N/A')}%<br>
            <b>Coverage:</b> {data.get('query_coverage', 'N/A')}%<br>
            <b>E-value:</b> {data.get('evalue', 'N/A')}<br>
            <h3>Alignment Visualization</h3>
            <div style="font-family:monospace; font-size:16px;">
            {alignment_html_doc}
            </div>
            {mutation_table_html}
            </body>
            </html>
            """
            
            st.download_button(
                label="⬇️ Download Alignment (.doc)",
                data=html_doc,
                file_name=f"alignment_{data.get('ref_id', 'reference')}.doc",
                mime="application/msword"
            )
    
    else:
        st.info("No mutation analysis results to display")

def generate_mutation_report(data: Dict) -> str:
    """Generate a comprehensive text report of mutation analysis"""
    report_lines = []
    report_lines.append("=" * 60)
    report_lines.append("MUTATION ANALYSIS REPORT")
    report_lines.append("=" * 60)
    report_lines.append("")
    
    # Header information
    if "user_id" in data:
        report_lines.append(f"Query Sequence ID: {data['user_id']}")
    if "ref_id" in data and "ref_title" in data:
        report_lines.append(f"Reference Sequence: {data['ref_id']}")
        report_lines.append(f"Description: {data['ref_title']}")
    report_lines.append("")
    
    # Alignment statistics
    if "alignment" in data:
        align = data["alignment"]
        report_lines.append("ALIGNMENT STATISTICS")
        report_lines.append("-" * 20)
        report_lines.append(f"Alignment Length: {align.get('length', 'N/A')} bp")
        report_lines.append(f"Identical Positions: {align.get('identity', 'N/A')}")
        report_lines.append(f"Percent Identity: {align.get('percent_identity', 0):.2f}%")
        report_lines.append(f"Query Coverage: {align.get('query_coverage', 0):.2f}%")
        report_lines.append(f"Alignment Score: {align.get('score', 'N/A')}")
        report_lines.append(f"E-value: {align.get('evalue', 'N/A')}")
        if "ref_start" in align and "ref_end" in align:
            report_lines.append(f"Reference Region: {align['ref_start']}-{align['ref_end']}")
        report_lines.append("")
    
    # Mutation summary
    mutations = data.get("mutations", [])
    report_lines.append(f"MUTATION SUMMARY")
    report_lines.append("-" * 15) 
    report_lines.append(f"Total Mutations Detected: {len(mutations)}")
    
    if mutations:
        report_lines.append("")
        report_lines.append("DETAILED MUTATION LIST")
        report_lines.append("-" * 22)
        report_lines.append("Position\tRef\tQuery\tMutation")
        report_lines.append("-" * 40)
        
        for mut in mutations:
            pos = mut.get('position', 'N/A')
            ref = mut.get('ref', '?')
            user = mut.get('user', '?')  
            notation = mut.get('mutation', f"{ref}>{user}")
            report_lines.append(f"{pos}\t{ref}\t{user}\t{notation}")
    else:
        report_lines.append("No mutations detected - perfect sequence match!")
    
    report_lines.append("")
    report_lines.append("=" * 60)
    report_lines.append("End of Report")
    report_lines.append("=" * 60)
    
    return "\n".join(report_lines)

def display_variant_results(data: Dict):
    """Display variant lookup results with same styling as mu_table.py"""
    if "variants" in data and data["variants"]:
        gene_symbol = data.get("gene_symbol", "Unknown Gene")
        st.subheader(f"📊 Genetic Variants for {gene_symbol} ({len(data['variants'])})")
        
        # Use the same table display function as mu_table.py
        display_variants_table(data["variants"])
        
        # Add download option
        if st.button("💾 Download Variant Table"):
            output = f"Variant Table for {gene_symbol}\n"
            output += "=" * 80 + "\n\n"
            
            table_data = []
            headers = ["Variant ID", "Position", "Change", "Description", "PubMed IDs", "rsID", "Link"]
            
            for variant in data["variants"]:
                pubmed_str = "-" if not variant.get('PubMed IDs') else ", ".join(variant['PubMed IDs'])
                row = [
                    variant.get('Variant ID', ''),
                    variant.get('Position', ''),
                    variant.get('Change', ''),
                    variant.get('Description', ''),
                    pubmed_str,
                    variant.get('rsID', ''),
                    variant.get('Link', '')
                ]
                table_data.append(row)
            
            output += tabulate(table_data, headers=headers, tablefmt="grid")
            
            st.download_button(
                label="📋 Download Table",
                data=output,
                file_name=f"variant_table_{gene_symbol}.txt",
                mime="text/plain"
            )
            
        # Information panel
        with st.expander("ℹ️ About Variant Tables"):
            st.markdown("""
            **Variant Types:**
            - **SNV**: Single Nucleotide Variant (e.g., A/G)
            - **DEL**: Deletion variant
            - **DELINS**: Deletion-insertion variant
            
            **Position Format:**
            - Format: chromosome:position
            - Example: 17:43126843 (chromosome 17, position 43126843)
            
            **Change Column:**
            - Shows reference/alternate alleles
            - Multiple alleles shown as A/C/G/T for complex variants
            
            **Links:**
            - Each RS ID links to dbSNP database
            - PubMed IDs link to research papers
            """)
    else:
        st.info("No variants found or error in lookup")

def display_literature_results(data: Dict):
    """Display literature search results with clickable PMID links"""
    if "papers" in data and data["papers"]:
        search_term = data.get("search_term", "")
        st.subheader(f"📚 Literature ({len(data['papers'])} papers)")
        
        if search_term:
            st.info(f"🔎 Search term: {search_term}")
        
        # Display papers in a table format similar to the screenshot
        st.markdown("### Research Papers Table")
        
        # Create table headers
        table_data = []
        headers = ["PMID", "Title", "Authors", "Journal (Year)", "Abstract"]
        
        for paper in data["papers"]:
            pmid = paper.get('pmid', 'N/A')
            title = paper.get('title', 'No title')[:100] + "..." if len(paper.get('title', '')) > 100 else paper.get('title', 'No title')
            authors = paper.get('authors', 'Not available')
            journal = paper.get('journal', 'Not available')
            abstract = (paper.get('abstract', 'No abstract available')[:200] + "...") if len(paper.get('abstract', '')) > 200 else paper.get('abstract', 'No abstract available')
            
            # Create clickable PMID link
            pmid_link = f"[{pmid}](https://pubmed.ncbi.nlm.nih.gov/{pmid}/)" if pmid != 'N/A' else pmid
            
            table_data.append([
                pmid_link,
                title,
                authors,
                journal,
                abstract
            ])
        
        # Display as markdown table
        table_md = tabulate(table_data, headers=headers, tablefmt="github")
        st.markdown(table_md, unsafe_allow_html=True)
        
        # Detailed expandable view for each paper
        st.markdown("### Detailed Paper Information")
        for i, paper in enumerate(data["papers"]):
            pmid = paper.get('pmid', 'N/A')
            pmid_link = f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/" if pmid != 'N/A' else "#"
            
            with st.expander(f"Paper {i+1}: {paper.get('title', 'No title')}"):
                col1, col2 = st.columns([1, 3])
                
                with col1:
                    st.write("**PMID:**")
                    if pmid != 'N/A':
                        st.markdown(f"[{pmid}]({pmid_link}) 🔗", unsafe_allow_html=True)
                    else:
                        st.write("Not available")
                    
                    st.write("**Authors:**")
                    st.write(paper.get('authors', 'Not available'))
                    
                    st.write("**Journal:**")
                    st.write(paper.get('journal', 'Not available'))
                
                with col2:
                    if paper.get('abstract'):
                        st.write("**Abstract:**")
                        st.write(paper['abstract'])
                    else:
                        st.write("**Abstract:** Not available")
        
        # Download option
        if st.button("💾 Download Literature Results"):
            output = f"Literature Search Results\n"
            output += f"Search term: {search_term}\n"
            output += "=" * 80 + "\n\n"
            
            for i, paper in enumerate(data["papers"]):
                pmid = paper.get('pmid', 'N/A')
                output += f"Paper {i+1}:\n"
                output += f"Title: {paper.get('title', 'No title')}\n"
                output += f"Authors: {paper.get('authors', 'Not available')}\n"
                output += f"Journal: {paper.get('journal', 'Not available')}\n"
                output += f"PMID: {pmid}\n"
                if pmid != 'N/A':
                    output += f"Link: https://pubmed.ncbi.nlm.nih.gov/{pmid}/\n"
                if paper.get('abstract'):
                    output += f"Abstract: {paper['abstract']}\n"
                output += "\n" + "-" * 60 + "\n\n"
            
            st.download_button(
                label="📋 Download Results",
                data=output,
                file_name=f"literature_search_results.txt",
                mime="text/plain"
            )
            
    else:
        st.info("No literature found")

def generate_emboss_alignment_html(data: Dict) -> str:
    """Generate EMBOSS-style HTML alignment visualization for Streamlit display (exactly like mutation.py emboss_style_html)"""
    if not ("alignment" in data and "ref_aligned" in data["alignment"] and "user_aligned" in data["alignment"]):
        return "<p>Alignment data not available</p>"
    
    aln_seqA = data["alignment"]["ref_aligned"]  # ref sequence 
    aln_seqB = data["alignment"]["user_aligned"]  # user sequence
    ref_id = data.get("ref_id", "Reference")
    user_id = data.get("user_id", "user_sequence")
    
    # Use exact parameters and labels like mutation.py
    name1 = f"Ref ({ref_id})"
    name2 = "User Seq"
    matrix = "EBLOSUM62"
    gap_open = 10.0
    gap_extend = 0.5
    label1 = name1
    label2 = name2
    
    # Calculate alignment statistics (exactly like mutation.py)
    def calculate_stats(aln_seqA, aln_seqB):
        length = len(aln_seqA)
        identity = sum(a == b and a != '-' for a, b in zip(aln_seqA, aln_seqB))
        gaps = sum(a == '-' or b == '-' for a, b in zip(aln_seqA, aln_seqB))
        similarity = identity  # Simplified similarity
        return length, identity, similarity, gaps
    
    length, identity, similarity, gaps = calculate_stats(aln_seqA, aln_seqB)
    
    # Start building HTML output (exactly matching mutation.py emboss_style_html)
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
    
    # Process alignment in blocks (exactly like mutation.py emboss_style_html)
    line_width = 60
    ref_pos = user_pos = 1
    label_width = max(len(label1), len(label2)) + 2
    
    for i in range(0, length, line_width):
        ref_block = aln_seqA[i:i+line_width]
        user_block = aln_seqB[i:i+line_width]
        match_line = ""
        
        # Generate colored sequences and match line (exactly matching mutation.py emboss_style_html)
        ref_colored = ""
        user_colored = ""
        
        for a, b in zip(ref_block, user_block):
            if a == b and a != '-':
                # Perfect match
                match_line += "|"
                ref_colored += a  # Plain text for matches
                user_colored += b  # Plain text for matches
            else:
                match_line += " "
                if a != '-' and b != '-':
                    # Substitution - exact colors from mutation.py
                    ref_colored += f"<span style='background-color:#f88'>{a}</span>"
                    user_colored += f"<span style='background-color:#8f8'>{b}</span>"
                elif a == '-' and b != '-':
                    # Insertion in user sequence
                    ref_colored += "-"  # Plain dash for gaps
                    user_colored += f"<span style='background-color:#8ff'>{b}</span>"
                elif a != '-' and b == '-':
                    # Deletion in user sequence  
                    ref_colored += f"<span style='background-color:#fbb'>{a}</span>"
                    user_colored += "-"  # Plain dash for gaps
        
        # Count residues for position tracking (exactly like mutation.py)
        def count_residues(seq):
            return sum(1 for c in seq if c != '-')
        
        ref_end = ref_pos + count_residues(ref_block) - 1
        user_end = user_pos + count_residues(user_block) - 1
        
        # Format exactly like mutation.py emboss_style_html
        html_output += (
            f"{label2:<{label_width}}{user_pos:>4} {user_colored} {user_end}\n"
            f"{'':<{label_width}}     {match_line}\n"
            f"{label1:<{label_width}}{ref_pos:>4} {ref_colored} {ref_end}\n\n"
        )
        
        ref_pos = ref_end + 1
        user_pos = user_end + 1
    
    html_output += "</pre>"
    return html_output

def generate_emboss_alignment_html_doc(data: Dict) -> str:
    """Generate EMBOSS-style HTML alignment for document download (exactly like mutation.py emboss_style_html_doc)"""
    if not ("alignment" in data and "ref_aligned" in data["alignment"] and "user_aligned" in data["alignment"]):
        return "<p>Alignment data not available</p>"
    
    aln_seqA = data["alignment"]["ref_aligned"]  # ref sequence 
    aln_seqB = data["alignment"]["user_aligned"]  # user sequence
    ref_id = data.get("ref_id", "Reference")
    user_id = data.get("user_id", "user_sequence")
    
    # Use exact parameters and labels like mutation.py
    name1 = f"Ref ({ref_id})"
    name2 = "User Seq"
    matrix = "EBLOSUM62"
    gap_open = 10.0
    gap_extend = 0.5
    label1 = name1
    label2 = name2
    
    # Calculate alignment statistics (exactly like mutation.py)
    def calculate_stats(aln_seqA, aln_seqB):
        length = len(aln_seqA)
        identity = sum(a == b and a != '-' for a, b in zip(aln_seqA, aln_seqB))
        gaps = sum(a == '-' or b == '-' for a, b in zip(aln_seqA, aln_seqB))
        similarity = identity  # Simplified similarity
        return length, identity, similarity, gaps
    
    length, identity, similarity, gaps = calculate_stats(aln_seqA, aln_seqB)
    
    # Start building HTML output (exactly matching mutation.py emboss_style_html_doc)
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
    
    # Process alignment in blocks (exactly like mutation.py emboss_style_html_doc)
    line_width = 60
    ref_pos = user_pos = 1
    label_width = max(len(label1), len(label2)) + 2
    
    for i in range(0, length, line_width):
        ref_block = aln_seqA[i:i+line_width]
        user_block = aln_seqB[i:i+line_width]
        
        # Generate colored sequences (exactly matching mutation.py emboss_style_html_doc)
        ref_colored = ""
        user_colored = ""
        
        for a, b in zip(ref_block, user_block):
            if a == b and a != '-':
                # Perfect match - no coloring (exactly like mutation.py)
                ref_colored += f"<span>{a}</span>"
                user_colored += f"<span>{b}</span>"
            elif a != '-' and b != '-':
                # Substitution - exact colors from mutation.py
                ref_colored += f"<span style='background-color:#f88'>{a}</span>"
                user_colored += f"<span style='background-color:#8f8'>{b}</span>"
            elif a == '-' and b != '-':
                # Insertion (in user) - exact colors from mutation.py
                ref_colored += "<span style='background-color:#e0e0e0'>-</span>"
                user_colored += f"<span style='background-color:#8ff'>{b}</span>"
            elif a != '-' and b == '-':
                # Deletion (in user) - exact colors from mutation.py
                ref_colored += f"<span style='background-color:#fbb'>{a}</span>"
                user_colored += "<span style='background-color:#e0e0e0'>-</span>"
        
        # Count residues for position tracking (exactly like mutation.py)
        def count_residues(seq):
            return sum(1 for c in seq if c != '-')
        
        ref_end = ref_pos + count_residues(ref_block) - 1
        user_end = user_pos + count_residues(user_block) - 1
        
        # Format exactly like mutation.py emboss_style_html_doc with <br> tags and bold labels
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

def main():
    """Main application"""
    init_session_state()
    
    st.title("🧬 Genomics Research Assistant")
    st.markdown("*AI-powered bioinformatics analysis and research support*")
    
    # Sidebar
    with st.sidebar:
        st.header("Navigation")
        
        # Quick file upload in sidebar
        st.subheader("📁 Quick Upload")
        sidebar_upload = st.file_uploader(
            "Upload FASTA", 
            type=['fasta', 'fa', 'fas', 'txt'],
            key="sidebar_upload"
        )
        
        if sidebar_upload is not None:
            sequence_from_sidebar = parse_fasta_file(sidebar_upload)
            if sequence_from_sidebar:
                st.session_state.current_sequence = sequence_from_sidebar
                st.success(f"✅ Loaded {len(sequence_from_sidebar)} chars")
        
        # Clear session
        if st.button("🗑️ Clear Session", use_container_width=True):
            for key in st.session_state.keys():
                del st.session_state[key]
            st.rerun()
        
        # Help section
        with st.expander("ℹ️ How to Use"):
            st.markdown("""
            **Chat Interface:**
            - Ask questions about genes, mutations, diseases
            - Request analysis by typing "analyze this sequence: ATCG..."
            - Get explanations of results
            
            **File Upload:**
            - Supports .fasta, .fa, .fas, .txt formats
            - Multiple sequences: choose which to analyze
            - Large sequences: no copy-paste limits
            
            **Analysis Tools:**
            - Upload FASTA files or paste sequences
            - Analyze mutations against reference
            - Look up genetic variants
            - Search scientific literature
            
            **FASTA Format Example:**
            ```
            >sequence_name_or_description
            ATCGATCGATCGATCG...
            ```
            
            **Features:**
            - Multi-agent AI coordination
            - Context-aware responses
            - Interactive visualizations
            - Research literature integration
            """)
        
        # Sample files section
        with st.expander("📄 Sample Data"):
            st.markdown("""
            **Test Sequences:**
            
            **DNA Sample:**
            ```
            >sample_dna
            ATCGATCGATCGATCGATCG
            ```
            
            **Protein Sample:**
            ```
            >sample_protein
            MKWVTFISLLLLFSSAYSRGVFRRDAHKSEVAHRFKDLGEENFKALVLIAFAQYLQQCPFEDHVKLVNEVTEFAKTCVADESAENCDKSLHTLFGDKLCTVATLRETYGEMADCCAKQEPERNECFLQHKDDNPNLPRLVRPEVDVMCTAFHDNEETFLKKYLYEIARRHPYFYAPELLFFAKRYKAAFTECCQAADKAACLLPKLDELRDEGKASSAKQRLKCASLQKFGERAFKAWAVARLSQRFPKAEFAEVSKLVTDLTKVHTECCHGDLLECADDRADLAKYICENQDSISSKLKECCEKPLLEKSHCIAEVENDEMPADLPSLAADFVESKDVCKNYAEAKDVFLGMFLYEYARRHPDYSVVLLLRLAKTYETTLEKCCAAADPHECYAKVFDEFKPLVEEPQNLIKQNCELFEQLGEYKFQNALLVRYTKKVPQVSTPTLVEVSRNLGKVGSKCCKHPEAKRMPCAEDYLSVVLNQLCVLHEKTPVSDRVTKCCTESLVNRRPCFSALEVDETYVPKEFNAETFTFHADICTLSEKERQIKKQTALVELVKHKPKATKEQLKAVMDDFAAFVEKCCKADDKETCFAEEGKKLVAASQAALGL
            ```
            
            Copy and save as .fasta file to test!
            """)
    
    # Main content area
    col1, col2 = st.columns([1, 1])
    
    with col1:
        display_chat_interface()
    
    with col2:
        display_analysis_panel()
    
    # Results section (full width)
    display_results()

if __name__ == "__main__":
    main()