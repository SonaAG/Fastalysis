import streamlit as st
import tempfile
import os
from controller import run_full_pipeline

st.set_page_config(page_title="Genomic Risk Assessment", page_icon="🧬")

st.title("🧬 Genomic Sequence Analysis & Risk Assessment")

# Sequence input
st.header("Sequence Input")
input_method = st.radio("Choose input method:", ["Enter Sequence", "Upload FASTA File"])

sequence = None
if input_method == "Enter Sequence":
    sequence = st.text_area("Enter DNA/RNA sequence:", height=150)
else:
    uploaded_file = st.file_uploader("Upload FASTA file", type=["fasta", "fa"])
    if uploaded_file:
        # Save uploaded file
        with tempfile.NamedTemporaryFile(delete=False, suffix=".fasta") as tmp:
            tmp.write(uploaded_file.getvalue())
            temp_fasta = tmp.name

if (sequence or 'temp_fasta' in locals()) and st.button("Run Analysis"):
    with st.spinner("Running genomic analysis..."):
        try:
            # Create temporary file if sequence was entered manually
            if sequence:
                temp_fasta = tempfile.mktemp(suffix=".fasta")
                with open(temp_fasta, "w") as f:
                    f.write(f">input_sequence\n{sequence}\n")
            
            # Run analysis
            results = run_full_pipeline(
                fasta_path=temp_fasta,
                num_blast_hits=5,
                max_variants=10,
                max_pubmed_results=5
            )
            
            # Display risk assessment
            if 'risk_assessment' in results:
                risk = results['risk_assessment']
                
                # Risk Level Header with colored background
                risk_level = risk['risk_level']
                if risk_level == "High":
                    st.error(f"🚨 Risk Level: {risk_level} ({risk['risk_score']}%)")
                elif risk_level == "Moderate":
                    st.warning(f"⚠️ Risk Level: {risk_level} ({risk['risk_score']}%)")
                else:
                    st.success(f"✅ Risk Level: {risk_level} ({risk['risk_score']}%)")
                
                # Display confidence score
                st.metric("Analysis Confidence", f"{risk['confidence']}%")
                
                # Risk factors
                st.subheader("Risk Factors Identified:")
                for factor in risk['risk_factors']:
                    st.warning(factor)
                
                # Recommendations
                st.subheader("Recommendations:")
                for rec in risk['recommendations']:
                    st.info(rec)
            
            # Display BLAST results
            if 'blast' in results and results['blast'].get('hits'):
                st.subheader("Top BLAST Matches")
                for hit in results['blast']['hits'][:3]:
                    st.write(f"🔹 **{hit.get('title', 'Unknown')}**")
                    st.write(f"   Identity: {hit.get('percent_identity', 0):.1f}%")
            
            # Display mutation summary
            if 'mutation_top1' in results:
                mutations = results['mutation_top1'].get('mutations', [])
                st.subheader(f"Mutations Detected: {len(mutations)}")
                
                if mutations:
                    st.write("Key Mutations:")
                    for mut in mutations[:5]:
                        st.code(f"Position {mut['position']}: {mut['ref']} → {mut['user']}")
            
            # Display full summary
            st.header("Analysis Summary")
            if 'summary' in results:
                st.markdown(results['summary'])
            
        except Exception as e:
            st.error(f"Error during analysis: {str(e)}")
        finally:
            # Clean up temporary file
            if 'temp_fasta' in locals():
                os.unlink(temp_fasta)

# Add helpful information
with st.expander("ℹ️ About this tool"):
    st.markdown("""
    This tool performs comprehensive genomic sequence analysis including:
    - BLAST sequence similarity search
    - Mutation detection and analysis
    - Risk assessment based on:
        - Sequence similarity to known pathogens
        - Presence of critical mutations
        - Literature-based evidence
    - Clinical recommendations
    
    The risk score is calculated using multiple factors and should be interpreted by qualified healthcare professionals.
    """)