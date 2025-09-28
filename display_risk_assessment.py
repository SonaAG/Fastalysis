"""
Functions to display risk assessment results in the Streamlit frontend
"""

import streamlit as st
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
from typing import Dict, List
from risk_assessment import InfectionRiskAssessment

def display_risk_assessment(blast_results: Dict = None, mutation_results: Dict = None, 
                       literature_data: Dict = None, risk_results: Dict = None) -> Dict:
    """
    Display comprehensive risk assessment for a biological sequence without showing BLAST or mutation analysis
    
    Args:
        blast_results: Dictionary containing BLAST results (used for calculation only)
        mutation_results: Dictionary containing mutation analysis results (used for calculation only)
        literature_data: Optional dictionary containing literature search results
        risk_results: Optional pre-calculated risk results
    
    Returns:
        Dictionary containing risk assessment results
    """
    # Use provided risk results if available, otherwise calculate them
    if risk_results is None and blast_results is not None and mutation_results is not None:
        # Create risk assessment object
        risk_assessor = InfectionRiskAssessment()
        
        # Calculate risk
        risk_results = risk_assessor.calculate_infection_risk(
            blast_results, 
            mutation_results,
            literature_data
        )
    elif risk_results is None:
        st.error("Cannot display risk assessment - missing required data.")
        st.info("Please run BLAST and Mutation Analysis first to generate risk assessment data.")
        return None
    
    # Display risk dashboard
    st.header("🛡️ Biosafety Risk Assessment")
    
    # Add guidance note
    st.info("📌 This risk assessment is calculated using sequence similarity, mutation analysis, and taxonomic classification. Use the buttons below to view detailed analysis results.")
    
    # Top level risk metrics
    col1, col2, col3 = st.columns(3)
    
    # Risk Score with gauge
    with col1:
        risk_score = risk_results.get("risk_score", 0)
        risk_level = risk_results.get("risk_level", "Unknown")
        risk_color = risk_results.get("color", "gray")
        
        # Create gauge chart
        fig = go.Figure(go.Indicator(
            mode="gauge+number",
            value=risk_score,
            domain={'x': [0, 1], 'y': [0, 1]},
            title={'text': "Risk Score", 'font': {'size': 24}},
            gauge={
                'axis': {'range': [0, 100], 'tickwidth': 1, 'tickcolor': "darkblue"},
                'bar': {'color': risk_color},
                'bgcolor': "white",
                'borderwidth': 2,
                'bordercolor': "gray",
                'steps': [
                    {'range': [0, 40], 'color': 'green'},
                    {'range': [40, 70], 'color': 'yellow'},
                    {'range': [70, 100], 'color': 'red'}
                ],
            }
        ))
        st.plotly_chart(fig, use_container_width=True)
    
    # Risk level and confidence
    with col2:
        # Create colored box based on risk level
        if risk_level == "High":
            st.error(f"### 🚨 {risk_level} Risk")
        elif risk_level == "Moderate":
            st.warning(f"### ⚠️ {risk_level} Risk")
        elif risk_level == "Low":
            st.success(f"### ✅ {risk_level} Risk")
        else:
            st.info(f"### ℹ️ {risk_level} Risk")
        
        confidence = risk_results.get("confidence", 0)
        st.metric("Confidence", f"{confidence}%")
        
        # Show taxonomic classification if available
        taxonomy = risk_results.get("taxonomy_classification", {})
        if taxonomy:
            st.write("#### Taxonomic Classification:")
            taxa_md = ""
            
            # Display known taxonomy levels
            for level in ["superkingdom", "family", "genus", "species"]:
                if taxonomy.get(level):
                    taxa_md += f"- **{level.title()}**: {taxonomy[level]}\n"
            
            # Display risk group if available
            if taxonomy.get("risk_group"):
                if taxonomy["risk_group"] == 4:
                    taxa_md += f"- **Risk Group**: 🚨 RG4 (highest)\n"
                elif taxonomy["risk_group"] == 3:
                    taxa_md += f"- **Risk Group**: ⚠️ RG3 (high)\n"
                elif taxonomy["risk_group"] == 2:
                    taxa_md += f"- **Risk Group**: ⚠️ RG2 (moderate)\n"
                else:
                    taxa_md += f"- **Risk Group**: ✅ RG1 (low/none)\n"
            
            st.markdown(taxa_md)
    
    # Risk factors and recommendations
    with col3:
        # Risk factors
        risk_factors = risk_results.get("risk_factors", [])
        if risk_factors:
            st.write("#### Risk Factors:")
            for factor in risk_factors[:5]:  # Show top 5 factors
                st.markdown(f"- {factor}")
            
            if len(risk_factors) > 5:
                with st.expander(f"Show {len(risk_factors)-5} more factors"):
                    for factor in risk_factors[5:]:
                        st.markdown(f"- {factor}")
    
    # Recommendations section
    st.subheader("📋 Safety Recommendations")
    recommendations = risk_results.get("recommendations", [])
    if recommendations:
        for rec in recommendations:
            st.markdown(f"- {rec}")
    else:
        st.info("No specific recommendations available")
    
    # Risk score breakdown
    st.subheader("📊 Risk Score Breakdown")
    
    # Get risk contributions
    risk_details = risk_results.get("analysis_details", {}).get("risk_contributions", {})
    
    if risk_details:
        # Create a DataFrame for the risk components
        risk_df = pd.DataFrame({
            "Component": list(risk_details.keys()),
            "Score": list(risk_details.values())
        })
        
        # Sort by score
        risk_df = risk_df.sort_values("Score", ascending=False)
        
        # Create bar chart
        fig = px.bar(
            risk_df, 
            x="Score", 
            y="Component",
            orientation="h",
            title="Risk Score Components",
            color="Score",
            color_continuous_scale=["green", "yellow", "red"],
            labels={"Score": "Risk Points", "Component": "Risk Factor"}
        )
        st.plotly_chart(fig, use_container_width=True)
    
    # Scientific basis expander
    with st.expander("🧪 Scientific Basis for Risk Assessment"):
        st.write("""
        This risk assessment is based on established biosafety guidelines and pathogen classification systems:
        
        1. **Risk Group Classification**: Based on NIH Guidelines (Appendix B) and CDC/WHO biosafety levels
        2. **Taxonomic Risk**: Based on known pathogenicity of taxonomic groups
        3. **Sequence Similarity**: Higher similarity to known pathogens increases risk (BLAST analysis)
        4. **Functional Impact**: Mutations in functional sites may affect pathogenicity (Mutation analysis)
        5. **Literature Evidence**: Scientific literature about pathogenicity
        
        The system quantifies these factors to provide an overall risk score and actionable recommendations.
        """)
        
        # Add citation
        st.markdown("#### References:")
        st.markdown("""
        - NIH Guidelines for Research Involving Recombinant or Synthetic Nucleic Acid Molecules
        - WHO Laboratory Biosafety Manual
        - CDC Biosafety in Microbiological and Biomedical Laboratories (BMBL)
        - American Biological Safety Association (ABSA) Risk Group Database
        """)
        
        st.info("""
        **Note on Data Sources**: This risk assessment incorporates data from:
        
        - **Sequence similarity analysis**: BLAST results provide taxonomic identification and closest matches
        - **Mutation analysis**: Functional impact assessment of variations
        - **Taxonomic data**: Risk group classification from established databases
        
        For detailed sequence analysis visualizations, use the dedicated BLAST and Mutation Analysis tabs.
        """)
    
    return risk_results

def display_mini_risk_assessment(risk_results: Dict):
    """
    Display a compact version of risk assessment for embedding in other displays
    
    Args:
        risk_results: Dictionary containing risk assessment results
    """
    risk_score = risk_results.get("risk_score", 0)
    risk_level = risk_results.get("risk_level", "Unknown")
    risk_color = risk_results.get("color", "gray")
    
    # Risk level with appropriate styling
    if risk_level == "High":
        st.error(f"### 🚨 {risk_level} Risk ({risk_score:.1f}/100)")
    elif risk_level == "Moderate":
        st.warning(f"### ⚠️ {risk_level} Risk ({risk_score:.1f}/100)")
    else:
        st.success(f"### ✅ {risk_level} Risk ({risk_score:.1f}/100)")
    
    # Top 3 risk factors
    risk_factors = risk_results.get("risk_factors", [])
    if risk_factors:
        st.write("**Key Risk Factors:**")
        for factor in risk_factors[:3]:
            st.markdown(f"- {factor}")
        
        if len(risk_factors) > 3:
            st.markdown(f"*Plus {len(risk_factors)-3} more factors*")
    
    # Top recommendation
    recommendations = risk_results.get("recommendations", [])
    if recommendations:
        st.write("**Key Recommendation:**")
        st.markdown(f"- {recommendations[0]}")
        
        if len(recommendations) > 1:
            st.markdown(f"*Plus {len(recommendations)-1} more recommendations*")