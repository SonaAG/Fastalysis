from typing import Dict, List
import numpy as np

def calculate_infection_risk(blast_results: Dict, mutation_results: Dict, literature_data: Dict = None) -> Dict:
    """
    Calculate infection risk and severity based on multiple factors:
    1. Sequence similarity to known pathogens
    2. Presence of critical mutations
    3. Associated clinical significance from literature
    """
    risk_score = 0
    confidence = 0
    risk_factors = []
    recommendations = []
    
    # Analyze BLAST hits
    pathogen_match = False
    highest_identity = 0
    matching_pathogen = ""
    
    if "hits" in blast_results:
        for hit in blast_results["hits"]:
            title = hit.get("title", "").lower()
            identity = float(hit.get("percent_identity", 0))
            
            # Check for pathogen keywords
            pathogen_keywords = ["virus", "pathogen", "infection", "disease", "bacterial", "fungal"]
            if any(keyword in title for keyword in pathogen_keywords):
                pathogen_match = True
                if identity > highest_identity:
                    highest_identity = identity
                    matching_pathogen = hit.get("title", "")
    
    # Score based on sequence similarity
    if pathogen_match:
        if highest_identity >= 95:
            risk_score += 40
            risk_factors.append(f"High sequence similarity ({highest_identity:.1f}%) to known pathogen: {matching_pathogen}")
            recommendations.append("Immediate clinical assessment recommended")
        elif highest_identity >= 85:
            risk_score += 30
            risk_factors.append(f"Moderate sequence similarity ({highest_identity:.1f}%) to known pathogen: {matching_pathogen}")
            recommendations.append("Further diagnostic testing recommended")
        elif highest_identity >= 70:
            risk_score += 20
            risk_factors.append(f"Low sequence similarity ({highest_identity:.1f}%) to known pathogen: {matching_pathogen}")
            recommendations.append("Monitor for symptoms and consider follow-up testing")
    
    # Analyze mutations
    if "mutations" in mutation_results:
        mutation_count = len(mutation_results["mutations"])
        critical_mutations = []
        
        # Define critical regions/positions (example)
        critical_positions = [10, 20, 30, 40, 50]  # Example positions
        
        for mutation in mutation_results["mutations"]:
            pos = mutation.get("position", 0)
            if pos in critical_positions:
                critical_mutations.append(mutation)
        
        if mutation_count > 10:
            risk_score += 20
            risk_factors.append(f"High number of mutations detected ({mutation_count})")
        elif mutation_count > 5:
            risk_score += 10
            risk_factors.append(f"Moderate number of mutations detected ({mutation_count})")
        
        if critical_mutations:
            risk_score += 20
            risk_factors.append(f"Mutations found in {len(critical_mutations)} critical positions")
            recommendations.append("Genetic counseling recommended due to mutations in critical regions")
    
    # Calculate confidence based on data quality
    confidence = min(100, max(0, (
        blast_results.get("hits", []).__len__() * 10 +  # More hits = more confidence
        (mutation_results.get("alignment", {}).get("percent_identity", 0)) +  # Higher identity = more confidence
        (50 if literature_data else 0)  # Literature data adds confidence
    )))
    
    # Determine risk level and color
    if risk_score >= 70:
        risk_level = "High"
        color = "red"
    elif risk_score >= 40:
        risk_level = "Moderate"
        color = "yellow"
    else:
        risk_level = "Low"
        color = "green"
    
    return {
        "risk_score": risk_score,
        "risk_level": risk_level,
        "color": color,
        "confidence": confidence,
        "risk_factors": risk_factors,
        "recommendations": recommendations,
        "analysis_details": {
            "sequence_match": {
                "pathogen_match": pathogen_match,
                "highest_identity": highest_identity,
                "matching_pathogen": matching_pathogen if pathogen_match else None
            },
            "mutations": {
                "total_count": len(mutation_results.get("mutations", [])),
                "critical_regions_affected": bool(critical_mutations) if "mutations" in mutation_results else False
            }
        }
    }