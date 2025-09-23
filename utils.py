def format_analysis_result(result: dict, query_type: str) -> str:
    """Format analysis results into a readable message"""
    if not result or "status" not in result:
        return "Analysis failed: No results returned"
    
    if result["status"] == "error":
        return f"Analysis failed: {result.get('message', 'Unknown error')}"
    
    if "blast" in query_type.lower():
        if "data" in result and result["data"]:
            hits = result["data"].get("hits", [])
            response = "🧬 **BLAST Search Results:**\n\n"
            for i, hit in enumerate(hits[:5], 1):
                response += f"{i}. **{hit.get('title', 'Unknown')}**\n"
                response += f"   - Score: {hit.get('score', 'N/A')}\n"
                response += f"   - E-value: {hit.get('evalue', 'N/A')}\n"
                response += f"   - Identity: {hit.get('identity', 'N/A')}%\n\n"
            return response
    
    elif "literature" in query_type.lower():
        if "data" in result and result["data"]:
            papers = result["data"].get("papers", [])
            response = "📚 **Literature Search Results:**\n\n"
            for i, paper in enumerate(papers[:5], 1):
                response += f"{i}. **{paper.get('title', 'Unknown')}**\n"
                response += f"   Authors: {paper.get('authors', 'N/A')}\n"
                response += f"   Journal: {paper.get('journal', 'N/A')}\n"
                response += f"   PMID: [{paper.get('pmid', 'N/A')}](https://pubmed.ncbi.nlm.nih.gov/{paper.get('pmid', '')})\n\n"
            return response
    
    else:  # General analysis
        if "data" in result and result["data"]:
            data = result["data"]
            response = "🔬 **Sequence Analysis Results:**\n\n"
            
            if "type" in data:
                response += f"**Sequence Type:** {data['type']}\n"
            if "length" in data:
                response += f"**Length:** {data['length']} bp/aa\n"
            if "features" in data:
                response += "\n**Detected Features:**\n"
                for feature in data["features"]:
                    response += f"- {feature}\n"
            if "blast_hits" in data:
                response += "\n**Top BLAST Hit:**\n"
                hit = data["blast_hits"][0]
                response += f"- **{hit.get('title', 'Unknown')}**\n"
                response += f"- Score: {hit.get('score', 'N/A')}\n"
                response += f"- Identity: {hit.get('identity', 'N/A')}%\n"
            
            return response
    
    return "No relevant results found in the analysis"