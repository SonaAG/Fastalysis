# 🧬 Genomics Research Assistant - Setup Guide

## System Architecture

```
Frontend (Streamlit)  →  FastAPI Gateway  →  AutoGen Agents  →  RAG System
                                         ↓
                                  Bioinformatics Tools
                                  (BLAST, Mutations, PubMed)
```

## Quick Start

### 1. Install Dependencies
```bash
cd d:\Fastalysis\app
pip install -r requirements.txt
```

### 2. Start FastAPI Backend
```bash
python main.py
# API will be available at http://localhost:8000
# Interactive docs at http://localhost:8000/docs
```

### 3. Start Enhanced Streamlit Frontend
```bash
streamlit run enhanced_streamlit_app.py
# App will be available at http://localhost:8501
```

## Configuration

### Chat Models Setup
Edit `genomics_agents.py` to configure your preferred models:

```python
llm_config = {
    "config_list": [
        {
            "model": "kimi-k2",  # Your model choice
            "api_key": "your-api-key-if-needed",
            "api_base": "your-model-endpoint",  # e.g., local or cloud
        }
    ],
    "temperature": 0.7,
}
```

### Environment Variables
Create `.env` file:
```
NCBI_EMAIL=your_email@domain.com
OPENAI_API_KEY=your_openai_key  # if using OpenAI models
HUGGINGFACE_TOKEN=your_hf_token  # if using HuggingFace models
```

## Features

### 🤖 AI Agent System
- **Coordinator Agent**: Understands user queries, orchestrates analysis
- **Analysis Agent**: Executes bioinformatics tools (BLAST, mutations)
- **Literature Agent**: Searches PubMed, summarizes research
- **Knowledge Agent**: Provides context using RAG system

### 🔬 Analysis Capabilities
- BLAST sequence similarity search
- Mutation analysis with visualization
- Genetic variant lookup
- Literature search and summarization
- Interactive result visualization

### 💬 Natural Language Interface
- "Analyze mutations in BRCA1 gene"
- "What do you know about this sequence: ATCGATCG..."
- "Find recent papers about TP53 mutations"
- "Explain this BLAST result"

## Sample Prompts & Questions

### 🧬 Sequence Analysis Prompts
**Basic Analysis:**
- "Analyze this sequence" 
- "What type of sequence is this?"
- "Identify this DNA/protein sequence"
- "Run BLAST search on my sequence"
- "Find similar sequences in the database"

**Detailed Analysis:**
- "Analyze mutations in this sequence compared to reference"
- "What mutations are present in this sequence?"
- "Compare this sequence to the reference genome"
- "Show me the alignment and mutations"

### 🔬 Gene-Specific Questions
**Gene Information:**
- "Tell me about the BRCA1 gene"
- "What is the function of TP53?"
- "What diseases are associated with CFTR mutations?"
- "Explain the role of APOE in Alzheimer's disease"

**Mutation Analysis:**
- "What are common BRCA1 mutations?"
- "Analyze mutations in the TP53 gene"
- "What are the clinical implications of this mutation?"
- "Is this mutation pathogenic or benign?"

### 📚 Literature & Research Questions
**Research Queries:**
- "Find recent papers about CRISPR gene editing"
- "What's new in cancer genomics research?"
- "Show me literature on BRCA1 mutations and breast cancer"
- "Find reviews about precision medicine in oncology"

**Clinical Context:**
- "What are the clinical guidelines for BRCA testing?"
- "How is this gene used in personalized medicine?"
- "What are the therapeutic targets for this pathway?"

### 🧪 Variant & Disease Questions
**Variant Analysis:**
- "Look up variants for the BRCA1 gene"
- "What are known pathogenic variants in TP53?"
- "Find ClinVar entries for this gene"
- "Show me population frequencies for these variants"

**Disease Associations:**
- "What diseases are linked to this gene?"
- "How does this mutation cause disease?"
- "What is the inheritance pattern?"
- "What are the risk factors?"

### 🔍 Technical Analysis Prompts
**BLAST & Similarity:**
- "Find homologous sequences"
- "What organism is this sequence from?"
- "Search for similar proteins/genes"
- "Compare this to known sequences"

**Functional Analysis:**
- "Predict the protein structure"
- "What domains are present in this protein?"
- "Analyze the evolutionary conservation"
- "What are the functional consequences?"

### 🤖 AI Assistant Capabilities
**Explanation & Education:**
- "Explain this result in simple terms"
- "What does this E-value mean?"
- "How do I interpret this alignment?"
- "Teach me about molecular genetics"

**Research Assistance:**
- "Help me design primers for PCR"
- "Suggest experiments to validate this finding"
- "What controls should I include?"
- "How do I confirm this mutation?"

### 📊 Data Interpretation
**Results Analysis:**
- "Interpret these BLAST results"
- "What do these statistics mean?"
- "Is this alignment significant?"
- "How reliable are these predictions?"

**Quality Control:**
- "Is this sequence quality good?"
- "Are there any problematic regions?"
- "What could cause this result?"
- "How can I improve my analysis?"

## Advanced Query Examples

### 🧬 Multi-Step Analysis
```
"I have a patient with breast cancer family history. 
Analyze this BRCA1 sequence for mutations, 
find literature on clinical significance, 
and suggest follow-up testing."
```

### 🔬 Research Workflow
```
"Compare this novel sequence to known oncogenes,
identify potential functional domains,
and find recent papers on similar discoveries."
```

### 📈 Clinical Decision Support
```
"This variant was found in genetic testing.
Explain its pathogenicity, clinical implications,
and current treatment guidelines."
```

## API Endpoints

### Core Analysis
- `POST /blast` - BLAST sequence search
- `POST /mutation` - Mutation analysis
- `POST /variants` - Variant lookup
- `POST /pubmed` - Literature search
- `POST /analyze` - Full pipeline analysis

### AI Interface
- `POST /chat` - Intelligent chat with context
- `GET /models` - Available AI models

## Advanced Usage

### N8N Integration
For workflow automation, create N8N workflows that call your FastAPI endpoints:

1. **Sequence Processing Workflow**:
   ```
   Trigger → HTTP Request (BLAST) → HTTP Request (Mutation) → Format Results
   ```

2. **Literature Review Workflow**:
   ```
   Gene List → For Each → PubMed Search → Summarize → Email Report
   ```

### RAG System Enhancement
Add custom knowledge:

```python
from genomics_rag import GenomicsRAG

rag = GenomicsRAG()
rag.add_pubmed_abstracts(your_abstracts)  # Add literature
rag.populate_custom_genes(your_gene_data)  # Add gene info
```

### Extending Agents
Add custom agents in `genomics_agents.py`:

```python
custom_agent = AssistantAgent(
    name="ProteinStructureAgent",
    system_message="Analyze protein structures and domains...",
    llm_config=llm_config
)
```

## Performance Optimization

### Caching
- Enable FastAPI caching for repeated BLAST queries
- Cache agent system initialization
- Use Redis for session management

### Scaling
- Deploy FastAPI with gunicorn/uvicorn workers
- Use Docker containers for deployment
- Consider cloud deployment (AWS/GCP/Azure)

## Troubleshooting

### Common Issues

1. **BLAST Timeout**: Increase timeout in BLAST functions
2. **Model API Errors**: Check your model endpoint configuration
3. **Memory Issues**: Reduce agent conversation history
4. **Database Errors**: Ensure ChromaDB permissions

### Logs
Check logs in:
- FastAPI: Console output
- Streamlit: Browser console and terminal
- ChromaDB: `./genomics_db/chroma.log`

## Next Steps

### Recommended Enhancements
1. **Add more specialized agents** (Protein structure, Phylogenetics)
2. **Implement user authentication** for session management
3. **Add result export** (PDF reports, CSV downloads)
4. **Create mobile-responsive UI** with React/Next.js
5. **Add real-time collaboration** features
6. **Integrate with lab instruments** and databases

### Production Deployment
1. **Containerization** with Docker
2. **Load balancing** with nginx
3. **Database backup** strategies
4. **Monitoring** with Prometheus/Grafana
5. **Security** (HTTPS, API keys, rate limiting)

## Support

For issues or enhancements:
1. Check existing analysis logs
2. Review API documentation at `/docs`
3. Test individual components separately
4. Consider fallback modes for critical functions

## 💡 Tips for Best Results

### 🎯 Getting Better Chat Responses
- **Be specific**: "Analyze BRCA1 mutations" vs "Tell me about genes"
- **Use context**: Upload sequence first, then ask questions
- **Ask follow-up questions**: "Explain that result in more detail"
- **Combine requests**: "Analyze this sequence and find related literature"

### 🧬 Sequence Analysis Tips
- **Valid formats**: FASTA, raw sequence (DNA/RNA: ATCG, Protein: standard amino acids)
- **Sequence length**: Longer sequences (>50bp) give better BLAST results
- **Quality check**: Remove ambiguous characters (N, X) if possible
- **Context matters**: Mention species, gene name, or source if known

### 📚 Literature Search Tips
- **Use gene symbols**: "BRCA1" instead of "breast cancer gene 1"
- **Add keywords**: "BRCA1 mutations breast cancer treatment"
- **Ask for summaries**: "Summarize recent CRISPR research"
- **Request specific info**: "What are clinical guidelines for BRCA testing?"

### 🔍 Interpretation Help
- **Ask for explanations**: "What does this E-value mean?"
- **Request context**: "Is this mutation clinically significant?"
- **Get recommendations**: "What experiments should I do next?"
- **Compare results**: "How does this compare to normal variants?"

### 🚀 Advanced Usage Examples
```
"I uploaded a BRCA1 sequence. Please analyze it for mutations, 
find recent literature on clinical significance, and explain 
what this means for patient care."
```

```
"This BLAST result shows 98% identity. What does this mean? 
Is this sequence from human or another organism?"
```

```
"Compare this protein sequence to known oncogenes and 
suggest functional studies to characterize it."
```