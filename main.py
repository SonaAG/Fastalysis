from fastapi import FastAPI, HTTPException, UploadFile, File
from fastapi.middleware.cors import CORSMiddleware
from pydantic import BaseModel
from typing import Optional, List
import tempfile
import os
import json
import traceback
import re
from dotenv import load_dotenv
load_dotenv('.env.local')

# Import controller functions directly - these should work
try:
    from controller import (
        run_blast_controller,
        run_mutation_controller,
        run_variant_table_controller,
        run_pubmed_controller,
        run_full_pipeline
    )
    print("✅ Successfully imported controller functions")
except Exception as e:
    print(f"❌ Error importing controller: {e}")

# Data models
class BlastRequest(BaseModel):
    # Create dummy functions as fallback
    def run_blast_controller(*args, **kwargs):
        return {"error": "Controller import failed"}
    def run_mutation_controller(*args, **kwargs):
        return {"error": "Controller import failed"}
    def run_variant_table_controller(*args, **kwargs):
        return {"error": "Controller import failed"}
    def run_pubmed_controller(*args, **kwargs):
        return {"error": "Controller import failed"}
    def run_full_pipeline(*args, **kwargs):
        return {"error": "Controller import failed"}

app = FastAPI(
    title="🧬 Fastalysis Genomics API",
    description="Comprehensive bioinformatics analysis API with BLAST, mutation analysis, variant lookup, and literature search",
    version="1.0.0"
)

# Enable CORS for frontend integration
app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

# Pydantic models for request/response
class BlastRequest(BaseModel):
    sequence: str
    num_hits: int = 5

class MutationRequest(BaseModel):
    sequence: str
    hit_index: int = 0

class VariantRequest(BaseModel):
    gene_or_accession: str
    by_accession: bool = True
    max_variants: int = 20

class PubMedRequest(BaseModel):
    accession: str
    additional_terms: Optional[str] = None
    max_results: int = 10

class FullPipelineRequest(BaseModel):
    sequence: str
    variant_gene_or_accession: Optional[str] = None
    pubmed_accession: Optional[str] = None
    additional_pubmed_terms: Optional[str] = None
    num_blast_hits: int = 5
    max_variants: int = 20
    max_pubmed_results: int = 10
    include_risk_assessment: bool = True  # Whether to include infection risk assessment

class ChatRequest(BaseModel):
    message: str
    model: str = "kimi-k2"
    context: Optional[dict] = None

# Utility function to create temporary FASTA file
def create_temp_fasta(sequence: str) -> str:
    """Create a temporary FASTA file from sequence string"""
    temp_file = tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False)
    temp_file.write(f">user_sequence\n{sequence}")
    temp_file.close()
    return temp_file.name

def analyze_user_intent(user_message: str, has_sequence: bool = False) -> dict:
    """
    Intelligently analyze user intent using context analysis instead of simple keywords.
    Returns a dictionary with intent classification and extracted parameters.
    """
    
    # Convert to lowercase for analysis
    message_lower = user_message.lower()
    
    # Initialize intent analysis result
    intent_result = {
        'primary_intent': 'general_question',
        'wants_literature': False,
        'wants_sequence_analysis': False,
        'wants_variants': False,
        'wants_mutations': False,
        'search_terms': [],
        'confidence': 0.5
    }
    
    # Context-aware analysis patterns
    literature_patterns = [
        r'find.*(?:literature|papers|research|articles|publications|studies)',
        r'search.*(?:literature|papers|research|articles|publications|studies)',
        r'(?:literature|papers|research|articles|publications|studies).*(?:about|on|for)',
        r'what.*(?:research|studies|papers).*(?:available|published|exist)',
        r'recent.*(?:publications|studies|research|papers)',
        r'literature review',
        r'bibliography',
        r'pubmed.*search'
    ]
    
    sequence_analysis_patterns = [
        r'analyze.*(?:sequence|dna|rna|protein|fasta)',
        r'(?:blast|alignment|homology|similarity).*(?:search|analysis)',
        r'identify.*(?:sequence|organism|species)',
        r'what.*(?:is this|organism|species|protein|gene)',
        r'compare.*(?:sequence|dna|rna|protein)',
        r'sequence.*(?:analysis|identification|comparison)',
        r'blast.*(?:search|analysis)',
        r'find.*(?:similar|homologous|matching).*sequences'
    ]
    
    variant_patterns = [
        r'(?:variants|mutations|snp|polymorphism).*(?:analysis|search|find)',
        r'find.*(?:variants|mutations|snp|polymorphisms)',
        r'genetic.*(?:variants|variations|mutations)',
        r'(?:clinical|pathogenic).*(?:variants|mutations)',
        r'variant.*(?:database|search|analysis)',
        r'mutation.*(?:analysis|screening|detection)',
        r'snp.*(?:analysis|detection|search)'
    ]
    
    # Check for literature intent
    for pattern in literature_patterns:
        if re.search(pattern, message_lower):
            intent_result['wants_literature'] = True
            intent_result['primary_intent'] = 'literature_search'
            intent_result['confidence'] = 0.8
            break
    
    # Check for sequence analysis intent
    for pattern in sequence_analysis_patterns:
        if re.search(pattern, message_lower):
            intent_result['wants_sequence_analysis'] = True
            if intent_result['primary_intent'] == 'general_question':
                intent_result['primary_intent'] = 'sequence_analysis'
                intent_result['confidence'] = 0.8
            break
    
    # Check for variant analysis intent
    for pattern in variant_patterns:
        if re.search(pattern, message_lower):
            intent_result['wants_variants'] = True
            if intent_result['primary_intent'] == 'general_question':
                intent_result['primary_intent'] = 'variant_analysis'
                intent_result['confidence'] = 0.8
            break
    
    # Extract potential search terms (gene names, proteins, diseases)
    # Look for capitalized words that might be gene names
    words = user_message.split()
    for i, word in enumerate(words):
        cleaned_word = re.sub(r'[^\w]', '', word)
        
        # Gene names are often capitalized (BRCA1, TP53, etc.)
        if cleaned_word.isupper() and len(cleaned_word) > 2:
            intent_result['search_terms'].append(cleaned_word)
        
        # Look for words after "for", "about", "on"
        if i > 0 and words[i-1].lower() in ['for', 'about', 'on', 'regarding', 'concerning']:
            if len(cleaned_word) > 2:
                intent_result['search_terms'].append(cleaned_word)
        
        # Look for quoted terms
        if word.startswith('"') and word.endswith('"'):
            clean_quoted = word.strip('"')
            if len(clean_quoted) > 2:
                intent_result['search_terms'].append(clean_quoted)
    
    # Context-based adjustments
    if has_sequence:
        if intent_result['wants_literature'] and not intent_result['wants_sequence_analysis']:
            # If user has sequence and wants literature, they probably want analysis first
            intent_result['wants_sequence_analysis'] = True
        
        if not any([intent_result['wants_literature'], intent_result['wants_sequence_analysis'], intent_result['wants_variants']]):
            # If user uploaded sequence but intent is unclear, assume they want analysis
            intent_result['wants_sequence_analysis'] = True
            intent_result['primary_intent'] = 'sequence_analysis'
            intent_result['confidence'] = 0.7
    
    # Detect mutation-specific requests
    if re.search(r'mutation.*(?:analysis|detection|comparison)', message_lower):
        intent_result['wants_mutations'] = True
    
    return intent_result

@app.get("/")
async def root():
    return {"message": "🧬 Fastalysis Genomics API is running!"}

@app.post("/blast")
async def blast_analysis(request: BlastRequest):
    """Run BLAST analysis on user sequence"""
    try:
        temp_fasta = create_temp_fasta(request.sequence)
        result = run_blast_controller(temp_fasta, num_hits=request.num_hits)
        os.unlink(temp_fasta)  # Clean up temp file
        return {"status": "success", "data": result}
    except Exception as e:
        # Return detailed error info for debugging
        import traceback
        error_details = traceback.format_exc()
        print(f"BLAST Error: {error_details}")  # Log to console
        return {
            "status": "error", 
            "message": str(e),
            "details": error_details
        }

@app.post("/mutation")
async def mutation_analysis(request: MutationRequest):
    """Run mutation analysis comparing user sequence with BLAST hits"""
    try:
        temp_fasta = create_temp_fasta(request.sequence)
        result = run_mutation_controller(temp_fasta, hit_index=request.hit_index)
        os.unlink(temp_fasta)
        return {"status": "success", "data": result}
    except Exception as e:
        import traceback
        error_details = traceback.format_exc()
        print(f"Mutation Error: {error_details}")
        return {
            "status": "error", 
            "message": str(e),
            "details": error_details
        }

@app.post("/variants")
async def variant_lookup(request: VariantRequest):
    """Look up genetic variants for a gene or accession"""
    try:
        result = run_variant_table_controller(
            request.gene_or_accession,
            by_accession=request.by_accession,
            max_variants=request.max_variants
        )
        return {"status": "success", "data": result}
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))

@app.post("/pubmed")
async def pubmed_search(request: PubMedRequest):
    """Search PubMed for relevant literature"""
    try:
        result = run_pubmed_controller(
            request.accession,
            additional_terms=request.additional_terms,
            max_results=request.max_results
        )
        return {"status": "success", "data": result}
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))

@app.post("/analyze")
async def full_analysis(request: FullPipelineRequest):
    """Run complete genomics analysis pipeline"""
    try:
        temp_fasta = create_temp_fasta(request.sequence)
        result = run_full_pipeline(
            temp_fasta,
            variant_accession_or_gene=request.variant_gene_or_accession,
            pubmed_accession=request.pubmed_accession,
            additional_pubmed_terms=request.additional_pubmed_terms,
            num_blast_hits=request.num_blast_hits,
            max_variants=request.max_variants,
            max_pubmed_results=request.max_pubmed_results
        )
        os.unlink(temp_fasta)
        return {"status": "success", "data": result}
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))

@app.post("/upload-fasta")
async def upload_fasta_file(file: UploadFile = File(...)):
    """Upload and parse FASTA file"""
    try:
        # Read file content
        content = await file.read()
        file_content = content.decode('utf-8')
        
        # Simple FASTA parsing
        lines = file_content.strip().split('\n')
        sequences = []
        current_header = ""
        current_sequence = ""
        
        for line in lines:
            line = line.strip()
            if line.startswith('>'):
                # Save previous sequence if exists
                if current_sequence:
                    sequences.append({
                        "header": current_header,
                        "sequence": current_sequence
                    })
                current_header = line[1:]  # Remove '>'
                current_sequence = ""
            elif line:
                current_sequence += line
        
        # Add last sequence
        if current_sequence:
            sequences.append({
                "header": current_header,
                "sequence": current_sequence
            })
        
        return {
            "status": "success",
            "filename": file.filename,
            "sequences_found": len(sequences),
            "sequences": sequences
        }
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))

@app.post("/blast-file")
async def blast_file_analysis(file: UploadFile = File(...), num_hits: int = 5):
    """Run BLAST analysis on uploaded FASTA file"""
    try:
        # Parse uploaded file
        content = await file.read()
        file_content = content.decode('utf-8')
        
        # Extract first sequence
        lines = file_content.strip().split('\n')
        sequence = ""
        for line in lines:
            if not line.startswith('>') and line.strip():
                sequence += line.strip()
        
        if not sequence:
            raise HTTPException(status_code=400, detail="No valid sequence found in file")
        
        # Create temporary file and run analysis
        temp_fasta = create_temp_fasta(sequence)
        result = run_blast_controller(temp_fasta, num_hits=num_hits)
        os.unlink(temp_fasta)
        
        return {"status": "success", "filename": file.filename, "data": result}
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))

@app.post("/chat")
async def chat_endpoint(request: ChatRequest):
    """Intelligent chat interface with sequence analysis capabilities"""
    try:
        user_message = request.message
        sequence_data = request.context.get('sequence') if request.context else None
        
        # Use intelligent intent analysis instead of simple keywords
        intent_analysis = analyze_user_intent(user_message, has_sequence=bool(sequence_data))
        
        # Debug logging
        print(f"DEBUG: User message: {user_message}")
        print(f"DEBUG: Has sequence: {bool(sequence_data)}")
        print(f"DEBUG: Intent analysis: {intent_analysis}")
        
        # Extract intents from analysis
        wants_analysis = intent_analysis['wants_sequence_analysis']
        wants_literature = intent_analysis['wants_literature'] 
        wants_variants = intent_analysis['wants_variants']
        wants_mutations = intent_analysis['wants_mutations']
        search_terms = intent_analysis['search_terms']
        primary_intent = intent_analysis['primary_intent']
        
        response = {
            "message": f"Processing your genomics query: {request.message}",
            "model": request.model,
            "status": "intelligent_mode",
            "intent_analysis": {
                "detected_intent": primary_intent,
                "confidence": intent_analysis['confidence'],
                "extracted_terms": search_terms
            }
        }
        
        # If user wants analysis and we have sequence data, perform it
        if (wants_analysis or wants_literature or wants_variants) and sequence_data:
            print(f"DEBUG: Starting analysis - wants_analysis: {wants_analysis}, wants_literature: {wants_literature}, wants_variants: {wants_variants}")
            try:
                # Create temporary file for analysis
                temp_fasta = create_temp_fasta(sequence_data)
                
                analysis_results = {}
                
                # Run BLAST analysis first to get gene/accession info
                blast_result = run_blast_controller(temp_fasta, num_hits=5)
                analysis_results['blast'] = blast_result
                
                # If user wants literature, search using BLAST results
                if wants_literature and 'hits' in blast_result and blast_result['hits']:
                    literature_results = []
                    
                    # Try to search literature using top BLAST hits
                    for i, hit in enumerate(blast_result['hits'][:3]):  # Use top 3 hits
                        accession = hit.get('accession', hit.get('Accession', ''))
                        title = hit.get('title', hit.get('Title', ''))
                        
                        if accession:
                            try:
                                # Search PubMed using accession
                                lit_result = run_pubmed_controller(accession, max_results=5)
                                if lit_result and 'papers' in lit_result and lit_result['papers']:
                                    literature_results.extend(lit_result['papers'][:3])  # Take top 3 papers per hit
                            except Exception as lit_error:
                                print(f"Literature search failed for {accession}: {lit_error}")
                        
                        # Also try searching with gene name if available in title
                        if 'gene' in title.lower() or 'protein' in title.lower():
                            # Extract potential gene names from title
                            words = title.split()
                            for word in words:
                                if len(word) > 2 and word.isupper():  # Likely gene symbol
                                    try:
                                        gene_lit_result = run_pubmed_controller(word, max_results=3)
                                        if gene_lit_result and 'papers' in gene_lit_result and gene_lit_result['papers']:
                                            literature_results.extend(gene_lit_result['papers'][:2])
                                            break
                                    except:
                                        continue
                    
                    # Remove duplicates based on PMID
                    seen_pmids = set()
                    unique_papers = []
                    for paper in literature_results:
                        pmid = paper.get('pmid', '')
                        if pmid and pmid not in seen_pmids:
                            seen_pmids.add(pmid)
                            unique_papers.append(paper)
                    
                    analysis_results['literature'] = {'papers': unique_papers[:10]}  # Limit to top 10
                
                # Run mutation analysis if requested or if doing sequence analysis
                if wants_mutations or wants_analysis or 'mutation' in user_message.lower():
                    try:
                        print(f"DEBUG: Running mutation analysis")
                        mutation_result = run_mutation_controller(temp_fasta, hit_index=0)
                        analysis_results['mutation'] = mutation_result
                        print(f"DEBUG: Mutation analysis completed: {mutation_result}")
                    except Exception as mut_error:
                        print(f"DEBUG: Mutation analysis failed: {mut_error}")
                        analysis_results['mutation'] = {'error': str(mut_error)}
                
                # Run variant analysis if requested
                if wants_variants and 'hits' in blast_result and blast_result['hits']:
                    # Try to get variants using first BLAST hit
                    first_hit = blast_result['hits'][0]
                    accession = first_hit.get('accession', first_hit.get('Accession', ''))
                    if accession:
                        try:
                            variant_result = run_variant_table_controller(accession, by_accession=True, max_variants=10)
                            analysis_results['variants'] = variant_result
                        except:
                            pass
                
                # Clean up temp file
                os.unlink(temp_fasta)
                
                # Format response with analysis results
                formatted_response = format_analysis_response(analysis_results, sequence_data, user_message)
                response = {
                    "message": formatted_response,
                    "model": request.model,
                    "status": "analysis_complete",
                    "analysis_data": analysis_results
                }
                
            except Exception as analysis_error:
                response = {
                    "message": f"I tried to analyze your sequence but encountered an error: {str(analysis_error)}. Please try using the analysis tools directly.",
                    "model": request.model,
                    "status": "analysis_failed",
                    "error": str(analysis_error)
                }
        
        # If user wants literature or variants without sequence, use extracted search terms
        elif wants_literature and not sequence_data:
            # Use intelligently extracted search terms instead of manual parsing
            if search_terms:
                try:
                    literature_results = []
                    for term in search_terms[:3]:  # Search top 3 terms
                        lit_result = run_pubmed_controller(term, max_results=5)
                        if lit_result and 'papers' in lit_result and lit_result['papers']:
                            literature_results.extend(lit_result['papers'])
                    
                    # Remove duplicates
                    seen_pmids = set()
                    unique_papers = []
                    for paper in literature_results:
                        pmid = paper.get('pmid', '')
                        if pmid and pmid not in seen_pmids:
                            seen_pmids.add(pmid)
                            unique_papers.append(paper)
                    
                    if unique_papers:
                        response_text = f"📚 **Literature Search Results for: {', '.join(search_terms)}**\n\n"
                        response_text += f"Found {len(unique_papers)} relevant papers:\n\n"
                        
                        for i, paper in enumerate(unique_papers[:10], 1):
                            title = paper.get('title', 'Untitled')
                            authors = paper.get('authors', 'Unknown authors')
                            year = paper.get('date', '')[:4] if paper.get('date') else 'Unknown'
                            pmid = paper.get('pmid', '')
                            
                            if len(title) > 100:
                                title = title[:97] + "..."
                            
                            response_text += f"**{i}. {title}** ({year})\n"
                            response_text += f"   Authors: {authors}\n"
                            if pmid:
                                response_text += f"   PMID: {pmid}\n"
                            response_text += "\n"
                        
                        response = {
                            "message": response_text,
                            "model": request.model,
                            "status": "literature_found",
                            "analysis_data": {"literature": {"papers": unique_papers[:10]}}
                        }
                    else:
                        response = {
                            "message": f"📚 I searched for literature about '{', '.join(search_terms)}' but couldn't find relevant papers. Try being more specific with gene names or topics.",
                            "model": request.model,
                            "status": "no_literature"
                        }
                
                except Exception as lit_error:
                    response = {
                        "message": f"❌ Literature search encountered an error: {str(lit_error)}. Please try again or contact support.",
                        "model": request.model,
                        "status": "literature_error",
                        "error": str(lit_error)
                    }
            else:
                response = {
                    "message": "📚 I'd be happy to search literature for you! Please specify what you're looking for (e.g., 'Find literature for BRCA1' or 'Search papers about COVID-19 mutations').",
                    "model": request.model,
                    "status": "needs_search_terms"
                }
        
        # Special case: If we have sequence data but no clear intent, assume analysis is wanted
        elif sequence_data and not (wants_analysis or wants_literature or wants_variants):
            print(f"DEBUG: Have sequence but no clear intent, assuming analysis is wanted")
            try:
                # Create temporary file for analysis
                temp_fasta = create_temp_fasta(sequence_data)
                
                analysis_results = {}
                
                # Run BLAST analysis by default
                blast_result = run_blast_controller(temp_fasta, num_hits=5)
                analysis_results['blast'] = blast_result
                
                # Clean up temp file
                os.unlink(temp_fasta)
                
                # Format response with analysis results
                formatted_response = format_analysis_response(analysis_results, sequence_data, user_message)
                response = {
                    "message": formatted_response,
                    "model": request.model,
                    "status": "auto_analysis_complete",
                    "analysis_data": analysis_results
                }
                
            except Exception as analysis_error:
                response = {
                    "message": f"I detected you have a sequence, but encountered an error analyzing it: {str(analysis_error)}. Please try using the analysis tools directly.",
                    "model": request.model,
                    "status": "auto_analysis_failed",
                    "error": str(analysis_error)
                }
        
        # Provide intelligent fallback based on detected intent
        else:
            if primary_intent == 'general_question' or not sequence_data:
                # For general questions without sequence data, try agent system first
                try:
                    from genomics_agents import create_genomics_agent_system
                    agent_system = create_genomics_agent_system()
                    result = await agent_system.process_query(
                        user_query=request.message,
                        sequence_data=sequence_data
                    )
                    response = result
                except Exception as agent_error:
                    print(f"Agent system failed, trying direct AI response: {agent_error}")
                    
                    # Try direct AI response using Groq API
                    try:
                        import requests
                        import json
                        
                        groq_api_key = os.getenv('GROQ_API_KEY')
                        
                        if not groq_api_key:
                            print("ERROR: GROQ_API_KEY not found in environment variables")
                            raise Exception("GROQ_API_KEY not found in environment variables")
                        
                        groq_response = requests.post(
                            "https://api.groq.com/openai/v1/chat/completions",
                            headers={
                                "Authorization": f"Bearer {groq_api_key}",
                                "Content-Type": "application/json"
                            },
                            json={
                                "model": request.model or "llama-3.1-8b-instant",
                                "messages": [
                                    {
                                        "role": "system", 
                                        "content": "You are a helpful genomics and bioinformatics assistant. Answer questions about genes, proteins, diseases, mutations, viruses, and biological processes. Be informative, accurate, and helpful."
                                    },
                                    {
                                        "role": "user",
                                        "content": user_message
                                    }
                                ],
                                "temperature": 0.7,
                                "max_tokens": 1024
                            },
                            timeout=30
                        )
                        
                        if groq_response.status_code == 200:
                            groq_data = groq_response.json()
                            ai_answer = groq_data["choices"][0]["message"]["content"]
                            
                            response = {
                                "message": ai_answer,
                                "model": request.model,
                                "status": "ai_response"
                            }
                        else:
                            print(f"DEBUG: Groq API response status: {groq_response.status_code}")
                            print(f"DEBUG: Groq API response: {groq_response.text}")
                            raise Exception(f"Groq API error: {groq_response.status_code} - {groq_response.text}")
                            
                    except Exception as groq_error:
                        print(f"Groq API failed: {groq_error}")
                        # Provide basic knowledge responses for common queries
                        user_lower = user_message.lower()
                        if 'brca1' in user_lower:
                            response = {
                                "message": """🧬 **BRCA1 (Breast Cancer 1 Gene)**

BRCA1 is a human tumor suppressor gene that produces a protein involved in DNA repair. Here are key facts:

🔬 **Function:**
• DNA damage repair and cell cycle regulation
• Maintains genomic stability
• Critical for homologous recombination repair

🏥 **Clinical Significance:**
• Mutations increase breast cancer risk by 55-65%
• Also increases ovarian cancer risk by 39-46%
• Associated with hereditary breast-ovarian cancer syndrome

🧪 **Molecular Details:**
• Located on chromosome 17 (17q21.31)
• Encodes 1,863 amino acid protein
• Contains RING finger and BRCT domains

💡 **Testing:** Genetic testing available for pathogenic variants. Preventive measures recommended for carriers.

Would you like me to search for recent literature about BRCA1 or analyze a BRCA1 sequence?""",
                                "model": request.model,
                                "status": "knowledge_response"
                            }
                        elif 'dengue' in user_lower:
                            response = {
                                "message": """🦟 **Dengue Virus**

Dengue virus is a mosquito-borne viral infection causing dengue fever. Key information:

🔬 **Virus Classification:**
• Family: Flaviviridae, Genus: Flavivirus  
• Four serotypes: DENV-1, DENV-2, DENV-3, DENV-4
• Single-stranded positive-sense RNA genome (~10.7 kb)

📊 **Epidemiology:**
• Endemic in tropical/subtropical regions
• ~390 million infections annually worldwide
• Transmitted by Aedes aegypti and Aedes albopictus mosquitoes

🏥 **Clinical Manifestations:**
• Dengue fever: fever, headache, muscle/joint pain
• Severe dengue: plasma leakage, bleeding, organ involvement
• Can be life-threatening if not properly managed

🧬 **Genome Structure:**
• Encodes 3 structural proteins (C, prM, E) 
• 7 non-structural proteins (NS1, NS2A, NS2B, NS3, NS4A, NS4B, NS5)
• Envelope (E) protein is main target for vaccines/therapeutics

Would you like me to analyze a dengue sequence or search for recent research?""",
                                "model": request.model,
                                "status": "knowledge_response"
                            }
                        else:
                            # Generic helpful response
                            response = {
                                "message": f"""I understand you asked about: "{request.message}"

I'm a genomics research assistant and can help with:
• **Genes & Proteins:** Information about specific genes like BRCA1, TP53, etc.
• **Diseases & Pathogens:** Cancer, viruses (like dengue), genetic disorders
• **Sequence Analysis:** Upload FASTA files for BLAST, mutation analysis
• **Literature Search:** Find recent papers on specific topics
• **Genetic Variants:** Look up known mutations and their significance

Please feel free to ask more specific questions or upload a sequence for analysis!""",
                                "model": request.model,
                                "status": "helpful_guidance"
                            }
                    else:
                        # Provide context-aware suggestions
                        suggestions = []
                        if search_terms:
                            suggestions.append(f"💡 I detected you mentioned: {', '.join(search_terms)}")
                            suggestions.append("Try uploading a sequence file and asking me to analyze it.")
                            suggestions.append("Or be more specific about what analysis you need.")
                        else:
                            suggestions.extend([
                                "💡 Here's what I can help you with:",
                                "• Answer questions about genes like BRCA1, TP53, etc.",
                                "• Upload a FASTA file and ask me to analyze it",
                                "• Search for literature about specific genes (e.g., 'Find BRCA1 papers')", 
                                "• Look up genetic variants for genes or sequences",
                                "• Run BLAST, mutation, or variant analysis"
                            ])
                        
                        response = {
                            "message": f"I understand you asked: '{request.message}'\n\n" + "\n".join(suggestions),
                            "model": request.model,
                            "status": "helpful_suggestions",
                            "detected_intent": primary_intent,
                            "confidence": intent_analysis['confidence']
                        }
            
            elif primary_intent in ['literature_search', 'variant_analysis', 'sequence_analysis']:
                # Provide specific guidance based on detected intent
                intent_guidance = {
                    'literature_search': "📚 I can help you search scientific literature! Try: 'Find literature for BRCA1' or upload a sequence and ask for related papers.",
                    'variant_analysis': "🧬 I can analyze genetic variants! Upload a sequence file or provide a gene name, then ask about variants.",
                    'sequence_analysis': "🔬 I can analyze sequences! Please upload a FASTA file and ask me to identify or analyze it."
                }
                
                response = {
                    "message": intent_guidance.get(primary_intent, "I'm ready to help with genomics analysis!"),
                    "model": request.model, 
                    "status": "intent_guidance",
                    "detected_intent": primary_intent,
                    "extracted_terms": search_terms
                }
        
        return {"status": "success", "response": response}
    except Exception as e:
        # Fallback response
        response = f"I'm here to help with genomics questions. You asked: '{request.message}'. Please try rephrasing or use the analysis tools."
        return {"status": "success", "response": response, "error": str(e)}

def format_analysis_response(analysis_results: dict, sequence: str, user_message: str = "") -> str:
    """Format analysis results into a human-readable response"""
    response_parts = []
    
    # Add sequence info
    seq_length = len(sequence)
    seq_type = "DNA/RNA" if set(sequence.upper()).issubset(set("ATCGNU")) else "protein"
    response_parts.append(f"🧬 **Sequence Analysis Complete**")
    response_parts.append(f"📏 Length: {seq_length} characters")
    response_parts.append(f"🔬 Type: {seq_type}")
    response_parts.append("")
    
    # Check if mutation analysis is the primary focus
    has_mutations = 'mutation' in analysis_results and analysis_results['mutation'].get('mutations')
    mutation_focus = 'mutation' in user_message.lower() or 'analyze' in user_message.lower()
    
    # Format mutation results first and prominently if available
    if 'mutation' in analysis_results:
        mutation_data = analysis_results['mutation']
        if mutation_data.get('error'):
            response_parts.append("🔬 **Mutation Analysis:** Error occurred during analysis")
            response_parts.append("")
        elif 'mutations' in mutation_data and mutation_data['mutations']:
            # Show reference information
            if 'ref_id' in mutation_data and 'ref_title' in mutation_data:
                response_parts.append(f"📍 **Reference Match:** {mutation_data['ref_id']}")
                response_parts.append(f"   {mutation_data['ref_title'][:80]}...")
                response_parts.append("")
            
            # Show alignment stats
            if 'alignment' in mutation_data:
                align = mutation_data['alignment']
                response_parts.append("📊 **Alignment Quality:**")
                response_parts.append(f"   Identity: {align.get('percent_identity', 0):.1f}% ({align.get('identity', 0)}/{align.get('length', 0)} positions)")
                response_parts.append(f"   Coverage: {align.get('query_coverage', 0):.1f}%")
                response_parts.append("")
            
            response_parts.append("� **Mutation Analysis:**")
            response_parts.append(f"Found {len(mutation_data['mutations'])} mutations compared to reference")
            
            # Show first few mutations as examples
            for i, mut in enumerate(mutation_data['mutations'][:5]):
                response_parts.append(f"   • Position {mut.get('position', 'N/A')}: {mut.get('ref', '?')} → {mut.get('user', '?')}")
            
            if len(mutation_data['mutations']) > 5:
                response_parts.append(f"   ... and {len(mutation_data['mutations']) - 5} more mutations")
            
            # Add reference to visualization
            response_parts.append("")
            response_parts.append("� **Rich Visualizations Available:**")
            response_parts.append("   • Interactive mutation scatter plots")
            response_parts.append("   • Mutation track diagrams") 
            response_parts.append("   • Detailed sequence alignment views")
            response_parts.append("   • Downloadable HTML reports")
            response_parts.append("")
        else:
            response_parts.append("✅ **Mutation Analysis:** Perfect match - no mutations detected!")
            if 'alignment' in mutation_data:
                align = mutation_data['alignment']
                response_parts.append(f"   Identity: {align.get('percent_identity', 0):.1f}% - Excellent sequence quality")
            response_parts.append("")
    
    # Only show BLAST results if mutation analysis is not the focus or failed
    if not mutation_focus and not has_mutations:
        if 'blast' in analysis_results:
            blast_data = analysis_results['blast']
            if 'hits' in blast_data and blast_data['hits']:
                response_parts.append("🎯 **Sequence Identification:**")
                response_parts.append(f"Top match from database search:")
                
                # Show only the top hit for simplicity
                hit = blast_data['hits'][0]
                accession = hit.get('accession', hit.get('Accession', 'N/A'))
                identity = hit.get('percent_identity', hit.get('Identity (%)', 0))
                title = hit.get('title', hit.get('Title', ''))[:100]
                response_parts.append(f"   **{accession}** - {identity}% identity")
                response_parts.append(f"   {title}...")
                response_parts.append("")
    
    # Format Literature results
    if 'literature' in analysis_results:
        lit_data = analysis_results['literature']
        if 'papers' in lit_data and lit_data['papers']:
            response_parts.append("📚 **Literature Search:**")
            response_parts.append(f"Found {len(lit_data['papers'])} relevant research papers:")
            
            for i, paper in enumerate(lit_data['papers'][:3]):  # Show top 3
                title = paper.get('title', 'Untitled')
                authors = paper.get('authors', 'Unknown authors')
                year = paper.get('date', '')[:4] if paper.get('date') else 'Unknown'
                pmid = paper.get('pmid', '')
                
                # Truncate long titles
                if len(title) > 80:
                    title = title[:77] + "..."
                
                response_parts.append(f"  {i+1}. **{title}** ({year})")
                response_parts.append(f"     Authors: {authors}")
                if pmid:
                    response_parts.append(f"     PMID: {pmid}")
            response_parts.append("")
    
    # Format Variant results
    if 'variants' in analysis_results:
        variant_data = analysis_results['variants']
        if 'variants' in variant_data and variant_data['variants']:
            response_parts.append("🧬 **Genetic Variants:**")
            response_parts.append(f"Found {len(variant_data['variants'])} known variants:")
            
            for i, variant in enumerate(variant_data['variants'][:3]):  # Show top 3
                var_id = variant.get('id', 'Unknown')
                var_type = variant.get('type', 'Unknown type')
                significance = variant.get('clinical_significance', 'Unknown')
                
                response_parts.append(f"  {i+1}. **{var_id}** - {var_type}")
                response_parts.append(f"     Clinical Significance: {significance}")
            response_parts.append("")
    
    response_parts.append("💡 View detailed results with interactive visualizations in the Analysis Results section below!")
    
    return "\n".join(response_parts)

@app.get("/models")
async def get_available_models():
    """Get list of available chat models"""
    models = [
        "kimi-k2",
        "meta-llama/llama-4-scout-17b-16e-instruct", 
        "llama-3.1-8b-instant",
        "deepseek-r1-distill-llama-70b",
        "llama-3.3-70b-versatile"
    ]
    return {"models": models}

if __name__ == "__main__":
    import uvicorn
    uvicorn.run(app, host="0.0.0.0", port=8000)