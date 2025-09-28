"""
Enhanced Streamlit Frontend for Genomics Research Assistant
Integrates with FastAPI backend and agent system
"""

# First load environment variables from .env.local
from env_loader import load_env_file, ensure_api_keys
print("Loading environment variables...")
loaded_env = load_env_file(".env.local")
print(f"Loaded {len(loaded_env)} environment variables")

# Now check if the API key is available
api_keys_available = ensure_api_keys()
if api_keys_available:
    print("Required API keys found!")
    print(f"GROQ_API_KEY present: {'GROQ_API_KEY' in loaded_env}")
else:
    print("WARNING: Required API keys are missing!")
    print("Please create a .env.local file with your API keys")
    print("Copy .env.local.template to .env.local and add your keys")
    print("Or use the Debug Tools in the sidebar to troubleshoot")
    
# Import debug tools
try:
    from app.debug_tools import display_debug_info, display_enhanced_debug_tools
except ImportError:
    # Define fallback functions if import fails
    def display_debug_info():
        if st.session_state.get('show_debug', False):
            st.sidebar.error("Debug tools module not found")
            
    def display_enhanced_debug_tools():
        if st.session_state.get('show_enhanced_debug', False):
            with st.sidebar.expander("🔍 Advanced Debug Tools"):
                st.subheader("RAG System")
                
                if st.button("Test RAG System"):
                    if rag_system is not None:
                        with st.spinner("Testing RAG system..."):
                            try:
                                test_result = rag_system.get_gene_context("CFTR")
                                st.text_area("Results", test_result, height=200)
                                if "Error" not in test_result:
                                    st.success("RAG system working correctly")
                                else:
                                    st.warning("RAG system returned an error")
                            except Exception as e:
                                st.error(f"RAG system test failed: {str(e)}")
                    else:
                        st.error("RAG system not initialized")
                
                if st.button("Reload Knowledge Base"):
                    if rag_system is not None:
                        with st.spinner("Reloading knowledge base..."):
                            try:
                                success = rag_system.reload_knowledge_base()
                                if success:
                                    st.success("Knowledge base reloaded successfully")
                                else:
                                    st.error("Failed to reload knowledge base")
                            except Exception as e:
                                st.error(f"Error reloading knowledge base: {str(e)}")
                    else:
                        st.error("RAG system not initialized")
                
                st.subheader("Collection Status")
                if rag_system is not None:
                    try:
                        collections = {
                            "genes_collection": rag_system.genes_collection.count(),
                            "diseases_collection": rag_system.diseases_collection.count(),
                            "mutations_collection": rag_system.mutations_collection.count(),
                            "literature_collection": rag_system.literature_collection.count()
                        }
                        for name, count in collections.items():
                            st.text(f"{name}: {count} items")
                    except Exception as e:
                        st.error(f"Error accessing collections: {str(e)}")
                else:
                    st.error("RAG system not initialized")

import streamlit as st
import requests
import json
import asyncio
import threading
import random
import time
from datetime import datetime, timedelta
import os
import socket
import time
from pathlib import Path
from typing import Dict, Optional, List, Any
import plotly.express as px
import pandas as pd
from Bio import Entrez, SeqIO
from Bio.Blast import NCBIWWW, NCBIXML
from xml.etree import ElementTree as ET
from tabulate import tabulate
import plotly.express as px
import plotly.graph_objects as go
# Import functions for sequence analysis
from app import detect_sequence_type, run_blast_top_hits
from risk_assessment import InfectionRiskAssessment
from display_risk_assessment import display_risk_assessment, display_mini_risk_assessment

# Import RAG and Agent systems - with fallbacks if imports fail
try:
    from genomics_rag import GenomicsRAG
    from genomics_agents import GenomicsAgentSystem
    from rag_agent_helper import use_rag_system_with_chat, use_agent_system_with_chat
    AGENT_IMPORTS_AVAILABLE = True
except ImportError as e:
    print(f"Warning: Advanced AI features unavailable - {str(e)}")
    # Create dummy classes/functions to avoid errors
    class GenomicsRAG:
        def __init__(self, *args, **kwargs):
            pass
        def initialize_knowledge_base(self):
            pass
        def query_knowledge(self, *args, **kwargs):
            return {"results": []}
        def get_mutation_context(self, *args, **kwargs):
            return ""
    
    class GenomicsAgentSystem:
        def __init__(self, *args, **kwargs):
            pass
        def register_bioinformatics_functions(self):
            pass
    
    def use_rag_system_with_chat(message, sequence, rag_system):
        return {"message": message, "context": ""}
    
    def use_agent_system_with_chat(message, sequence, agent_system):
        return {"status": "error", "response": "Agent system not available"}
    
    AGENT_IMPORTS_AVAILABLE = False

# Configuration
API_BASE_URL = "http://localhost:8000"
Entrez.email = "your.email@example.com"  # Replace with your real email

# Initialize RAG and Agent systems
rag_system = None
agent_system = None

def initialize_ai_systems():
    """Initialize RAG and Agent systems if not already initialized"""
    global rag_system, agent_system
    
    # Check if imports are available
    if not AGENT_IMPORTS_AVAILABLE:
        st.warning("Advanced AI features unavailable due to missing dependencies.")
        return
    
    # Initialize RAG system
    try:
        if rag_system is None:
            with st.spinner("Initializing knowledge base..."):
                # Initialize RAG system with timeout protection
                try:
                    print("Creating GenomicsRAG instance...")
                    
                    # First, check if we have a working database path
                    import os
                    import glob
                    
                    # Try to find the most recent genomics_db directory
                    db_dirs = glob.glob("./genomics_db*")
                    if not db_dirs:
                        print("No genomics_db directories found, using default")
                        db_path = "./genomics_db"
                    else:
                        # Sort by modification time, newest first
                        db_dirs.sort(key=os.path.getmtime, reverse=True)
                        db_path = db_dirs[0]
                        print(f"Using most recent database: {db_path}")
                    
                    # Initialize with a 30-second timeout
                    import threading
                    import time
                    
                    initialization_complete = False
                    temp_rag_system = None
                    
                    def initialize_rag():
                        nonlocal initialization_complete, temp_rag_system
                        try:
                            # Try to initialize the RAG system
                            temp_rag = GenomicsRAG(db_path=db_path)
                            # Set flag to indicate initialization is done
                            temp_rag_system = temp_rag
                            initialization_complete = True
                        except Exception as e:
                            print(f"Error in initialization thread: {str(e)}")
                            initialization_complete = True  # Mark as complete even if it failed
                    
                    # Start initialization in separate thread
                    init_thread = threading.Thread(target=initialize_rag)
                    init_thread.daemon = True
                    init_thread.start()
                    
                    # Wait for initialization with timeout
                    timeout_seconds = 30
                    start_time = time.time()
                    
                    while not initialization_complete and (time.time() - start_time) < timeout_seconds:
                        time.sleep(0.5)
                        # Show progress
                        if (int(time.time() - start_time) % 5) == 0:
                            print(f"Waiting for RAG initialization... ({int(time.time() - start_time)}s)")
                    
                    if not initialization_complete:
                        print(f"RAG initialization timed out after {timeout_seconds} seconds")
                        raise TimeoutError(f"RAG initialization timed out after {timeout_seconds} seconds")
                    
                    if temp_rag_system is None:
                        raise Exception("RAG initialization failed")
                        
                    # Now assign to the global variable
                    rag_system = temp_rag_system
                        
                    print("RAG system initialized, now initializing knowledge base...")
                    
                    # Initialize knowledge base with timeout protection
                    initialization_complete = False
                    
                    def init_knowledge():
                        nonlocal initialization_complete
                        try:
                            rag_system.initialize_knowledge_base()
                            initialization_complete = True
                        except Exception as e:
                            print(f"Error initializing knowledge base: {str(e)}")
                            initialization_complete = True  # Mark as complete even if it failed
                    
                    # Start knowledge base initialization in separate thread
                    init_thread = threading.Thread(target=init_knowledge)
                    init_thread.daemon = True
                    init_thread.start()
                    
                    # Wait for knowledge base initialization with timeout
                    start_time = time.time()
                    
                    while not initialization_complete and (time.time() - start_time) < timeout_seconds:
                        time.sleep(0.5)
                        # Show progress
                        if (int(time.time() - start_time) % 5) == 0:
                            print(f"Waiting for knowledge base initialization... ({int(time.time() - start_time)}s)")
                    
                    if not initialization_complete:
                        print(f"Knowledge base initialization timed out after {timeout_seconds} seconds")
                        # Continue anyway, might have partial functionality
                    
                    # Test RAG system by querying for a known gene
                    try:
                        test_result = rag_system.get_gene_context("CFTR")
                        if test_result and "Error" not in test_result:
                            print("RAG system initialized and tested successfully")
                        else:
                            print("RAG system initialized but test query failed")
                    except Exception as e:
                        print(f"Test query failed: {str(e)}")
                        
                except Exception as inner_e:
                    import traceback
                    print(f"Error initializing RAG system: {str(inner_e)}")
                    print(traceback.format_exc())
                    st.error(f"Failed to initialize knowledge base: {str(inner_e)}")
                    # Create a minimal RAG system for fallback - with simple functionality
                    try:
                        print("Creating minimal fallback RAG system...")
                        
                        # Create a minimal class to handle RAG functionality
                        class MinimalRAG:
                            def get_gene_context(self, gene):
                                return f"Gene information for {gene} is not available (RAG system is in minimal mode)"
                                
                            def get_mutation_context(self, mutation):
                                return f"Mutation information for {mutation} is not available (RAG system is in minimal mode)"
                        
                        rag_system = MinimalRAG()
                        print("Minimal RAG system created")
                    except:
                        rag_system = None
    except Exception as e:
        import traceback
        print(f"Error in RAG initialization outer block: {str(e)}")
        print(traceback.format_exc())
        st.warning(f"Warning: Knowledge base initialization issue: {str(e)}")
        rag_system = None
        # We'll continue without RAG if there's an error
    
    # Try to initialize agent system
    try:
        if agent_system is None:
            with st.spinner("Setting up AI assistants..."):
                # Check for AutoGen dependencies
                try:
                    # Verify autogen import
                    from autogen import AssistantAgent
                except ImportError:
                    st.warning("Warning: AutoGen not available. Agent system will be disabled.")
                    return
                
                # Use our custom Groq compatibility module
                from groq_compat import create_groq_config, verify_groq_api_key
                
                # First check if API key is available and valid
                if verify_groq_api_key():
                    # Create a configuration using the helper function
                    llm_config = create_groq_config("llama-3.1-8b-instant")
                    
                    # Show message to user
                    if llm_config:
                        st.info("Using Groq API as the backend for AI assistants")
                    else:
                        st.warning("Issue creating Groq config - agent system may not work")
                        llm_config = None
                else:
                    st.warning("Groq API key not found or invalid - agent system disabled")
                    # Set to None to skip initialization
                    llm_config = None
                
                # Check if we have a valid llm_config with proper API key
                if not llm_config or not llm_config.get("config_list") or not llm_config["config_list"][0].get("api_key"):
                    st.warning("Cannot initialize agent system - missing API key configuration")
                    agent_system = None
                else:
                    # Create and initialize the agent system
                    try:
                        # Pass the selected model to the agent system
                        agent_system = GenomicsAgentSystem(llm_config)
                        agent_system.register_bioinformatics_functions()
                        st.success("Agent system initialized successfully")
                    except Exception as e:
                        st.warning(f"Agent system initialization warning: {str(e)}")
                        agent_system = None
    except Exception as e:
        st.warning(f"Warning: Agent system initialization failed: {str(e)}")
        # Continue without agent system

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
    if "rag_initialized" not in st.session_state:
        st.session_state.rag_initialized = False
    if "blast_results" not in st.session_state:
        st.session_state.blast_results = []
    if "show_debug" not in st.session_state:
        st.session_state.show_debug = False
    if "show_enhanced_debug" not in st.session_state:
        st.session_state.show_enhanced_debug = False
    if "rate_limit_errors" not in st.session_state:
        st.session_state.rate_limit_errors = {}  # Track rate limit errors by model
    if "api_calls" not in st.session_state:
        st.session_state.api_calls = []  # Track recent API calls for throttling
    if "show_detailed_risk" not in st.session_state:
        st.session_state.show_detailed_risk = False


def debug_environment():
    """Debug function to check environment variables and configuration"""
    # Show environment variables status
    st.subheader("🔧 Environment Variables")
    
    # Check .env.local file
    env_file = Path(".env.local")
    if env_file.exists():
        st.success(f"✅ .env.local file found at {env_file.absolute()}")
        
        # Try to read the file
        try:
            with open(env_file, "r") as f:
                content = f.read()
                # Mask API keys for security
                masked_content = re.sub(r'(API_KEY=")([^"]{5})([^"]+)(")', r'\1\2***\4', content)
                masked_content = re.sub(r'(API_KEY=)([^"][^\n]{5})([^\n]+)', r'\1\2***', masked_content)
                st.code(masked_content)
        except Exception as e:
            st.error(f"❌ Error reading .env.local file: {str(e)}")
    else:
        st.error(f"❌ .env.local file not found at {env_file.absolute()}")
    
    # Check if specific API keys are in environment
    required_keys = ["GROQ_API_KEY", "HUGGINGFACE_API_KEY"]
    for key in required_keys:
        value = os.environ.get(key, "")
        if value:
            masked_value = value[:5] + "***" + value[-4:] if len(value) > 9 else "***"
            st.success(f"✅ {key} found in environment: {masked_value}")
        else:
            st.error(f"❌ {key} not found in environment variables")
    
    # Show current directory and path information
    st.subheader("📁 File System Info")
    st.info(f"Current working directory: {os.getcwd()}")
    
    # Check AutoGen installation
    st.subheader("🤖 AutoGen Status")
    try:
        import autogen
        st.success(f"✅ AutoGen installed: version {autogen.__version__}")
    except ImportError:
        st.error("❌ AutoGen not installed")
    except Exception as e:
        st.error(f"❌ Error checking AutoGen: {str(e)}")

def perform_sequence_analysis(sequence: str) -> Dict:
    """Perform comprehensive sequence analysis using FastAPI backend"""
    try:
        # First run BLAST analysis
        blast_response = requests.post(
            "http://localhost:8000/blast",
            json={"sequence": sequence, "num_hits": 5}
        )
        
        # Get mutation analysis
        mutation_response = requests.post(
            "http://localhost:8000/mutation",
            json={"sequence": sequence}
        )
        
        # Search literature (optional)
        literature_data = None
        try:
            # Try to get top hits to use for literature search
            if blast_response.ok:
                blast_data = blast_response.json()
                top_hits = blast_data.get("data", {}).get("hits", [])
                if top_hits and len(top_hits) > 0:
                    # Use top hit title for literature search
                    search_term = top_hits[0].get("title", "").split("[")[0].strip()
                    if search_term:
                        literature_response = requests.post(
                            "http://localhost:8000/literature",
                            json={"query": search_term, "max_results": 5}
                        )
                        if literature_response.ok:
                            literature_data = literature_response.json().get("data", {})
        except Exception as e:
            print(f"Literature search failed: {str(e)}")
            # Continue without literature data
        
        # Combine results
        results = {
            "blast": blast_response.json() if blast_response.ok else {"status": "error"},
            "mutation": mutation_response.json() if mutation_response.ok else {"status": "error"}
        }
        
        # Add risk assessment
        try:
            blast_data = results["blast"].get("data", {})
            mutation_data = results["mutation"].get("data", {})
            
            # Create risk assessment object
            risk_assessor = InfectionRiskAssessment()
            
            # Calculate risk
            risk_results = risk_assessor.calculate_infection_risk(
                blast_data, 
                mutation_data,
                literature_data
            )
            
            # Add to results
            results["risk_assessment"] = {
                "status": "success",
                "data": risk_results
            }
        except Exception as e:
            print(f"Risk assessment failed: {str(e)}")
            results["risk_assessment"] = {"status": "error", "message": str(e)}
        
        # Format the results into a readable message
        blast_hits = results["blast"].get("data", {}).get("hits", [])
        mutations = results["mutation"].get("data", {}).get("mutations", [])
        
        response = "🧬 **Sequence Analysis Results**\n\n"
        
        # Add BLAST results
        response += "**BLAST Results:**\n"
        if blast_hits:
            for i, hit in enumerate(blast_hits[:3], 1):
                response += f"{i}. **{hit.get('title', 'Unknown Sequence')}**\n"
                response += f"   - Identity: {hit.get('identity', 'N/A')}%\n"
                response += f"   - E-value: {hit.get('evalue', 'N/A')}\n"
                response += f"   - Accession: {hit.get('accession', 'N/A')}\n\n"
        else:
            response += "No significant BLAST hits found\n\n"
        
        # Add mutation results if available
        if mutations:
            response += "**Mutation Analysis:**\n"
            for mut in mutations[:5]:
                response += f"- Position {mut['position']}: {mut['reference']} → {mut['variant']}\n"
            response += "\n"
        
        return {
            "status": "success",
            "response": response
        }
    except Exception as e:
        return {
            "status": "error",
            "response": f"Analysis failed: {str(e)}"
        }

def get_fallback_model(current_model: str) -> str:
    """Get a fallback model when the current model is rate limited"""
    # Define model tiers (from least to most likely to be rate limited)
    model_tiers = {
        "tier1": ["claude-3-haiku-20240307"],  # Usually less rate limited (different API)
        "tier2": ["llama-3.1-8b-instant"],   # Small models, less likely to be rate limited
        "tier3": ["kimi-k2", "deepseek-r1-distill-llama-70b"],  # Medium models
        "tier4": ["meta-llama/llama-4-scout-17b-16e-instruct", "llama-3.3-70b-versatile"]  # Large models
    }
    
    # Find current model tier
    current_tier = None
    for tier, models in model_tiers.items():
        if current_model in models:
            current_tier = tier
            break
    
    # If model not found in tiers or is in tier1 already, use a tier1 model
    if not current_tier or current_tier == "tier1":
        return random.choice(model_tiers["tier1"])
    
    # Otherwise choose a model from a lower tier (less likely to be rate limited)
    if current_tier == "tier4":
        return random.choice(model_tiers["tier3"] + model_tiers["tier2"])
    elif current_tier == "tier3":
        return random.choice(model_tiers["tier2"] + model_tiers["tier1"])
    else:  # tier2
        return random.choice(model_tiers["tier1"])

def call_chat_api(message: str, model: str = "llama-3.1-8b-instant", sequence: str = None) -> Dict:
    """Make API calls to Node.js chat backend or sequence analysis based on message
    
    Enhanced with RAG and Agent systems if available
    Includes error handling and rate limit management with fallback models
    """
    # Track API call start time to calculate latency
    start_time = time.time()
    
    try:
        # Check for analysis commands
        if sequence and message in ["ANALYZE_SEQUENCE", "RUN_BLAST", "FIND_LITERATURE"]:
            if message == "ANALYZE_SEQUENCE":
                return perform_sequence_analysis(sequence)
            elif message == "RUN_BLAST":
                response = requests.post(
                    "http://localhost:8000/blast",
                    json={"sequence": sequence, "num_hits": 5}
                )
            elif message == "FIND_LITERATURE":
                response = requests.post(
                    "http://localhost:8000/pubmed",
                    json={"sequence": sequence, "max_results": 10}
                )
                
            if response.ok:
                result = response.json()
                return {
                    "status": "success",
                    "response": format_analysis_result(result, message)
                }
        
        # First try to use the agent system if available for normal queries
        global agent_system, rag_system
        
        # Initialize systems if not already done
        if not rag_system or not agent_system:
            initialize_ai_systems()
        
        # Try to use RAG system to enhance the response with relevant context
        rag_enhanced = None
        if rag_system:
            try:
                st.info("Enhancing response with genomics knowledge base...")
                rag_enhanced = use_rag_system_with_chat(message, sequence, rag_system)
                enhanced_message = rag_enhanced.get("message", message)
                context = rag_enhanced.get("context", "")
                
                # If we have context, include it in our log
                if context:
                    print(f"Using RAG context: {context[:100]}...")
                    st.success("Found relevant genomics information to enhance the response.")
                    
                # If context was found, update the message with the context
                message_with_context = message
                if context:
                    message_with_context = f"""User query: {message}
                    
Additional context for answering:
{context}

Please answer the user's query using this additional context when relevant."""
                    message = message_with_context
            except Exception as e:
                print(f"RAG enhancement error (continuing without it): {str(e)}")
                st.warning(f"Knowledge base enhancement encountered an issue: {str(e)}")
        
        # Try to use agent system if available
        if agent_system:
            try:
                st.info("Consulting specialized AI assistants for your query...")
                agent_result = use_agent_system_with_chat(message, sequence, agent_system)
                if agent_result.get("status") == "success":
                    st.success("Response generated by AI assistant team.")
                    return {"status": "success", "response": agent_result.get("response", "")}
                else:
                    st.warning("AI assistants couldn't process this query, falling back to general response.")
            except Exception as e:
                print(f"Agent system error (falling back to chat API): {str(e)}")
                st.warning("AI assistant team unavailable, using standard chat model instead.")
        
        # Fall back to Node.js chat API with rate limit handling
        try:
            # First check if we've had recent rate limit errors with this model
            if hasattr(st.session_state, "rate_limit_errors") and model in st.session_state.rate_limit_errors:
                last_error_time, count = st.session_state.rate_limit_errors[model]
                time_since_error = time.time() - last_error_time
                
                # If we've had multiple rate limit errors recently, enforce a cooldown period
                if count > 2 and time_since_error < 60:  # 60 second cooldown after 3 errors
                    wait_time = 60 - time_since_error
                    st.warning(f"⏳ Rate limit cooldown: Please wait {int(wait_time)} more seconds before trying again.")
                    return {
                        "status": "error", 
                        "message": f"Rate limit cooldown: Please wait {int(wait_time)} more seconds before trying again.",
                        "rate_limited": True
                    }
            
            # Add API key and environment info to help diagnose issues
            api_keys_info = {
                "GROQ_API_KEY": "✅ Present" if os.environ.get("GROQ_API_KEY") else "❌ Missing",
                "HUGGINGFACE_API_KEY": "✅ Present" if os.environ.get("HUGGINGFACE_API_KEY") else "❌ Missing"
            }
            
            # Initialize rate limit tracking if not exists
            if not hasattr(st.session_state, "rate_limited_models"):
                st.session_state.rate_limited_models = {}
            
            # Check if current model was recently rate limited
            current_time = datetime.now()
            cooldown_minutes = 10  # Time to wait before trying a rate-limited model again
            
            if model in st.session_state.rate_limited_models:
                limited_time = st.session_state.rate_limited_models[model]
                time_since_limit = (current_time - limited_time).total_seconds() / 60
                
                # If model was recently rate limited, use a fallback model
                if time_since_limit < cooldown_minutes:
                    original_model = model
                    model = get_fallback_model(model)
                    st.warning(f"Model {original_model} was recently rate limited. Using {model} instead.")
                    print(f"Switching from {original_model} to {model} due to recent rate limit")
            
            # Make the API call with timeout and retry logic
            max_retries = 2
            retry_delay = 5  # Start with 5 seconds
            
            for attempt in range(max_retries):
                try:
                    # Add timestamp to help with rate limit debugging
                    timestamp = current_time.strftime("%Y-%m-%d %H:%M:%S")
                    
                    st.info(f"Requesting response from {model} (attempt {attempt+1}/{max_retries})...")
                    
                    response = requests.post(
                        "http://localhost:3000/chat",
                        json={
                            "message": message, 
                            "model": model,
                            "sequence": sequence,  # Pass sequence data to the API
                            "timestamp": timestamp  # For tracking requests
                        },
                        timeout=60  # Extended timeout (60 seconds)
                    )
                    
                    # Calculate response latency
                    latency = time.time() - start_time
                    
                    if response.ok:
                        # Update rate limit tracking - this model worked!
                        if model in st.session_state.rate_limited_models:
                            del st.session_state.rate_limited_models[model]
                        
                        return {
                            "status": "success", 
                            "response": response.json().get("response", ""),
                            "latency": f"{latency:.2f}s",
                            "model_used": model
                        }
                    
                    # If we got a rate limit error, mark this model as rate limited
                    if response.status_code == 429:
                        st.session_state.rate_limited_models[model] = current_time
                        
                        # Try a fallback model if not the last attempt
                        if attempt < max_retries - 1:
                            fallback_model = get_fallback_model(model)
                            st.warning(f"Rate limit detected. Switching to {fallback_model}...")
                            model = fallback_model
                            time.sleep(2)  # Short delay before trying new model
                            continue
                    else:
                        # For non-rate limit errors, just break
                        break
                        
                except (requests.exceptions.ConnectionError, requests.exceptions.Timeout) as e:
                    if attempt < max_retries - 1:
                        st.info(f"Connection error. Retrying in {retry_delay} seconds... (Attempt {attempt+1}/{max_retries})")
                        time.sleep(retry_delay)
                        retry_delay *= 2  # Exponential backoff
                    else:
                        return {
                            "status": "error",
                            "message": f"Connection error after {max_retries} attempts: {str(e)}"
                        }
            
            # If we get here, all retries failed or non-429 error occurred
            if response.status_code == 429:
                # Handle rate limit error - store it for future reference
                now = datetime.now()
                st.session_state.rate_limited_models[model] = now
                
                # Suggest available models
                available_models = ["llama-3.1-8b-instant", "kimi-k2", "claude-3-haiku-20240307"]
                suggested_model = random.choice(available_models)
                
                return {
                    "status": "error",
                    "message": f"Error: Unable to fetch response due to rate limiting. Please try again in a few minutes or switch to another model like {suggested_model}.",
                    "rate_limited": True
                }
                
                # This code is unreachable since we already returned
                # Keep for reference
                """
                # Update rate limit tracking
                now = time.time()
                if model in st.session_state.rate_limit_errors:
                    _, count = st.session_state.rate_limit_errors[model]
                    st.session_state.rate_limit_errors[model] = (now, count + 1)
                else:
                    st.session_state.rate_limit_errors[model] = (now, 1)
                """
            else:
                # For other errors, return appropriate message
                if response.status_code == 429:
                    return {
                        "status": "error",
                        "message": f"Error: Unable to fetch response from Groq API. API Keys configured:\n- GROQ_API_KEY: {api_keys_info['GROQ_API_KEY']}\n- HUGGINGFACE_API_KEY: {api_keys_info['HUGGINGFACE_API_KEY']}\n\nGroq API Error: Request failed with status code 429 (rate limit exceeded).\n\nTry again in a few minutes or switch to a different model.",
                        "rate_limited": True
                    }
                else:
                    return {
                        "status": "error", 
                        "message": f"Chat API error: {response.status_code} - {response.text}\n\nAPI Keys configured:\n- GROQ_API_KEY: {api_keys_info['GROQ_API_KEY']}\n- HUGGINGFACE_API_KEY: {api_keys_info['HUGGINGFACE_API_KEY']}"
                    }
        except requests.exceptions.Timeout:
            return {
                "status": "error",
                "message": "Request timed out. The server may be busy or experiencing issues."
            }
        except requests.exceptions.ConnectionError:
            return {
                "status": "error",
                "message": "Connection error. Make sure the chat server is running on port 3000."
            }
    except Exception as e:
        return {"status": "error", "message": f"Error: {str(e)}"}

def call_api(endpoint: str, data: Dict) -> Dict:
    """Make API calls to FastAPI backend for non-chat endpoints"""
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
        
        # Ensure we return dictionaries with consistent keys
        for i, record in enumerate(records.get("PubmedArticle", [])):
            try:
                article = record.get("MedlineCitation", {}).get("Article", {})
                
                # Extract title
                title = article.get("ArticleTitle", "Untitled")
                
                # Extract authors
                author_list = article.get("AuthorList", [])
                authors = []
                for author in author_list:
                    if "LastName" in author and "ForeName" in author:
                        authors.append(f"{author['LastName']} {author['ForeName'][0]}")
                authors_str = ", ".join(authors) if authors else "Unknown"
                
                # Extract journal and year
                journal_info = article.get("Journal", {})
                journal_name = journal_info.get("Title", "Unknown Journal")
                
                pub_date = None
                if "PubDate" in journal_info.get("JournalIssue", {}):
                    pub_date = journal_info["JournalIssue"]["PubDate"]
                    
                year = pub_date.get("Year", "") if pub_date else ""
                
                # Extract abstract
                abstract_text = ""
                if "Abstract" in article and "AbstractText" in article["Abstract"]:
                    abstract_text = " ".join([str(text) for text in article["Abstract"]["AbstractText"]])
                
                papers.append({
                    "pmid": pmids[i],
                    "title": title,
                    "authors": authors_str,
                    "journal": journal_name,
                    "year": year,
                    "abstract": abstract_text
                })
            except Exception as e:
                # If parsing a specific record fails, add a placeholder
                papers.append({
                    "pmid": pmids[i] if i < len(pmids) else "Unknown",
                    "title": "Error parsing publication details",
                    "authors": "Unknown",
                    "journal": "Unknown",
                    "year": "",
                    "abstract": ""
                })
                continue
        
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
    # Define available models
    available_models = [
        "llama-3.1-8b-instant",  # Default and reliable model
        "llama-3.3-70b-versatile",
        "deepseek-r1-distill-llama-70b",
        "meta-llama/llama-4-scout-17b-16e-instruct",
        "microsoft/BioGPT-Large"  # Microsoft's BioGPT model for biomedical text
    ]
    
    # Model selection
    selected_model = st.selectbox(
        "Choose AI Model", 
        available_models,
        index=0,  # Default to llama-3.1-8b-instant
        help="Select the AI model for conversations"
    )
    
    # Model-specific information can be added here if needed
    # (BioGPT model option has been removed)
    
    # Chat input with file context
    if st.session_state.current_sequence:
        seq_preview = st.session_state.current_sequence[:50] + "..." if len(st.session_state.current_sequence) > 50 else st.session_state.current_sequence
        st.info(f"📄 Current sequence available: {seq_preview} ({len(st.session_state.current_sequence)} chars)")
    
    user_input = st.chat_input("Ask me about genomics, genes, mutations, or request analysis...")
    
    # Quick action buttons for sequence analysis
    if st.session_state.current_sequence:
        col1, col2, col3 = st.columns(3)
        with col1:
            if st.button("� Analyze Sequence", key="chat_analyze", help="Run comprehensive analysis including BLAST search and mutation detection"):
                with st.spinner("🔬 Analyzing sequence..."):
                    try:
                        seq = st.session_state.current_sequence
                        
                        # Create a user-friendly analysis summary
                        analysis_summary = "📊 **Sequence Analysis Summary**\n\n"
                        
                        # Basic sequence stats
                        length = len(seq)
                        gc_content = (seq.upper().count('G') + seq.upper().count('C')) / length * 100
                        
                        analysis_summary += f"I analyzed your sequence of {length} nucleotides.\n\n"
                        analysis_summary += f"**Key characteristics:**\n"
                        analysis_summary += f"• Length: {length} nucleotides\n"
                        analysis_summary += f"• GC content: {gc_content:.1f}%\n"
                        
                        # Add a note about what was found
                        if gc_content > 60:
                            analysis_summary += f"\nThe high GC content ({gc_content:.1f}%) suggests this sequence may come from a GC-rich organism or region.\n"
                        elif gc_content < 35:
                            analysis_summary += f"\nThe low GC content ({gc_content:.1f}%) suggests this sequence may come from an AT-rich organism or region.\n"
                        
                        analysis_summary += "\nTo learn more about this sequence, try using the BLAST Search button to find similar sequences, or ask me specific questions about genes, mutations, or other features you're interested in."
                        
                        # Add results to chat
                        st.session_state.messages.append({"role": "assistant", "content": analysis_summary})
                    except Exception as e:
                        error_message = f"❌ Error during sequence analysis: {str(e)}"
                        st.session_state.messages.append({"role": "assistant", "content": error_message})
                    result = perform_sequence_analysis(st.session_state.current_sequence)
                    if result["status"] == "success":
                        st.session_state.messages.append({"role": "assistant", "content": result["response"]})
                    else:
                        st.error(result["response"])
                st.rerun()
                
        with col2:
            if st.button("🔍 BLAST Search", key="chat_blast", help="Search NCBI database for similar sequences"):
                with st.spinner("Running enhanced BLAST search with knowledge retrieval..."):
                    try:
                        # Get the sequence
                        seq = st.session_state.current_sequence
                        
                        # Define local detect_sequence_type function
                        def local_detect_sequence_type(seq):
                            protein_set = set("ACDEFGHIKLMNPQRSTVWY")
                            nucleotide_set = set("ACGTU")
                            seq_set = set(seq.upper())
                            if seq_set.issubset(nucleotide_set) or len(seq_set & nucleotide_set) / len(seq_set) > 0.95:
                                return "nucleotide"
                            elif seq_set.issubset(protein_set):
                                return "protein"
                            else:
                                return "nucleotide"  # Default to nucleotide
                        
                        # Determine sequence type
                        seq_type = local_detect_sequence_type(seq)
                        
                        # Use our robust BLAST search implementation with fallbacks
                        from robust_blast import RobustBLASTSearch
                        
                        # Create robust BLAST searcher with multiple fallback mechanisms
                        blast_searcher = RobustBLASTSearch(email=Entrez.email)
                        
                        # Run the search with retry and fallback mechanisms
                        with st.status("Running BLAST search with multiple fallback mechanisms..."):
                            st.write("Searching sequence databases...")
                            blast_result = blast_searcher.run_blast_search(
                                sequence=seq, 
                                max_retries=2,
                                hitlist_size=5
                            )
                            st.write(f"Search completed via {blast_searcher.search_method_used} method")
                        
                        # Process and save results for RAG
                        blast_results = blast_result["data"]["hits"]
                        
                        # Store results for future use
                        st.session_state.blast_results = blast_results
                        
                        # Show whether we got actual BLAST results or fallback results
                        if blast_searcher.search_method_used in ["direct_ncbi", "local_api"]:
                            st.success("Successfully retrieved BLAST results")
                        elif blast_searcher.search_method_used == "alternative_provider":
                            st.info("Retrieved BLAST results from alternative provider")
                        else:
                            st.warning("Network issues prevented BLAST search. Using offline analysis.")
                        
                        # Use RAG to enhance results with contextual knowledge
                        blast_summary = "🧬 **BLAST Search Results:**\n\n"
                        if blast_results:
                            # Extract organism name from first hit
                            first_hit_title = blast_results[0]["title"]
                            organism_name = ""
                            
                            # Try to extract a clean organism name
                            import re
                            match = re.search(r'\[(.*?)\]', first_hit_title)
                            if match:
                                organism_name = match.group(1)
                            else:
                                # Alternative approach - get the most relevant part
                                parts = first_hit_title.split('|')
                                if len(parts) >= 3:
                                    organism_name = parts[-1].split(',')[0].strip()
                                else:
                                    words = first_hit_title.split()
                                    if len(words) >= 3:
                                        organism_name = ' '.join(words[1:4])
                            
                            # Get the top match identity
                            # Use the percent_identity that's already calculated in blast_results
                            top_identity = blast_results[0]["percent_identity"]
                                
                            # Human-friendly summary
                            if top_identity > 95:
                                match_quality = "very strong"
                            elif top_identity > 85:
                                match_quality = "strong"
                            elif top_identity > 70:
                                match_quality = "moderate"
                            else:
                                match_quality = "weak"
                            
                            # Create a conversational response
                            if "virus" in first_hit_title.lower() or "viral" in first_hit_title.lower():
                                blast_summary += f"I found a **{match_quality} match** ({top_identity:.1f}%) to sequences from the **{organism_name or 'virus'}**.\n\n"
                                
                                if top_identity > 90:
                                    blast_summary += f"This suggests your sequence is very likely from this virus or a closely related strain. "
                                    blast_summary += f"The extremely high similarity ({top_identity:.1f}%) indicates this is almost certainly the correct identification.\n\n"
                                elif top_identity > 80:
                                    blast_summary += f"This suggests your sequence is probably from this virus or a related strain. "
                                    blast_summary += f"The high similarity ({top_identity:.1f}%) makes this a reliable match.\n\n"
                                else:
                                    blast_summary += f"There is some similarity to this virus, but the match is not definitive. "
                                    blast_summary += f"The moderate similarity ({top_identity:.1f}%) suggests a possible relationship.\n\n"
                            else:
                                blast_summary += f"I found a **{match_quality} match** ({top_identity:.1f}%) to sequences from **{organism_name or 'the organism'}**.\n\n"
                            
                            # Add secondary matches if they exist and are different
                            if len(blast_results) > 1:
                                other_organisms = set()
                                for hit in blast_results[1:3]:  # Just check next 2
                                    match = re.search(r'\[(.*?)\]', hit["title"])
                                    if match and match.group(1) not in other_organisms:
                                        other_organisms.add(match.group(1))
                                
                                if other_organisms:
                                    blast_summary += "Other similar sequences were found in:\n"
                                    for org in other_organisms:
                                        blast_summary += f"• {org}\n"
                            
                            # Get additional context from RAG system
                            gene_context = ""
                            try:
                                # First ensure RAG system is initialized
                                if rag_system is None:
                                    # Try to initialize it
                                    initialize_ai_systems()
                                
                                # Check again if initialization succeeded
                                if rag_system is not None:
                                    if organism_name:
                                        # Use RAG to get additional knowledge context
                                        gene_context = rag_system.get_gene_context(organism_name)
                                    else:
                                        # Try using the title
                                        gene_context = rag_system.get_gene_context(first_hit_title[:50])
                                else:
                                    gene_context = "Additional context not available - knowledge base not initialized."
                            except Exception as e:
                                gene_context = f"Additional context unavailable - {str(e)}"
                            
                            # Add RAG-enhanced implications
                            blast_summary += "\n**What this means:**\n"
                            if "dengue" in first_hit_title.lower() or (organism_name and "dengue" in organism_name.lower()):
                                blast_summary += "This sequence appears to be from Dengue virus, which causes dengue fever. "
                                blast_summary += "Dengue is a mosquito-borne viral disease that has rapidly spread in all regions in recent years."
                            elif "coronavirus" in first_hit_title.lower() or "sars" in first_hit_title.lower():
                                blast_summary += "This sequence appears to be from a coronavirus. Further analysis would be needed to determine the exact strain and potential implications."
                            else:
                                blast_summary += "This sequence has significant similarity to known genomic sequences. Further analysis is recommended to understand its biological significance."
                            
                            # Add RAG knowledge enhancement
                            if gene_context and len(gene_context) > 20 and not gene_context.startswith("Additional context"):
                                blast_summary += "\n\n**Additional Knowledge Context:**\n"
                                blast_summary += gene_context
                            else:
                                # Provide hardcoded context for common viruses if RAG system didn't return anything
                                if "dengue" in first_hit_title.lower() or (organism_name and "dengue" in organism_name.lower()):
                                    blast_summary += "\n\n**Additional Knowledge Context:**\n"
                                    blast_summary += "Dengue virus is a single-stranded RNA virus that causes dengue fever. " 
                                    blast_summary += "It's primarily transmitted by Aedes mosquitoes, especially Aedes aegypti. " 
                                    blast_summary += "There are four serotypes (DENV-1, DENV-2, DENV-3, and DENV-4). "
                                    blast_summary += "Infection with one serotype provides lifelong immunity against that specific serotype, but only partial protection against others. "
                                    blast_summary += "Secondary infection with a different serotype increases the risk of developing severe dengue. "
                                    blast_summary += "The virus genome encodes for three structural proteins (capsid, membrane, envelope) and seven non-structural proteins."
                                # Only mention context issues if not related to initialization and we didn't provide hardcoded context
                                elif gene_context and not gene_context.startswith("Additional context not available"):
                                    blast_summary += "\n\n**Additional Knowledge Context:** Not available for this sequence."
                        elif blast_searcher.search_method_used == "fallback_analysis":
                            # This is our fallback analysis with no BLAST hits
                            blast_summary += "⚠️ **Network Connectivity Issues**\n\n"
                            blast_summary += "I couldn't connect to the BLAST databases. Here's a basic analysis of your sequence:\n\n"
                            
                            # Get sequence info from the fallback results
                            if blast_results and "sequence_length" in blast_results[0]:
                                seq_len = blast_results[0]["sequence_length"]
                                gc_content = blast_results[0].get("gc_content")
                                
                                blast_summary += f"• Sequence length: {seq_len} base pairs\n"
                                if gc_content:
                                    blast_summary += f"• GC content: {gc_content}%\n"
                                
                                # Add suggestions
                                blast_summary += "\nSuggestions:\n"
                                blast_summary += "• Try again when your internet connection is more stable\n"
                                blast_summary += "• You can still perform local analysis of this sequence\n"
                                blast_summary += "• Consider using a smaller section of the sequence if it's very large\n"
                        else:
                            blast_summary += "❌ I couldn't find any significant matches for this sequence in the database. This might indicate a novel sequence or potential sequencing errors."
                        
                        # Add results to chat
                        st.session_state.messages.append({"role": "assistant", "content": blast_summary})
                    except Exception as e:
                        error_message = f"❌ Error during BLAST search: {str(e)}"
                        st.session_state.messages.append({"role": "assistant", "content": error_message})
                
        with col3:
            if st.button("📚 Find Literature", key="chat_literature", help="Search for relevant research papers"):
                with st.spinner("Searching scientific literature with knowledge enhancement..."):
                    try:
                        # Get the sequence
                        seq = st.session_state.current_sequence
                        
                        # Try to determine what to search for
                        search_term = None
                        
                        # If we have BLAST results with high identity, use those
                        if "blast_results" in st.session_state and st.session_state.get("blast_results"):
                            for hit in st.session_state.blast_results:
                                if hit.get("percent_identity", 0) > 90:
                                    # Extract search term from hit title
                                    title = hit.get("title", "")
                                    import re
                                    match = re.search(r'\[(.*?)\]', title)
                                    if match:
                                        search_term = match.group(1)
                                    else:
                                        words = title.split()
                                        if len(words) > 2:
                                            search_term = " ".join(words[1:3])
                                    break
                        
                        # Try using RAG to determine best search term if no BLAST results
                        if not search_term and "blast_results" in st.session_state and len(st.session_state.blast_results) > 0:
                            # Use the RAG system to help determine the most relevant search term
                            try:
                                hit_title = st.session_state.blast_results[0].get("title", "")
                                # Try to query for relevant gene context
                                context = rag_system.query_knowledge(hit_title, "genes_info")
                                if context and 'results' in context and context['results']:
                                    # Extract potential search terms from context
                                    for result in context['results']:
                                        if "gene" in result.lower():
                                            gene_match = re.search(r'gene\s+([A-Za-z0-9]+)', result, re.IGNORECASE)
                                            if gene_match:
                                                search_term = gene_match.group(1)
                                                break
                            except Exception as e:
                                # Fall back to standard search terms
                                pass
                        
                        if not search_term:
                            # Try DNA analysis
                            if len(seq) > 500:
                                search_term = "genomic sequence analysis"
                            else:
                                search_term = "DNA sequence analysis"
                        
                        # Search PubMed using our existing function
                        pmids = search_pubmed(search_term, max_results=5)
                        
                        # Format results in a conversational style with RAG enhancement
                        lit_summary = f"📚 **Literature Search Results for '{search_term}'**\n\n"
                        
                        # Get additional knowledge context from RAG
                        knowledge_context = ""
                        try:
                            # Try to get relevant literature context from our RAG system
                            context_results = rag_system.query_knowledge(search_term, "literature", n_results=2)
                            if context_results and 'results' in context_results and context_results['results']:
                                knowledge_context = "\n".join(context_results['results'])
                        except Exception as e:
                            # Continue without RAG enhancement if there's an error
                            pass
                        
                        if pmids:
                            # Fetch paper details
                            paper_details = fetch_pubmed_details(pmids)
                            
                            # Add RAG-enhanced introduction if available
                            if knowledge_context:
                                lit_summary += f"**Research Context:**\n{knowledge_context}\n\n"
                            
                            lit_summary += "I found these relevant scientific publications:\n\n"
                            
                            for i, paper in enumerate(paper_details, 1):
                                title = paper.get("title", "Untitled")
                                authors = paper.get("authors", "Unknown")
                                if isinstance(authors, str) and len(authors) > 50:
                                    authors = authors[:50] + "..."
                                journal = paper.get("journal", "")
                                year = paper.get("year", "")
                                
                                lit_summary += f"**{i}. {title}**\n"
                                lit_summary += f"   Authors: {authors}\n"
                                lit_summary += f"   Published in: {journal} ({year})\n"
                                
                                # Try to get paper summary from RAG if possible
                                try:
                                    abstract = paper.get("abstract", "")
                                    if abstract and len(abstract) > 100:
                                        # Include a brief insight about this paper from RAG if possible
                                        paper_context = rag_system.query_knowledge(abstract[:200], "literature")
                                        if paper_context and 'results' in paper_context and paper_context['results']:
                                            lit_summary += f"   **Key insight:** {paper_context['results'][0]}\n"
                                except Exception:
                                    pass
                                
                                lit_summary += "\n"
                                
                            lit_summary += "These papers might provide scientific context for your sequence. Would you like me to summarize any of these papers in more detail?"
                        else:
                            if knowledge_context:
                                lit_summary += f"While I couldn't find specific scientific publications directly related to your sequence, here's some relevant information:\n\n{knowledge_context}\n\n"
                            else:
                                lit_summary += "I couldn't find any scientific publications directly related to your sequence. Try modifying your search terms or ask me to search for a specific topic."
                        
                        # Add to chat
                        st.session_state.messages.append({"role": "assistant", "content": lit_summary})
                    except Exception as e:
                        error_message = f"❌ Error searching literature: {str(e)}"
                        st.session_state.messages.append({"role": "assistant", "content": error_message})
    
    if user_input:
        # Add user message
        st.session_state.messages.append({"role": "user", "content": user_input})
        
        # Get current sequence if available
        sequence = st.session_state.current_sequence if "current_sequence" in st.session_state else None
        
        # Call Node.js chat API
        with st.spinner("🧠 AI is thinking..."):
            response = call_chat_api(user_input, selected_model, sequence)
        
        if response["status"] == "success":
            # Add AI response to chat with model info if available
            model_used = response.get("model_used", selected_model)
            response_content = response["response"]
            
            # Add a small indicator of which model was used (only for non-default models)
            if model_used != selected_model:
                response_content += f"\n\n<small>*Response generated using {model_used} due to rate limits on {selected_model}*</small>"
            
            st.session_state.messages.append({
                "role": "assistant",
                "content": response_content
            })
        else:
            error_msg = response.get('message', 'Unknown error')
            st.error(f"Chat error: {error_msg}")
            
            # Special handling for rate limit errors
            if response.get('rate_limited', False):
                st.warning("⚠️ Rate limit exceeded. Please try one of these options:")
                
                # Show model switching options
                available_models = ["llama-3.1-8b-instant", "kimi-k2", "claude-3-haiku-20240307"]
                col1, col2, col3 = st.columns(3)
                
                with col1:
                    if st.button(f"Try {available_models[0]}", key="try_model1"):
                        st.session_state.selected_model = available_models[0]
                        st.experimental_rerun()
                        
                with col2:
                    if st.button(f"Try {available_models[1]}", key="try_model2"):
                        st.session_state.selected_model = available_models[1]
                        st.experimental_rerun()
                        
                with col3:
                    if st.button("Wait 2 minutes", key="wait_option"):
                        # Add a helpful message
                        st.info("Waiting for rate limits to reset... You can ask a different question while waiting.")
                        
                # Show rate limit guide if it exists
                rate_limit_guide = Path("RATE_LIMIT_GUIDE.md")
                if rate_limit_guide.exists():
                    with st.expander("📋 How to fix rate limit issues"):
                        with open(rate_limit_guide, "r") as f:
                            st.markdown(f.read())
    
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
                    
                    # Run risk assessment automatically if BLAST data is available
                    if result.get("status") == "success" and "blast" in st.session_state.analysis_results:
                        with st.spinner("Running risk assessment..."):
                            try:
                                blast_data = st.session_state.analysis_results["blast"].get("data", {})
                                mutation_data = result.get("data", {})
                                
                                # Create risk assessment object
                                risk_assessor = InfectionRiskAssessment()
                                
                                # Calculate risk
                                risk_results = risk_assessor.calculate_infection_risk(
                                    blast_data, 
                                    mutation_data,
                                    None  # No literature data
                                )
                                
                                # Add to results
                                st.session_state.analysis_results["risk_assessment"] = {
                                    "status": "success",
                                    "data": risk_results
                                }
                            except Exception as e:
                                st.warning(f"Risk assessment failed: {str(e)}")
                                st.session_state.analysis_results["risk_assessment"] = {
                                    "status": "error", 
                                    "message": str(e)
                                }
        
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
            # Risk Assessment button
            st.markdown("**Risk Assessment**")
            
            if st.button("🛡️ Risk Assessment", use_container_width=True):
                with st.spinner("Running comprehensive risk assessment..."):
                    # Check if we have BLAST and mutation data
                    has_blast = "blast" in st.session_state.analysis_results and \
                                st.session_state.analysis_results["blast"].get("status") == "success"
                    has_mutation = "mutation" in st.session_state.analysis_results and \
                                  st.session_state.analysis_results["mutation"].get("status") == "success"
                    
                    # If missing data, run the necessary analysis
                    if not has_blast:
                        with st.spinner("Running BLAST analysis first..."):
                            blast_result = call_api("blast", {
                                "sequence": st.session_state.current_sequence,
                                "num_hits": 5
                            })
                            st.session_state.analysis_results["blast"] = blast_result
                            has_blast = blast_result.get("status") == "success"
                    
                    if not has_mutation:
                        with st.spinner("Running mutation analysis first..."):
                            mutation_result = call_api("mutation", {
                                "sequence": st.session_state.current_sequence,
                                "hit_index": 0
                            })
                            st.session_state.analysis_results["mutation"] = mutation_result
                            has_mutation = mutation_result.get("status") == "success"
                    
                    # Run risk assessment if we have the needed data
                    if has_blast and has_mutation:
                        try:
                            blast_data = st.session_state.analysis_results["blast"].get("data", {})
                            mutation_data = st.session_state.analysis_results["mutation"].get("data", {})
                            
                            # Get literature data if available
                            literature_data = None
                            if "literature" in st.session_state.analysis_results and \
                               st.session_state.analysis_results["literature"].get("status") == "success":
                                literature_data = st.session_state.analysis_results["literature"].get("data", {})
                            
                            # Create risk assessment object
                            risk_assessor = InfectionRiskAssessment()
                            
                            # Calculate risk
                            risk_results = risk_assessor.calculate_infection_risk(
                                blast_data, 
                                mutation_data,
                                literature_data
                            )
                            
                            # Add to results
                            st.session_state.analysis_results["risk_assessment"] = {
                                "status": "success",
                                "data": risk_results
                            }
                            
                            st.success("Risk assessment completed! View results in the Analysis Results section.")
                        except Exception as e:
                            st.error(f"Risk assessment failed: {str(e)}")
                            st.session_state.analysis_results["risk_assessment"] = {
                                "status": "error", 
                                "message": str(e)
                            }
                    else:
                        st.error("Cannot run risk assessment - BLAST and mutation analysis required.")
            
            # Add a divider
            st.markdown("---")
            
            # Literature Search (keep this too)
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
        
        # Add risk assessment to the beginning if it exists
        if "risk_assessment" in st.session_state.analysis_results:
            # Move risk assessment to the front
            result_types = ["risk_assessment"] + [rt for rt in result_types if rt != "risk_assessment"]
        
        # Custom tab labels with icons
        tab_labels = {
            "blast": "🧬 BLAST",
            "mutation": "🔬 Mutations",
            "variants": "📊 Variants",
            "literature": "📚 Literature",
            "risk_assessment": "🛡️ Risk Assessment"
        }
        
        # Set default active tab if not already set
        if "active_tab" not in st.session_state:
            st.session_state["active_tab"] = result_types[0] if result_types else None
        
        # Create tab index mapping for programmatic access
        tab_index_map = {rt: i for i, rt in enumerate(result_types)}
        
        # Create the tabs with selected tab based on session state
        tab_index = tab_index_map.get(st.session_state["active_tab"], 0) if st.session_state["active_tab"] in tab_index_map else 0
        tabs = st.tabs([tab_labels.get(rt, f"🔍 {rt.title()}") for rt in result_types])
        
        for i, (result_type, tab) in enumerate(zip(result_types, tabs)):
            with tab:
                # Update active tab when this tab is selected
                if i == tab_index:
                    st.session_state["active_tab"] = result_type
                
                result_data = st.session_state.analysis_results[result_type]
                
                if result_data.get("status") == "success":
                    data = result_data["data"]
                    
                    if result_type == "blast":
                        display_blast_results(data)
                    elif result_type == "mutation":
                        # Add risk assessment data to mutation results if available
                        if "risk_assessment" in st.session_state.analysis_results and \
                           st.session_state.analysis_results["risk_assessment"].get("status") == "success":
                            data["risk_assessment_data"] = st.session_state.analysis_results["risk_assessment"]["data"]
                        display_mutation_results(data)
                    elif result_type == "variants":
                        display_variant_results(data)
                    elif result_type == "literature":
                        display_literature_results(data)
                    elif result_type == "risk_assessment":
                        display_risk_assessment_results(data)
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
    """Display mutation analysis results with rich visualizations"""
    if "mutations" in data:
        mutations = data.get("mutations", [])
        num_mutations = len(mutations)
        
        # Basic mutation statistics
        st.subheader("🧬 Mutation Statistics")
        col1, col2 = st.columns(2)
        
        with col1:
            st.metric("Total Mutations", num_mutations)
        with col2:
            if "alignment" in data:
                align_data = data["alignment"]
                identity_val = align_data.get('percent_identity', 0)
                st.metric("Sequence Identity", f"{identity_val:.1f}%")
        
        # Main mutation analysis dashboard
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

def display_risk_assessment_results(data: Dict):
    """Display comprehensive risk assessment results without showing BLAST or mutation analysis"""
    if data:
        # Use the display_risk_assessment function from our module, passing only the risk results
        display_risk_assessment(None, None, None, data)
        
        # Add source information
        st.info("ℹ️ Risk assessment is calculated using sequence similarity, taxonomic information, mutation analysis, and literature data when available.")
    else:
        st.warning("No risk assessment data available")

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
    
    # Initialize AI systems on first load
    if not st.session_state.get("rag_initialized", False):
        initialize_ai_systems()
        st.session_state.rag_initialized = True
    
    st.title("🧬 Genomics Research Assistant")
    st.markdown("*AI-powered bioinformatics analysis and research support with RAG and Multi-Agent systems*")
    
    # Sidebar
    with st.sidebar:
        st.header("Navigation")
        
        # Debug tools in expander
        with st.expander("🛠️ Debug Tools"):
            if st.button("Check Environment Variables"):
                st.session_state.show_debug = True
            if st.button("Advanced Troubleshooting"):
                st.session_state.show_debug = True
                st.session_state.show_enhanced_debug = True
            
            # Show help link
            debug_guide_path = Path("DEBUG_GUIDE.md")
            if debug_guide_path.exists():
                if st.button("📋 View Debug Guide"):
                    with open(debug_guide_path, "r") as f:
                        st.markdown(f.read())
    
    # Display debug information if requested
    display_debug_info()
    if st.session_state.get("show_enhanced_debug", False):
        display_enhanced_debug_tools()
        
    # Show RAG system status
    if rag_system is None:
        st.sidebar.error("⚠️ RAG system unavailable. Some features like gene information retrieval will not work correctly.")
        
        # Show help button for fixing RAG
        with st.sidebar.expander("🔧 Fix RAG System"):
            st.markdown("""
            ### Steps to fix the RAG system:
            
            1. **Close any running Streamlit applications**
            2. Open a command prompt and run:
               ```
               cd d:\\Fastalysis\\app
               python create_fresh_rag_db.py
               ```
            3. Wait for the script to complete successfully
            4. Update the `db_path` in `enhanced_streamlit_app.py` to the new path
            5. Restart the application
            """)
            if st.button("Create Fresh RAG DB (may take a minute)"):
                try:
                    import subprocess
                    import sys
                    
                    with st.spinner("Creating fresh RAG database... (this may take a minute)"):
                        # Run the script to create a fresh RAG database
                        process = subprocess.Popen(
                            [sys.executable, "create_fresh_rag_db.py"], 
                            stdout=subprocess.PIPE,
                            stderr=subprocess.PIPE,
                            text=True
                        )
                        stdout, stderr = process.communicate(timeout=180)  # 3 minutes timeout
                        
                        if process.returncode == 0:
                            st.success("Fresh RAG database created successfully! Please restart the application.")
                            st.code(stdout)
                        else:
                            st.error(f"Failed to create fresh RAG database: {stderr}")
                except Exception as e:
                    st.error(f"Failed to create fresh RAG database: {str(e)}")
    
    elif not hasattr(rag_system, 'genes_collection'):
        # Minimal RAG system or incompatible RAG implementation
        st.sidebar.warning("⚠️ RAG system is in minimal mode. Gene and mutation information will be limited.")
        
    elif hasattr(rag_system, 'genes_collection') and rag_system.genes_collection.count() == 0:
        st.sidebar.warning("⚠️ Knowledge base is empty. Use Debug Tools to reload if needed.")
                
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