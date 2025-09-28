"""
Genomics Research Agent System using AutoGen
Multi-agent coordination for intelligent bioinformatics analysis
"""

from typing import Dict, List, Optional, Any
import json
import asyncio
import os
from datetime import datetime

# AutoGen imports
from autogen import AssistantAgent, UserProxyAgent, GroupChat, GroupChatManager
from autogen.agentchat.contrib.retrieve_assistant_agent import RetrieveAssistantAgent
from autogen.agentchat.contrib.retrieve_user_proxy_agent import RetrieveUserProxyAgent

# Your existing modules
from controller import (
    run_blast_controller,
    run_mutation_controller, 
    run_variant_table_controller,
    run_pubmed_controller,
    run_full_pipeline
)

class GenomicsAgentSystem:
    def __init__(self, chat_model_config: Dict):
        # Import our custom Groq compatibility module if needed
        try:
            from groq_compat import create_groq_config
            # If config looks like a placeholder, recreate it with Groq compatibility
            if (not chat_model_config.get("config_list") or 
                chat_model_config["config_list"][0].get("api_key") == "placeholder"):
                self.llm_config = create_groq_config()
            else:
                self.llm_config = chat_model_config
        except ImportError:
            # Fall back to provided config
            self.llm_config = chat_model_config
            
        # Initialize agents
        self.setup_agents()
        
    def setup_agents(self):
        """Initialize all specialized agents"""
        # Check if we have a valid config
        if (not self.llm_config or not self.llm_config.get("config_list") or 
            not self.llm_config["config_list"][0].get("api_key")):
            raise ValueError("Missing valid LLM configuration with API key")
        
        # 1. Chat Coordinator Agent - Main interface
        self.coordinator = AssistantAgent(
            name="GenomicsCoordinator",
            system_message="""You are a Genomics Research Coordinator. Your role is to:
            1. Understand user queries about genomics, mutations, genes, and sequences
            2. Determine which analysis tools are needed
            3. Coordinate with specialist agents
            4. Provide human-friendly explanations of results
            
            Available tools: BLAST search, mutation analysis, variant lookup, literature search.
            Always explain results in researcher-friendly language.""",
            llm_config=self.llm_config
        )
        
        # 2. Analysis Agent - Handles bioinformatics operations  
        self.analyzer = AssistantAgent(
            name="BioinformaticsAnalyzer", 
            system_message="""You are a Bioinformatics Analysis Specialist. You:
            1. Execute BLAST searches, mutation analysis, variant lookups
            2. Interpret sequence data and alignment results
            3. Provide technical details about genomic findings
            4. Recommend follow-up analyses
            
            You have access to comprehensive bioinformatics tools.""",
            llm_config=self.llm_config
        )
        
        # 3. Literature Agent - Handles PubMed and knowledge
        self.literature_agent = AssistantAgent(
            name="LiteratureSpecialist",
            system_message="""You are a Scientific Literature Specialist. You:
            1. Search PubMed for relevant research papers
            2. Summarize key findings from genomics literature
            3. Connect research findings to user queries
            4. Provide context about gene functions, diseases, mutations
            
            You excel at explaining complex research in accessible terms.""",
            llm_config=self.llm_config
        )
        
        # 4. RAG Agent - Knowledge retrieval (will be enhanced)
        self.knowledge_agent = RetrieveAssistantAgent(
            name="GenomicsKnowledge",
            system_message="""You are a Genomics Knowledge Base. You provide:
            1. Background information about genes, proteins, pathways
            2. Disease associations and clinical significance
            3. Evolutionary and functional context
            4. Current research trends and discoveries""",
            llm_config=self.llm_config,
            # Will be configured with vector database
        )
        
        # 5. User Proxy - Executes functions
        self.user_proxy = UserProxyAgent(
            name="user_proxy",
            human_input_mode="NEVER",
            max_consecutive_auto_reply=10,
            code_execution_config={"work_dir": ".", "use_docker": False},
        )
        
    async def process_query(self, user_query: str, sequence_data: Optional[str] = None) -> Dict:
        """Main entry point for processing user queries with robust error handling"""
        
        try:
            # Try to import our rate limiting utility
            try:
                from app.api_utils import with_rate_limit
                use_rate_limiting = True
            except ImportError:
                use_rate_limiting = False
            
            # Create group chat for agent coordination
            groupchat = GroupChat(
                agents=[self.coordinator, self.analyzer, self.literature_agent, self.user_proxy],
                messages=[],
                max_round=10
            )
            manager = GroupChatManager(groupchat=groupchat, llm_config=self.llm_config)
            
            # Enhanced prompt with context
            enhanced_query = f"""
            User Query: {user_query}
            Sequence Data: {sequence_data if sequence_data else 'None provided'}
            
            Please coordinate to provide a comprehensive response including:
            1. Understanding of the user's request
            2. Appropriate bioinformatics analysis if needed  
            3. Literature context and background
            4. Clear, researcher-friendly explanation
            """
            
            # Start the conversation with rate limiting if available
            if use_rate_limiting:
                response = with_rate_limit(
                    self.user_proxy.initiate_chat,
                    manager, 
                    message=enhanced_query
                )
            else:
                response = self.user_proxy.initiate_chat(manager, message=enhanced_query)
            
            return {
                "status": "success",
                "timestamp": datetime.now().isoformat(),
                "query": user_query,
                "response": response,
                "agents_involved": ["coordinator", "analyzer", "literature_agent"]
            }
        
        except Exception as e:
            # Detailed error reporting
            error_msg = str(e)
            error_type = type(e).__name__
            
            # Check if this is a rate limit error
            is_rate_limit = False
            if '429' in error_msg or 'rate limit' in error_msg.lower() or 'too many requests' in error_msg.lower():
                is_rate_limit = True
                
            # Prepare user-friendly error message
            if is_rate_limit:
                friendly_msg = ("The AI service is currently experiencing high demand. "
                              "Please wait a minute and try again.")
            else:
                friendly_msg = "There was an issue processing your request with the AI assistants."
            
            return {
                "status": "error",
                "error_type": error_type,
                "error": error_msg,
                "user_message": friendly_msg,
                "timestamp": datetime.now().isoformat(),
                "query": user_query,
                "retry_suggested": is_rate_limit
            }

    def register_bioinformatics_functions(self):
        """Register your controller functions with agents"""
        
        @self.user_proxy.register_for_execution()
        @self.analyzer.register_for_llm(description="Run BLAST search on sequence")
        def blast_search(sequence: str, num_hits: int = 5) -> Dict:
            """Execute BLAST search"""
            import tempfile
            import os
            temp_file = tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False)
            temp_file.write(f">user_sequence\n{sequence}")
            temp_file.close()
            
            result = run_blast_controller(temp_file.name, num_hits=num_hits)
            os.unlink(temp_file.name)
            return result
            
        @self.user_proxy.register_for_execution()
        @self.analyzer.register_for_llm(description="Analyze mutations between sequences")  
        def mutation_analysis(sequence: str, hit_index: int = 0) -> Dict:
            """Execute mutation analysis"""
            import tempfile
            import os
            temp_file = tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False)
            temp_file.write(f">user_sequence\n{sequence}")
            temp_file.close()
            
            result = run_mutation_controller(temp_file.name, hit_index=hit_index)
            os.unlink(temp_file.name)
            return result
            
        @self.user_proxy.register_for_execution()
        @self.literature_agent.register_for_llm(description="Search PubMed literature")
        def search_literature(gene_or_accession: str, additional_terms: str = None) -> Dict:
            """Search scientific literature"""
            return run_pubmed_controller(gene_or_accession, additional_terms=additional_terms)
        
        @self.user_proxy.register_for_execution() 
        @self.analyzer.register_for_llm(description="Look up genetic variants")
        def lookup_variants(gene_or_accession: str, by_accession: bool = True) -> Dict:
            """Look up genetic variants"""
            return run_variant_table_controller(gene_or_accession, by_accession=by_accession)

# Example usage configuration
def create_genomics_agent_system():
    """Factory function to create the agent system"""
    
    # Use our custom Groq compatibility module if available
    try:
        from groq_compat import create_groq_config
        llm_config = create_groq_config("llama-3.1-8b-instant")
        
        # Check if config is empty (no API key found)
        if not llm_config:
            # Fall back to default config
            llm_config = {
                "config_list": [
                    {
                        "model": "llama-3.1-8b-instant", 
                        "api_key": os.environ.get("GROQ_API_KEY", "missing-key"),
                        "base_url": "https://api.groq.com/openai/v1",
                    }
                ],
                "temperature": 0.7,
                "timeout": 120,
            }
    except ImportError:
        # Fall back to default config
        llm_config = {
            "config_list": [
                {
                    "model": "llama-3.1-8b-instant",
                    "api_key": os.environ.get("GROQ_API_KEY", "missing-key"),
                    "base_url": "https://api.groq.com/openai/v1",
                }
            ],
            "temperature": 0.7,
            "timeout": 120,
        }
    
    agent_system = GenomicsAgentSystem(llm_config)
    agent_system.register_bioinformatics_functions()
    
    return agent_system

# Test function
async def test_agent_system():
    """Test the agent system with a sample query"""
    agent_system = create_genomics_agent_system()
    
    result = await agent_system.process_query(
        user_query="I want to analyze mutations in the BRCA1 gene. Can you help me understand what mutations are commonly found?",
        sequence_data="ATCGATCGATCG"  # sample sequence
    )
    
    print(json.dumps(result, indent=2))

if __name__ == "__main__":
    asyncio.run(test_agent_system())