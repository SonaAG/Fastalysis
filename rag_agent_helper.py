"""
Helper functions to add RAG and Agent systems to the enhanced Streamlit app
"""

import asyncio
import concurrent.futures
from concurrent.futures import ThreadPoolExecutor

def use_rag_system_with_chat(message, sequence, rag_system):
    """Process a user message using RAG to enhance the response"""
    try:
        # Get context from RAG system
        context = ""
        
        # Try to get relevant context based on the message
        context_results = rag_system.query_knowledge(message, "genes_info", n_results=2)
        if context_results and 'results' in context_results and context_results['results']:
            context += "\n".join(context_results['results'])
        
        # If we have a sequence, try to get context for it too
        if sequence and len(sequence) > 20:
            # Use first part of sequence for context lookup
            seq_context = rag_system.get_mutation_context(sequence[:50])
            if seq_context:
                context += "\n" + seq_context
        
        # Return the enhanced context
        return {
            "message": message,
            "context": context
        }
    except Exception as e:
        print(f"RAG processing error: {str(e)}")
        return {"message": message, "context": ""}

def use_agent_system_with_chat(message, sequence, agent_system):
    """Process a user message using the agent system"""
    if not agent_system:
        return {"status": "error", "response": "Agent system not initialized"}
    
    try:
        # First check if the agent system has the process_query method
        if not hasattr(agent_system, 'process_query'):
            print("Agent system missing 'process_query' method")
            return {"status": "error", "response": "Agent system missing required functionality"}
        
        # Check if the method is callable
        if not callable(getattr(agent_system, 'process_query')):
            print("Agent system's 'process_query' is not callable")
            return {"status": "error", "response": "Agent system method not callable"}
        
        # Run the agent system asynchronously
        def run_async(coro):
            loop = asyncio.new_event_loop()
            asyncio.set_event_loop(loop)
            try:
                return loop.run_until_complete(coro)
            finally:
                loop.close()
        
        # Execute the agent query in a separate thread
        with ThreadPoolExecutor() as executor:
            try:
                # Get coroutine from process_query
                coro = agent_system.process_query(message, sequence)
                
                # Submit to executor
                future = executor.submit(run_async, coro)
                
                # Get result with timeout
                agent_response = future.result(timeout=30)  # 30 second timeout
                
                return {"status": "success", "response": agent_response}
            except concurrent.futures.TimeoutError:
                print("Agent system query timed out after 30 seconds")
                return {"status": "error", "response": "The AI assistant team took too long to respond. Please try again with a simpler query."}
            except AttributeError as ae:
                print(f"Agent system attribute error: {str(ae)}")
                return {"status": "error", "response": "The AI assistant system is not properly configured."}
            except TypeError as te:
                print(f"Agent system type error: {str(te)}")
                return {"status": "error", "response": "The AI assistant system received an invalid input type."}
    
    except Exception as e:
        print(f"Agent system error: {str(e)}")
        return {"status": "error", "response": str(e)}