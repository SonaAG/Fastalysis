"""
Debug tools for Fastalysis application
Provides utilities for troubleshooting environment variables and dependencies
"""

import os
import streamlit as st
import sys
from typing import Dict, List, Optional

def display_debug_info():
    """Display debug information for troubleshooting."""
    if st.session_state.get('show_debug', False):
        st.sidebar.markdown("### 🔍 Debug Information")
        
        # Display environment variables (safely)
        st.sidebar.markdown("#### Environment Variables")
        debug_env_vars = {
            "GROQ_API_KEY": f"{'✓ Present' if os.environ.get('GROQ_API_KEY') else '❌ Missing'}",
            "BIOGPT_API_KEY": f"{'✓ Present' if os.environ.get('BIOGPT_API_KEY') else '❌ Missing'}",
            "OPENAI_API_KEY": f"{'✓ Present' if os.environ.get('OPENAI_API_KEY') else '❌ Missing'}"
        }
        for key, value in debug_env_vars.items():
            st.sidebar.text(f"{key}: {value}")
        
        # Check if Groq key is properly formatted
        if os.environ.get('GROQ_API_KEY'):
            try:
                from app.groq_compat import verify_groq_key
                is_valid = verify_groq_key(os.environ.get('GROQ_API_KEY', ''))
                st.sidebar.text(f"GROQ Key Format Valid: {'✓ Yes' if is_valid else '❌ No'}")
            except ImportError:
                st.sidebar.text("Unable to verify GROQ key format (groq_compat not found)")
        
        # Display Python environment details
        st.sidebar.markdown("#### Python Environment")
        st.sidebar.text(f"Python Version: {sys.version.split()[0]}")
        
        # Display critical imports status
        st.sidebar.markdown("#### Dependencies")
        dependencies = {
            "autogen": check_dependency("autogen"),
            "chromadb": check_dependency("chromadb"),
            "langchain": check_dependency("langchain"),
            "biopython": check_dependency("Bio"),
            "streamlit": "✓ Available"
        }
        for dep, status in dependencies.items():
            st.sidebar.text(f"{dep}: {status}")
        
        # Check local files
        st.sidebar.markdown("#### Critical Files")
        files = {
            ".env.local": "✓ Found" if os.path.exists(".env.local") else "❌ Missing",
            "env_loader.py": "✓ Found" if os.path.exists("env_loader.py") else "❌ Missing",
            "groq_compat.py": "✓ Found" if os.path.exists("app/groq_compat.py") else "❌ Missing"
        }
        for file, status in files.items():
            st.sidebar.text(f"{file}: {status}")
        
        # Add button to hide debug info
        if st.sidebar.button("Hide Debug Info"):
            st.session_state.show_debug = False
            st.rerun()

def check_dependency(module_name: str) -> str:
    """Check if a dependency is installed and return status"""
    try:
        __import__(module_name)
        return "✓ Available"
    except ImportError:
        return "❌ Missing"

def check_env_variable_sources() -> Dict[str, List[str]]:
    """Check all possible sources of environment variables"""
    sources = {
        "os.environ": [],
        ".env file": [],
        ".env.local file": [],
    }
    
    # Check os.environ
    for key in ["GROQ_API_KEY", "BIOGPT_API_KEY", "OPENAI_API_KEY"]:
        if key in os.environ:
            sources["os.environ"].append(key)
    
    # Check .env file
    if os.path.exists(".env"):
        with open(".env", "r") as f:
            content = f.read()
            for key in ["GROQ_API_KEY", "BIOGPT_API_KEY", "OPENAI_API_KEY"]:
                if key in content:
                    sources[".env file"].append(key)
    
    # Check .env.local file
    if os.path.exists(".env.local"):
        with open(".env.local", "r") as f:
            content = f.read()
            for key in ["GROQ_API_KEY", "BIOGPT_API_KEY", "OPENAI_API_KEY"]:
                if key in content:
                    sources[".env.local file"].append(key)
    
    return sources

def test_groq_api_connection() -> Dict[str, str]:
    """Test connection to Groq API"""
    import requests
    
    result = {
        "status": "unknown",
        "message": "",
        "model_availability": {}
    }
    
    api_key = os.environ.get("GROQ_API_KEY", "")
    if not api_key:
        result["status"] = "failed"
        result["message"] = "No API key found"
        return result
    
    # Test models endpoint
    try:
        response = requests.get(
            "https://api.groq.com/v1/models", 
            headers={"Authorization": f"Bearer {api_key}"}
        )
        
        if response.status_code == 200:
            result["status"] = "success"
            result["message"] = "Successfully connected to Groq API"
            
            # Check specific models
            models = response.json()
            for model in models.get("data", []):
                model_id = model.get("id")
                if model_id:
                    result["model_availability"][model_id] = "Available"
        else:
            result["status"] = "failed"
            result["message"] = f"API error: {response.status_code} - {response.text}"
    except Exception as e:
        result["status"] = "failed"
        result["message"] = f"Connection error: {str(e)}"
    
    return result

def display_enhanced_debug_tools():
    """Display advanced debugging tools"""
    st.sidebar.markdown("#### Advanced Troubleshooting")
    
    # Test connectivity to Groq API
    if st.sidebar.button("Test Groq API Connection"):
        with st.sidebar:
            with st.spinner("Testing connection..."):
                result = test_groq_api_connection()
                
                if result["status"] == "success":
                    st.success(result["message"])
                    st.write("Available models:")
                    for model, status in result["model_availability"].items():
                        st.text(f"- {model}: {status}")
                else:
                    st.error(result["message"])
    
    # Check variable sources
    if st.sidebar.button("Check Env Variable Sources"):
        sources = check_env_variable_sources()
        st.sidebar.markdown("#### Environment Variable Sources:")
        for source, vars in sources.items():
            if vars:
                st.sidebar.text(f"{source}: {', '.join(vars)}")
            else:
                st.sidebar.text(f"{source}: None found")