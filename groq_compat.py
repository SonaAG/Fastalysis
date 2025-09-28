"""
Utility module for customizing AutoGen to work with Groq API
With rate limiting and retry mechanisms for handling API quotas
"""

import os
import sys
import time
import random
from typing import Dict, Any, List, Optional
import logging

# Configure logging
logging.basicConfig(level=logging.INFO,
                    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s')
logger = logging.getLogger("groq_compat")

# Import rate limiting utilities
try:
    from app.api_utils import with_rate_limit
    RATE_LIMIT_AVAILABLE = True
except ImportError:
    logger.warning("Rate limiting utilities not available")
    RATE_LIMIT_AVAILABLE = False
    
    # Define a passthrough function if api_utils isn't available
    def with_rate_limit(func, *args, **kwargs):
        return func(*args, **kwargs)

def create_groq_config(model: str = "llama-3.1-8b-instant") -> Dict[str, Any]:
    """
    Create a configuration for AutoGen that works with Groq API
    
    Args:
        model: The model name to use from Groq
        
    Returns:
        AutoGen-compatible configuration using Groq API
    """
    # Get the GROQ API key from environment
    groq_api_key = os.environ.get("GROQ_API_KEY", "")
    
    if not groq_api_key:
        logger.warning("GROQ_API_KEY not found in environment variables")
        print("WARNING: GROQ_API_KEY not found in environment variables - AI assistants won't work")
        return {}
    
    # Define valid Groq models
    valid_models = [
        "llama-3.1-8b-instant",
        "llama-3.3-70b-versatile",
        "deepseek-r1-distill-llama-70b",
        "meta-llama/llama-4-scout-17b-16e-instruct"
    ]
    
    # Use default if specified model is not valid
    if model not in valid_models:
        logging.warning(f"Specified model '{model}' not in list of valid Groq models. Using default.")
        model = "llama-3.1-8b-instant"
    
    # Create configuration
    config = {
        "config_list": [
            {
                "model": model,
                "api_key": groq_api_key,
                "base_url": "https://api.groq.com/openai/v1"
            }
        ],
        "temperature": 0.7,
        "timeout": 120,
    }
    
    return config


def verify_groq_api_key() -> bool:
    """
    Verify if GROQ_API_KEY is present in environment variables
    
    Returns:
        bool: True if key is present, False otherwise
    """
    groq_api_key = os.environ.get("GROQ_API_KEY", "")
    
    print(f"Verifying GROQ_API_KEY...")
    
    if not groq_api_key:
        logger.warning("GROQ_API_KEY not found in environment variables")
        print("GROQ_API_KEY not found in environment variables!")
        return False
        
    # Check if key looks valid (just basic format check)
    if groq_api_key.startswith("gsk_") and len(groq_api_key) > 20:
        logger.info("GROQ_API_KEY found and appears valid")
        print(f"GROQ_API_KEY found and appears valid: {groq_api_key[:5]}...{groq_api_key[-4:]}")
        return True
    else:
        logger.warning("GROQ_API_KEY found but format appears invalid")
        print(f"GROQ_API_KEY found but format appears invalid: {groq_api_key[:5]}...") 
        return True  # Still return True even if format doesn't match expected pattern

def verify_groq_key(api_key: str) -> bool:
    """
    Verify if a provided Groq API key appears to be valid
    
    Args:
        api_key: The API key to verify
    
    Returns:
        bool: True if key format appears valid, False otherwise
    """
    if not api_key:
        logger.warning("Empty API key provided to verification function")
        return False
        
    # Check if key matches expected format
    if api_key.startswith("gsk_") and len(api_key) > 20:
        return True
    else:
        logger.warning(f"API key format appears invalid: {api_key[:5]}...")
        return False

class GroqCompatibilityWrapper:
    """
    Wrapper to make Groq API work with AutoGen by intercepting key validation
    Includes rate limiting and retry logic to handle quotas
    """
    
    def __init__(self, api_key: str, base_url: str = "https://api.groq.com/openai/v1"):
        """Initialize wrapper with API key and base URL"""
        self.api_key = api_key
        self.base_url = base_url
        
        # Track request timing for rate limiting
        self.last_request_time = 0
        self.min_request_interval = 1.0  # seconds between requests
        
        # Log status
        if not api_key:
            logger.warning("Empty API key provided to GroqCompatibilityWrapper")
    
    def _respect_rate_limit(self):
        """Ensure we don't exceed rate limits by adding delays between requests"""
        current_time = time.time()
        elapsed = current_time - self.last_request_time
        
        # If we've made a request too recently, wait
        if elapsed < self.min_request_interval:
            sleep_time = self.min_request_interval - elapsed
            logger.info(f"Rate limiting: Waiting {sleep_time:.2f}s between requests")
            time.sleep(sleep_time)
        
        # Update the last request time
        self.last_request_time = time.time()
        
    def create_client(self, **kwargs):
        """Create a client that works with Groq API"""
        try:
            import openai
            
            # Create a standard client first
            standard_client = openai.OpenAI(
                api_key=self.api_key,
                base_url=self.base_url
            )
            
            if RATE_LIMIT_AVAILABLE:
                # Wrap the client's completions.create method with rate limiting
                original_create = standard_client.chat.completions.create
                
                def rate_limited_create(*args, **kwargs):
                    self._respect_rate_limit()
                    return with_rate_limit(original_create, *args, **kwargs)
                
                # Replace the method with our rate-limited version
                standard_client.chat.completions.create = rate_limited_create
                logger.info("Applied rate limiting wrapper to Groq API client")
            
            return standard_client
        except ImportError:
            logging.error("OpenAI package not available")
            return None
            
    def get_config(self, model: str = "llama-3.1-8b-instant") -> Dict[str, Any]:
        """Get AutoGen-compatible configuration"""
        return {
            "config_list": [
                {
                    "model": model,
                    "api_key": self.api_key,
                    "base_url": self.base_url
                }
            ],
            "temperature": 0.7,
            "timeout": 120,
        }