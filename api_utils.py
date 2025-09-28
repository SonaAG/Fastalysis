"""
Rate limiting and retry utilities for API calls
Helps manage API request quotas and handles rate limit errors gracefully
"""

import time
import random
from typing import Callable, Any, Dict, Optional, TypeVar
import logging

# Configure logging
logging.basicConfig(level=logging.INFO,
                    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s')
logger = logging.getLogger("api_utils")

# Type variable for generic function return type
T = TypeVar('T')

class RateLimitHandler:
    """Handler for API rate limits with exponential backoff"""
    
    def __init__(self, 
                 base_delay: float = 1.0,
                 max_delay: float = 60.0, 
                 max_retries: int = 5,
                 jitter: bool = True):
        """
        Initialize rate limit handler
        
        Args:
            base_delay: Initial delay in seconds
            max_delay: Maximum delay in seconds
            max_retries: Maximum number of retry attempts
            jitter: Whether to add randomness to delay
        """
        self.base_delay = base_delay
        self.max_delay = max_delay
        self.max_retries = max_retries
        self.jitter = jitter
        
        # Track request timestamps by endpoint to manage rate limits
        self.request_history = {}
        
    def with_retry(self, func: Callable[..., T], *args, **kwargs) -> T:
        """
        Execute a function with retry logic for rate limits
        
        Args:
            func: Function to execute
            *args: Positional arguments for the function
            **kwargs: Keyword arguments for the function
            
        Returns:
            Result from the function call
            
        Raises:
            Exception: If all retries fail
        """
        last_exception = None
        
        # Try the function with retries
        for attempt in range(self.max_retries):
            try:
                return func(*args, **kwargs)
            except Exception as e:
                last_exception = e
                
                # Check if this is a rate limit error
                is_rate_limit = self._is_rate_limit_error(e)
                
                if is_rate_limit:
                    # Calculate delay with exponential backoff
                    delay = min(self.base_delay * (2 ** attempt), self.max_delay)
                    
                    # Add jitter if enabled
                    if self.jitter:
                        delay = delay * (0.5 + random.random())
                    
                    logger.warning(f"Rate limit hit. Retrying in {delay:.2f} seconds. Attempt {attempt+1}/{self.max_retries}")
                    time.sleep(delay)
                else:
                    # Not a rate limit error, don't retry
                    raise
        
        # If we get here, all retries failed
        logger.error(f"All {self.max_retries} retry attempts failed")
        raise last_exception
    
    def _is_rate_limit_error(self, exception: Exception) -> bool:
        """
        Determine if an exception is due to rate limiting
        
        Args:
            exception: The exception to check
            
        Returns:
            True if rate limit error, False otherwise
        """
        # Check for common rate limit error patterns
        error_str = str(exception).lower()
        
        # OpenAI/Groq style status code check
        if '429' in error_str or 'too many requests' in error_str:
            return True
        
        # More specific OpenAI/Groq messages
        if 'rate limit' in error_str or 'quota exceeded' in error_str:
            return True
            
        # Check for requests library status codes
        if hasattr(exception, 'status_code') and getattr(exception, 'status_code') == 429:
            return True
            
        return False

# Global rate limit handler instance
global_rate_limiter = RateLimitHandler()

def with_rate_limit(func: Callable[..., T], *args, **kwargs) -> T:
    """
    Convenience wrapper to use global rate limiter
    
    Args:
        func: Function to execute with rate limiting
        *args: Positional arguments for the function
        **kwargs: Keyword arguments for the function
        
    Returns:
        Result from the function call
    """
    return global_rate_limiter.with_retry(func, *args, **kwargs)