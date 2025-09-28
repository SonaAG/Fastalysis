"""
Environment variable loader utility for Fastalysis
"""

import os
import re
from pathlib import Path
from typing import Dict, Optional

def load_env_file(env_file: str = ".env.local") -> Dict[str, str]:
    """
    Load environment variables from a file
    
    Args:
        env_file: Path to the .env file
        
    Returns:
        Dictionary of loaded environment variables
    """
    loaded_vars = {}
    
    # Try to locate the .env file
    env_path = Path(env_file)
    if not env_path.is_absolute():
        # Check current directory
        if Path(env_file).exists():
            env_path = Path(env_file)
        # Check app directory
        elif Path("app") / env_file.exists():
            env_path = Path("app") / env_file
        # Check parent directory
        elif Path("..") / env_file.exists():
            env_path = Path("..") / env_file
        else:
            print(f"Warning: Could not find {env_file}")
            return loaded_vars
    
    try:
        with open(env_path, "r") as f:
            for line in f:
                line = line.strip()
                
                # Skip empty lines or comments
                if not line or line.startswith("#"):
                    continue
                    
                # Match VAR=VALUE or VAR="VALUE" format
                match = re.match(r'^([A-Za-z0-9_]+)=(?:"([^"]*)"|(.*))$', line)
                if match:
                    key = match.group(1)
                    # Group 2 is for quoted values, group 3 for unquoted
                    value = match.group(2) if match.group(2) is not None else match.group(3)
                    
                    # Set the environment variable
                    os.environ[key] = value
                    loaded_vars[key] = value
                    print(f"Loaded environment variable: {key}")
    
    except Exception as e:
        print(f"Error loading environment variables from {env_file}: {str(e)}")
    
    return loaded_vars

def ensure_api_keys() -> bool:
    """
    Ensure that required API keys are loaded
    
    Returns:
        True if all required keys are loaded, False otherwise
    """
    # First try to load from .env.local
    load_env_file(".env.local")
    
    # Check for required API keys
    required_keys = ["GROQ_API_KEY"]
    missing_keys = []
    
    for key in required_keys:
        if not os.environ.get(key):
            missing_keys.append(key)
    
    if missing_keys:
        print(f"Warning: Missing required API keys: {', '.join(missing_keys)}")
        return False
        
    return True

if __name__ == "__main__":
    # Test the environment loader
    loaded = load_env_file()
    print(f"Loaded {len(loaded)} environment variables")
    
    # Check if required keys are present
    if ensure_api_keys():
        print("All required API keys are present")
    else:
        print("Some required API keys are missing")