"""
API Quota Monitor and Troubleshooting Tool
Helps check API status and diagnose rate limiting issues
"""

import os
import sys
import time
import json
import argparse
import requests
from datetime import datetime, timedelta
import logging
from pathlib import Path

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
    handlers=[
        logging.StreamHandler(),
        logging.FileHandler("api_quota_monitor.log")
    ]
)
logger = logging.getLogger("api_quota")

def load_env_from_dotenv():
    """Load environment variables from .env.local file"""
    # First check current directory, then parent directories
    dotenv_paths = [".env.local", "../.env.local", "../../.env.local"]
    
    # Add the app directory relative path
    script_dir = Path(__file__).parent
    root_dir = script_dir.parent
    dotenv_paths.extend([
        str(script_dir / ".env.local"),
        str(root_dir / ".env.local")
    ])
    
    # Try to load from each potential path
    for env_path in dotenv_paths:
        env_file = Path(env_path)
        if env_file.exists():
            print(f"Found environment file at {env_path}")
            loaded_vars = 0
            
            try:
                with open(env_file, "r") as f:
                    for line in f:
                        line = line.strip()
                        if line and not line.startswith("#"):
                            try:
                                key, value = line.split("=", 1)
                                key = key.strip()
                                value = value.strip().strip("'\"")
                                os.environ[key] = value
                                loaded_vars += 1
                                
                                # Don't print actual API keys for security
                                if "api_key" in key.lower() or "secret" in key.lower():
                                    logger.info(f"Loaded {key} (value hidden)")
                                else:
                                    logger.info(f"Loaded {key}={value}")
                            except ValueError:
                                # Skip lines that don't have key=value format
                                continue
                
                print(f"✅ Loaded {loaded_vars} environment variables from {env_path}")
                return True
            except Exception as e:
                print(f"❌ Error loading {env_path}: {str(e)}")
                logger.error(f"Failed to load environment file {env_path}: {str(e)}")
    
    # Check if .env.local.template exists but .env.local doesn't
    template_paths = [".env.local.template", "../.env.local.template", 
                     str(script_dir / ".env.local.template"),
                     str(root_dir / ".env.local.template")]
    
    for template_path in template_paths:
        if Path(template_path).exists():
            print("❌ No .env.local file found, but .env.local.template exists")
            print("ℹ️  Please copy .env.local.template to .env.local and add your API keys")
            return False
    
    print("❌ No .env.local file found")
    print("ℹ️  Please create a .env.local file with your API keys (GROQ_API_KEY=your_key)")
    return False

def check_groq_api_status():
    """Check Groq API status and quota information"""
    
    # First try to load environment variables from .env.local
    load_env_from_dotenv()
    
    api_key = os.environ.get("GROQ_API_KEY")
    if not api_key:
        logger.error("GROQ_API_KEY not found in environment variables")
        print("❌ GROQ_API_KEY not found. Please set this environment variable.")
        return False
    
    # Check if the key looks valid
    if not api_key.startswith("gsk_"):
        logger.warning(f"API key doesn't match expected format (should start with 'gsk_')")
        print("⚠️ API key doesn't match expected format (should start with 'gsk_')")
    
    # Test the models endpoint which is typically rate-limited less aggressively
    try:
        print("\nTesting connection to Groq API...")
        response = requests.get(
            "https://api.groq.com/openai/v1/models",  # Updated to correct Groq API endpoint path
            headers={"Authorization": f"Bearer {api_key}"}
        )
        
        if response.status_code == 200:
            models = response.json().get("data", [])
            print(f"✅ Successfully connected to Groq API")
            print(f"Available models: {len(models)}")
            
            # Show a sample of available models
            if models:
                print("\nSample of available models:")
                for model in models[:3]:  # Show first 3 models
                    print(f"- {model.get('id')}")
                if len(models) > 3:
                    print(f"- ... and {len(models) - 3} more")
            
            return True
        elif response.status_code == 429:
            print("❌ Rate limit exceeded! You are currently being rate limited.")
            print("Recommendation: Wait a few minutes before trying again.")
            
            # Check headers for rate limit info
            if 'retry-after' in response.headers:
                retry_after = int(response.headers['retry-after'])
                print(f"API suggests waiting {retry_after} seconds before retrying")
            
            return False
        elif response.status_code == 401:
            print("❌ Authentication failed. Your API key may be invalid.")
            return False
        else:
            print(f"❌ API request failed with status {response.status_code}")
            print(f"Response: {response.text}")
            return False
            
    except Exception as e:
        logger.error(f"Error checking Groq API: {str(e)}")
        print(f"❌ Error checking API: {str(e)}")
        return False

def test_api_rate_limits(model="llama-3.1-8b-instant", requests_count=3):
    """Test API rate limits by making multiple simple requests"""
    
    # Make sure environment variables are loaded
    load_env_from_dotenv()
    
    api_key = os.environ.get("GROQ_API_KEY")
    if not api_key:
        print("❌ GROQ_API_KEY not found. Please set this environment variable.")
        return
    
    # Allow user to choose test mode
    print("\n=== API Rate Limit Testing ===")
    print("1. Standard test (spaced requests)")
    print("2. Stress test (rapid requests)")
    test_mode = 1  # Default to standard test
    
    print(f"\nRunning test mode {test_mode} with {requests_count} requests using model '{model}'")
    delay = 0.1 if test_mode == 2 else 1.0  # Short delay for stress test
    
    successful = 0
    rate_limited = 0
    
    for i in range(requests_count):
        try:
            print(f"Request {i+1}/{requests_count}... ", end="", flush=True)
            
            # Make a minimal request
            response = requests.post(
                "https://api.groq.com/openai/v1/chat/completions",  # Updated to correct Groq API endpoint path
                headers={
                    "Authorization": f"Bearer {api_key}",
                    "Content-Type": "application/json"
                },
                json={
                    "model": model,
                    "messages": [{"role": "user", "content": "Hello"}],
                    "max_tokens": 5  # Keep response tiny
                }
            )
            
            if response.status_code == 200:
                print("✅ Success")
                successful += 1
            elif response.status_code == 429:
                print("❌ Rate limited")
                rate_limited += 1
                
                # Check headers for rate limit info
                if 'retry-after' in response.headers:
                    retry_after = int(response.headers['retry-after'])
                    print(f"   API suggests waiting {retry_after} seconds before retrying")
                    time.sleep(min(retry_after, 5))  # Wait up to 5 seconds
                else:
                    # Standard backoff
                    time.sleep(2)
            else:
                print(f"❓ Unexpected response: {response.status_code}")
                print(f"   {response.text}")
                
        except Exception as e:
            print(f"❌ Error: {str(e)}")
        
        # Small delay between requests to avoid hammering the API
        if i < requests_count - 1:
            time.sleep(delay)
    
    print(f"\nResults: {successful} successful, {rate_limited} rate-limited")
    
    if rate_limited > 0:
        print("\nRecommendation: Your application is hitting rate limits.")
        print("1. Implement exponential backoff retry logic")
        print("2. Add delay between requests (at least 1 second)")
        print("3. Consider upgrading your API tier if available")
        print("\nDetails:")
        print(f"- Model tested: {model}")
        print(f"- Delay between requests: {delay} seconds")
        print(f"- Failure rate: {rate_limited/requests_count*100:.1f}%")
        print("\nSee RATE_LIMIT_GUIDE.md for more troubleshooting tips")
    else:
        print("\nYour API quota appears to be working correctly!")
        print("However, when used in the full application, you may still hit limits due to:")
        print("1. Multiple concurrent requests from different components")
        print("2. Larger context sizes in real usage")
        print("3. Longer sessions with more frequent requests")

def get_api_service_status():
    """Check the service status page if available"""
    try:
        # Note: Groq doesn't have a public status page API endpoint
        # This is a placeholder - consider replacing with the actual status page check
        print("\nChecking Groq service status...")
        print("Note: Automated status checks not available.")
        print("Please check https://status.groq.com/ manually for service status.")
    except Exception as e:
        print(f"Error checking service status: {str(e)}")

def main():
    parser = argparse.ArgumentParser(description="API Quota Monitor and Troubleshooting")
    parser.add_argument("--check", action="store_true", help="Check API status")
    parser.add_argument("--test", action="store_true", help="Test API rate limits")
    parser.add_argument("--count", type=int, default=3, help="Number of test requests")
    parser.add_argument("--status", action="store_true", help="Check service status page")
    parser.add_argument("--env", action="store_true", help="Show environment variables")
    
    args = parser.parse_args()
    
    # Header
    print("\n" + "="*50)
    print("API QUOTA MONITOR AND TROUBLESHOOTING TOOL")
    print("="*50 + "\n")
    
    # Always load environment variables first
    load_env_from_dotenv()
    
    if args.env:
        print("\n=== Environment Variables ===")
        # Only show if the API key exists, not its value for security
        groq_key = os.environ.get("GROQ_API_KEY", "")
        print(f"GROQ_API_KEY: {'✅ Present' if groq_key else '❌ Missing'}")
        if groq_key:
            # Show only first and last 4 chars for verification without revealing the whole key
            masked_key = f"{groq_key[:4]}...{groq_key[-4:]}" if len(groq_key) > 8 else "***"
            print(f"  Key format check: {'✅ Valid' if groq_key.startswith('gsk_') else '⚠️ Unexpected format'}")
            print(f"  Key preview: {masked_key}")
        
        # Check other common API keys
        biogpt_key = os.environ.get("BIOGPT_API_KEY", "")
        print(f"BIOGPT_API_KEY: {'✅ Present' if biogpt_key else '❌ Missing'}")
        
        openai_key = os.environ.get("OPENAI_API_KEY", "")
        print(f"OPENAI_API_KEY: {'✅ Present' if openai_key else '❌ Missing'}")
    
    if args.check or not (args.test or args.status or args.env):
        check_groq_api_status()
    
    if args.test:
        test_api_rate_limits(requests_count=args.count)
    
    if args.status:
        get_api_service_status()
    
    # Footer
    print("\n" + "="*50)
    print(f"Report generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print("="*50 + "\n")

if __name__ == "__main__":
    main()