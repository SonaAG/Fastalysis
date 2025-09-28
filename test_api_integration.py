"""
Test script to verify API integration in Fastalysis
Run this before starting the main application to check for issues
"""

import os
import sys
import time
import json
import requests
from pathlib import Path
import argparse
import logging

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger("api_test")

def load_env_variables():
    """Load environment variables from .env.local"""
    env_file = Path(".env.local")
    if not env_file.exists():
        print("❌ .env.local file not found")
        print("Please create this file with your API keys")
        return False
    
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
                    except ValueError:
                        continue
        print("✅ Loaded environment variables from .env.local")
        return True
    except Exception as e:
        print(f"❌ Error loading environment variables: {str(e)}")
        return False

def test_groq_api():
    """Test connection to Groq API"""
    api_key = os.environ.get("GROQ_API_KEY")
    if not api_key:
        print("❌ GROQ_API_KEY not found")
        return False
    
    print("\nTesting Groq API connection...")
    try:
        # Test models endpoint
        response = requests.get(
            "https://api.groq.com/openai/v1/models",
            headers={"Authorization": f"Bearer {api_key}"}
        )
        
        if response.status_code == 200:
            print("✅ Successfully connected to Groq API")
            models = response.json().get("data", [])
            print(f"  Found {len(models)} available models")
            return True
        else:
            print(f"❌ API error: {response.status_code} - {response.text}")
            return False
    except Exception as e:
        print(f"❌ Connection error: {str(e)}")
        return False

def test_chat_completion():
    """Test chat completion API"""
    api_key = os.environ.get("GROQ_API_KEY")
    if not api_key:
        return False
    
    print("\nTesting chat completion API...")
    try:
        response = requests.post(
            "https://api.groq.com/openai/v1/chat/completions",
            headers={
                "Authorization": f"Bearer {api_key}",
                "Content-Type": "application/json"
            },
            json={
                "model": "llama-3.1-8b-instant",
                "messages": [
                    {"role": "system", "content": "You are a helpful assistant."},
                    {"role": "user", "content": "What is genomics?"}
                ],
                "max_tokens": 100
            }
        )
        
        if response.status_code == 200:
            print("✅ Successfully received chat completion")
            response_text = response.json()["choices"][0]["message"]["content"]
            print(f"  Response preview: {response_text[:50]}...")
            return True
        elif response.status_code == 429:
            print("❌ Rate limit exceeded (429)")
            print("  Wait a few minutes and try again")
            return False
        else:
            print(f"❌ API error: {response.status_code}")
            print(f"  Response: {response.text}")
            return False
    except Exception as e:
        print(f"❌ Error: {str(e)}")
        return False

def check_required_files():
    """Check if all required files exist"""
    print("\nChecking required files...")
    required_files = [
        "enhanced_streamlit_app.py",
        "groq_compat.py",
        "api_utils.py"
    ]
    
    all_found = True
    for file in required_files:
        if Path(file).exists() or Path(f"app/{file}").exists():
            print(f"✅ Found {file}")
        else:
            print(f"❌ Missing {file}")
            all_found = False
    
    return all_found

def simulate_agent_system():
    """Simulate the agent system calls to test API integration"""
    print("\nSimulating agent system API calls...")
    
    api_key = os.environ.get("GROQ_API_KEY")
    if not api_key:
        return False
    
    # Test parameters that match what the agent system would use
    try:
        response = requests.post(
            "https://api.groq.com/openai/v1/chat/completions",
            headers={
                "Authorization": f"Bearer {api_key}",
                "Content-Type": "application/json"
            },
            json={
                "model": "llama-3.1-8b-instant",
                "messages": [
                    {"role": "system", "content": "You are a genomics research coordinator."},
                    {"role": "user", "content": "What are mutations?"}
                ],
                "temperature": 0.7,
                "max_tokens": 500
            }
        )
        
        if response.status_code == 200:
            print("✅ Agent system simulation successful")
            return True
        elif response.status_code == 429:
            print("❌ Rate limit exceeded in agent system simulation")
            return False
        else:
            print(f"❌ Error: {response.status_code} - {response.text}")
            return False
    except Exception as e:
        print(f"❌ Error: {str(e)}")
        return False

def main():
    parser = argparse.ArgumentParser(description="Test API integration for Fastalysis")
    parser.add_argument("--complete", action="store_true", help="Run all tests")
    parser.add_argument("--connection", action="store_true", help="Test API connection only")
    parser.add_argument("--chat", action="store_true", help="Test chat completion only")
    parser.add_argument("--agent", action="store_true", help="Test agent system simulation")
    parser.add_argument("--files", action="store_true", help="Check required files")
    
    args = parser.parse_args()
    
    # Header
    print("\n" + "="*60)
    print("FASTALYSIS API INTEGRATION TEST")
    print("="*60)
    
    # Load environment variables
    load_env_variables()
    
    # Determine which tests to run
    run_all = args.complete or not any([args.connection, args.chat, args.agent, args.files])
    
    # Run tests
    results = {}
    
    if run_all or args.connection:
        results["api_connection"] = test_groq_api()
    
    if run_all or args.chat:
        results["chat_completion"] = test_chat_completion()
    
    if run_all or args.files:
        results["files_check"] = check_required_files()
    
    if run_all or args.agent:
        results["agent_system"] = simulate_agent_system()
    
    # Summary
    print("\n" + "="*60)
    print("TEST SUMMARY")
    print("="*60)
    
    for test, result in results.items():
        print(f"{test}: {'✅ PASSED' if result else '❌ FAILED'}")
    
    all_passed = all(results.values())
    print("\nOverall result:", "✅ ALL TESTS PASSED" if all_passed else "❌ SOME TESTS FAILED")
    
    if not all_passed:
        print("\nTroubleshooting tips:")
        print("1. Check your API keys in .env.local")
        print("2. Run 'python api_quota_monitor.py --check' to verify API status")
        print("3. Wait a few minutes if you're experiencing rate limits")
        print("4. See RATE_LIMIT_GUIDE.md for more troubleshooting tips")
    
    print("\n" + "="*60)

if __name__ == "__main__":
    main()