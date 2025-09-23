#!/usr/bin/env python3
"""
Simple test script to check if the Groq API key is working
"""

import os
import requests
from dotenv import load_dotenv

# Load environment variables from .env.local
load_dotenv('.env.local')

def test_groq_api():
    print("=== Groq API Test ===")
    
    # Get the API key
    api_key = os.getenv('GROQ_API_KEY')
    
    print(f"API Key found: {'Yes' if api_key else 'No'}")
    if api_key:
        print(f"API Key format: {api_key[:10]}..." if len(api_key) > 10 else api_key)
        print(f"API Key length: {len(api_key)}")
        print(f"Starts with 'gsk_': {api_key.startswith('gsk_')}")
    else:
        print("ERROR: No API key found!")
        print("Available environment variables with 'GROQ':", [k for k in os.environ.keys() if 'GROQ' in k.upper()])
        return
    
    # Test API call
    print("\n=== Making API Call ===")
    
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
                    {"role": "user", "content": "Hello, please just say 'API test successful'"}
                ],
                "max_tokens": 50
            },
            timeout=10
        )
        
        print(f"Response Status: {response.status_code}")
        print(f"Response Headers: {dict(response.headers)}")
        
        if response.status_code == 200:
            data = response.json()
            message = data.get('choices', [{}])[0].get('message', {}).get('content', 'No content')
            print(f"✅ SUCCESS: {message}")
        else:
            print(f"❌ FAILED: {response.status_code}")
            print(f"Response: {response.text}")
            
            # Parse error details if JSON
            try:
                error_data = response.json()
                if 'error' in error_data:
                    print(f"Error Type: {error_data['error'].get('type', 'Unknown')}")
                    print(f"Error Message: {error_data['error'].get('message', 'No message')}")
            except:
                print("Could not parse error response as JSON")
                
    except Exception as e:
        print(f"❌ EXCEPTION: {e}")

if __name__ == "__main__":
    test_groq_api()