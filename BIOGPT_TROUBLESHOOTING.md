# BioGPT Troubleshooting Guide

## Issue: Unable to fetch response from BioGPT API

If you're seeing the error message:

```
Error: Unable to fetch response from BioGPT API. Please check if the HUGGINGFACE_API_KEY is set.
```

This indicates that the application is unable to communicate with the Hugging Face API to use the BioGPT model. Follow these steps to resolve the issue:

## Solution: Set up the HUGGINGFACE_API_KEY

### Method 1: Using the Setup Script (Recommended)

1. Open a terminal/command prompt
2. Navigate to the Fastalysis app directory:
   ```
   cd d:\Fastalysis\app
   ```
3. Run the setup script:
   ```
   node setup_biogpt.js
   ```
4. Follow the on-screen instructions to enter your Hugging Face API key

### Method 2: Manual Setup

1. **Get a Hugging Face API Key**:
   - Go to [Hugging Face](https://huggingface.co/) and create an account if you don't have one
   - Navigate to your profile → Settings → Access Tokens
   - Click "New token"
   - Give it a name (e.g., "Fastalysis")
   - Select "Read" role is sufficient
   - Click "Generate a token"
   - Copy the generated API key

2. **Add the API Key to your environment**:
   - In the Fastalysis root folder, create or edit a file named `.env.local`
   - Add the following line:
     ```
     HUGGINGFACE_API_KEY=your_api_key_here
     ```
   - Save the file

   - Alternatively, create the file in the app directory:
     ```
     d:\Fastalysis\app\.env.local
     ```

3. **Verify the API Key is Set**:
   - After saving the file, verify that it contains the correct key
   - The file should have the format:
     ```
     HUGGINGFACE_API_KEY=hf_xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
     ```

### After Setting the API Key

1. **Restart the Node.js server**:
   - If the server is running, stop it (Ctrl+C)
   - Restart it with:
     ```
     node app/server.js
     ```

2. **Restart the Streamlit application**:
   - If Streamlit is running, stop it (Ctrl+C)
   - Restart it with:
     ```
     cd app
     streamlit run enhanced_streamlit_app.py
     ```

3. **Test BioGPT**:
   - In the chat interface, select "microsoft/BioGPT-Large" from the model dropdown
   - Try asking a genomics-related question
   - The response should now come from BioGPT instead of showing the error

## Checking API Key Status

To check if your API key is properly recognized by the system:

1. In the Streamlit application, click on "Debug Tools" in the sidebar
2. Click "Check Environment Variables"
3. Look for the HUGGINGFACE_API_KEY entry
4. It should show "✅ Present" if the key is properly set

## Troubleshooting Common Issues

### API Key Not Recognized

If your key is in the .env.local file but the system doesn't recognize it:

1. Ensure there are no spaces around the equals sign
2. Make sure the file is named exactly `.env.local` (not `.env.local.txt`)
3. Check that the key format is correct (should start with "hf_")
4. Try moving the .env.local file to both the project root and the app directory

### Permission Issues

If you can't save the .env.local file:

1. Run your editor as administrator
2. Or manually create the file in a location where you have write permissions, then move it to the correct location

### Server Issues

If restarting the server doesn't help:

1. Check the server logs for any specific error messages
2. Ensure no other processes are using the same port
3. Try rebooting your system to clear any lingering processes

### Need Help?

If you continue to experience issues after following these steps, please check the additional documentation in the README_BIOGPT.md file or contact support.