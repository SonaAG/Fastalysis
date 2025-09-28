# Using Microsoft BioGPT with Fastalysis

This document explains how to set up and use Microsoft's BioGPT model with the Fastalysis application.

## What is BioGPT?

BioGPT is a large language model developed by Microsoft specifically for biomedical text generation. It was trained on a vast corpus of biomedical literature, making it particularly suited for:

- Genomics research questions
- Understanding mutations and their effects
- Disease associations with genes
- Molecular pathways and interactions
- Interpreting biomedical research findings

## Setup Instructions

To use BioGPT with Fastalysis, follow these steps:

1. **Get a Hugging Face API Key**:
   - Go to [Hugging Face](https://huggingface.co/) and create an account if you don't have one
   - Navigate to your profile settings and create an API key
   - Copy the generated API key

2. **Add the API Key to your environment**:
   - In the root folder of Fastalysis, create or edit the `.env.local` file
   - Add the following line:
     ```
     HUGGINGFACE_API_KEY=your_api_key_here
     ```
   - Save the file

3. **Restart the Node.js server**:
   - Stop any running instance of the Node.js server
   - Restart it with:
     ```
     node app/server.js
     ```

## Using BioGPT in the Application

1. Launch the Streamlit application
2. In the chat interface, select "microsoft/BioGPT-Large" from the model dropdown
3. Your queries will now be processed by the BioGPT model

## Best Practices for BioGPT Queries

For optimal results with BioGPT:

- Be specific with biomedical terminology
- Provide context about the gene, disease, or mutation you're inquiring about
- Frame questions in a scientific manner
- Consider longer, more detailed prompts for complex topics

## Limitations

- BioGPT may have knowledge cutoff from when it was trained
- For very recent genomic discoveries, other models might be more up-to-date
- Processing time may be slower than some other models
- Maximum context length may be more limited than larger models

## Troubleshooting

If you encounter errors:
- Verify your Hugging Face API key is correctly set in `.env.local`
- Check the terminal logs for any API-related errors
- Ensure you have internet connectivity for API calls
- If queries time out, try breaking them into smaller, more focused questions