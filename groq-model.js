// groq-model.js
const axios = require("axios");

const MODEL_MAP = {
  "kimi-k2": "moonshotai/kimi-k2-instruct",
  "meta-llama/llama-4-scout-17b-16e-instruct": "meta-llama/llama-4-scout-17b-16e-instruct",
  "llama-3.1-8b-instant": "llama-3.1-8b-instant",
  "deepseek-r1-distill-llama-70b": "deepseek-r1-distill-llama-70b",
  "llama-3.3-70b-versatile": "llama-3.3-70b-versatile",
  "claude-3-haiku-20240307": "claude-3-haiku-20240307", // Alternative to BioGPT
};

async function generate(message, modelID, apiKey) {
  const modelName = MODEL_MAP[modelID] || MODEL_MAP["kimi-k2"];
  
  // Use Groq API for all models
  // BioGPT has been removed
  return await generateWithGroq(message, modelName, apiKey);
}

async function generateWithGroq(message, modelName, apiKey) {
  // Implement exponential backoff retry logic
  const maxRetries = 3;
  let retryDelay = 1000; // Start with 1 second delay
  
  for (let attempt = 0; attempt < maxRetries; attempt++) {
    try {
      console.log(`Attempt ${attempt + 1}/${maxRetries} for model ${modelName}`);
      
      const response = await axios.post(
        "https://api.groq.com/openai/v1/chat/completions",
        {
          model: modelName,
          messages: [{ role: "user", content: message }],
        },
        {
          headers: {
            "Authorization": `Bearer ${apiKey}`,
            "Content-Type": "application/json",
          },
          timeout: 30000, // 30 second timeout
        }
      );
      
      // If successful, return the content
      return response.data.choices[0].message.content;
      
    } catch (err) {
      // Log detailed error info
      console.error(`Groq API Error (Attempt ${attempt + 1}/${maxRetries}):`, err.message);
      
      // Check if this is a rate limit error (HTTP 429)
      if (err.response && err.response.status === 429) {
        // Get retry delay from response headers if available
        const retryAfter = err.response.headers['retry-after'];
        const waitTime = retryAfter ? parseInt(retryAfter) * 1000 : retryDelay;
        
        console.log(`Rate limited. Waiting ${waitTime/1000} seconds before retry.`);
        
        // Don't retry on last attempt
        if (attempt < maxRetries - 1) {
          await new Promise(resolve => setTimeout(resolve, waitTime));
          retryDelay *= 2; // Exponential backoff
          continue; // Try again
        }
        
        return "Error: Unable to fetch response from Groq API. Rate limit exceeded. Please try again in a few minutes or switch to a different model.";
      }
      
      // For other errors or final attempt, return error message
      return `Error: Unable to fetch response from Groq API. ${err.message}`;
    }
  }
  
  return "Error: Maximum retry attempts exceeded. Please try again later.";
}

module.exports = { model: { generate } };
