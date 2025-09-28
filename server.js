// server.js
require('dotenv').config({ path: '.env.local' });
const express = require('express');
const bodyParser = require('body-parser');
// Import your model logic here (using API keys from environment)
const { model } = require('./groq-model');

// Log available API keys (without showing the actual keys)
console.log("API Keys configured:");
console.log("- GROQ_API_KEY:", process.env.GROQ_API_KEY ? "✅ Present" : "❌ Missing");
console.log("- HUGGINGFACE_API_KEY:", process.env.HUGGINGFACE_API_KEY ? "✅ Present" : "❌ Missing");

const app = express();
app.use(bodyParser.json());

app.post('/chat', async (req, res) => {
  try {
    const { message, model: modelID, sequence } = req.body;
    console.log(`Received chat request: model=${modelID}, has_sequence=${!!sequence}`);
    
    // Add sequence context if available
    let enhancedMessage = message;
    
    if (sequence) {
      enhancedMessage = `User query: ${message}\n\nCurrent sequence in context: ${sequence.substring(0, 100)}${sequence.length > 100 ? '...' : ''}\n\nPlease consider the above sequence when answering the query.`;
      console.log("Added sequence context to query");
    }
    
    // Use GROQ_API_KEY from process.env.GROQ_API_KEY
    const response = await model.generate(enhancedMessage, modelID, process.env.GROQ_API_KEY);
    res.json({ response });
  } catch (error) {
    console.error("Error in chat endpoint:", error);
    res.status(500).json({ 
      error: "An error occurred processing your request",
      message: error.message
    });
  }
});

app.listen(3000, () => console.log('Chat API running on port 3000'));