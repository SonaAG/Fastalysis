// server.js
require('dotenv').config({ path: '.env.local' });
const express = require('express');
const bodyParser = require('body-parser');
// Import your model logic here (using GROQ_API_KEY)
const { model } = require('./groq-model');

const app = express();
app.use(bodyParser.json());

app.post('/chat', async (req, res) => {
  const { message, model: modelID } = req.body;
  // Use GROQ_API_KEY from process.env.GROQ_API_KEY
  // Generate response using your model logic
  const response = await model.generate(message, modelID, process.env.GROQ_API_KEY);
  res.json({ response });
});

app.listen(3000, () => console.log('Chat API running on port 3000'));