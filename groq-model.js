// groq-model.js
const axios = require("axios");

const MODEL_MAP = {
  "kimi-k2": "moonshotai/kimi-k2-instruct",
  "meta-llama/llama-4-scout-17b-16e-instruct": "meta-llama/llama-4-scout-17b-16e-instruct",
  "llama-3.1-8b-instant": "llama-3.1-8b-instant",
  "deepseek-r1-distill-llama-70b": "deepseek-r1-distill-llama-70b",
  "llama-3.3-70b-versatile": "llama-3.3-70b-versatile",
};

async function generate(message, modelID, apiKey) {
  const modelName = MODEL_MAP[modelID] || MODEL_MAP["kimi-k2"];
  try {
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
      }
    );
    return response.data.choices[0].message.content;
  } catch (err) {
    return "Error: Unable to fetch response from Groq API.";
  }
}

module.exports = { model: { generate } };
