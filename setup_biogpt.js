/**
 * Setup BioGPT API Key
 * This script helps set up the HUGGINGFACE_API_KEY in your .env.local file
 */

const fs = require('fs');
const path = require('path');
const readline = require('readline');

const rl = readline.createInterface({
  input: process.stdin,
  output: process.stdout
});

console.log('=====================================');
console.log('BioGPT Setup for Fastalysis');
console.log('=====================================');
console.log('\nThis script will help you set up your HUGGINGFACE_API_KEY for BioGPT integration.\n');

// Check if .env.local exists
const envLocalPath = path.join(__dirname, '..', '.env.local');
const appEnvLocalPath = path.join(__dirname, '.env.local');
let existingContent = '';
let targetPath = '';

if (fs.existsSync(envLocalPath)) {
  existingContent = fs.readFileSync(envLocalPath, 'utf8');
  targetPath = envLocalPath;
  console.log('Found .env.local file in the project root directory.');
} else if (fs.existsSync(appEnvLocalPath)) {
  existingContent = fs.readFileSync(appEnvLocalPath, 'utf8');
  targetPath = appEnvLocalPath;
  console.log('Found .env.local file in the app directory.');
} else {
  console.log('No existing .env.local file found. Will create one in the project root.');
  targetPath = envLocalPath;
}

// Check if key already exists
const keyRegex = /HUGGINGFACE_API_KEY=([^\n]*)/;
const keyMatch = existingContent.match(keyRegex);

if (keyMatch) {
  console.log('\nExisting HUGGINGFACE_API_KEY found: ' + keyMatch[1]);
  rl.question('\nWould you like to replace it? (y/n): ', (answer) => {
    if (answer.toLowerCase() === 'y') {
      promptForNewKey();
    } else {
      console.log('\nKeeping existing API key. No changes made.');
      rl.close();
    }
  });
} else {
  console.log('\nNo existing HUGGINGFACE_API_KEY found in .env.local.');
  promptForNewKey();
}

function promptForNewKey() {
  console.log('\n1. Go to https://huggingface.co/ and create an account if you don\'t have one');
  console.log('2. Navigate to your profile settings and create an API key');
  console.log('3. Copy the generated API key and paste it below\n');

  rl.question('Enter your Hugging Face API Key: ', (apiKey) => {
    if (!apiKey.trim()) {
      console.log('\nAPI key cannot be empty. Please try again.');
      promptForNewKey();
      return;
    }

    try {
      let newContent = existingContent;
      if (keyMatch) {
        // Replace existing key
        newContent = existingContent.replace(keyRegex, `HUGGINGFACE_API_KEY=${apiKey}`);
      } else {
        // Add new key
        if (existingContent && !existingContent.endsWith('\n')) {
          newContent += '\n';
        }
        newContent += `HUGGINGFACE_API_KEY=${apiKey}\n`;
      }

      fs.writeFileSync(targetPath, newContent);
      console.log(`\n✅ Successfully ${keyMatch ? 'updated' : 'added'} HUGGINGFACE_API_KEY in ${targetPath}`);
      console.log('\nNext steps:');
      console.log('1. Restart your Node.js server with: node app/server.js');
      console.log('2. Restart your Streamlit application');
      console.log('3. Select "microsoft/BioGPT-Large" from the model dropdown in the chat interface');
      rl.close();
    } catch (error) {
      console.error('\n❌ Error writing to .env.local file:', error.message);
      console.log('Please make sure you have write permissions for this directory.');
      rl.close();
    }
  });
}