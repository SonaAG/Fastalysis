# 🧬 Fastalysis - Quick Start Guide

## 🚀 How to Run Your Genomics Research Assistant

### Option 1: Automatic Startup (Recommended)

**For Windows:**
```bash
# Navigate to your app directory
cd d:\Fastalysis\app

# Run the startup script
start_app.bat
```

**For Linux/Mac:**
```bash
# Navigate to your app directory  
cd /path/to/Fastalysis/app

# Make script executable
chmod +x start_app.sh

# Run the startup script
./start_app.sh
```

### Option 2: Manual Step-by-Step

1. **Open Terminal/Command Prompt**
   ```bash
   cd d:\Fastalysis\app
   ```

2. **Install Dependencies** (first time only)
   ```bash
   pip install -r requirements.txt
   ```

3. **Start FastAPI Backend** (Terminal 1)
   ```bash
   python main.py
   ```
   - Should show: "Uvicorn running on http://0.0.0.0:8000"
   - Keep this terminal open

4. **Start Streamlit Frontend** (Terminal 2)
   ```bash
   streamlit run enhanced_streamlit_app.py
   ```
   - Should show: "You can now view your Streamlit app in your browser"
   - App opens at: http://localhost:8501

## 🌐 Access Points

- **Main App**: http://localhost:8501 (Streamlit UI)
- **API Docs**: http://localhost:8000/docs (FastAPI documentation)
- **API Base**: http://localhost:8000 (Backend endpoints)

## 📋 Prerequisites Check

### Required Software:
- Python 3.8 or higher
- pip (Python package manager)

### Check Your Setup:
```bash
# Check Python version
python --version

# Check pip
pip --version

# Check if files exist
dir main.py enhanced_streamlit_app.py controller.py  # Windows
ls main.py enhanced_streamlit_app.py controller.py   # Linux/Mac
```

## 🔧 Troubleshooting

### Common Issues:

1. **"Python not found"**
   - Install Python from https://python.org
   - Make sure Python is in your PATH

2. **"Module not found" errors**
   ```bash
   pip install --upgrade pip
   pip install -r requirements.txt
   ```

3. **Port already in use**
   - FastAPI (8000): `netstat -an | findstr :8000`
   - Streamlit (8501): `netstat -an | findstr :8501`
   - Kill processes or use different ports

4. **Permission errors**
   - Windows: Run Command Prompt as Administrator
   - Linux/Mac: Use `sudo` if needed

### Manual Package Installation:
```bash
# Core packages
pip install fastapi uvicorn streamlit
pip install biopython requests plotly pandas

# AI/ML packages  
pip install autogen-agentchat chromadb sentence-transformers

# If you get dependency conflicts:
pip install --upgrade pip setuptools wheel
```

## 🎯 First Time Setup

1. **Run the app** using Option 1 or 2 above
2. **Open browser** to http://localhost:8501
3. **Test with sample data**:
   - Upload a test FASTA file or paste sequence
   - Try chat: "What is BRCA1?"
   - Click analysis buttons

## 📱 Using the App

### Chat Interface:
- Ask: "Analyze mutations in BRCA1"
- Ask: "What does this sequence do: ATCGATCG..."
- Ask: "Find papers about TP53 mutations"

### Analysis Tools:
- Upload FASTA files (supports .fasta, .fa, .fas, .txt)
- Run BLAST searches
- Analyze mutations
- Look up genetic variants
- Search literature

### File Upload:
- Drag & drop FASTA files
- Supports multiple sequences
- Auto-detects DNA vs Protein
- No size limits for sequences

## 🆘 Getting Help

If you encounter issues:

1. **Check the terminal output** for error messages
2. **Verify all files are present**:
   ```
   main.py
   enhanced_streamlit_app.py
   controller.py
   genomics_agents.py
   genomics_rag.py
   requirements.txt
   ```
3. **Check Python version**: Must be 3.8+
4. **Try manual installation** of problem packages

## 🔄 Stopping the App

- **Windows**: Close both terminal windows or press Ctrl+C
- **Linux/Mac**: Press Ctrl+C in the terminal running the script

## 🚀 Next Steps

Once running successfully:
1. Test with your own FASTA files
2. Configure your preferred AI models
3. Add custom genomics knowledge to RAG system
4. Explore API endpoints at http://localhost:8000/docs