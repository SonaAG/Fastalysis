#!/bin/bash

echo "🧬 Starting Fastalysis Genomics Research Assistant..."
echo ""

# Check if we're in the right directory
if [ ! -f "main.py" ]; then
    echo "❌ Error: main.py not found. Make sure you're in the app directory."
    echo "Current directory: $(pwd)"
    echo "Expected files: main.py, enhanced_streamlit_app.py, controller.py"
    exit 1
fi

echo "📋 Checking Python installation..."
if ! command -v python3 &> /dev/null; then
    echo "❌ Python3 not found. Please install Python 3.8+ first."
    exit 1
else
    echo "✅ Python found: $(python3 --version)"
fi

echo ""
echo "📦 Installing dependencies..."
python3 -m pip install -r requirements.txt
if [ $? -ne 0 ]; then
    echo "⚠️ Some packages might have failed to install. Continuing anyway..."
fi

echo ""
echo "🚀 Starting FastAPI backend..."
python3 main.py &
FASTAPI_PID=$!

echo "⏳ Waiting for FastAPI to start (5 seconds)..."
sleep 5

echo ""
echo "🎨 Starting Streamlit frontend..."
echo "📱 Your app will open at: http://localhost:8501"
echo "🔧 API documentation at: http://localhost:8000/docs"
echo ""
echo "🛑 To stop: Press Ctrl+C"
echo ""

# Function to cleanup background processes
cleanup() {
    echo ""
    echo "🛑 Stopping services..."
    kill $FASTAPI_PID 2>/dev/null
    exit 0
}

# Set trap to cleanup on script exit
trap cleanup EXIT INT TERM

streamlit run enhanced_streamlit_app.py