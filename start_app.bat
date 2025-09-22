@echo off
echo 🧬 Starting Fastalysis Genomics Research Assistant...
echo.

REM Check if we're in the right directory
if not exist "main.py" (
    echo ❌ Error: main.py not found. Make sure you're in the app directory.
    echo Current directory: %CD%
    echo Expected files: main.py, enhanced_streamlit_app.py, controller.py
    pause
    exit /b 1
)

echo 📋 Checking Python installation...
python --version >nul 2>&1
if errorlevel 1 (
    echo ❌ Python not found. Please install Python 3.8+ first.
    pause
    exit /b 1
) else (
    echo ✅ Python found
)

echo.
echo 📦 Installing dependencies...
pip install -r requirements.txt
if errorlevel 1 (
    echo ⚠️ Some packages might have failed to install. Continuing anyway...
)

echo.
echo 🚀 Starting FastAPI backend (this will run in background)...
start "FastAPI Backend" cmd /k "echo 🔧 FastAPI Server Starting... && python main.py"

echo ⏳ Waiting for FastAPI to start (5 seconds)...
timeout /t 5 /nobreak >nul

echo.
echo 🎨 Starting Streamlit frontend...
echo 📱 Your app will open at: http://localhost:8501
echo 🔧 API documentation at: http://localhost:8000/docs
echo.
echo 🛑 To stop: Close both terminal windows or press Ctrl+C
echo.

streamlit run enhanced_streamlit_app.py

pause