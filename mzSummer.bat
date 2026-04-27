@echo off
REM Check if Python is installed
where python >nul 2>&1
IF ERRORLEVEL 1 (
    echo Python is not installed. Please install Python and try again.
    pause
    exit /b
)

REM Check if pip is installed
python -m pip --version >nul 2>&1
IF ERRORLEVEL 1 (
    echo pip is not installed. Installing pip...
    python -m ensurepip --upgrade
)

REM Install required dependencies
echo Installing dependencies from requirements.txt...
python -m pip install --upgrade pip
python -m pip install -r "%~dp0requirements.txt"

REM Run the Python script
echo Running the program...
python "%~dp0GUI3.py"

pause
