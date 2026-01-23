@echo off
setlocal

pushd "%~dp0"

echo Checking for Python installation...
python --version >nul 2>&1
if %errorlevel% neq 0 (
    echo.
    echo ERROR: Python is not installed or not on PATH.
    echo Please install Python 3.9+ from https://www.python.org/downloads/
    exit /b 1
)

echo Python found: 
python --version

echo.
echo Creating virtual environment in .venv ...

python -m venv .venv
if %errorlevel% neq 0 (
    echo.
    echo ERROR: Failed to create virtual environment.
    exit /b 1
)

echo.
echo Activating virtual environment...
call .venv\Scripts\activate

if %errorlevel% neq 0 (
    echo.
    echo ERROR: Failed to activate virtual environment.
    exit /b 1
)

echo.
if exist requirements.txt (
    echo Installing dependencies from requirements.txt...
    pip install --upgrade pip
    pip install -r requirements.txt
) else (
    echo No requirements.txt found. Skipping dependency installation.
)

call .venv\Scripts\deactivate.bat

echo.
echo Setup complete.
echo.

popd