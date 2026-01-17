#!/bin/bash
# MAGVID Simulator Setup Script

echo "Setting up MAGVID Electromagnetic Particle Simulator..."

# Check if Python is available
if ! command -v python3 &> /dev/null; then
    echo "Python 3 is required but not installed. Please install Python 3.7+ first."
    exit 1
fi

# Create virtual environment if it doesn't exist
if [ ! -d "magvid_env" ]; then
    echo "Creating virtual environment..."
    python3 -m venv magvid_env
fi

# Activate virtual environment
echo "Activating virtual environment..."
source magvid_env/bin/activate

# Upgrade pip
echo "Upgrading pip..."
pip install --upgrade pip

# Install requirements
echo "Installing requirements..."
pip install -r requirements.txt

echo "Installation complete!"
echo ""
echo "To run the simulator:"
echo "  source magvid_env/bin/activate"
echo "  python magvid_simulator.py"
echo ""
echo "To deactivate the environment when done:"
echo "  deactivate"