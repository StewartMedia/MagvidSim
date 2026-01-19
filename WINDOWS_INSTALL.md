# Windows Installation Guide

This guide addresses Windows-specific installation issues, particularly with the CADQuery dependencies required for the EM field solver.

## Quick Start (Particle Simulator Only)

If you only need the particle trajectory simulator and don't need the EM field solver, you can skip CADQuery:

```powershell
# Create virtual environment
python -m venv magvid_env
.\magvid_env\Scripts\activate

# Install core dependencies (excluding CADQuery)
pip install numpy scipy matplotlib jupyter ipywidgets casadi
```

Then run:
```powershell
python simulator.py
```

## Full Installation (Including EM Field Solver)

The EM field solver requires CADQuery, which has complex native dependencies on Windows.

### Option 1: Conda (Recommended)

Conda provides pre-compiled binaries that avoid DLL issues:

```powershell
# Install Miniconda from: https://docs.conda.io/en/latest/miniconda.html

# Create environment
conda create -n magvid_env python=3.9
conda activate magvid_env

# Install CADQuery from conda-forge
conda install -c conda-forge -c cadquery cadquery

# Install remaining dependencies
pip install pyvista numpy scipy matplotlib jupyter ipywidgets casadi trimesh
```

### Option 2: Fix Visual C++ Runtime

If you prefer pip installation, the OCP DLL error usually means missing Visual C++ runtime libraries:

1. **Download and install Visual C++ Redistributable**:
   - 64-bit: https://aka.ms/vs/17/release/vc_redist.x64.exe
   - Run the installer and restart your computer

2. **Install Python packages**:
   ```powershell
   python -m venv magvid_env
   .\magvid_env\Scripts\activate
   pip install -r requirements.txt
   ```

3. **Verify installation**:
   ```powershell
   python -c "import cadquery as cq; print('CADQuery OK')"
   ```

### Option 3: Docker (Most Reliable)

Docker provides a consistent Linux environment on Windows:

1. **Install Docker Desktop**:
   - Download from: https://www.docker.com/products/docker-desktop

2. **Create Dockerfile** (in project root):
   ```dockerfile
   FROM python:3.9-slim

   WORKDIR /app

   # Install system dependencies
   RUN apt-get update && apt-get install -y \
       libgl1-mesa-glx \
       libglib2.0-0 \
       libsm6 \
       libxext6 \
       libxrender-dev \
       && rm -rf /var/lib/apt/lists/*

   # Copy requirements and install Python dependencies
   COPY requirements.txt .
   RUN pip install --no-cache-dir -r requirements.txt

   # Copy project files
   COPY . .

   CMD ["python", "em_solver.py"]
   ```

3. **Build and run**:
   ```powershell
   docker build -t magvidsim .
   docker run -v ${PWD}:/app magvidsim
   ```

## Common Issues

### Issue 1: "DLL load failed while importing OCP"

**Cause**: Missing Visual C++ runtime or incompatible OpenCascade binaries

**Solutions**:
1. Install VC++ Redistributable (see Option 2 above)
2. Use Conda instead of pip (see Option 1 above)
3. Use Docker (see Option 3 above)

### Issue 2: "pyvista" or "cadquery" not found

**Cause**: Out-of-date requirements.txt

**Solution**: The current `requirements.txt` includes both packages. Ensure you're using the latest version from the repository:
```powershell
git pull origin main
pip install -r requirements.txt
```

### Issue 3: Python version mismatch

**Symptom**: Various import errors or compatibility warnings

**Solution**: This project requires Python 3.9-3.11. Check your version:
```powershell
python --version
```

If needed, install the correct version or use Conda to create an environment with the right Python version.

### Issue 4: Interactive visualization doesn't work

**Symptom**: PyVista viewer fails to open or crashes

**Solution**: The EM solver automatically saves to HTML format for browser viewing:
```powershell
python em_solver.py
# Opens: magvid_field_visualization.html
```

Open the generated HTML file in your browser (Chrome, Firefox, Edge all work).

## Testing Your Installation

### Test 1: Particle Simulator
```powershell
python simulator.py
```
Should display 3D particle trajectories.

### Test 2: EM Field Solver (Basic)
```powershell
python -c "import cadquery as cq; import pyvista as pv; print('Dependencies OK')"
```
Should print "Dependencies OK" without errors.

### Test 3: Full EM Solver
```powershell
python em_solver.py
```
Should generate `magvid_field_visualization.html` with 3D field visualization.

## Python Version Compatibility

| Python Version | Status | Notes |
|----------------|--------|-------|
| 3.8 | Not tested | May work but not recommended |
| 3.9 | ✅ Recommended | Best tested compatibility |
| 3.10 | ✅ Supported | Known to work |
| 3.11 | ✅ Supported | Known to work |
| 3.12 | ✅ Supported | Tested on macOS |
| 3.13 | ✅ Works | Confirmed working on Windows (requires VTK/PyVista reinstall, see below) |

### Python 3.13 Specific Issues

If using Python 3.13, you may need to reinstall VTK and PyVista packages:

```powershell
pip uninstall vtk pyvista trame trame-vtk trame-vuetify trame-client trame-server pooch scooby -y
pip cache purge
pip install --upgrade --force-reinstall vtk pyvista[all]
```

**Note**: Some packages call nonexistent OS library APIs in Python 3.13. The reinstall process ensures compatible binary versions are installed.

## Getting Help

If you encounter issues not covered here:

1. Check that all dependencies are installed:
   ```powershell
   pip list | findstr "cadquery pyvista numpy"
   ```

2. Verify Python version compatibility:
   ```powershell
   python --version
   ```

3. Try the Conda installation method (most reliable on Windows)

4. Report the issue with full error output and system information:
   - Python version: `python --version`
   - OS version: `winver`
   - Installation method: pip/conda/docker
   - Full error traceback
