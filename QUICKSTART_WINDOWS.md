# Quick Start for Windows Users

**Problem**: Getting "DLL load failed while importing OCP" error on Windows?

**Quick Solution**: Use Conda instead of pip.

## Fastest Working Method (5 minutes)

1. **Install Miniconda** (if not already installed):
   - Download: https://docs.conda.io/en/latest/miniconda.html
   - Run installer, accept defaults

2. **Open Anaconda Prompt** (not regular command prompt)

3. **Run these commands**:
   ```bash
   # Navigate to project folder
   cd path\to\magvidsim

   # Create environment with Python 3.9
   conda create -n magvid_env python=3.9

   # Activate environment
   conda activate magvid_env

   # Install CADQuery from conda-forge (this avoids the DLL issue)
   conda install -c conda-forge -c cadquery cadquery

   # Install remaining dependencies
   pip install pyvista numpy scipy matplotlib jupyter ipywidgets casadi trimesh
   ```

4. **Test it**:
   ```bash
   python em_solver.py
   ```

   This should generate `magvid_field_visualization.html` - open it in your browser!

## Why This Works

The DLL error happens because:
- `cadquery-ocp` package contains compiled C++ code (OpenCascade bindings)
- On Windows, it needs specific Visual C++ runtime DLLs
- Pip's wheel might not include all required DLLs
- Conda includes pre-compiled binaries with all dependencies bundled

## Alternative: Just Fix the DLL Issue

If you prefer to stick with pip:

1. **Install Visual C++ Redistributable**:
   - Download: https://aka.ms/vs/17/release/vc_redist.x64.exe
   - Install and restart your computer

2. **Try again**:
   ```bash
   python em_solver.py
   ```

## Python 3.13 Users

If you're on Python 3.13, you may need to reinstall VTK/PyVista after installing requirements:

```bash
pip uninstall vtk pyvista trame trame-vtk trame-vuetify trame-client trame-server pooch scooby -y
pip cache purge
pip install --upgrade --force-reinstall vtk pyvista[all]
```

This fixes issues where packages try to use nonexistent OS library APIs.

## Still Having Issues?

See `WINDOWS_INSTALL.md` for comprehensive troubleshooting, Docker setup, and more.

## Just Want the Particle Simulator?

The particle trajectory simulator doesn't need CADQuery:

```bash
# Regular Python environment works fine
python -m venv magvid_env
.\magvid_env\Scripts\activate
pip install numpy scipy matplotlib jupyter ipywidgets casadi
python simulator.py
```
