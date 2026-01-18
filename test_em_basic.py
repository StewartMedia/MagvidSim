#!/usr/bin/env python3
"""Quick test of electromagnetic solver components."""

import numpy as np
import cadquery as cq
import pyvista as pv

print("="*60)
print("Testing EM Solver Components")
print("="*60)

# Test 1: Load STEP file
print("\n1. Loading STEP file...")
try:
    shape = cq.importers.importStep('out.STEP')
    print("   ✓ STEP file loaded successfully")
except Exception as e:
    print(f"   ✗ Error: {e}")
    exit(1)

# Test 2: Export to STL for tessellation
print("\n2. Tessellating geometry...")
try:
    import tempfile
    import os
    with tempfile.NamedTemporaryFile(suffix='.stl', delete=False) as tmp:
        tmp_name = tmp.name

    cq.exporters.export(shape, tmp_name)

    import trimesh
    mesh = trimesh.load(tmp_name)
    os.unlink(tmp_name)

    print(f"   ✓ Mesh: {len(mesh.vertices)} vertices, {len(mesh.faces)} faces")
    print(f"   ✓ Bounds: {mesh.bounds}")
except Exception as e:
    print(f"   ✗ Error: {e}")
    import traceback
    traceback.print_exc()
    exit(1)

# Test 3: Simple field calculation
print("\n3. Testing field calculation...")
vertices = mesh.vertices * 0.001  # Convert mm to m
faces = mesh.faces

# Calculate field at one test point
test_point = np.array([0.05, 0.0, 0.0])  # 5cm from origin

print(f"   Test point: {test_point}")
print(f"   Calculating B-field using Biot-Savart...")

# Simplified calculation - just test the math
MU_0 = 4 * np.pi * 1e-7
current = -500.0  # Amperes
B = np.zeros(3)

# Sample a few faces
for i, face in enumerate(faces[:100]):  # Just first 100 faces for speed
    v0, v1, v2 = vertices[face]
    center = (v0 + v1 + v2) / 3

    edge1 = v1 - v0
    edge2 = v2 - v0
    normal = np.cross(edge1, edge2)
    area = 0.5 * np.linalg.norm(normal)

    r_vec = test_point - center
    r = np.linalg.norm(r_vec)

    if r > 1e-9:
        J_mag = current / (mesh.area * 0.001**2)  # Current density
        J = J_mag * normal / (2 * area) if area > 0 else np.zeros(3)
        dB = (MU_0 / (4 * np.pi)) * np.cross(J, r_vec) / (r**3) * area
        B += dB

print(f"   ✓ B-field at test point: {B}")
print(f"   ✓ |B| = {np.linalg.norm(B):.6e} T")

# Test 4: PyVista visualization
print("\n4. Testing PyVista...")
try:
    faces_pv = np.hstack([[3] + list(f) for f in faces[:1000]])  # Subset for speed
    pv_mesh = pv.PolyData(vertices, faces_pv)

    plotter = pv.Plotter(off_screen=True)
    plotter.add_mesh(pv_mesh, color='lightgray')

    # Try to export HTML
    print("   Exporting to HTML...")
    plotter.export_html('test_visualization.html')
    print("   ✓ Exported to test_visualization.html")

except Exception as e:
    print(f"   ✗ Error: {e}")
    import traceback
    traceback.print_exc()

print("\n" + "="*60)
print("All tests completed!")
print("="*60)
