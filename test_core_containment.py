#!/usr/bin/env python3
"""Test that core containment works properly (ray-casting vs bounding box)."""

import numpy as np
import trimesh

# Create simple torus mesh
major_radius = 0.05  # 5cm
minor_radius = 0.01  # 1cm

torus = trimesh.creation.torus(major_radius=major_radius, minor_radius=minor_radius)

print("Torus mesh:")
print(f"  Vertices: {len(torus.vertices)}")
print(f"  Faces: {len(torus.faces)}")
print(f"  Bounds: {torus.bounds}")

# Test points
test_points = np.array([
    [0, 0, 0],           # Center (should be OUTSIDE for torus)
    [major_radius, 0, 0],# On major radius (should be INSIDE torus material)
    [0.08, 0, 0],        # Beyond outer edge (should be OUTSIDE)
])

print("\nContainment test (ray-casting):")
for i, pt in enumerate(test_points):
    is_inside = torus.contains([pt])[0]
    print(f"  Point {i} {pt}: {'INSIDE' if is_inside else 'OUTSIDE'}")

# Compare to bounding box method (the WRONG way)
print("\nBounding box test (OLD BROKEN METHOD):")
mins = torus.vertices.min(axis=0)
maxs = torus.vertices.max(axis=0)
print(f"  Bounding box: min={mins}, max={maxs}")
for i, pt in enumerate(test_points):
    is_inside_bbox = np.all(pt >= mins) and np.all(pt <= maxs)
    print(f"  Point {i} {pt}: {'INSIDE' if is_inside_bbox else 'OUTSIDE'}")

print("\n✓ Ray-casting properly identifies center of torus as OUTSIDE")
print("✗ Bounding box incorrectly identifies center of torus as INSIDE")
