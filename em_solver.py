#!/usr/bin/env python3
"""
Complete Electromagnetic Field Solver for MAGVID
=================================================

Loads STEP AP214 files and computes magnetostatic fields.
"""

import numpy as np
import pyvista as pv
import cadquery as cq
from dataclasses import dataclass
from typing import List, Tuple
import re


MU_0 = 4 * np.pi * 1e-7  # Permeability of free space [H/m]


@dataclass
class CurrentSheet:
    """Current-carrying surface from STEP file."""
    name: str
    vertices: np.ndarray  # Nx3 array
    faces: np.ndarray     # Mx3 array of vertex indices
    current: float  # Amperes

    @property
    def area(self):
        """Total surface area."""
        areas = []
        for face in self.faces:
            v0, v1, v2 = self.vertices[face]
            areas.append(0.5 * np.linalg.norm(np.cross(v1 - v0, v2 - v0)))
        return np.sum(areas)


@dataclass
class MagneticCore:
    """High-permeability magnetic core."""
    vertices: np.ndarray
    faces: np.ndarray
    mu_r: float

    def extract_surface_region(self, winding_bounds: np.ndarray, margin: float = 0.001) -> tuple:
        """Extract core surface faces within winding bounding box.

        Args:
            winding_bounds: 2x3 array [[x_min, y_min, z_min], [x_max, y_max, z_max]]
            margin: Margin to add to bounds (meters)

        Returns:
            (vertices, faces) - the subset of core surface in the bounded region
        """
        # Expand bounds by margin
        bounds_min = winding_bounds[0] - margin
        bounds_max = winding_bounds[1] + margin

        # Vectorized computation of face centroids
        # Get all vertices for all faces at once
        v0 = self.vertices[self.faces[:, 0]]
        v1 = self.vertices[self.faces[:, 1]]
        v2 = self.vertices[self.faces[:, 2]]

        # Compute all centroids at once
        centroids = (v0 + v1 + v2) / 3.0

        # Check which centroids are within bounds (vectorized)
        in_bounds = np.all((centroids >= bounds_min) & (centroids <= bounds_max), axis=1)

        # Select faces where centroid is in bounds
        faces_in_region = self.faces[in_bounds]

        print(f"  Extracted {len(faces_in_region)}/{len(self.faces)} core faces "
              f"in winding region", flush=True)

        return self.vertices, faces_in_region

    def _ray_triangle_intersection(self, ray_origin, ray_dir, v0, v1, v2):
        """Möller–Trumbore ray-triangle intersection algorithm."""
        EPSILON = 1e-9

        edge1 = v1 - v0
        edge2 = v2 - v0
        h = np.cross(ray_dir, edge2)
        a = np.dot(edge1, h)

        if abs(a) < EPSILON:
            return False

        f = 1.0 / a
        s = ray_origin - v0
        u = f * np.dot(s, h)

        if u < 0.0 or u > 1.0:
            return False

        q = np.cross(s, edge1)
        v = f * np.dot(ray_dir, q)

        if v < 0.0 or u + v > 1.0:
            return False

        t = f * np.dot(edge2, q)
        return t > EPSILON

    def contains_point(self, point: np.ndarray) -> bool:
        """Check if point is inside core using ray-casting (odd/even rule)."""
        # Cast a ray from point in +X direction and count triangle intersections
        ray_origin = point
        ray_dir = np.array([1.0, 0.0, 0.0])

        intersection_count = 0
        for face in self.faces:
            v0, v1, v2 = self.vertices[face]
            if self._ray_triangle_intersection(ray_origin, ray_dir, v0, v1, v2):
                intersection_count += 1

        # Odd = inside, Even = outside
        return intersection_count % 2 == 1

    def contains_points(self, points: np.ndarray) -> np.ndarray:
        """Vectorized version for multiple points."""
        return np.array([self.contains_point(pt) for pt in points])


class STEPLoader:
    """Load geometry from STEP AP214 files using CADQuery."""

    def __init__(self, filepath: str):
        self.filepath = filepath
        self.shape = None
        self.winding_layers = {}
        self.step_entities = {}

    def load(self):
        """Load STEP file."""
        print(f"Loading STEP file: {self.filepath}")
        try:
            self.shape = cq.importers.importStep(self.filepath)
            print(f"  Successfully loaded STEP file")
            self._parse_step_entities()
            return True
        except Exception as e:
            print(f"  Error loading STEP: {e}")
            return False

    def _parse_step_entities(self):
        """Parse STEP file to identify shells and layer assignments."""
        with open(self.filepath, 'r') as f:
            content = f.read()

        # Find CLOSED_SHELL (solid core) and OPEN_SHELL (winding surfaces)
        closed_shell_pattern = r'#(\d+)\s*=\s*CLOSED_SHELL\s*\(\s*\'[^\']*\'\s*,\s*\(([^)]+)\)'
        open_shell_pattern = r'#(\d+)\s*=\s*OPEN_SHELL\s*\(\s*\'[^\']*\'\s*,\s*\(([^)]+)\)'

        self.step_entities['closed_shells'] = []
        self.step_entities['open_shells'] = []

        # Find closed shells (core)
        for match in re.finditer(closed_shell_pattern, content):
            shell_id = int(match.group(1))
            # Parse face references, handling spaces and # symbols
            face_refs_str = match.group(2).replace(' ', '')  # Remove spaces
            face_refs = [int(x.strip('#')) for x in face_refs_str.split(',') if x.strip()]
            self.step_entities['closed_shells'].append((shell_id, face_refs))
            print(f"  Found CLOSED_SHELL #{shell_id} with {len(face_refs)} faces (core)")

        # Find open shells (windings)
        for match in re.finditer(open_shell_pattern, content):
            shell_id = int(match.group(1))
            # Parse face references, handling spaces and # symbols
            face_refs_str = match.group(2).replace(' ', '')  # Remove spaces
            face_refs = [int(x.strip('#')) for x in face_refs_str.split(',') if x.strip()]
            self.step_entities['open_shells'].append((shell_id, face_refs))
            print(f"  Found OPEN_SHELL #{shell_id} with {len(face_refs)} faces (winding)")

        # Find layer assignments for open shells
        layer_pattern = r'PRESENTATION_LAYER_ASSIGNMENT\s*\(\s*\'(Winding\d+)\'\s*,\s*\'[^\']*\'\s*,\s*\(\s*#(\d+)\s*\)'
        self.step_entities['layer_assignments'] = {}

        for match in re.finditer(layer_pattern, content):
            winding_name = match.group(1)
            styled_item_id = int(match.group(2))
            self.step_entities['layer_assignments'][styled_item_id] = winding_name
            print(f"  Layer assignment: {winding_name} -> STYLED_ITEM #{styled_item_id}")

    def extract_winding_layers_from_step_text(self) -> dict:
        """Parse STEP file text to find winding layer assignments."""
        windings = {}

        with open(self.filepath, 'r') as f:
            content = f.read()

        # Find PRESENTATION_LAYER_ASSIGNMENT entities
        pattern = r'#(\d+)\s*=\s*PRESENTATION_LAYER_ASSIGNMENT\s*\(\s*\'(Winding\d+)\'[^;]+;'
        matches = re.findall(pattern, content)

        for entity_id, winding_name in matches:
            if winding_name not in windings:
                windings[winding_name] = []
            windings[winding_name].append(int(entity_id))

        print(f"  Found {len(windings)} winding layers: {list(windings.keys())}")
        return windings

    def tessellate_shape(self, shell_type='all') -> Tuple[np.ndarray, np.ndarray]:
        """
        Convert CADQuery shape to triangular mesh.

        Args:
            shell_type: 'all' (default), 'closed' (core only), or 'open' (windings only)
        """
        try:
            import tempfile
            import trimesh
            import os

            with tempfile.NamedTemporaryFile(suffix='.stl', delete=False) as tmp:
                tmp_name = tmp.name

            cq.exporters.export(self.shape, tmp_name)
            mesh = trimesh.load(tmp_name)
            os.unlink(tmp_name)

            # Convert mm to m
            mesh.vertices *= 0.001

            return mesh.vertices, mesh.faces

        except Exception as e:
            print(f"  Warning: Tessellation error: {e}")
            return np.array([]), np.array([])

    def separate_components(self) -> dict:
        """
        Separate mesh into connected components (core and windings).

        Returns dict with 'core' and 'winding' meshes.
        """
        try:
            import tempfile
            import trimesh
            import os

            with tempfile.NamedTemporaryFile(suffix='.stl', delete=False) as tmp:
                tmp_name = tmp.name

            cq.exporters.export(self.shape, tmp_name)
            mesh = trimesh.load(tmp_name)
            os.unlink(tmp_name)

            # Convert mm to m
            mesh.vertices *= 0.001

            # Split into connected components
            components = mesh.split(only_watertight=False)
            print(f"\n  Mesh split into {len(components)} connected components:")

            result = {}

            for i, comp in enumerate(components):
                n_faces = len(comp.faces)
                is_watertight = comp.is_watertight
                volume = comp.volume if is_watertight else 0

                print(f"    Component {i}: {n_faces} faces, "
                      f"watertight={is_watertight}, volume={volume:.6f} m³")

                # Heuristic: largest watertight component is the core
                # Smaller non-watertight components are winding surfaces
                if is_watertight and volume > 0:
                    if 'core' not in result or volume > result['core'].volume:
                        result['core'] = comp
                else:
                    # Winding surfaces (open shells)
                    if 'winding' not in result:
                        result['winding'] = comp
                    else:
                        # Merge winding components
                        result['winding'] = trimesh.util.concatenate([result['winding'], comp])

            # Ensure we have both components
            if 'core' not in result:
                print("    WARNING: No watertight core found!")
                # Use largest component as core
                result['core'] = max(components, key=lambda c: len(c.faces))

            if 'winding' not in result:
                print("    WARNING: No winding surfaces found!")
                # Use smallest component as winding
                result['winding'] = min(components, key=lambda c: len(c.faces))

            return result

        except Exception as e:
            print(f"  Error separating components: {e}")
            return {}


class MagnetostaticSolver:
    """Solve magnetostatic field from current sheets and cores."""

    def __init__(self, cores: List[MagneticCore], windings: List[CurrentSheet]):
        self.cores = cores
        self.windings = windings

    def is_in_core(self, point: np.ndarray) -> bool:
        """Check if point is inside any core."""
        for core in self.cores:
            if core.contains_point(point):
                return True
        return False

    def biot_savart_sheet(self, point: np.ndarray, sheet: CurrentSheet) -> np.ndarray:
        """
        Compute B-field from current sheet using Biot-Savart law.

        Treats sheet as collection of triangular current elements.
        """
        B = np.zeros(3)

        for face_idx, face in enumerate(sheet.faces):
            # Triangle vertices
            v0, v1, v2 = sheet.vertices[face]

            # Triangle properties
            center = (v0 + v1 + v2) / 3
            edge1 = v1 - v0
            edge2 = v2 - v0
            normal = np.cross(edge1, edge2)
            area = 0.5 * np.linalg.norm(normal)
            normal = normal / (2 * area) if area > 0 else np.zeros(3)

            # Current density (uniform over sheet)
            J_magnitude = sheet.current / sheet.area
            J = J_magnitude * normal

            # Distance vector
            r_vec = point - center
            r = np.linalg.norm(r_vec)

            if r < 1e-9:  # Avoid singularity
                continue

            # Biot-Savart: dB = (μ₀/4π) * (J × r) / r³ * dA
            dB = (MU_0 / (4 * np.pi)) * np.cross(J, r_vec) / (r**3) * area
            B += dB

        return B

    def calculate_field(self, point: np.ndarray) -> np.ndarray:
        """Calculate total B-field at point."""
        # Don't calculate inside cores
        if self.is_in_core(point):
            return np.array([np.nan, np.nan, np.nan])

        B_total = np.zeros(3)

        # Sum contributions from all windings
        for winding in self.windings:
            B_winding = self.biot_savart_sheet(point, winding)

            # Apply permeability enhancement (simplified model)
            # In presence of high-μ cores, field is enhanced
            mu_factor = 1.0
            if len(self.cores) > 0:
                # Rough approximation: fields near cores are amplified
                mu_factor = 10.0  # Empirical factor

            B_total += B_winding * mu_factor

        return B_total

    def calculate_grid(self, bounds: np.ndarray, resolution: int = 15):
        """Calculate field on 3D grid."""
        x = np.linspace(bounds[0, 0], bounds[0, 1], resolution)
        y = np.linspace(bounds[1, 0], bounds[1, 1], resolution)
        z = np.linspace(bounds[2, 0], bounds[2, 1], resolution)

        points = []
        for xi in x:
            for yi in y:
                for zi in z:
                    points.append([xi, yi, zi])

        points = np.array(points)
        B_vectors = np.zeros_like(points)

        # Pre-filter points that are inside cores (vectorized)
        inside_mask = np.zeros(len(points), dtype=bool)
        for core in self.cores:
            core_inside = core.contains_points(points)
            print(f"  Core detection: {core_inside.sum()}/{len(points)} points inside core")
            inside_mask |= core_inside

        points_outside = ~inside_mask
        n_outside = points_outside.sum()

        print(f"\nCalculating B-field at {n_outside}/{len(points)} points (excluding core interior)...")

        # Only calculate field at points outside cores
        calc_idx = 0
        for i in range(len(points)):
            if inside_mask[i]:
                B_vectors[i] = np.array([np.nan, np.nan, np.nan])
            else:
                if calc_idx % 100 == 0:
                    print(f"  Progress: {calc_idx}/{n_outside} ({100*calc_idx/n_outside:.1f}%)")
                B_vectors[i] = self.calculate_field(points[i])
                calc_idx += 1

        return points, B_vectors


class FieldVisualizer:
    """Interactive 3D visualization with PyVista."""

    def __init__(self):
        self.plotter = pv.Plotter()
        self.plotter.set_background('white')

    def add_core(self, core: MagneticCore):
        """Add core geometry."""
        faces_pv = np.hstack([[3] + list(f) for f in core.faces])
        mesh = pv.PolyData(core.vertices, faces_pv)
        print(f"  Adding core: {len(core.vertices)} vertices, {len(core.faces)} faces")
        self.plotter.add_mesh(mesh, color='gray', opacity=0.6,
                            label=f'Core (μr={core.mu_r})', show_edges=True)

    def add_winding(self, winding: CurrentSheet):
        """Add winding geometry."""
        faces_pv = np.hstack([[3] + list(f) for f in winding.faces])
        mesh = pv.PolyData(winding.vertices, faces_pv)
        self.plotter.add_mesh(mesh, color='orange', opacity=0.7,
                            label=f'{winding.name} ({winding.current}A)')

    def add_field_magnitude(self, points: np.ndarray, B_vectors: np.ndarray):
        """Add field vectors with arrows colored by magnitude."""
        B_mag = np.linalg.norm(B_vectors, axis=1)

        # Filter out NaN (inside cores)
        valid = ~np.isnan(B_mag)
        valid_points = points[valid]
        valid_vectors = B_vectors[valid]
        valid_mag = B_mag[valid]

        if len(valid_points) == 0:
            print("Warning: No valid field points")
            return

        # Normalize vectors for consistent arrow size
        # Scale by log of magnitude to handle large range
        scale_factor = 0.0001  # meters (0.1mm arrows)
        normalized = valid_vectors / (valid_mag[:, np.newaxis] + 1e-10)
        scaled_vectors = normalized * scale_factor * np.log10(valid_mag[:, np.newaxis] / valid_mag.min() + 1)

        # Create arrow glyphs
        cloud = pv.PolyData(valid_points)
        cloud['vectors'] = scaled_vectors
        cloud['|B|'] = valid_mag

        # Add arrows
        arrows = cloud.glyph(orient='vectors', scale=False, factor=1.0)
        self.plotter.add_mesh(arrows, scalars='|B|', cmap='viridis',
                            scalar_bar_args={'title': '|B| (Tesla)'},
                            label='B-field')

    def add_streamlines(self, points: np.ndarray, B_vectors: np.ndarray,
                       n_lines: int = 30):
        """Add magnetic field streamlines."""
        valid = ~np.isnan(B_vectors[:, 0])
        valid_points = points[valid]
        valid_vectors = B_vectors[valid]

        if len(valid_points) < 10:
            print("Warning: Not enough points for streamlines")
            return

        # Create unstructured grid
        cloud = pv.PolyData(valid_points)
        cloud['B'] = valid_vectors

        # Generate seed points for streamlines
        bounds = valid_points.min(axis=0), valid_points.max(axis=0)
        seed_points = []
        for _ in range(n_lines):
            pt = np.random.uniform(bounds[0], bounds[1])
            seed_points.append(pt)
        seed_points = np.array(seed_points)

        # Add seeds
        seeds = pv.PolyData(seed_points)

        # Create streamlines
        try:
            streamlines = cloud.streamlines_from_source(seeds, vectors='B',
                                                       max_steps=500,
                                                       integration_direction='both')
            self.plotter.add_mesh(streamlines, color='black', line_width=2,
                                label='Field lines')
        except Exception as e:
            print(f"Warning: Could not generate streamlines: {e}")

    def show(self, save_html: str = None):
        """Show interactive visualization."""
        self.plotter.add_legend()
        self.plotter.show_axes()
        self.plotter.enable_trackball_style()

        if save_html:
            # Save as interactive HTML file
            print(f"\nSaving interactive visualization to: {save_html}")
            print("  Open this file in a web browser for full 3D interaction")
            print("  Controls: Left-click drag=Rotate, Right-click drag=Pan, Scroll=Zoom")
            self.plotter.export_html(save_html)
        else:
            # Try to show interactively
            print("\nOpening interactive 3D viewer...")
            print("  Left mouse: Rotate")
            print("  Middle mouse: Pan")
            print("  Right mouse / scroll: Zoom")
            try:
                self.plotter.show()
            except Exception as e:
                print(f"  Could not open interactive window: {e}")
                print("  Run with save_html parameter to export to browser-viewable file")


def main():
    import sys
    print("="*70, flush=True)
    print("MAGVID Electromagnetic Field Solver", flush=True)
    print("="*70, flush=True)
    sys.stdout.flush()

    # Configuration
    step_file = "out.STEP"
    winding_currents = {"Winding0": -500.0}
    core_mu_r = 5000.0
    field_resolution = 8  # Grid resolution (8³ = 512 points)

    print(f"\nConfiguration:", flush=True)
    print(f"  STEP file: {step_file}", flush=True)
    print(f"  Core μr: {core_mu_r}", flush=True)
    print(f"  Winding currents: {winding_currents}", flush=True)
    print(f"  Grid resolution: {field_resolution}³ points", flush=True)
    sys.stdout.flush()

    # Load STEP file
    print("\n" + "-"*70, flush=True)
    sys.stdout.flush()
    loader = STEPLoader(step_file)
    if not loader.load():
        print("ERROR: Failed to load STEP file", flush=True)
        return

    # Extract and separate geometry into core and windings
    print("\nSeparating geometry into components...", flush=True)
    sys.stdout.flush()
    components = loader.separate_components()

    if not components:
        print("ERROR: Failed to separate geometry", flush=True)
        return

    # Extract core mesh
    print("Extracting core mesh...", flush=True)
    sys.stdout.flush()
    core_mesh = components['core']
    core_vertices = core_mesh.vertices
    core_faces = core_mesh.faces

    print(f"\nCore geometry: {len(core_vertices)} vertices, {len(core_faces)} faces", flush=True)
    sys.stdout.flush()

    # Create core with FULL mesh (needed for accurate point-in-core detection)
    # IMPORTANT: Don't subsample core faces - ray-casting needs watertight mesh!
    core = MagneticCore(vertices=core_vertices, faces=core_faces, mu_r=core_mu_r)

    # Create subsampled version for visualization only
    core_faces_viz = core_faces
    if len(core_faces) > 500:
        print(f"  Creating subsampled version for visualization...")
        target_faces = min(300, len(core_faces))
        step = max(1, len(core_faces) // target_faces)
        core_faces_viz = core_faces[::step]
        print(f"  Visualization mesh: {len(core_vertices)} vertices, {len(core_faces_viz)} faces")

    # Extract winding mesh to get bounds
    winding_mesh = components['winding']
    winding_bounds = winding_mesh.bounds

    print(f"\nWinding bounds:")
    print(f"  X: [{winding_bounds[0,0]:.4f}, {winding_bounds[1,0]:.4f}]")
    print(f"  Y: [{winding_bounds[0,1]:.4f}, {winding_bounds[1,1]:.4f}]")
    print(f"  Z: [{winding_bounds[0,2]:.4f}, {winding_bounds[1,2]:.4f}]")

    # Extract core surface in winding region
    print("\nExtracting core surface in winding region...")
    winding_vertices, winding_faces = core.extract_surface_region(
        winding_bounds=winding_bounds,
        margin=0.002  # 2mm margin
    )

    # Subsample winding faces for speed (field calc is O(n_points × n_faces))
    if len(winding_faces) > 100:
        print(f"  Subsampling winding surface for practical computation time...", flush=True)
        sys.stdout.flush()
        # Target ~50-100 faces for reasonable speed
        target_faces = min(100, len(winding_faces))
        step = max(1, len(winding_faces) // target_faces)
        winding_faces = winding_faces[::step]
        print(f"  Subsampled to: {len(winding_faces)} faces", flush=True)
        sys.stdout.flush()

    # Create winding current sheets using core surface
    print("\nCreating winding current sheets from core surface...", flush=True)
    sys.stdout.flush()
    windings = []
    for winding_name, current in winding_currents.items():
        print(f"  Creating {winding_name}...", flush=True)
        sys.stdout.flush()
        winding = CurrentSheet(
            name=winding_name,
            vertices=winding_vertices,
            faces=winding_faces,
            current=current
        )
        windings.append(winding)
        print(f"  {winding_name}: {current}A, area={winding.area:.6f} m²", flush=True)
        sys.stdout.flush()

    # Calculate bounds for field grid (use core bounds with small margin)
    print("\nCalculating field grid bounds...", flush=True)
    sys.stdout.flush()
    mins = core_vertices.min(axis=0)
    maxs = core_vertices.max(axis=0)
    margin = 0.015  # 1.5cm margin for tighter visualization focus
    bounds = np.array([[mins[0]-margin, maxs[0]+margin],
                      [mins[1]-margin, maxs[1]+margin],
                      [mins[2]-margin, maxs[2]+margin]])

    print(f"\nField calculation bounds:", flush=True)
    print(f"  X: [{bounds[0,0]:.3f}, {bounds[0,1]:.3f}] m", flush=True)
    print(f"  Y: [{bounds[1,0]:.3f}, {bounds[1,1]:.3f}] m", flush=True)
    print(f"  Z: [{bounds[2,0]:.3f}, {bounds[2,1]:.3f}] m", flush=True)
    sys.stdout.flush()

    # Solve fields
    print("\n" + "-"*70, flush=True)
    print("Creating solver...", flush=True)
    sys.stdout.flush()
    solver = MagnetostaticSolver(cores=[core], windings=windings)
    print("Calculating field grid...", flush=True)
    sys.stdout.flush()
    points, B_vectors = solver.calculate_grid(bounds, resolution=field_resolution)

    # Statistics
    B_mag = np.linalg.norm(B_vectors, axis=1)
    valid_B = B_mag[~np.isnan(B_mag)]
    if len(valid_B) > 0:
        print(f"\nField statistics:")
        print(f"  |B| max: {valid_B.max():.6e} T")
        print(f"  |B| mean: {valid_B.mean():.6e} T")
        print(f"  |B| min: {valid_B.min():.6e} T")

    # Visualize
    print("\n" + "-"*70)
    print("Creating visualization...")

    # Create visualization-only core with subsampled mesh
    core_viz = MagneticCore(vertices=core_vertices, faces=core_faces_viz, mu_r=core_mu_r)

    viz = FieldVisualizer()
    viz.add_core(core_viz)
    for winding in windings:
        viz.add_winding(winding)
    viz.add_field_magnitude(points, B_vectors)
    viz.add_streamlines(points, B_vectors, n_lines=25)

    # Save as interactive HTML (works in headless mode)
    viz.show(save_html="magvid_field_visualization.html")


if __name__ == "__main__":
    main()
