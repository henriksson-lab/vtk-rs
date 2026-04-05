#!/usr/bin/env python3
"""Generate VTK C++ reference data for validation tests.

Requires: pip install vtk

Outputs JSON metadata files to tests/vtk_validation/reference/
Each file contains point counts, cell counts, bounds, scalar ranges, etc.
that Rust tests compare against.
"""

import vtk
import json
import os

REF_DIR = os.path.join(os.path.dirname(__file__), "reference")
os.makedirs(REF_DIR, exist_ok=True)

def save(name, data):
    path = os.path.join(REF_DIR, f"{name}.json")
    with open(path, "w") as f:
        json.dump(data, f, indent=2)
    print(f"  {name}: {data.get('num_points', '?')} pts, {data.get('num_cells', '?')} cells")

def mesh_info(pd):
    b = pd.GetBounds()
    info = {
        "num_points": pd.GetNumberOfPoints(),
        "num_cells": pd.GetNumberOfCells(),
        "bounds": [b[0], b[1], b[2], b[3], b[4], b[5]],
    }
    if pd.GetPointData().GetNormals():
        info["has_normals"] = True
        info["num_normal_tuples"] = pd.GetPointData().GetNormals().GetNumberOfTuples()
    if pd.GetPointData().GetScalars():
        r = pd.GetPointData().GetScalars().GetRange()
        info["scalar_range"] = [r[0], r[1]]
        info["scalar_name"] = pd.GetPointData().GetScalars().GetName()
    return info

def make_sphere(theta=32, phi=32):
    s = vtk.vtkSphereSource()
    s.SetThetaResolution(theta)
    s.SetPhiResolution(phi)
    s.SetRadius(0.5)
    s.Update()
    return s.GetOutput()

# ========== SOURCES ==========
print("=== Sources ===")

# Sphere
pd = make_sphere()
save("source_sphere_32x32", mesh_info(pd))

# Cube
s = vtk.vtkCubeSource(); s.Update()
save("source_cube", mesh_info(s.GetOutput()))

# Cone
s = vtk.vtkConeSource(); s.SetResolution(32); s.Update()
save("source_cone_32", mesh_info(s.GetOutput()))

# Cylinder
s = vtk.vtkCylinderSource(); s.SetResolution(32); s.Update()
save("source_cylinder_32", mesh_info(s.GetOutput()))

# Plane
s = vtk.vtkPlaneSource(); s.SetXResolution(10); s.SetYResolution(10); s.Update()
save("source_plane_10x10", mesh_info(s.GetOutput()))

# Arrow
s = vtk.vtkArrowSource(); s.Update()
save("source_arrow", mesh_info(s.GetOutput()))

# ========== FILTERS ==========
print("\n=== Filters ===")
sphere = make_sphere()

# Normals
f = vtk.vtkPolyDataNormals(); f.SetInputData(sphere); f.ComputePointNormalsOn(); f.Update()
save("filter_normals", mesh_info(f.GetOutput()))

# Triangulate
f = vtk.vtkTriangleFilter(); f.SetInputData(sphere); f.Update()
info = mesh_info(f.GetOutput())
info["all_triangles"] = True
save("filter_triangulate", info)

# Clean
duped = vtk.vtkAppendPolyData()
duped.AddInputData(sphere); duped.AddInputData(sphere); duped.Update()
f = vtk.vtkCleanPolyData(); f.SetInputData(duped.GetOutput()); f.Update()
info = mesh_info(f.GetOutput())
info["input_points"] = duped.GetOutput().GetNumberOfPoints()
save("filter_clean", info)

# Decimate
tri = vtk.vtkTriangleFilter(); tri.SetInputData(sphere); tri.Update()
f = vtk.vtkDecimatePro(); f.SetInputData(tri.GetOutput()); f.SetTargetReduction(0.5); f.Update()
info = mesh_info(f.GetOutput())
info["target_reduction"] = 0.5
info["input_cells"] = tri.GetOutput().GetNumberOfCells()
save("filter_decimate_50", info)

# Smooth
f = vtk.vtkSmoothPolyDataFilter(); f.SetInputData(sphere)
f.SetNumberOfIterations(20); f.SetRelaxationFactor(0.5); f.Update()
info = mesh_info(f.GetOutput())
info["iterations"] = 20
save("filter_smooth_20", info)

# Connectivity
# Two separate spheres
s2 = vtk.vtkSphereSource(); s2.SetCenter(3,0,0); s2.SetThetaResolution(8); s2.SetPhiResolution(8); s2.Update()
s1 = vtk.vtkSphereSource(); s1.SetThetaResolution(8); s1.SetPhiResolution(8); s1.Update()
app = vtk.vtkAppendPolyData(); app.AddInputData(s1.GetOutput()); app.AddInputData(s2.GetOutput()); app.Update()
f = vtk.vtkConnectivityFilter(); f.SetInputData(app.GetOutput()); f.SetExtractionModeToAllRegions()
f.ColorRegionsOn(); f.Update()
info = mesh_info(f.GetOutput())
info["num_regions"] = f.GetNumberOfExtractedRegions()
save("filter_connectivity", info)

# Feature edges
f = vtk.vtkFeatureEdges(); f.SetInputData(sphere)
f.BoundaryEdgesOn(); f.FeatureEdgesOff(); f.ManifoldEdgesOff(); f.NonManifoldEdgesOff()
f.Update()
info = mesh_info(f.GetOutput())
info["boundary_edges"] = f.GetOutput().GetNumberOfCells()
save("filter_feature_edges", info)

# Clip by plane
plane = vtk.vtkPlane(); plane.SetOrigin(0,0,0); plane.SetNormal(1,0,0)
f = vtk.vtkClipPolyData(); f.SetInputData(tri.GetOutput()); f.SetClipFunction(plane); f.Update()
info = mesh_info(f.GetOutput())
info["clip_normal"] = [1,0,0]
save("filter_clip_plane", info)

# Elevation
f = vtk.vtkElevationFilter(); f.SetInputData(sphere)
f.SetLowPoint(0,0,-0.5); f.SetHighPoint(0,0,0.5); f.Update()
info = mesh_info(f.GetOutput())
save("filter_elevation", info)

# Threshold (on elevation)
elev = f.GetOutput()
t = vtk.vtkThreshold(); t.SetInputData(elev)
t.SetLowerThreshold(0.3); t.SetUpperThreshold(0.7)
t.SetThresholdFunction(vtk.vtkThreshold.THRESHOLD_BETWEEN)
t.SetInputArrayToProcess(0, 0, 0, 0, "Elevation")
t.Update()
info = {"num_cells": t.GetOutput().GetNumberOfCells(), "num_points": t.GetOutput().GetNumberOfPoints()}
info["threshold_range"] = [0.3, 0.7]
save("filter_threshold", info)

# Marching cubes
img = vtk.vtkImageData(); img.SetDimensions(64, 64, 64); img.SetSpacing(1,1,1)
img.AllocateScalars(vtk.VTK_DOUBLE, 1)
arr = img.GetPointData().GetScalars()
for k in range(64):
    for j in range(64):
        for i in range(64):
            arr.SetValue(k*64*64+j*64+i, (i-32)**2 + (j-32)**2 + (k-32)**2)
mc = vtk.vtkMarchingCubes(); mc.SetInputData(img); mc.SetValue(0, 400); mc.Update()
info = mesh_info(mc.GetOutput())
info["isovalue"] = 400
save("filter_marching_cubes_64", info)

# Cell centers
f = vtk.vtkCellCenters(); f.SetInputData(sphere); f.Update()
info = {"num_points": f.GetOutput().GetNumberOfPoints()}
save("filter_cell_centers", info)

# Mass properties
tri2 = vtk.vtkTriangleFilter(); tri2.SetInputData(sphere); tri2.Update()
f = vtk.vtkMassProperties(); f.SetInputData(tri2.GetOutput())
info = {"surface_area": f.GetSurfaceArea(), "volume": f.GetVolume()}
save("filter_mass_properties", info)

# Curvatures
f = vtk.vtkCurvatures(); f.SetInputData(tri2.GetOutput()); f.SetCurvatureTypeToMean(); f.Update()
info = mesh_info(f.GetOutput())
info["curvature_type"] = "mean"
save("filter_curvatures_mean", info)

# Shrink
f = vtk.vtkShrinkFilter(); f.SetInputData(sphere); f.SetShrinkFactor(0.5); f.Update()
info = mesh_info(f.GetOutput())
info["shrink_factor"] = 0.5
save("filter_shrink", info)

# Reverse sense
f = vtk.vtkReverseSense(); f.SetInputData(tri2.GetOutput()); f.ReverseNormalsOn(); f.ReverseCellsOn(); f.Update()
info = mesh_info(f.GetOutput())
info["reversed"] = True
save("filter_reverse_sense", info)

# Extract edges
f = vtk.vtkExtractEdges(); f.SetInputData(sphere); f.Update()
info = {"num_edges": f.GetOutput().GetNumberOfCells(), "num_points": f.GetOutput().GetNumberOfPoints()}
save("filter_extract_edges", info)

# Delaunay 2D
pts = vtk.vtkPoints()
import random; random.seed(42)
for _ in range(100):
    pts.InsertNextPoint(random.uniform(-1,1), random.uniform(-1,1), 0)
pd2d = vtk.vtkPolyData(); pd2d.SetPoints(pts)
f = vtk.vtkDelaunay2D(); f.SetInputData(pd2d); f.Update()
info = mesh_info(f.GetOutput())
save("filter_delaunay_2d", info)

# ========== I/O ==========
print("\n=== I/O ===")
import tempfile

for fmt, writer_cls, reader_cls, ext in [
    ("vtk_legacy", vtk.vtkPolyDataWriter, vtk.vtkPolyDataReader, ".vtk"),
    ("stl", vtk.vtkSTLWriter, vtk.vtkSTLReader, ".stl"),
    ("obj", vtk.vtkOBJWriter, vtk.vtkOBJReader, ".obj"),
    ("ply", vtk.vtkPLYWriter, vtk.vtkPLYReader, ".ply"),
]:
    path = os.path.join(tempfile.gettempdir(), f"vtk_ref{ext}")
    w = writer_cls(); w.SetFileName(path); w.SetInputData(tri2.GetOutput()); w.Write()
    r = reader_cls(); r.SetFileName(path); r.Update()
    out = r.GetOutput()
    info = mesh_info(out)
    info["format"] = fmt
    save(f"io_{fmt}", info)

print(f"\nDone. {len(os.listdir(REF_DIR))} reference files in {REF_DIR}")
