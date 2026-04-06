#!/usr/bin/env python3
"""Generate VTK C++ reference data — batch 3."""
import json, os, tempfile
import vtk

OUT = "tests/vtk_validation/reference"

def save(name, data):
    with open(f"{OUT}/{name}.json", "w") as f:
        json.dump(data, f, indent=2)
    print(f"  {name}")

def sphere(res=32):
    s = vtk.vtkSphereSource()
    s.SetThetaResolution(res); s.SetPhiResolution(res); s.Update()
    return s.GetOutput()

def tri_sphere(res=32):
    t = vtk.vtkTriangleFilter(); t.SetInputData(sphere(res)); t.Update()
    return t.GetOutput()

def basic(pd):
    return {"num_points": pd.GetNumberOfPoints(), "num_cells": pd.GetNumberOfCells()}

# --- Calculator ---
print("calculator...")
pd = sphere(32)
calc = vtk.vtkArrayCalculator(); calc.SetInputData(pd)
calc.AddCoordinateScalarVariable("coordsX", 0)
calc.AddCoordinateScalarVariable("coordsY", 1)
calc.AddCoordinateScalarVariable("coordsZ", 2)
calc.SetFunction("coordsX*coordsX + coordsY*coordsY + coordsZ*coordsZ")
calc.SetResultArrayName("Result"); calc.Update()
r = calc.GetOutput()
save("filter_calculator", {"num_points": r.GetNumberOfPoints(), "has_result": r.GetPointData().GetArray("Result") is not None})

# --- Cell to point data ---
print("cell_to_point_data...")
pd = sphere(32)
elev = vtk.vtkElevationFilter(); elev.SetInputData(pd)
elev.SetLowPoint(0,0,-0.5); elev.SetHighPoint(0,0,0.5); elev.Update()
p2c = vtk.vtkPointDataToCellData(); p2c.SetInputData(elev.GetOutput()); p2c.Update()
c2p = vtk.vtkCellDataToPointData(); c2p.SetInputData(p2c.GetOutput()); c2p.Update()
save("filter_cell_to_point_data", basic(c2p.GetOutput()))

# --- Point to cell data ---
print("point_to_cell_data...")
save("filter_point_to_cell_data", basic(p2c.GetOutput()))

# --- Clip scalar ---
print("clip_scalar...")
pd_elev = elev.GetOutput()
clip = vtk.vtkClipPolyData(); clip.SetInputData(pd_elev)
clip.SetValue(0.5); clip.Update()
save("filter_clip_scalar", basic(clip.GetOutput()))

# --- Curvatures large ---
print("curvatures_large...")
pd = tri_sphere(128)
curv = vtk.vtkCurvatures(); curv.SetInputData(pd); curv.SetCurvatureTypeToMean(); curv.Update()
save("filter_curvatures_large", basic(curv.GetOutput()))

# --- Dihedral angles ---
print("dihedral_angles...")
# No direct VTK filter, store basic info
pd = tri_sphere(32)
save("filter_dihedral_angles", {"num_cells": pd.GetNumberOfCells()})

# --- Distance to origin ---
print("distance_to_origin...")
save("filter_distance_to_origin", {"num_points": sphere(32).GetNumberOfPoints()})

# --- Extract cells half ---
print("extract_cells_half...")
pd = tri_sphere(32)
nc = pd.GetNumberOfCells()
ids = vtk.vtkIdTypeArray()
for i in range(0, nc, 2):
    ids.InsertNextValue(i)
sel = vtk.vtkSelectionNode(); sel.SetFieldType(vtk.vtkSelectionNode.CELL); sel.SetContentType(vtk.vtkSelectionNode.INDICES)
sel.SetSelectionList(ids)
selection = vtk.vtkSelection(); selection.AddNode(sel)
ext = vtk.vtkExtractSelection(); ext.SetInputData(0, pd); ext.SetInputData(1, selection); ext.Update()
surf = vtk.vtkDataSetSurfaceFilter(); surf.SetInputData(ext.GetOutput()); surf.Update()
save("filter_extract_cells_half", basic(surf.GetOutput()))

# --- Hedgehog ---
print("hedgehog...")
pd = sphere(32)
nf = vtk.vtkPolyDataNormals(); nf.SetInputData(pd); nf.Update()
hh = vtk.vtkHedgeHog(); hh.SetInputData(nf.GetOutput()); hh.SetScaleFactor(0.1); hh.Update()
save("filter_hedgehog", {"num_points": hh.GetOutput().GetNumberOfPoints(), "num_lines": hh.GetOutput().GetNumberOfLines()})

# --- Offset surface ---
print("offset_surface...")
pd = tri_sphere(32)
save("filter_offset_surface", {"num_points": pd.GetNumberOfPoints()})

# --- Plane 32 ---
print("plane_32...")
p = vtk.vtkPlaneSource(); p.SetResolution(32,32); p.Update()
save("source_plane_32x32", basic(p.GetOutput()))

# --- Poly data distance ---
print("poly_data_distance...")
s1 = tri_sphere(32)
s2 = vtk.vtkSphereSource(); s2.SetCenter(0.1,0,0); s2.SetThetaResolution(32); s2.SetPhiResolution(32); s2.Update()
t2 = vtk.vtkTriangleFilter(); t2.SetInputData(s2.GetOutput()); t2.Update()
dd = vtk.vtkDistancePolyDataFilter(); dd.SetInputData(0, s1); dd.SetInputData(1, t2.GetOutput()); dd.Update()
save("filter_poly_data_distance", basic(dd.GetOutput()))

# --- Reflect ---
print("reflect...")
pd = tri_sphere(32)
ref = vtk.vtkReflectionFilter(); ref.SetInputData(pd); ref.SetPlaneToX(); ref.CopyInputOn(); ref.Update()
r = ref.GetOutput()
save("filter_reflect", {"num_points": r.GetNumberOfPoints(), "num_cells": r.GetNumberOfCells()})

# --- Ribbon ---
print("ribbon...")
line = vtk.vtkLineSource(); line.SetPoint1(0,0,0); line.SetPoint2(1,0,0); line.SetResolution(20); line.Update()
rib = vtk.vtkRibbonFilter(); rib.SetInputData(line.GetOutput()); rib.SetWidth(0.1); rib.Update()
save("filter_ribbon", basic(rib.GetOutput()))

# --- Ruled surface ---
print("ruled_surface...")
l1 = vtk.vtkLineSource(); l1.SetPoint1(0,0,0); l1.SetPoint2(1,0,0); l1.SetResolution(10); l1.Update()
l2 = vtk.vtkLineSource(); l2.SetPoint1(0,1,0); l2.SetPoint2(1,1,0.5); l2.SetResolution(10); l2.Update()
app = vtk.vtkAppendPolyData(); app.AddInputData(l1.GetOutput()); app.AddInputData(l2.GetOutput()); app.Update()
rs = vtk.vtkRuledSurfaceFilter(); rs.SetInputData(app.GetOutput()); rs.Update()
save("filter_ruled_surface", basic(rs.GetOutput()))

# --- Silhouette ---
print("silhouette...")
pd = tri_sphere(32)
cam = vtk.vtkCamera(); cam.SetPosition(0,0,5)
sil = vtk.vtkPolyDataSilhouette(); sil.SetInputData(pd); sil.SetCamera(cam); sil.Update()
save("filter_silhouette", {"num_lines": sil.GetOutput().GetNumberOfLines()})

# --- Sphere 64 ---
print("sphere_64...")
s = vtk.vtkSphereSource(); s.SetThetaResolution(64); s.SetPhiResolution(64); s.Update()
save("source_sphere_64x64", basic(s.GetOutput()))

# --- Topology ---
print("topology_analysis...")
pd = tri_sphere(32)
fe = vtk.vtkFeatureEdges(); fe.SetInputData(pd)
fe.BoundaryEdgesOn(); fe.FeatureEdgesOff(); fe.ManifoldEdgesOff(); fe.NonManifoldEdgesOff(); fe.Update()
save("filter_topology_analysis", {"num_boundary_edges": fe.GetOutput().GetNumberOfLines(), "is_closed": fe.GetOutput().GetNumberOfLines() == 0})

# --- VTK legacy roundtrip ---
print("vtk_roundtrip...")
pd = tri_sphere(32)
path = os.path.join(tempfile.gettempdir(), "test_rt.vtk")
w = vtk.vtkPolyDataWriter(); w.SetFileName(path); w.SetInputData(pd); w.Write()
rd = vtk.vtkPolyDataReader(); rd.SetFileName(path); rd.Update()
save("io_vtk_roundtrip", basic(rd.GetOutput()))

# --- VTK large ---
print("vtk_large...")
pd = tri_sphere(128)
path = os.path.join(tempfile.gettempdir(), "test_large.vtk")
w = vtk.vtkPolyDataWriter(); w.SetFileName(path); w.SetInputData(pd); w.Write()
rd = vtk.vtkPolyDataReader(); rd.SetFileName(path); rd.Update()
save("io_vtk_large", basic(rd.GetOutput()))

# --- Transform large ---
print("transform_large...")
pd = sphere(128)
t = vtk.vtkTransformPolyDataFilter()
tr = vtk.vtkTransform(); tr.RotateZ(45); tr.Scale(2,2,2)
t.SetTransform(tr); t.SetInputData(pd); t.Update()
save("filter_transform_large", basic(t.GetOutput()))

print(f"Done.")
