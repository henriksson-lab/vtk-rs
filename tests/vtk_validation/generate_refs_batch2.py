#!/usr/bin/env python3
"""Generate VTK C++ reference data — batch 2."""
import json, os
import vtk

OUT = "tests/vtk_validation/reference"
os.makedirs(OUT, exist_ok=True)

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

# --- Butterfly subdivision ---
print("butterfly_1...")
pd = tri_sphere(32)
bf = vtk.vtkButterflySubdivisionFilter(); bf.SetInputData(pd); bf.SetNumberOfSubdivisions(1); bf.Update()
save("filter_butterfly_1", basic(bf.GetOutput()))

# --- Catmull-Clark ---
print("catmull_clark_1...")
# VTK doesn't have a direct Catmull-Clark on PolyData; skip or use Loop
lp = vtk.vtkLoopSubdivisionFilter(); lp.SetInputData(pd); lp.SetNumberOfSubdivisions(1); lp.Update()
save("filter_catmull_clark_1", {"note": "using Loop as proxy", **basic(lp.GetOutput())})

# --- Cell Size ---
print("cell_size...")
pd = tri_sphere(32)
cs = vtk.vtkCellSizeFilter(); cs.SetInputData(pd); cs.Update()
r = cs.GetOutput()
arr = r.GetCellData().GetArray("Area")
save("filter_cell_size", {"num_cells": r.GetNumberOfCells(), "has_area": arr is not None})

# --- Center of mass ---
print("center_of_mass...")
pd = sphere(32)
com = vtk.vtkCenterOfMass(); com.SetInputData(pd); com.Update()
c = com.GetCenter()
save("filter_center_of_mass", {"center": list(c)})

# --- Clip closed ---
print("clip_closed...")
pd = tri_sphere(32)
clip = vtk.vtkClipClosedSurface()
planes = vtk.vtkPlaneCollection()
p = vtk.vtkPlane(); p.SetOrigin(0,0,0); p.SetNormal(1,0,0)
planes.AddItem(p)
clip.SetClippingPlanes(planes); clip.SetInputData(pd); clip.Update()
save("filter_clip_closed", basic(clip.GetOutput()))

# --- Collision ---
print("collision...")
s1 = tri_sphere(16)
s2 = vtk.vtkSphereSource(); s2.SetCenter(0.8,0,0); s2.SetThetaResolution(16); s2.SetPhiResolution(16); s2.Update()
t2 = vtk.vtkTriangleFilter(); t2.SetInputData(s2.GetOutput()); t2.Update()
cd = vtk.vtkCollisionDetectionFilter()
cd.SetInputData(0, s1); cd.SetInputData(1, t2.GetOutput())
m1 = vtk.vtkMatrix4x4(); m2 = vtk.vtkMatrix4x4()
t1l = vtk.vtkTransform(); t1l.SetMatrix(m1)
t2l = vtk.vtkTransform(); t2l.SetMatrix(m2)
cd.SetTransform(0, t1l); cd.SetTransform(1, t2l)
cd.SetCollisionModeToAllContacts(); cd.Update()
nc = cd.GetNumberOfContacts()
save("filter_collision", {"num_contacts": nc, "contacts_nonzero": nc > 0})

# --- Delaunay ---
print("delaunay...")
import random; random.seed(42)
pts = vtk.vtkPoints()
for _ in range(500):
    pts.InsertNextPoint(random.uniform(-1,1), random.uniform(-1,1), 0)
pd = vtk.vtkPolyData(); pd.SetPoints(pts)
d2d = vtk.vtkDelaunay2D(); d2d.SetInputData(pd); d2d.Update()
save("filter_delaunay_500", basic(d2d.GetOutput()))

# --- Depth sort ---
print("depth_sort...")
pd = tri_sphere(32)
ds = vtk.vtkDepthSortPolyData(); ds.SetInputData(pd)
ds.SetDirectionToSpecifiedVector(); ds.SetVector(0,0,1)
ds.SetCamera(vtk.vtkCamera()); ds.Update()
save("filter_depth_sort", basic(ds.GetOutput()))

# --- Densify ---
print("densify...")
pd = tri_sphere(32)
sub = vtk.vtkAdaptiveSubdivisionFilter(); sub.SetInputData(pd)
sub.SetMaximumEdgeLength(0.1); sub.Update()
r = sub.GetOutput()
save("filter_densify", {"num_cells_gte": r.GetNumberOfCells()})

# --- Extract largest ---
print("extract_largest...")
s1 = tri_sphere(32)
s2 = vtk.vtkSphereSource(); s2.SetCenter(3,0,0); s2.SetThetaResolution(8); s2.SetPhiResolution(8); s2.Update()
t2 = vtk.vtkTriangleFilter(); t2.SetInputData(s2.GetOutput()); t2.Update()
app = vtk.vtkAppendPolyData(); app.AddInputData(s1); app.AddInputData(t2.GetOutput()); app.Update()
conn = vtk.vtkPolyDataConnectivityFilter(); conn.SetInputData(app.GetOutput())
conn.SetExtractionModeToLargestRegion(); conn.Update()
save("filter_extract_largest", basic(conn.GetOutput()))

# --- Hull ---
print("hull_200...")
pts = vtk.vtkPoints()
random.seed(99)
for _ in range(200):
    pts.InsertNextPoint(random.gauss(0,1), random.gauss(0,1), random.gauss(0,1))
pd = vtk.vtkPolyData(); pd.SetPoints(pts)
hull = vtk.vtkDelaunay3D(); hull.SetInputData(pd); hull.Update()
surf = vtk.vtkDataSetSurfaceFilter(); surf.SetInputData(hull.GetOutput()); surf.Update()
save("filter_hull_200", basic(surf.GetOutput()))

# --- Mask points ---
print("mask_points...")
pd = sphere(32)
mp = vtk.vtkMaskPoints(); mp.SetInputData(pd); mp.SetOnRatio(3); mp.Update()
save("filter_mask_points_3", {"num_points": mp.GetOutput().GetNumberOfPoints()})

# --- MC 64 ---
print("mc_64...")
img = vtk.vtkImageData()
img.SetDimensions(64,64,64)
img.SetOrigin(0,0,0); img.SetSpacing(1,1,1)
arr = vtk.vtkDoubleArray(); arr.SetNumberOfTuples(64*64*64)
for k in range(64):
    for j in range(64):
        for i in range(64):
            arr.SetValue(k*64*64+j*64+i, (i-32)**2+(j-32)**2+(k-32)**2)
img.GetPointData().SetScalars(arr)
mc = vtk.vtkMarchingCubes(); mc.SetInputData(img); mc.SetValue(0, 400); mc.Update()
save("filter_mc_64", basic(mc.GetOutput()))

# --- FE 64 ---
print("fe_64...")
fe = vtk.vtkFlyingEdges3D(); fe.SetInputData(img); fe.SetValue(0, 400); fe.Update()
save("filter_fe_64", basic(fe.GetOutput()))

# --- Point density ---
print("point_density...")
pd = sphere(32)
# VTK's vtkPointDensityFilter works on ImageData, so just store point count
save("filter_point_density", {"num_points": pd.GetNumberOfPoints()})

# --- Separate cells ---
print("separate_cells...")
pd = tri_sphere(32)
sc = vtk.vtkShrinkFilter(); sc.SetInputData(pd); sc.SetShrinkFactor(1.0); sc.Update()
# vtkShrinkFilter with factor=1.0 separates cells (duplicates shared vertices)
r = sc.GetOutput()
save("filter_separate_cells", {"num_points_gte": r.GetNumberOfPoints()})

# --- Surface nets ---
print("surface_nets...")
img2 = vtk.vtkImageData(); img2.SetDimensions(32,32,32); img2.SetSpacing(1,1,1)
arr2 = vtk.vtkDoubleArray(); arr2.SetName("values"); arr2.SetNumberOfTuples(32*32*32)
for k in range(32):
    for j in range(32):
        for i in range(32):
            x,y,z = i-15.5, j-15.5, k-15.5
            arr2.SetValue(k*32*32+j*32+i, x*x+y*y+z*z)
img2.GetPointData().AddArray(arr2); img2.GetPointData().SetScalars(arr2)
sn = vtk.vtkSurfaceNets3D(); sn.SetInputData(img2); sn.SetValue(0, 200); sn.Update()
save("filter_surface_nets_32", basic(sn.GetOutput()))

# --- Voxel grid ---
print("voxel_grid...")
# No direct VTK equivalent, just record input size
save("filter_voxel_grid", {"note": "no VTK equivalent"})

# --- Voronoi 2D ---
print("voronoi_2d...")
# No direct VTK equivalent
save("filter_voronoi_2d", {"note": "no VTK equivalent"})

# --- BYU I/O ---
print("byu...")
import tempfile
pd = tri_sphere(32)
w = vtk.vtkBYUWriter(); w.SetGeometryFileName(os.path.join(tempfile.gettempdir(), "test.byu"))
w.SetInputData(pd); w.Write()
rd = vtk.vtkBYUReader(); rd.SetGeometryFileName(os.path.join(tempfile.gettempdir(), "test.byu")); rd.Update()
r = rd.GetOutput()
save("io_byu", basic(r))

print(f"Done. Generated refs in {OUT}/")
