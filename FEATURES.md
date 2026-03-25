# FEATURES.md

Feature tracking for vtk-rs — a pure Rust reimplementation of VTK.

Last updated: 2026-03-25 | Tests: 174 | Clippy: clean

---

## Data Model

### Implemented

- [x] `DataArray<T>` — generic contiguous tuple array (f32, f64, i8–u64)
- [x] `AnyDataArray` — type-erased enum over all `DataArray<T>` variants
- [x] `CellArray` — offsets + connectivity (mirrors `vtkCellArray`)
- [x] `Points<T>` — 3-component point coordinates (default f64)
- [x] `FieldData` — named collection of `AnyDataArray`
- [x] `DataSetAttributes` — active scalars/vectors/normals/tcoords/tensors
- [x] `DataObject` trait, `DataSet` trait
- [x] `PolyData` — 4 cell arrays (verts/lines/polys/strips) + point/cell data
- [x] `ImageData` — regular grid with extent/spacing/origin
- [x] `UnstructuredGrid` — mixed-cell mesh with explicit connectivity + cell types
- [x] `BoundingBox`
- [x] `Scalar` trait + `ScalarType` enum (10 types: f32, f64, i8–u64)
- [x] `CellType` enum — all VTK linear + quadratic cell types
- [x] `VtkError` — Io, Parse, InvalidData, Unsupported, IndexOutOfBounds

### Not Yet Implemented

- [ ] `StructuredGrid` — curvilinear grid with explicit points
- [ ] `RectilinearGrid` — axis-aligned grid with per-axis coordinate arrays
- [ ] `Table` — columnar data (rows × named columns)
- [ ] `Graph` / `Tree` — graph data structures
- [ ] `MultiBlockDataSet` / `PartitionedDataSet` — composite datasets
- [ ] `HyperTreeGrid` — AMR-style hierarchical grids
- [ ] `ExplicitStructuredGrid`
- [ ] `Selection` / `SelectionNode`
- [ ] Implicit functions (`Plane`, `Sphere`, `Box` as evaluatable objects)
- [ ] Spatial locators (`KdTree`, `OctreePointLocator`, `CellLocator`)
- [ ] Higher-order / Lagrange / Bezier cell support
- [ ] `Molecule` data type

---

## Geometry Sources (11 / ~40 in VTK)

### Implemented

- [x] `sphere` — UV sphere with configurable resolution
- [x] `cube` — 6-face cube with per-face normals (24 vertices)
- [x] `cone` — configurable direction, resolution, optional cap
- [x] `cylinder` — Y-axis aligned, optional caps
- [x] `plane` — parallelogram with subdivisions and texture coordinates
- [x] `arrow` — composite cylinder shaft + cone tip
- [x] `disk` — flat disk or annulus with configurable inner/outer radius and resolution
- [x] `line` — line segment with configurable resolution (polyline)
- [x] `point_source` — random point cloud within a sphere (rejection sampling)
- [x] `regular_polygon` — regular N-gon with optional outline-only mode
- [x] `arc` — circular arc as a polyline with configurable angles and normal

### Not Yet Implemented

- [ ] `superquadric` — superellipsoid source
- [ ] `platonic_solid` — tetrahedron, octahedron, icosahedron, dodecahedron
- [ ] `frustum` — truncated pyramid
- [ ] `parametric_function` — evaluate parametric surfaces
- [ ] `text` — 2D/3D text geometry

---

## Processing Filters (34 / ~300+ in VTK)

### Implemented

| Filter | Description |
|--------|-------------|
| `normals` | Smooth vertex normals (Newell's method + averaging) |
| `triangulate` | Quads/polygons/strips → triangles (fan triangulation) |
| `append` | Merge multiple PolyData with index renumbering |
| `clean` | Merge duplicate points (spatial hash), remove degenerate cells |
| `transform` | Apply 4×4 matrix to points and normals |
| `marching_cubes` | Isosurface extraction from scalar field on ImageData |
| `clip` | Clip mesh by plane, splitting crossing triangles |
| `slice` | Plane–mesh intersection → line segments |
| `contour` | Contour lines at scalar isovalues (2D marching) |
| `elevation` | Scalar from projection along an axis |
| `threshold` | Extract cells by scalar value range |
| `decimate` | Quadric error metric mesh simplification |
| `smooth` | Laplacian smoothing with boundary preservation |
| `subdivide` | Loop subdivision (each triangle → 4) |
| `warp` | Displace by scalar (along normals) or vector field |
| `connectivity` | Connected components via union-find |
| `extract_surface` | Boundary surface of UnstructuredGrid |
| `feature_edges` | Boundary, feature, manifold, non-manifold edge extraction |
| `extract_edges` | Extract all unique edges as line segments |
| `reflect` | Mirror across coordinate plane with winding reversal |
| `shrink` | Shrink cells toward centroids |
| `orient` | Consistent polygon winding via BFS edge traversal |
| `cell_centers` | Generate vertex points at cell centroids |
| `cell_data_to_point_data` | Average cell values to shared points |
| `point_data_to_cell_data` | Average point values per cell |
| `generate_ids` | Generate sequential PointIds and/or CellIds arrays |
| `mask_points` | Subsample points (every Nth or random ratio) |
| `curvatures` | Discrete Gaussian and mean curvature (angle deficit + cotangent weights) |
| `gradient` | Least-squares gradient of scalar field over one-ring neighborhoods |
| `mass_properties` | Surface area, volume, and centroid of closed triangle meshes |
| `densify` | Subdivide polygon edges exceeding a max length |
| `hull` | 3D convex hull via incremental algorithm |
| `texture_map` | Generate texture coords by plane projection or spherical mapping |
| `glyph` | Place scaled mesh copies at input points |
| `tube` | Generate tubes with caps around line cells |

### High-Priority — Not Yet Implemented

- [ ] `delaunay_2d` — 2D Delaunay triangulation
- [ ] `delaunay_3d` — 3D Delaunay tetrahedralization
- [ ] `boolean` — Boolean operations on PolyData (union, intersection, difference)
- [ ] `probe` — Interpolate source data at probe point locations
- [ ] `stream_tracer` — Streamline integration through vector fields
- [ ] `windowed_sinc_smooth` — Windowed sinc smoothing (better than Laplacian)
- [ ] `distance_poly_data` — Signed distance between two meshes
- [ ] `spline` — Fit splines through points
- [ ] `flying_edges_3d` — Parallelizable marching cubes variant
- [ ] `clip_data_set` — Clip UnstructuredGrid (not just PolyData)

### Lower-Priority — Not Yet Implemented

- [ ] `resample_with_dataset` — Resample one dataset onto another
- [ ] `integrate_attributes` — Integrate field data over cells
- [ ] `texture_map_to_plane` / `texture_map_to_sphere` — Texture coordinate generation
- [ ] `voxel_modeller` — Convert PolyData to ImageData
- [ ] `implicit_modeller` — Distance field from PolyData
- [ ] `sample_function` — Evaluate implicit function on grid
- [ ] `temporal_interpolator` — Interpolate between time steps
- [ ] `tensor_glyph` — Visualize tensor fields with ellipsoids

---

## I/O Formats (8 / ~60+ in VTK)

### Implemented

| Format | Crate | Read | Write | Notes |
|--------|-------|:----:|:-----:|-------|
| VTK Legacy `.vtk` | vtk-io-legacy | yes | yes | PolyData + ImageData + UnstructuredGrid, ASCII & binary |
| VTK XML `.vtp` | vtk-io-xml | yes | yes | PolyData, ASCII only |
| VTK XML `.vtu` | vtk-io-xml | yes | yes | UnstructuredGrid, ASCII only |
| VTK XML `.vti` | vtk-io-xml | yes | yes | ImageData, ASCII only |
| STL `.stl` | vtk-io-stl | yes | yes | ASCII & binary |
| Wavefront `.obj` | vtk-io-obj | yes | yes | Vertices, normals, texture coords, faces |
| Stanford PLY `.ply` | vtk-io-ply | yes | yes | ASCII & binary little-endian |

### High-Priority — Not Yet Implemented

- [ ] VTK XML binary/appended — Binary and appended-data modes for `.vtp`/`.vtu`/`.vti`
- [ ] VTK XML `.vtr` — RectilinearGrid XML format
- [ ] VTK XML `.vts` — StructuredGrid XML format
- [ ] VTK XML `.vtm` — MultiBlock XML format
- [ ] VTP binary/appended — Binary and appended-data modes for `.vtp`
- [ ] glTF `.gltf` / `.glb` — For interchange with 3D tools and web viewers

### Lower-Priority — Not Yet Implemented

- [ ] EnSight — EnSight Gold format
- [ ] Exodus — Exodus II (HDF5-based FEM format)
- [ ] CGNS — CFD General Notation System
- [ ] NetCDF — Network Common Data Format
- [ ] XDMF — eXtensible Data Model and Format
- [ ] LSDyna — LS-DYNA result files
- [ ] Alembic — Alembic interchange format
- [ ] OpenVDB — Sparse volumetric data
- [ ] USD — Universal Scene Description

---

## Rendering

### Implemented

- [x] `Camera` — position, focal point, view up, FOV, clip planes, orbit, dolly
- [x] `Scene` / `Actor` — scene graph with actors holding PolyData + coloring
- [x] `Renderer` trait — backend-agnostic rendering abstraction
- [x] `ColorMap` — jet, viridis, cool-to-warm, grayscale with linear interpolation
- [x] `Coloring::Solid` / `Coloring::ScalarMap` — per-actor coloring modes
- [x] wgpu backend — Phong shading, smooth normals, scalar color mapping
- [x] Mouse interaction — orbit (left drag), zoom (scroll)
- [x] Depth buffer

### Not Yet Implemented

- [ ] Wireframe rendering mode
- [ ] Points rendering mode
- [ ] Edge overlay (surface + edges)
- [ ] Transparency / alpha blending
- [ ] Multiple lights / light types (point, spot, ambient)
- [ ] Texture mapping (2D textures on surfaces)
- [ ] Volume rendering (ray casting on ImageData)
- [ ] Axes / orientation widget
- [ ] Scalar bar (color legend)
- [ ] Text / label rendering (2D overlay)
- [ ] Screenshot / offscreen rendering
- [ ] Picking (point/cell selection via mouse click)
- [ ] Anti-aliasing (MSAA)
- [ ] PBR materials (physically-based rendering)
- [ ] Silhouette / outline rendering
- [ ] LOD (level of detail) for large meshes
- [ ] Instanced rendering for glyphs

---

## Infrastructure

### Implemented

- [x] Workspace with 12 crates
- [x] 4 interactive examples (triangle, shapes, isosurface, scalar_colors)
- [x] 174 unit tests
- [x] Clippy-clean (`-D warnings`)

### Not Yet Implemented

- [ ] Pipeline system (lazy evaluation, caching, automatic updates)
- [ ] Parallel filter execution (rayon integration)
- [ ] WASM / web target support
- [ ] Python bindings (PyO3)
- [ ] Benchmarks
- [ ] CI (GitHub Actions)
- [ ] Published to crates.io
- [ ] Documentation (rustdoc with examples)
- [ ] Property / fuzz testing
