# FEATURES.md

Feature tracking for vtk-rs — a pure Rust reimplementation of VTK.

Last updated: 2026-03-26 | Tests: 296 | Clippy: clean

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
- [x] `RectilinearGrid` — axis-aligned grid with non-uniform per-axis coordinate arrays
- [x] `StructuredGrid` — curvilinear grid with explicit points and structured i×j×k topology
- [x] `BoundingBox`
- [x] `Scalar` trait + `ScalarType` enum (10 types: f32, f64, i8–u64)
- [x] `CellType` enum — all VTK linear + quadratic cell types
- [x] `VtkError` — Io, Parse, InvalidData, Unsupported, IndexOutOfBounds

### Not Yet Implemented

- [x] `Table` — columnar data (rows × named columns)
- [ ] `Graph` / `Tree` — graph data structures
- [x] `MultiBlockDataSet` — composite dataset with named heterogeneous blocks
- [ ] `HyperTreeGrid` — AMR-style hierarchical grids
- [ ] `ExplicitStructuredGrid`
- [x] `Selection` / `SelectionNode` — content types: Indices, Thresholds, GlobalIds, Frustum, Blocks
- [x] Implicit functions (`ImplicitPlane`, `ImplicitSphere`, `ImplicitBox` + `ImplicitFunction` trait)
- [x] `KdTree` — k-d tree for nearest-neighbor, k-NN, and radius queries
- [x] `OctreePointLocator` — octree for nearest-neighbor and radius queries
- [ ] Spatial locators (`CellLocator`)
- [ ] Higher-order / Lagrange / Bezier cell support
- [ ] `Molecule` data type

---

## Geometry Sources (17 / ~40 in VTK)

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
- [x] `superquadric` — superellipsoid with configurable roundness exponents
- [x] `platonic_solid` — tetrahedron, octahedron, icosahedron, dodecahedron
- [x] `frustum` — truncated cone/pyramid with configurable radii and caps

- [x] `parametric_function` — evaluate parametric surfaces (+ torus, Klein bottle helpers)
- [x] `bounding_box_source` — wireframe bounding box from bounds or PolyData
- [x] `axes` — XYZ axis triad with optional arrowhead cones

### Not Yet Implemented

- [ ] `text` — 2D/3D text geometry

---

## Processing Filters (65 / ~300+ in VTK)

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
| `delaunay_2d` | 2D Delaunay triangulation (Bowyer-Watson algorithm) |
| `spline` | Catmull-Rom spline interpolation along polylines |
| `clip_data_set` | Clip UnstructuredGrid by plane (whole-cell removal) |
| `windowed_sinc_smooth` | Low-pass windowed sinc smoothing (Taubin method) |
| `probe` | Interpolate source data at probe points (nearest-neighbor) |
| `stream_tracer` | RK4 streamline integration through vector fields |
| `voxel_modeller` | Convert PolyData surface to binary voxel ImageData |
| `sample_function` | Evaluate scalar function on ImageData grid |
| `integrate_attributes` | Integrate point data over surface (area-weighted) |
| `distance_poly_data` | Minimum distance from target points to source surface |
| `implicit_modeller` | Distance field from PolyData surface on ImageData grid |
| `tensor_glyph` | Place tensor-transformed ellipsoid glyphs at input points |
| `resample` | Resample source PolyData onto target ImageData grid |
| `calculator` | Expression-based scalar/vector field computation |
| `select_enclosed_points` | Ray-casting inside/outside classification |
| `extract_cells` | Extract cells by index or predicate |
| `icp` | Iterative Closest Point rigid registration (SVD-based) |
| `extract_points` | Extract points by index or scalar range |
| `signed_distance` | Signed distance field from closed surface |
| `cell_size` | Compute area/length of polygon/line cells |
| `extrude` | Linear extrusion of 2D geometry along a direction |
| `cell_quality` | Aspect ratio, min/max angle, area metrics per cell |
| `reverse_sense` | Flip polygon winding and normals |
| `strip` | Convert triangles to triangle strips (greedy) |
| `interpolate` | Inverse distance weighting interpolation |
| `fill_holes` | Close open boundary loops with fan triangulation |
| `center_of_mass` | Point and area-weighted center of mass |
| `project_points` | Project points onto nearest surface location |
| `subdivide_midpoint` | Midpoint subdivision (triangles → 4, quads → 4) |
| `random_attributes` | Generate random scalar/vector point and cell data |
| `image_to_poly_data` | Convert ImageData surface to PolyData quads |

### High-Priority — Not Yet Implemented

- [ ] `delaunay_3d` — 3D Delaunay tetrahedralization
- [ ] `boolean` — Boolean operations on PolyData (union, intersection, difference)
- [ ] `flying_edges_3d` — Parallelizable marching cubes variant

### Lower-Priority — Not Yet Implemented

- [ ] `temporal_interpolator` — Interpolate between time steps

---

## I/O Formats (11 / ~60+ in VTK)

### Implemented

| Format | Crate | Read | Write | Notes |
|--------|-------|:----:|:-----:|-------|
| VTK Legacy `.vtk` | vtk-io-legacy | yes | yes | PolyData + ImageData + UnstructuredGrid, ASCII & binary |
| VTK XML `.vtp` | vtk-io-xml | yes | yes | PolyData, ASCII only |
| VTK XML `.vtu` | vtk-io-xml | yes | yes | UnstructuredGrid, ASCII only |
| VTK XML `.vti` | vtk-io-xml | yes | yes | ImageData, ASCII only |
| VTK XML `.vtr` | vtk-io-xml | yes | yes | RectilinearGrid, ASCII only |
| VTK XML `.vts` | vtk-io-xml | yes | yes | StructuredGrid, ASCII only |
| STL `.stl` | vtk-io-stl | yes | yes | ASCII & binary |
| Wavefront `.obj` | vtk-io-obj | yes | yes | Vertices, normals, texture coords, faces |
| Stanford PLY `.ply` | vtk-io-ply | yes | yes | ASCII & binary little-endian |
| VTK XML `.vtm` | vtk-io-xml | | yes | MultiBlock index file (write only) |

### High-Priority — Not Yet Implemented

- [ ] VTK XML binary/appended — Binary and appended-data modes for `.vtp`/`.vtu`/`.vti`
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
- [x] 296 unit tests
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
