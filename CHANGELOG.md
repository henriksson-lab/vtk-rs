# Changelog

All notable changes to vtk-rs are documented here.

## [Unreleased]

### Data Model
- 18+ data types: PolyData, ImageData, UnstructuredGrid, RectilinearGrid, StructuredGrid, MultiBlockDataSet, Table, Graph, Tree, ExplicitStructuredGrid, HyperTreeGrid, Molecule, TemporalDataSet
- Full trait coverage: Display, PartialEq, Clone, Debug, Index, FromIterator, Extend, From, Send+Sync
- Prelude modules for quick imports
- Builder patterns on all major types
- DataArray: from_fn, from_slice, filled, map, scale, normalize, magnitude, extract_component, concat, compose, statistics, iter_tuples
- Points: from_fn, to_vec, as_flat_slice, centroid, transform, Iterator/IntoIterator/ExactSizeIterator
- PolyData: from_triangles/quads/lines/vertices/polyline/points/polygons/xyz_arrays, push_point/triangle/quad/line, append, to_table, add_scalars/vectors, get_scalars/vectors, summary, approx_eq, num_edges/triangles, reverse_cells
- CellArray: from_triangles/quads, cell_size/sizes, max_cell_size, is_homogeneous, PartialEq, IntoIterator
- Table: CSV read/write, filter_rows, sort_by_column, select_rows, value_f64, describe_column, correlation, linear_regression

### Filters
- 575+ filters, 40 sources
- Pipeline system with lazy evaluation, caching, invalidation, builder API
- Parallel filters (rayon): compute_normals_par, elevation_par, smooth_par, marching_cubes_par
- Topology analysis: manifold check, Euler characteristic, genus, boundary detection
- Dataset conversions: PolyData ↔ UnstructuredGrid, ImageData/RectilinearGrid/StructuredGrid → PolyData
- Selection extraction, merge, seed strategies, extract by cell type
- Unified I/O: auto-format read/write by file extension

### Rendering
- Blinn-Phong + Cook-Torrance PBR shading
- 4x MSAA anti-aliasing
- Per-actor model matrix (position, scale, visibility)
- Transparency with opaque-first rendering
- Edge overlay, wireframe, points representation modes
- Flat shading (dpdx/dpdy), backface culling
- 6 GPU clip planes
- Distance fog (linear, exponential, exponential²)
- Skybox gradient (solid, 2-color, 3-stop) — GPU rendered
- GPU volume rendering (3D texture, transfer function LUT, ray marching)
- Silhouette edge extraction
- 15 color map presets + by_name() lookup
- Scalar bar with bitmap font text labels
- Axes orientation widget
- 3D annotations (labels, rulers, protractors) — GPU rendered
- CPU + GPU picking
- LOD (level of detail)
- Instanced glyphs
- Camera: orbit, dolly, pan, unproject, look_at, standard views, turntable/zoom animation
- Offscreen rendering (render_to_image)
- Screenshot: save_ppm, save_bmp, save_tga
- Scene: builders, summary, print_info, find_actor_at, reset_camera, JSON export
- Measurement tools: surface area, edge lengths, point distance, angle, triangle area
- Animation: keyframe tracks, easing functions, turntable, zoom
- Transfer function editor
- Viewport system (split-screen)
- Shadow/bloom/stereo config (ready for GPU implementation)

### I/O
- VTK Legacy (.vtk): ASCII + binary, PolyData + ImageData + UnstructuredGrid
- VTK XML (.vtp/.vtu/.vti/.vtr/.vts/.vtm): ASCII + binary read/write
- STL: ASCII + binary
- OBJ: read/write
- PLY: ASCII + binary
- glTF (.glb): read/write
- EnSight Gold: ASCII read/write
- XDMF: inline XML writer
- LS-DYNA: keyword reader
- CSV: Table read/write

### Infrastructure
- 14 crates, 102K+ lines
- 2444+ tests (unit, integration, doctest, proptest, Send+Sync)
- 10 examples
- Python bindings (PyO3)
- WASM support (wasm32-unknown-unknown)
- Benchmarks
- Property-based testing (proptest)
- Crate-level documentation with examples
