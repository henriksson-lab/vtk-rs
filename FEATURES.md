# FEATURES.md

Feature tracking for vtk-rs — a pure Rust reimplementation of VTK.

Last updated: 2026-04-04 | ~296K lines Rust | Tests: ~9100 | Clippy: clean | 5282 source files across 48 crates | Coverage: ~25% of VTK C++ (1M lines)

---

## Data Model

### Implemented

- [x] `DataArray<T>` — `from_fn()`, `filled()`, `iter_tuples()`, `map()`, `scale()`, `normalize()`, `magnitude()`, `extract_component()`, `concat()`, `compose()`
- [x] `AnyDataArray` — type-erased enum over all `DataArray<T>` variants + `statistics()`, `range()`
- [x] `CellArray` — offsets + connectivity + `cell_size()`, `cell_sizes()`, `is_homogeneous()`, `IntoIterator`
- [x] `Points<T>` — `from_fn()`, `to_vec()`, `as_flat_slice()`, `centroid()`, `transform()`, `Iterator`/`IntoIterator`
- [x] `FieldData` — named collection + `has_array()`, `names()`, `clear()`
- [x] `DataSetAttributes` — active attributes + `has_array()`, `array_names()`, `remove_array()`, `iter()`
- [x] `DataObject` trait, `DataSet` trait
- [x] `PolyData` — `from_triangles/quads/lines/vertices/polyline/xyz_arrays`, `append()`, `to_table()`, `summary()`, `approx_eq()`, `num_edges()`, builder
- [x] `ImageData` — `from_function()`, `with_spacing/origin/point_array` builders
- [x] `UnstructuredGrid` — mixed-cell mesh + `from_tetrahedra/hexahedra` + builder pattern
- [x] `RectilinearGrid` — axis-aligned grid + `uniform()` + builder pattern
- [x] `StructuredGrid` — curvilinear grid + `from_function()`, `uniform()`, builder
- [x] `BoundingBox` — size, volume, contains, intersects, union, intersection, pad, from_corners
- [x] `Scalar` trait + `ScalarType` enum (10 types: f32, f64, i8–u64)
- [x] `CellType` enum — all VTK linear + quadratic cell types
- [x] `VtkError` — Io, Parse, InvalidData, Unsupported, IndexOutOfBounds, DimensionMismatch, EmptyData + helpers

- [x] `Table` — columnar data + `filter_rows()`, `sort_by_column()`, `select_rows()`, `value_f64()`, builder
- [x] `Graph` / `Tree` — directed/undirected graph with vertex/edge data, tree with parent/children/depth
- [x] `MultiBlockDataSet` — composite dataset + typed getters, `flatten()`, `remove_by_name()`, builder
- [x] `HyperTreeGrid` — AMR octree/quadtree grid with recursive subdivision and per-cell data
- [x] `ExplicitStructuredGrid` — structured grid with per-cell blanking/visibility
- [x] `Selection` / `SelectionNode` — Indices, Thresholds, GlobalIds, Frustum, Blocks + `extract_selection()` filter
- [x] Implicit functions (`ImplicitPlane`, `ImplicitSphere`, `ImplicitBox` + `ImplicitFunction` trait)
- [x] `KdTree` — k-d tree for nearest-neighbor, k-NN, and radius queries
- [x] `OctreePointLocator` — octree for nearest-neighbor and radius queries
- [x] `CellLocator` — BVH-based cell locator for closest-cell and radius queries
- [x] Higher-order / Lagrange / Bezier cell support — cell types (68-81), Lagrange/Bernstein basis evaluation, curve/quad tessellation
- [x] `Molecule` — atoms (atomic number, position) + bonds (connectivity, order), CPK colors, element symbols

### Not Yet Implemented

- [ ] `AMRDataSet` — multi-level adaptive mesh refinement dataset (VTK's `vtkOverlappingAMR`)
- [ ] `PartitionedDataSet` / `PartitionedDataSetCollection` — distributed data decomposition
- [ ] `CellGrid` — discontinuous Galerkin / high-order cell representation (VTK 9.3+ feature)
- [ ] `GenericDataSet` — generic cell interface for arbitrary cell types
- [ ] `TemporalDataSetCache` — in-memory cache for temporal datasets

---

## Geometry Sources (62 / ~45 in VTK)

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
- [x] `torus` — torus with configurable ring and cross-section radius/resolution
- [x] `helix` — helical spiral polyline with configurable turns, height, and radius
- [x] `ellipsoid` — UV ellipsoid with three semi-axis radii and smooth normals
- [x] `spring` — helical tube (coil) with configurable coil/tube radius and turns
- [x] `capsule` — cylinder with hemispherical caps, configurable radius/length
- [x] `geodesic_sphere` — icosphere by recursive subdivision of icosahedron
- [x] `grid` — rectangular grid of quads with texture coordinates and normals

- [x] `text_3d` — 3D text geometry with built-in vector font, optional extrusion
- [x] `wavelet` — analytic wavelet scalar field on ImageData for testing
- [x] `circle` — closed circle polyline with configurable normal direction

- [x] `mobius` — Möbius strip with half-twist
- [x] `star` — star polygon with configurable inner/outer radii
- [x] `noise_field` — 3D value-noise scalar field on ImageData
- [x] `ring` — thick ring (torus with small tube radius)
- [x] `klein_bottle` — Klein bottle immersion surface
- [x] `trefoil_knot` — trefoil knot tube geometry
- [x] `cross` — 3D cross / plus-sign geometry
- [x] `boy_surface` — Boy's surface immersion of the projective plane
- [x] `spiral` — spiral curve polyline
- [x] `icosphere` — icosphere by recursive subdivision
- [x] `mobius_strip` — Möbius strip surface variant
- [x] `gear` — gear/cog wheel polygon
- [x] `grid_2d` — 2D grid with configurable resolution
- [x] `earth` — Earth globe sphere with latitude/longitude arrays and equirectangular UV
- [x] `sector` — Sector (pie-slice) or annular sector in XY plane
- [x] `arrow_2d` — 2D arrow polygon glyph for vector field visualization
- [x] `prism` — Triangular/hexagonal prism with caps and sides
- [x] `hemisphere` — Hemisphere with configurable resolution and optional cap
- [x] `annulus` — Flat ring (annulus) in XY plane
- [x] `stadium` — Rounded rectangle (discorectangle) in XY plane
- [x] `lens` — Biconvex lens shape from intersecting sphere caps
- [x] `superellipse` — 2D superellipse + 3D superellipsoid with configurable exponents
- [x] `pyramid` — Square-base pyramid + truncated pyramid (frustum)
- [x] `polygon_extrude` — Extrude 2D polygon to 3D solid along Z or arbitrary direction
- [x] `wave_surface` — Sine wave + ripple (concentric wave) parametric surfaces
- [x] `terrain` — Procedural terrain generation (fBm value noise + ridged multifractal)
- [x] `cloth` — Cloth simulation mesh with mass/fixed arrays + spring-mass simulation step
- [x] `arrow_3d` — 3D arrow from point A to B with shaft cylinder + cone tip
- [x] `heightfield` — Mesh from 2D function f(x,y)→z or flat Z-value array
- [x] `coordinate_axes` — RGB-colored XYZ axes with optional arrowhead cones
- [x] `grid_3d` — 3D structured grid surface (boundary quads) and wireframe
- [x] `torus_knot` — (p,q) torus knot tube geometry with configurable radii
- [x] `spiral_surface` — Logarithmic spiral surface + spiral curve
- [x] `catenary` — Catenary curve, catenoid minimal surface, bridge cable
- [x] `seashell` — Parametric seashell/nautilus with spiral growth model

### Not Yet Implemented

---

## Processing Filters (~4093 modules / ~938 in VTK C++)

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
| `flying_edges` | Flying Edges 3D — efficient scanline marching cubes |
| `hausdorff` | Hausdorff and mean distance between point sets |
| `color_transfer` | Piecewise-linear color transfer function with presets |
| `extrude_normals` | Extrude surface along vertex normals |
| `cell_normals` | Per-cell face normals as cell data |
| `merge_points` | Tolerance-based coincident point merging |
| `bounding_box_filter` | Per-cell bounding boxes as cell data |
| `data_array_math` | Array arithmetic: add, subtract, multiply, scale, magnitude |
| `edge_lengths` | Min/max/mean edge lengths per cell |
| `extract_region` | Extract sub-region of ImageData by index extent |
| `warp_implicit` | Warp points along implicit function gradient |
| `delaunay_3d` | 3D Delaunay tetrahedralization (Bowyer-Watson) |
| `boolean` | Boolean operations on PolyData (union, intersection, difference) |
| `ribbon` | Create oriented ribbon strips from polyline cells |
| `rotation_extrude` | Rotational extrusion (surface of revolution) around Y axis |
| `geometry_filter` | Extract boundary surface from UnstructuredGrid or ImageData |
| `banded_contour` | Filled contour bands from scalar field with BandIndex cell data |
| `ruled_surface` | Create ruled surface between two polylines |
| `multi_threshold` | Threshold by multiple scalar intervals with ThresholdId |
| `data_set_triangulate` | Decompose hex/wedge/pyramid cells to tetrahedra |
| `point_sampler` | Area-proportional point sampling on triangle mesh surface |
| `outline` | Wireframe bounding box outline (12 line segments) |
| `histogram` | Scalar histogram computation returning Table with bin data |
| `attribute_smooth` | Iterative Laplacian smoothing of point data attributes |
| `dicer` | Spatial grid-based region decomposition with RegionId |
| `intersection_poly_data` | Compute intersection lines between two triangle meshes |
| `quadric_clustering` | Fast O(n) mesh simplification by spatial vertex clustering |
| `adaptive_subdivide` | Subdivide triangles exceeding max edge length |
| `poly_data_to_image_data` | Convert surface mesh to signed distance field on ImageData |
| `normal_estimation` | PCA-based normal estimation for unstructured point clouds |
| `array_rename` | Rename point or cell data arrays |
| `vertex_glue` | Grid-snapping vertex merging with tolerance |
| `iso_volume` | Extract surface between two scalar isovalues |
| `surface_nets` | Smooth quad-dominant isosurface extraction from ImageData |
| `voronoi_2d` | 2D Voronoi diagram from Delaunay dual with SiteId cell data |
| `median_smooth` | Median filtering of scalar data (robust to outliers) |
| `extract_component` | Extract connected component by index (largest=0) |
| `pass_arrays` | Select/remove specific data arrays by name |
| `winding_number` | Generalized winding number for inside/outside classification |
| `mark_boundary` / `extract_boundary` | Identify and extract boundary edges |
| `scalar_range` / `normalize_array` / `rescale_array` | Scalar range queries and normalization |
| `point_merge` | Merge multiple PolyData with k-d tree deduplication |
| `rectilinear_to_poly_data` | Convert RectilinearGrid surface to PolyData quads |
| `structured_to_poly_data` | Convert StructuredGrid surface to PolyData quads |
| `offset_surface` | Offset mesh along vertex normals by distance |
| `table_to_poly_data` | Convert Table columns (X,Y,Z) to PolyData points |
| `image_threshold` / `image_binary_threshold` | Threshold/mask ImageData scalar fields |
| `image_gaussian_smooth` | Separable 3D Gaussian blur on ImageData |
| `image_gradient` | Central-difference gradient + magnitude on ImageData |
| `image_dilate` / `image_erode` | Morphological dilation and erosion on ImageData |
| `programmable_filter` | Apply user-defined closure to point positions/scalars |
| `field_data_to_attribute` | Copy arrays between field data and point/cell data |
| `image_resample` | Trilinear resampling of ImageData to new resolution |
| `image_laplacian` | Discrete Laplacian (7-point stencil) on ImageData |
| `aggregate` / `quick_stats` | Statistical aggregation (min/max/mean/variance/std) of scalar arrays |
| `subdivide_butterfly` | Modified butterfly subdivision for smoother meshes |
| `poly_data_distance` / `distance_stats` | K-d tree point-to-point distance + symmetric stats |
| `image_flip` | Flip ImageData along X/Y/Z axes |
| `cell_data_to_point_data` / `point_data_to_cell_data` | Averaging attribute interpolation |
| `compute_area` / `cell_areas` / `compute_length` | Surface area and line length computation |
| `image_add` / `image_subtract` / `image_multiply` / `image_scale` | Element-wise ImageData math |
| `deep_copy` / `copy_attributes` | Deep copy with optional data stripping |
| `extract_surface_by_normal` | Extract faces by normal direction alignment |
| `select_by_scalar` | Select points/cells by scalar value range |
| `image_crop` | Crop ImageData to sub-extent by index ranges |
| `mesh_quality` / `mesh_quality_arrays` | Mesh quality stats + per-cell aspect ratio/area |
| `resample_to_image` | Resample PolyData scalars onto ImageData grid via k-d tree |
| `add_coordinate_arrays` / `distance_to_point` | Add X/Y/Z arrays and point-to-reference distance |
| `clip_by_scalar` | Clip mesh by scalar isovalue with triangle splitting |
| `image_distance_transform` | Chamfer distance transform on binary ImageData |
| `close_holes` | Close boundary loops by centroid fan triangulation |
| `image_statistics` / `image_histogram` | ImageData scalar statistics and histograms |
| `sample_along_line` | Sample scalar values along a line via k-d tree |
| `collapse_edges` | Collapse short edges with midpoint merging |
| `image_pad` | Pad ImageData with constant values on each side |
| `unstructured_to_poly_data` | Convert UnstructuredGrid cells to PolyData faces |
| `point_cloud_density` | Local point density via k-d tree radius search |
| `remesh` | Isotropic remeshing via iterative edge splitting |
| `image_sobel` | Sobel edge detection on ImageData |
| `poly_line_to_strip` | Join line segments into connected polylines |
| `centroid_filter` / `weighted_centroid` | Cell centroids and weighted center of mass |
| `image_median` | Median filtering on ImageData (robust to noise) |
| `triangulate_strips` / `triangles_to_strips` | Convert between triangle strips and triangles |
| `image_connected_components` | Flood-fill labeling of connected regions in ImageData |
| `temporal_interpolator` | Linear interpolation between two PolyData (positions + scalars) |
| `separate_cells` | Duplicate shared vertices so each cell is independent |
| `image_normalize` / `image_invert` | Normalize ImageData to [0,1] or invert values |
| `poly_data_summary` / `is_manifold` / `is_triangle_mesh` | Mesh topology queries |
| `image_and` / `image_or` / `image_xor` / `image_not` | Logical operations on binary ImageData |
| `dual_mesh` | Dual mesh from triangle mesh (centroids become vertices) |
| `spherical_coordinates` / `cylindrical_coordinates` | Add coordinate system arrays |
| `image_histogram_equalize` | Histogram equalization for contrast enhancement |
| `extract_largest` / `extract_n_largest` | Extract N largest connected components |
| `auto_orient_normals` | Orient polygon normals outward from mesh centroid |
| `image_convolve_3x3` / `image_sharpen` | Custom 3x3 convolution and unsharp masking |
| `merge_coplanar` | Group coplanar adjacent faces by normal similarity |
| `split_by_array` | Split PolyData into parts by cell data array values |
| `image_variance` | Local variance computation for texture analysis |
| `point_to_vertex` / `point_to_poly_vertex` | Add vertex cells to point clouds |
| `image_downsample` | Block-averaging downsampling of ImageData |
| `dihedral_angles` | Compute dihedral angles between adjacent triangles |
| `sample_on_sphere` | Sample scalar field onto spherical probe via k-d tree |
| `translate` / `scale` / `rotate` | Direct geometric transformations (Rodrigues) |
| `poly_data_to_table` | Convert PolyData coordinates + scalars to Table |
| `image_abs` / `image_sqrt` / `image_log` / `image_exp` / `image_clamp` | Unary math on ImageData |
| `bounding_box_overlap` / `bounding_boxes_overlap` | Bounding box intersection tests |
| `image_window_level` | Window/level contrast adjustment for medical imaging |
| `compute_tangents` | Tangent vectors from texture coords or edge direction |
| `solid_color` / `color_by_height` / `color_by_cell_index` | Vertex and cell coloring |
| `image_bilateral` | Edge-preserving bilateral filtering on ImageData |
| `barycentric_coordinates` | Barycentric coords for probe points on triangle mesh |
| `catmull_clark` | Catmull-Clark subdivision for quad/triangle meshes |
| `image_power_spectrum` | DFT-based power spectrum of 1D ImageData signals |
| `snap_to_grid` / `quantize_points` | Snap points to regular grid or quantize levels |
| `bounding_box_mesh` / `oriented_bounding_box` | AABB mesh and PCA-based OBB computation |
| `edge_collapse_quadric` | Garland-Heckbert quadric error metric mesh simplification |
| `otsu_threshold` / `image_otsu` | Otsu's automatic threshold + binary segmentation |
| `vertex_valence` / `valence_histogram` | Vertex valence computation and histogram |
| `laplacian_deform` | Handle-based Laplacian mesh deformation |
| `euler_characteristic` / `genus` / `count_boundary_loops` | Topological invariants |
| `star` source | Star polygon with inner/outer radii |
| `taubin_smooth` | Volume-preserving Taubin smoothing (alternating lambda/mu) |
| `image_min` / `image_max` / `image_diff` | Element-wise min/max/abs-diff of ImageData |
| `remove_duplicate_cells` / `remove_degenerate_cells` / `remove_unused_points` | Mesh repair utilities |
| `cross_section_areas` | Approximate cross-section areas along an axis |
| `array_add` / `array_subtract` / `array_multiply` / `array_divide` | Point data array arithmetic |
| `mirror` | Mirror mesh across coordinate plane with optional merge |
| `image_hessian_det` | Hessian determinant for blob/saddle detection |
| `center_mesh` / `normalize_mesh` / `align_to_z` | Mesh centering, normalization, PCA alignment |
| `sample_along_polyline` | Sample scalar field along an existing polyline |
| `curvature_flow` | Mean curvature flow smoothing |
| `image_local_maxima` / `image_local_minima` | Local extrema detection in ImageData |
| `measure_mesh` | Comprehensive mesh measurements (area, edges, cell types) |
| `image_percentile` / `image_clip_percentile` | Percentile queries and percentile-based clipping |
| `geodesic_distance` | Dijkstra-based geodesic distance from seed points |
| `image_divergence` / `image_curl_magnitude` | Divergence and curl of vector fields on ImageData |
| `translate_to_match` / `registration_error` | Translation-only point set registration + error |
| `cell_data_to_face_varying` | Convert cell data to per-face-vertex point data |
| `image_island_remove` | Remove small connected components from binary ImageData |
| `skeleton_2d` | 2D medial axis via Voronoi dual of Delaunay |
| `image_compose_vector` / `image_decompose_vector` / `image_vector_magnitude` | Vector field composition |
| `subdivide_sqrt3` | Sqrt(3) subdivision scheme |
| `symmetry_score` / `best_symmetry_plane` | Bilateral symmetry detection |
| `image_sphere_roi` / `image_box_roi` | Region-of-interest masking |
| `ring` source | Thick ring geometry |
| `ray_cast` / `ray_cast_all` | Möller–Trumbore ray-mesh intersection |
| `image_convolve_separable` / `box_kernel` / `gaussian_kernel` | Separable 3D convolution |
| `geodesic_path` / `geodesic_path_length` | Shortest path on mesh via Dijkstra |
| `vertex_cluster` | Farthest-point-sampling vertex clustering |
| `heat_diffusion` | Heat equation simulation on mesh connectivity |
| `image_adaptive_threshold` / `image_binarize` | Adaptive and global binarization |
| `local_feature_size` / `average_spacing` | Point cloud feature size and spacing |
| `planar_parameterize` / `cylindrical_parameterize` | UV parameterization methods |
| `pick_closest_point` / `pick_points_in_sphere` / `pick_closest_cell` | Point/cell picking |
| `sdf_from_oriented_points` | SDF from oriented point cloud via k-d tree |
| `edge_flip_delaunay` | Delaunay-criterion edge flipping for quality improvement |
| `image_watershed` | Marker-based watershed segmentation |
| `find_correspondences` / `transfer_attribute` | Point set correspondence and attribute transfer |
| `voxel_downsample` / `random_downsample` | Point cloud downsampling methods |
| `texture_atlas` | Per-triangle texture atlas UV parameterization |
| `image_slic` | SLIC superpixel segmentation |
| `signed_volume` / `surface_area` / `compactness` | Closed mesh volume, area, and sphericity |
| `heat_kernel_signature` | Heat Kernel Signature shape descriptor |
| `image_moments` | Spatial moments (center of mass, variance) of ImageData |
| `extract_boundary_loops` / `num_boundary_loops` | Extract and count boundary loops |
| `decimate_flat` | Angle-based decimation of coplanar regions |
| `principal_curvatures` | Principal curvatures (K1, K2) and shape index |
| `image_cumulative_sum_x/y/z` / `image_integral` | Cumulative sums and integral image |
| `face_coloring` / `chromatic_number` | Greedy face graph coloring |
| `extract_one_ring` / `extract_n_ring` | Vertex neighborhood extraction |
| `spin_image_density` | Spin image point density descriptor |
| `image_local_contrast` | Local contrast enhancement (CLAHE-like normalization) |
| `sample_surface_uniform` | Uniform stratified surface sampling |
| `remove_high_aspect_ratio` / `extract_high_aspect_ratio` | Sliver triangle filtering |
| `coarsen` | K-d tree based mesh coarsening by vertex merging |
| `image_morphological_gradient` / `opening` / `closing` | Morphological gradient and open/close |
| `is_convex` / `convexity_defect` | Mesh convexity testing |
| `angle_defect` / `gauss_bonnet_check` | Discrete angle defect and Gauss-Bonnet verification |
| `laplacian_coordinates` / `laplacian_magnitude` | Laplacian coordinates and detail magnitude |
| `image_multi_threshold` / `image_multi_otsu` | Multi-level thresholding and segmentation |
| `boundary_distance` | BFS hop distance from boundary vertices |
| `flatten_z` / `project_to_plane` / `scale_z` | Mesh flattening and Z scaling |
| `point_cloud_to_mesh` | K-NN triangle reconstruction from point cloud |
| `image_colorize` / `image_apply_lut` | Scalar-to-RGB colormapping and lookup table |
| `random_perturb` / `gaussian_perturb` | Add random/Gaussian noise to vertex positions |
| `dual_graph` | Face adjacency dual graph as lines between centroids |
| `extract_sharp_edges` / `mark_sharp_vertices` | Sharp edge detection by dihedral angle |
| `image_tile` | Tile (repeat) ImageData in X and Y |
| `explode` | Explode cells outward from mesh centroid |
| `twist` / `bend` | Axis-based twist and bend deformations |
| `simplify_flat_vertices` | Curvature-based flat vertex removal |
| `image_blend` / `image_weighted_blend` | Alpha and weighted blending of ImageData |
| `convex_hull_2d` / `point_in_convex_hull_2d` | 2D convex hull (Andrew's monotone chain) |
| `turning_angles` | Signed turning angles at polygon vertices |
| `cotangent_smooth` | Cotangent-weighted Laplacian smoothing |
| `image_thin` | Zhang-Suen morphological thinning (skeletonization) |
| `fiedler_vector` / `spectral_partition` | Spectral mesh partitioning via Fiedler vector |
| `bvh_depth` | BVH spatial subdivision depth per cell |
| `ambient_occlusion` | Hemisphere ray-casting ambient occlusion |
| `image_profile_row` / `column` / `diagonal` | 1D profile extraction from ImageData |
| `split_sharp_edges` | Duplicate vertices at sharp edges for flat shading |
| `topology_check` | Comprehensive mesh topology validation report |
| `directed_hausdorff_colored` | Per-point Hausdorff error with max/mean/rms stats |
| `image_stack` / `image_extract_slice` | Stack 2D slices to 3D and extract slices |
| `area_weighted_normals` / `angle_weighted_normals` | Advanced vertex normal methods |
| `inertia_tensor` / `principal_axes` | Moment of inertia and principal axis computation |
| `detect_holes` | Detect boundary holes with perimeter measurement |
| `anisotropic_diffusion` | Perona-Malik edge-preserving anisotropic diffusion |
| `select_patch` | Select face patch within geodesic radius of seed |
| `random_sample_cells` / `every_nth_cell` | Random and systematic cell sampling |
| `cage_deform` | Inverse-distance weighted cage deformation |
| `image_extract_component` / `image_merge_components` | Multi-component array manipulation |
| `boundary_fill` | Exponential decay flood-fill from mesh boundary |
| `normal_deviation` | Vertex normal vs face normal deviation detection |
| `remove_small_components` | Remove small connected components by face count |
| `image_gradient_direction` / `image_gradient_orientation` | Gradient direction/orientation |
| `unwrap_to_uv` / `uv_distortion` | UV unwrap visualization and distortion metrics |
| `wave_step` | Wave equation simulation on mesh connectivity |
| `equalize_edge_lengths` / `edge_length_stats` | Edge length equalization and statistics |
| `image_apply_mask` / `image_create_mask` | Mask application and creation |
| `thicken` | Surface-to-solid thickening along normals |
| `remove_small_faces` / `remove_large_faces` | Area-based face filtering |
| `poisson_disk_sample` | Poisson-disk subsampling for well-spaced points |
| `register_translation_2d` | 2D image translation registration via cross-correlation |
| `discrete_mean_curvature` / `curvature_classification` | Mean curvature and convex/concave/flat labeling |
| `vertex_color_transfer` / `transfer_all_point_data` | K-d tree based attribute transfer between meshes |
| `mesh_voronoi_partition` / `cvt_iterate` | Voronoi partition on mesh + centroidal Voronoi tessellation |
| `image_entropy` | Local Shannon entropy for texture analysis |
| `paint_sphere` / `paint_sphere_smooth` | Interactive sphere-brush scalar painting |
| `subdivide_long_edges` | One-pass long-edge subdivision with re-triangulation |
| `label_statistics` / `label_counts` | Per-label region statistics for segmented images |
| `merge_close_vertices` | K-d tree based efficient vertex merging |
| `orient_faces_consistent` | BFS-based consistent face winding orientation |
| `harmonic_solve` | Harmonic field interpolation with boundary conditions |
| `image_structure_tensor` | Structure tensor eigenvalues for corner/edge detection |
| `intrinsic_distance` / `geodesic_distance_between` | Weighted geodesic distance |
| `edge_graph` / `vertex_degree` / `edge_count` | Edge graph extraction and analysis |
| `bilateral_mesh_smooth` | Edge-preserving bilateral mesh smoothing |
| `image_dot_product` / `image_cross_product` / `image_scale_vector` | Vector field math |
| `spring_relaxation` | Spring-mass physics relaxation |
| `point_set_intersection` / `difference` / `union` | Point set boolean operations |
| `subdivide_quads` | Quad-to-4-quads subdivision with edge and center midpoints |
| `image_granulometry` / `erosion_thickness` | Granulometry and erosion-based thickness |
| `check_normal_consistency` / `count_inconsistent_edges` | Face winding consistency check |
| `weighted_smooth` | Scalar-weighted Laplacian smoothing |
| `fix_non_manifold` | Fix non-manifold edges by vertex duplication |
| `image_harris_corners` / `harris_corner_points` | Harris corner detection |
| `surface_integral` / `surface_flux` | Scalar integral and vector flux over surface |
| `region_grow` | Scalar-constrained region growing from seed vertices |
| `geodesic_voronoi` | Geodesic Voronoi partition via edge-weighted Dijkstra |
| `image_lbp` / `lbp_uniformity` | Local Binary Pattern texture descriptor |
| `split_vertex` | Split vertex into per-face-group copies |
| `detect_outliers` / `remove_outliers` | Statistical k-NN outlier detection and removal |
| `compute_frame_field` | Tangent/bitangent/normal orthonormal frame per vertex |
| `image_label_boundary` / `label_contact_area` | Label region boundaries and inter-label contact |
| `compute_vertex_rings` / `detect_irregular_vertices` | Vertex ring analysis and irregular detection |
| `select_faces_in_box` / `select_faces_in_sphere` | Spatial face selection by centroid |
| `shape_diameter_function` | Inward ray-casting thickness estimation |
| `image_histogram_match` | Histogram specification (match to reference) |
| `tutte_parameterize` | Tutte embedding for disk-topology meshes |
| `cell_to_point_area_weighted` | Area-weighted cell-to-point data interpolation |
| `adjacency_matrix` / `betweenness_centrality` | Graph algorithms on mesh connectivity |
| `image_gabor` | Gabor filter for oriented texture analysis |
| `conformal_factor` / `angle_distortion` | Mesh deformation quality metrics |
| `tri_to_quad` | Greedy triangle-pair to quad conversion |
| `minimum_spanning_tree` / `mst_weight` | Kruskal MST on mesh edge graph |
| `image_largest_component` / `image_component_sizes` | Binary image largest component extraction |
| `equalize_areas` / `area_variance` | Triangle area equalization |
| `scalar_field_critical_points` / `count_critical_points` | Min/max/saddle detection on scalar mesh |
| `mesh_level_set` | Isocontour extraction on triangle mesh scalar field |
| `image_resize_nearest` / `image_resize_by_factor` | Nearest-neighbor image resizing |
| `label_face_groups` | Edge-connected face component labeling |
| `circumradius` | Triangle circumradius and circumradius/inradius ratio |
| `spectral_descriptor` | Multi-scale heat diffusion shape descriptor |
| `image_template_match` | NCC-based template matching |
| `midpoint_refine` | Pure midpoint subdivision (tri→4, quad→4) |
| `bilateral_scalar_smooth` | Edge-preserving bilateral scalar smoothing on mesh |
| `mesh_saliency` | Multi-scale curvature saliency detection |
| `image_chamfer_distance` | 3-4-5 Borgefors chamfer distance transform |
| `has_self_intersection` / `count_self_intersections` | Triangle mesh self-intersection detection |
| `vertex_importance` | Multi-factor vertex importance for simplification |
| `displacement_field` / `cell_strain` | Mesh deformation displacement and strain metrics |
| `image_max_pool` / `image_min_pool` / `image_avg_pool` | Pooling operations on ImageData |
| `orient_point_cloud_normals` | BFS-based consistent normal orientation for point clouds |
| `cell_adjacency_list` / `face_graph_diameter` | Cell adjacency analysis |
| `heat_method_distance` | Heat method approximate geodesic distance |
| `image_pyramid` / `laplacian_pyramid` | Gaussian and Laplacian image pyramids |
| `point_pair_features` | PPF-based shape descriptor for point clouds |
| `displace_by_scalar` / `displace_by_vector` | Scalar/vector displacement of mesh vertices |
| `project_to_axis_plane` / `project_to_sphere` / `project_to_cylinder` | Geometric projections |
| `image_canny` | Approximate Canny edge detection with hysteresis |
| `extract_ridge_valley_lines` | Ridge and valley line extraction from dihedral angles |
| `topological_distance` / `eccentricity` | Hop-count distance and graph eccentricity |
| `dirichlet_energy` / `willmore_energy` / `smoothness_ratio` | Mesh energy metrics |
| `region_properties` | Per-label region props (area, centroid, bounding box, mean) |
| `aspect_ratio_histogram` / `area_histogram` | Mesh quality distribution histograms |
| `watertight_check` | Comprehensive watertight mesh validation |
| `geodesic_iso_contours` | Iso-distance contour extraction from geodesic field |
| `image_ssim` | Structural Similarity Index (SSIM) between images |
| `surface_area_density` / `area_statistics` | Per-vertex Voronoi area and statistics |
| `closest_pair` / `farthest_pair` / `point_set_diameter` | Extremal point pair queries |
| `tangential_smooth` | Tangent-plane-only Laplacian smoothing |
| `image_psnr` / `image_mae` | PSNR and MAE image quality metrics |
| `face_cluster_by_normal` | Normal-based face clustering via region growing |
| `depth_sort` / `depth_sort_front_to_back` | Depth-ordered face sorting for rendering |
| `arap_deform` | As-Rigid-As-Possible mesh deformation |
| `image_knn_classify` | K-nearest-neighbor scalar classification |
| `mesh_morph` / `mesh_morph_sequence` | Mesh morphing with easing functions |
| `array_statistics` | Comprehensive scalar statistics (mean/median/std/skew/IQR) |
| `select_cells_by_predicate` / `select_cells_any_vertex` | Predicate-based cell selection |
| `image_correlation` / `image_covariance` | Pearson correlation and covariance between images |
| `random_walk_distribution` / `random_walk_path` | Random walk simulation on mesh graph |
| `pca_align` | PCA-based point cloud alignment |
| `curvature_histogram` / `curvature_percentiles` | Curvature distribution analysis |
| `extract_slice_along_axis` / `max_intensity_projection` | Axis-aligned slice and MIP |
| `vertex_mask_to_cell_mask` / `cell_mask_to_vertex_mask` | Mask conversion between vertex/cell |
| `multi_contour_on_mesh` / `contour_length` | Multi-isocontour extraction with length |
| `procrustes_align` / `procrustes_distance` | Procrustes shape alignment and distance |
| `image_gamma` / `image_sigmoid` | Gamma correction and sigmoid contrast |
| `subdivide_ternary` | 1-to-9 ternary triangle subdivision |
| `scalar_gradient_on_mesh` | Per-face gradient averaged to vertices |
| `quality_improve` | Edge-length-weighted mesh quality improvement |
| `image_percentile_filter` | Generalized percentile filter (min/median/max) |
| `face_normal_angle` / `face_normal_dot` | Face normal vs direction angle/dot product |
| `neighborhood_stats` | Per-vertex min/max/mean/range of scalar neighborhood |
| `difference_of_gaussians_mesh` / `detect_scale_space_features` | DoG and scale-space features |
| `mesh_stats_report` / `is_all_triangles` / `is_all_quads` | Mesh info and type queries |
| `uv_checkerboard` / `uv_seam_length` | UV quality visualization and seam analysis |
| `image_dog` / `detect_blobs` | Difference of Gaussians blob detection on ImageData |
| `point_in_mesh` / `points_in_mesh` / `count_points_inside` | Ray-casting point-in-mesh test |
| `connected_threshold` | Connected threshold region growing segmentation |
| `convex_decompose` | Approximate convex decomposition by patch growing |
| `wireframe` / `boundary_wireframe` / `internal_wireframe` | Edge extraction variants |
| `laplacian_eigenvalues` / `spectral_gap` | Laplacian spectral analysis |
| `component_volumes_3d` / `filter_by_volume` | 3D connected component volume analysis |
| `kmeans_cluster` / `kmeans_inertia` | K-means vertex clustering with inertia metric |
| `signed_angle_field` | Direction field via tangent-plane angle projection |
| `collapse_short_edges_kdtree` | K-d tree accelerated short edge collapse |
| `image_non_local_means` | Non-local means denoising via patch similarity |
| `vertex_star_info` / `count_boundary_vertices` | Vertex star degree and boundary classification |
| `scatter_plot_3d` / `histogram_2d` | Data visualization: scatter plots and 2D histograms |
| `geodesic_farthest_point_sampling` / `coverage_radius` | Geodesic FPS and coverage metric |
| `image_quantize` / `image_dither_quantize` | Scalar quantization and dithered quantization |
| `vertices_only` / `cell_centroids_as_points` / `edge_midpoints` | Point extraction variants |
| `curvature_color` | Signed-curvature diverging RGB colormap |
| `fix_t_junctions` | T-junction repair by edge splitting |
| `image_local_range` / `image_local_std` | Local range and standard deviation filters |
| `heat_trace_descriptor` | Multi-scale heat trace shape descriptor |
| `triangles_to_greedy_strips` / `count_strips` | Greedy triangle strip generation |
| `dual_contour_mesh` | Topological dual mesh with actual polygon cells |
| `image_rle_encode` / `count_unique_values` | Run-length encoding and value analysis |
| `convex_layers_2d` / `num_convex_layers` | Convex onion peeling layers |
| `planar_cross_section` / `cross_section_area` | Plane-mesh intersection and area |
| `half_edge_valence` / `detect_non_manifold_vertices` | Half-edge analysis and non-manifold detection |
| `integral_image_2d` / `rect_sum` | Summed area table and O(1) rectangular queries |
| `spin_axis` | Rotational symmetry axis detection via PCA |
| `scatter_from_centroid` / `jitter_vertices` / `shrink_mesh` | Vertex displacement utilities |
| `aspect_ratio_report` | Comprehensive aspect ratio quality analysis |
| `image_cosine_similarity` / `image_euclidean_distance` / `image_l1_distance` | Image similarity metrics |
| `total_absolute_curvature` / `mean_total_curvature` | Curvature integral metrics |
| `dihedral_angle_stats` / `smoothness_index` / `flatness_score` | Dihedral angle analysis |
| `vertex_angle_extremes` / `mesh_angle_range` | Per-vertex min/max angle and global angle range |
| `recompute_normals` / `flip_normals` | Area-weighted normal recomputation and flipping |
| `image_local_mean` / `image_z_score` | Local mean and z-score outlier detection on ImageData |
| `cell_from_point_stats` | Per-cell min/max/mean/range from point data |
| `compute_edge_valence` / `edge_valence_stats` | Per-vertex edge valence and mesh-wide stats |
| `extract_boundary_edges` / `count_boundary_edges` | Boundary edge extraction and counting |
| `compute_face_areas` / `face_area_stats` | Per-face area computation and statistics |
| `image_gradient_magnitude` | Gradient magnitude via central differences on ImageData |
| `compute_vertex_neighbor_counts` / `n_ring_neighborhood` | Vertex neighbor counting and N-ring extraction |
| `compute_dihedral_angles` | Dihedral angle per interior edge as line PolyData |
| `image_laplacian` (filters) | Discrete Laplacian of scalar field on ImageData |
| `compute_centroid` / `compute_area_weighted_centroid` | Mesh centroid (simple and area-weighted) |
| `compute_aspect_ratio` | Per-triangle aspect ratio (circumradius / 2*inradius) |
| `non_manifold_edges` / `count_non_manifold_edges` | Non-manifold edge detection |
| `image_threshold_binary` / `image_threshold_binary_inverse` | Binary thresholding on ImageData |
| `smooth_laplacian_simple` | Iterative Laplacian smoothing with boundary preservation |
| `compute_mean_curvature` / `mean_curvature_stats` | Mean curvature via cotangent Laplacian |
| `check_normal_consistency` / `add_consistency_cell_data` | Face normal orientation consistency check |
| `compute_distance_to_surface` | Euclidean distance from ImageData voxels to PolyData |
| `subdivide_loop` | Loop subdivision (1:4 triangle split) |
| `compute_signed_volume` / `compute_volume` | Closed mesh volume via divergence theorem |
| `edge_length_stats` / `extract_edges_with_length` | Per-edge length computation and statistics |
| `binary_dilate` / `binary_erode` | Binary morphology on ImageData (3x3x3 structuring element) |
| `compute_gaussian_curvature` | Gaussian curvature via angle defect |
| `count_components` | Connected component counting with optional labeling |
| `point_to_cell_average` | Average point data to cells |
| `resample_nearest` | Nearest-neighbor ImageData resampling |
| `decimate_vertex_cluster` | Vertex clustering decimation |
| `compute_face_normals` | Per-face normal computation as cell data |
| `cell_to_point_average` | Average cell data to points |
| `mask_by_scalar_range` | Scalar-range masking on ImageData |
| `scale_uniform` / `scale_nonuniform` | Mesh scaling around center point |
| `translate` / `translate_to_origin` | Mesh translation and centering |
| `rotate_axis_angle` | Rodrigues' rotation around axis through center |
| `extract_slice_z` / `extract_slice_y` / `extract_slice_x` | 2D slice extraction from 3D ImageData |
| `remove_unused_points` | Remove unreferenced points and reindex cells |
| `barycentric_subdivide` | Barycentric subdivision (1:6 triangle split) |
| `select_faces_by_normal` | Select faces by normal direction alignment |
| `remap_range` | Linear scalar range remapping on ImageData |
| `compute_surface_area` | Total surface area computation |
| `clip_by_plane` | Plane clipping with triangle splitting |
| `mesh_distance_stats` | Symmetric Hausdorff-like distance between meshes |
| `image_gradient_vector` | 3-component gradient vector on ImageData |
| `merge_poly_data` | Concatenate two PolyData meshes |
| `ray_intersect_mesh` | Möller-Trumbore ray-mesh intersection |
| `points_to_vertices` / `vertices_to_points` | Point cloud ↔ vertex PolyData conversion |
| `max_projection_z` / `min_projection_z` | Maximum/minimum intensity projection |
| `triangulate_quads` | Quad-to-triangle conversion (shortest diagonal) |
| `merge_coplanar_faces` | Coplanar triangle merging |
| `color_edges_by_angle` | Edge coloring by dihedral angle |
| `compute_vertex_normals_weighted` | Vertex normals with area/angle/uniform weighting |
| `extract_largest_component` | Extract largest connected component |
| `compute_displacement` | Per-vertex displacement between two meshes |
| `extract_component` (image) | Extract single component from multi-component array |
| `midpoint_split` | Edge midpoint subdivision (1:4 triangle split) |
| `compute_signed_volume` / `compute_volume` | Mesh volume via divergence theorem |
| `compute_gaussian_curvature` | Gaussian curvature via angle defect |
| `count_components` | Connected component counting and labeling |
| `decimate_vertex_cluster` | Vertex clustering decimation |

#### Mesh Boolean & Clipping (additional)
| `poly_data_boolean_2d` | 2D polygon boolean operations |
| `mesh_boolean_exact` | Exact boolean operations on PolyData |
| `mesh_boolean_difference` | Boolean difference between meshes |
| `mesh_boolean_union_approx` | Approximate boolean union |
| `mesh_boolean_cells` | Cell-level boolean operations |
| `mesh_boolean_point_set` | Point set boolean operations (intersection/difference/union) |
| `mesh_clip_by_plane` | Plane clipping with triangle splitting |
| `mesh_cylinder_clip` | Clip mesh by cylinder volume |
| `mesh_sphere_clip` | Clip mesh by sphere volume |

#### Mesh Smoothing (additional)
| `mesh_smooth_taubin` | Volume-preserving Taubin smoothing |
| `mesh_smoothing_cotangent` | Cotangent-weighted Laplacian smoothing |
| `mesh_smooth_laplacian_simple` | Iterative Laplacian smoothing with boundary preservation |
| `mesh_smooth_cotangent` | Cotangent Laplacian mesh smoothing |
| `mesh_smooth_bilaplacian` | Bi-Laplacian (4th order) mesh smoothing |
| `mesh_smooth_hc` | HC (Humphrey's Classes) Laplacian smoothing |
| `mesh_bilateral_smooth` | Edge-preserving bilateral mesh smoothing |
| `mesh_weighted_smooth` | Scalar-weighted Laplacian smoothing |
| `mesh_tangent_smooth` | Tangent-plane-only Laplacian smoothing |
| `mesh_scalar_smoothing` | Scalar field Laplacian smoothing on mesh |

#### Mesh Subdivision (additional)
| `mesh_subdivision_sqrt3` | Sqrt(3) subdivision scheme |
| `mesh_subdivision_ternary` | 1-to-9 ternary triangle subdivision |
| `mesh_subdivide_quad` | Quad-to-4-quads subdivision |
| `mesh_subdivide_loop` | Loop subdivision (1:4 triangle split) |
| `mesh_subdivide_edge` | One-pass long-edge subdivision |
| `mesh_midpoint_refine` | Pure midpoint subdivision (tri→4, quad→4) |
| `mesh_barycentric_subdivide` | Barycentric subdivision (1:6 triangle split) |
| `mesh_edge_midpoint_split` | Edge midpoint subdivision (1:4 triangle split) |

#### Mesh Decimation & Simplification (additional)
| `mesh_decimate_angle` | Angle-based decimation of coplanar regions |
| `mesh_decimate_vertex_clustering` | Vertex clustering decimation |
| `mesh_decimate_quadric_error` | Garland-Heckbert quadric error metric decimation |
| `mesh_simplify_vertex_removal` | Curvature-based flat vertex removal |
| `mesh_coarsen` | K-d tree based mesh coarsening |
| `mesh_collapse_short_edges` | K-d tree accelerated short edge collapse |
| `mesh_aspect_ratio_improve` | Edge-length-weighted mesh quality improvement |

#### Mesh Curvature (additional)
| `mesh_curvature_tensor` | Curvature tensor per vertex |
| `mesh_curvature_mean_simple` | Mean curvature via cotangent Laplacian |
| `mesh_vertex_curvature_gaussian` | Gaussian curvature via angle defect |
| `mesh_principal_curvatures` | Principal curvatures (K1, K2) and shape index |
| `mesh_curvature_line` | Curvature line extraction |
| `mesh_curvature_histogram` | Curvature distribution analysis |
| `mesh_curvature_color` | Signed-curvature diverging RGB colormap |
| `mesh_curvature_color_map` | Curvature-to-color mapping |
| `mesh_curvature_map_to_color` | Map curvature values to vertex colors |
| `mesh_surface_curvature_integral` | Total absolute/mean curvature integrals |

#### Mesh Analysis & Queries (additional)
| `mesh_info` | Mesh topology and geometry summary statistics |
| `mesh_measure` | Comprehensive mesh measurements (area, edges, cell types) |
| `mesh_export_stats` | Mesh info and type queries |
| `mesh_aspect_ratio_stats` | Comprehensive aspect ratio quality analysis |
| `mesh_aspect_ratio_compute` | Per-triangle aspect ratio computation |
| `mesh_face_area_stats` | Per-face area computation and statistics |
| `mesh_face_data_stats` | Per-cell min/max/mean/range from point data |
| `mesh_edge_angle_stats` | Dihedral angle analysis and statistics |
| `mesh_edge_valence_stats` | Per-vertex edge valence and mesh-wide stats |
| `mesh_edge_length_stats` | Per-edge length computation and statistics |
| `mesh_vertex_neighbors` | Vertex neighbor counting and N-ring extraction |
| `mesh_vertex_star` | Vertex star degree and boundary classification |
| `mesh_dihedral_angle_array` | Dihedral angle per interior edge as line PolyData |
| `mesh_boundary_edges` | Boundary edge extraction and counting |
| `mesh_face_pair_angle` | Face normal vs direction angle/dot product |
| `mesh_face_normal_array` | Per-face normal computation as cell data |
| `mesh_face_centroid_array` | Per-face centroid computation as cell data |
| `mesh_face_centroid_filter` | Cell centroids and weighted center of mass |
| `mesh_centroid` | Mesh centroid (simple and area-weighted) |
| `mesh_non_manifold_edges` | Non-manifold edge detection |
| `mesh_point_data_range` | Point data scalar range queries |
| `mesh_face_skewness` | Per-cell face skewness metrics |
| `mesh_volume_compute` | Closed mesh volume via divergence theorem |
| `mesh_surface_area_compute` | Total surface area computation |
| `mesh_surface_area_density` | Per-vertex Voronoi area and statistics |
| `mesh_connected_components_count` | Connected component counting and labeling |
| `mesh_connected_face_groups` | Edge-connected face component labeling |
| `mesh_quality_histogram` | Aspect ratio and area distribution histograms |
| `mesh_circumradius` | Triangle circumradius and circumradius/inradius ratio |
| `mesh_face_normal_consistency` | Face normal orientation consistency check |
| `mesh_normal_consistency_check` | Face winding consistency check |
| `mesh_half_edge` | Half-edge valence and non-manifold detection |
| `mesh_cell_connectivity_matrix` | Cell adjacency analysis |
| `mesh_ring_buffer` | Vertex ring analysis and irregular detection |
| `mesh_topology_check` | Comprehensive mesh topology validation report |
| `mesh_topology_distance` | Hop-count distance and graph eccentricity |
| `mesh_topology_genus` | Euler characteristic and genus computation |
| `mesh_watertight_check` | Comprehensive watertight mesh validation |
| `mesh_convex_check` | Mesh convexity testing |
| `mesh_self_intersection` | Triangle mesh self-intersection detection |
| `mesh_neighborhood_stats` | Per-vertex min/max/mean/range of scalar neighborhood |
| `mesh_array_statistics` | Comprehensive scalar statistics (mean/median/std/skew/IQR) |

#### Mesh Topology & Repair (additional)
| `mesh_repair` | Remove duplicate/degenerate cells and unused points |
| `mesh_topology_repair` | Fix non-manifold edges by vertex duplication |
| `mesh_topology_simplify` | Topology simplification operations |
| `mesh_topological_noise` | Remove small connected components by face count |
| `mesh_remove_unused_points` | Remove unreferenced points and reindex cells |
| `mesh_orient_faces` | BFS-based consistent face winding orientation |
| `mesh_invert_faces` | Flip polygon winding and normals |
| `mesh_weld_vertices` | Tolerance-based vertex welding |
| `mesh_vertex_merge_by_distance` | K-d tree based efficient vertex merging |
| `mesh_duplicate_faces` | Detect and remove duplicate faces |
| `mesh_edge_swap` | Delaunay-criterion edge swapping for quality improvement |
| `mesh_edge_split` | Duplicate vertices at sharp edges for flat shading |
| `mesh_vertex_split` | Split vertex into per-face-group copies |

#### Mesh Geodesic & Distance (additional)
| `mesh_geodesic_path` | Shortest path on mesh via Dijkstra |
| `mesh_geodesic_distance_dijkstra` | Dijkstra-based geodesic distance from seed points |
| `mesh_geodesic_distance_heat` | Heat method approximate geodesic distance |
| `mesh_geodesic_voronoi` | Geodesic Voronoi partition via edge-weighted Dijkstra |
| `mesh_geodesic_farthest_point` | Geodesic farthest point sampling and coverage |
| `mesh_geodesic_heat_kernel` | Multi-scale heat diffusion shape descriptor |
| `mesh_heat_geodesic_iso` | Iso-distance contour extraction from geodesic field |
| `mesh_heat_method_distance` | Heat method approximate geodesic distance |
| `mesh_distance_field` | Signed distance field from mesh surface |
| `mesh_signed_distance_field` | Signed distance from closed surface on grid |
| `mesh_distance_between_meshes` | Symmetric Hausdorff-like distance between meshes |
| `mesh_intrinsic_distance` | Weighted geodesic distance |
| `mesh_closest_pair` | Closest and farthest point pair queries |
| `mesh_closest_point_on_surface` | Project points onto nearest surface location |

#### Mesh Deformation (additional)
| `mesh_twist` | Axis-based twist deformation |
| `mesh_bend` | Axis-based bend deformation |
| `mesh_taper` | Taper deformation along axis |
| `mesh_wave_deform` | Wave-based deformation |
| `mesh_explode` | Explode cells outward from mesh centroid |
| `mesh_flatten` | Mesh flattening and Z scaling |
| `mesh_project_to_plane` | Project mesh onto a plane |
| `mesh_random_perturb` | Add random/Gaussian noise to vertex positions |
| `mesh_noise_displacement` | Noise-based vertex displacement |
| `mesh_vertex_displacement` | Scalar/vector displacement of mesh vertices |
| `mesh_vertex_displacement_array` | Per-vertex displacement between two meshes |
| `mesh_vertex_scatter` | Scatter cells outward from centroid |
| `mesh_cage_deform` | Inverse-distance weighted cage deformation |
| `mesh_laplacian_deform_simple` | Handle-based Laplacian mesh deformation |
| `mesh_arap_deform` | As-Rigid-As-Possible mesh deformation |
| `mesh_spring_relaxation` | Spring-mass physics relaxation |
| `mesh_wave_simulation` | Wave equation simulation on mesh connectivity |
| `mesh_heat_equation` | Heat equation simulation on mesh connectivity |

#### Mesh Spatial & Transform (additional)
| `mesh_scale` | Mesh scaling around center point |
| `mesh_translate` | Mesh translation and centering |
| `mesh_rotate` | Rodrigues' rotation around axis through center |
| `mesh_center` | Mesh centering, normalization, PCA alignment |
| `mesh_merge` | Concatenate two PolyData meshes |
| `mesh_offset` | Offset mesh along vertex normals by distance |
| `mesh_thicken` | Surface-to-solid thickening along normals |
| `mesh_extrude_along_normal` | Extrude surface along vertex normals |
| `mesh_mirror` | Mirror mesh across coordinate plane |
| `mesh_axis_projection` | Geometric projections (axis plane, sphere, cylinder) |
| `mesh_axis_aligned_slice` | Axis-aligned mesh slicing |
| `mesh_pca_axes` | PCA-based principal axis computation |
| `mesh_obb` | Oriented bounding box computation |
| `mesh_spin_axis` | Rotational symmetry axis detection via PCA |

#### Mesh Parameterization & UV (additional)
| `mesh_parameterize` | Planar and cylindrical UV parameterization |
| `mesh_parameterize_tutte` | Tutte embedding for disk-topology meshes |
| `mesh_texture_atlas` | Per-triangle texture atlas UV parameterization |
| `mesh_unwrap` | UV unwrap visualization and distortion metrics |
| `mesh_uv_checker` | UV quality visualization and seam analysis |
| `mesh_uv_from_projection` | UV coordinates from planar/spherical projection |
| `mesh_conformal_factor` | Mesh deformation quality metrics (angle distortion) |

#### Mesh Descriptors & Features (additional)
| `mesh_heat_kernel` | Heat Kernel Signature shape descriptor |
| `mesh_spectral_descriptor` | Multi-scale heat diffusion shape descriptor |
| `mesh_point_pair_features` | PPF-based shape descriptor for point clouds |
| `mesh_shape_diameter` | Inward ray-casting thickness estimation |
| `mesh_saliency` | Multi-scale curvature saliency detection |
| `mesh_feature_size` | Point cloud feature size and spacing |
| `mesh_feature_line` | Ridge and valley line extraction |
| `mesh_feature_angle_detect` | Feature angle detection by dihedral angle |
| `mesh_extract_sharp_features` | Sharp edge detection by dihedral angle |
| `mesh_scale_space` | Difference of Gaussians and scale-space features |
| `mesh_vertex_importance` | Multi-factor vertex importance for simplification |

#### Mesh Coloring & Attributes (additional)
| `mesh_face_color` | Greedy face graph coloring |
| `mesh_vertex_color_transfer` | K-d tree based attribute transfer between meshes |
| `mesh_vertex_color_from_position` | Color vertices by position |
| `mesh_edge_color` | Edge coloring by dihedral angle |
| `mesh_point_data_interpolate` | Point data interpolation between meshes |
| `mesh_point_data_gradient` | Per-face gradient averaged to vertices |
| `mesh_scalar_gradient_on_mesh` | Scalar gradient computation on mesh |
| `mesh_point_to_cell_data` | Average point data to cells |
| `mesh_cell_to_point_data` | Average cell data to points |
| `mesh_cell_to_point_interpolate` | Area-weighted cell-to-point interpolation |
| `mesh_vertex_to_cell` | Vertex mask to cell mask conversion |
| `mesh_face_normal_angle` | Face normal vs direction angle/dot product |
| `mesh_select_by_normal` | Select faces by normal direction alignment |
| `mesh_select_cells_by_area` | Area-based cell selection |

#### Mesh Point Cloud (additional)
| `mesh_point_cloud_to_mesh` | K-NN triangle reconstruction from point cloud |
| `mesh_point_cloud_to_vertices` | Point cloud to vertex PolyData conversion |
| `mesh_point_cloud_from_mesh` | Extract point cloud from mesh |
| `mesh_point_cloud_normal_estimation` | PCA-based normal estimation for point clouds |
| `mesh_point_cloud_normals_orient` | BFS-based consistent normal orientation |
| `mesh_point_cloud_outlier` | Statistical k-NN outlier detection and removal |
| `mesh_point_cloud_align` | PCA-based point cloud alignment |
| `mesh_vertex_cluster` | Farthest-point-sampling vertex clustering |
| `mesh_vertex_cluster_kmeans` | K-means vertex clustering with inertia metric |
| `mesh_vertex_cluster_by_normal` | Normal-based vertex clustering |

#### Mesh Partitioning & Segmentation (additional)
| `mesh_face_cluster` | Normal-based face clustering via region growing |
| `mesh_face_group_by_normal` | Group faces by normal similarity |
| `mesh_region_grow` | Scalar-constrained region growing |
| `mesh_voronoi_cells` | Voronoi partition on mesh + CVT iteration |
| `mesh_split_by_connectivity` | Split mesh by connected components |
| `mesh_convex_decompose` | Approximate convex decomposition |
| `mesh_convex_layers` | Convex onion peeling layers (2D) |
| `mesh_convex_hull_2d` | 2D convex hull (Andrew's monotone chain) |
| `mesh_convex_hull_3d` | 3D convex hull computation |

#### Mesh Ray Casting & Picking (additional)
| `mesh_ray_cast` | Möller–Trumbore ray-mesh intersection |
| `mesh_ray_triangle_intersect` | Single ray-triangle intersection test |
| `mesh_pick` | Point/cell picking via closest-point queries |
| `mesh_point_in_mesh` | Ray-casting point-in-mesh test |
| `mesh_visibility` | Hemisphere ray-casting ambient occlusion |
| `mesh_ambient_occlusion` | Per-vertex ambient occlusion |

#### Mesh Graph & Spectral (additional)
| `mesh_edge_graph` | Edge graph extraction and analysis |
| `mesh_abstract_graph` | Graph algorithms on mesh connectivity |
| `mesh_laplacian_matrix` | Laplacian matrix computation |
| `mesh_laplacian_eigenmaps` | Spectral mesh analysis via Laplacian eigenvalues |
| `mesh_laplacian_spectrum` | Laplacian spectral analysis and spectral gap |
| `mesh_laplacian_energy` | Dirichlet/Willmore energy and smoothness metrics |
| `mesh_minimum_spanning_tree` | Kruskal MST on mesh edge graph |
| `mesh_random_walk` | Random walk simulation on mesh graph |
| `mesh_harmonic_map` | Harmonic field interpolation with boundary conditions |
| `mesh_scalar_field_topology` | Min/max/saddle detection on scalar mesh |
| `mesh_signed_angle_field` | Direction field via tangent-plane angle projection |

#### Mesh Morphing & Registration (additional)
| `mesh_morph` | Mesh morphing with easing functions |
| `mesh_procrustes` | Procrustes shape alignment and distance |
| `mesh_correspondence` | Point set correspondence and attribute transfer |
| `mesh_symmetry` | Bilateral symmetry detection |
| `mesh_symmetry_plane` | Best symmetry plane computation |
| `mesh_spherical_harmonic_decompose` | Spherical harmonic decomposition |

#### Mesh Extraction & Selection (additional)
| `mesh_boundary_loop` | Extract and count boundary loops |
| `mesh_extract_boundary_loop` | Extract boundary loops as polylines |
| `mesh_extract_largest_component` | Extract largest connected component |
| `mesh_extract_sharp_features` | Extract sharp edges and features |
| `mesh_one_ring` | Vertex neighborhood extraction |
| `mesh_patch_select` | Select face patch within geodesic radius |
| `mesh_random_sample_cells` | Random and systematic cell sampling |
| `mesh_sample_surface` | Uniform stratified surface sampling |
| `mesh_poisson_disk` | Poisson-disk subsampling for well-spaced points |
| `mesh_poisson_sample` | Poisson-disk surface sampling |
| `mesh_wireframe_extract` | Edge extraction variants (boundary/internal) |

#### Mesh Miscellaneous (additional)
| `mesh_depth_sort` | Depth-ordered face sorting for rendering |
| `mesh_silhouette` | Silhouette edge extraction |
| `mesh_border_distance` | BFS hop distance from boundary vertices |
| `mesh_boundary_fill` | Exponential decay flood-fill from boundary |
| `mesh_vertex_normal_deviation` | Vertex normal vs face normal deviation |
| `mesh_angle_defect` | Discrete angle defect and Gauss-Bonnet verification |
| `mesh_surface_integral` | Scalar integral and vector flux over surface |
| `mesh_frame_field` | Tangent/bitangent/normal frame per vertex |
| `mesh_tensor_field` | Tensor field computation on mesh |
| `mesh_point_normals_recompute` | Area-weighted normal recomputation and flipping |
| `mesh_area_weighted_normals` | Advanced area/angle-weighted vertex normals |
| `mesh_vertex_normal_weight` | Vertex normals with configurable weighting |
| `mesh_inertia` | Moment of inertia and principal axes |
| `mesh_hole_detect` | Detect boundary holes with perimeter measurement |
| `mesh_skeleton` | 2D medial axis via Voronoi dual of Delaunay |
| `mesh_skeleton_extract` | Skeleton curve extraction |
| `mesh_voxelize` | Convert mesh surface to voxel ImageData |
| `mesh_resample_on_grid` | Resample mesh data onto regular grid |
| `mesh_iso_line` | Isocontour extraction on triangle mesh |
| `mesh_contour_on_mesh` | Multi-isocontour extraction with length |
| `mesh_planar_section` | Plane-mesh intersection and area |
| `mesh_level_set` | Level set isocontour extraction |
| `mesh_thinning` | Medial axis thinning |
| `mesh_marching_squares` | 2D marching squares contour extraction |
| `mesh_delaunay_triangulate_2d` | 2D Delaunay triangulation variant |
| `mesh_voronoi_diagram_2d` | 2D Voronoi diagram variant |
| `mesh_quad_remesh` | Greedy triangle-pair to quad conversion |
| `mesh_remesh_isotropic` | Isotropic remeshing via iterative edge operations |
| `mesh_dual_contour` | Topological dual mesh with polygon cells |
| `mesh_inscribed_sphere` | Inscribed sphere computation |
| `mesh_circumscribed_sphere` | Circumscribed sphere computation |
| `mesh_sdf_from_points` | SDF from oriented point cloud |
| `mesh_area_equalize` | Triangle area equalization |
| `mesh_edge_length_equalize` | Edge length equalization |
| `mesh_paint` | Interactive sphere-brush scalar painting |
| `mesh_array_calculator` | Expression-based field computation on mesh |
| `mesh_scatter_plot` | 3D scatter plot and 2D histogram data visualization |
| `mesh_spin` | Spin image point density descriptor |
| `mesh_aabb_tree` | BVH spatial subdivision |
| `mesh_hausdorff_directed` | Per-point Hausdorff error with max/mean/rms stats |
| `mesh_triangulate_quads` | Quad-to-triangle conversion (shortest diagonal) |
| `mesh_coplanar_merge` | Coplanar triangle merging |
| `mesh_triangle_strip_convert` | Convert between triangle strips and triangles |
| `mesh_aspect_ratio_filter` | Sliver triangle filtering by aspect ratio |
| `mesh_face_area_filter` | Area-based face filtering |
| `mesh_edge_flip` | Delaunay-criterion edge flipping |
| `mesh_vertex_to_polydata` | Point extraction variants |

#### Image Processing (additional)
| `attribute_convert` | Attribute type conversion |
| `decimate_boundary` | Boundary-preserving decimation |
| `image_clamp_range` | Clamp ImageData values to range |
| `image_min_max` | Element-wise min/max/abs-diff of ImageData |
| `image_fft` | DFT-based power spectrum of ImageData |
| `image_bilateral_3d` | 3D bilateral filtering |
| `image_bilateral_filter` | Edge-preserving bilateral filtering |
| `image_roi` | Region-of-interest masking (sphere/box) |
| `image_binarize` | Adaptive and global binarization |
| `image_threshold_multi` | Multi-level thresholding and segmentation |
| `image_threshold_binary` | Binary thresholding on ImageData |
| `image_local_extrema` | Local extrema detection in ImageData |
| `image_local_extrema_detect` | Local maxima/minima point extraction |
| `image_local_contrast` | Local contrast enhancement (CLAHE-like) |
| `image_local_range` | Local range and standard deviation filters |
| `image_morphological_gradient` | Morphological gradient |
| `image_morphological_skeleton` | Morphological skeletonization |
| `image_opening_closing` | Morphological opening and closing |
| `image_top_hat` | Top-hat morphological transform |
| `image_erode_binary` | Binary morphological erosion |
| `image_erode_dilate_binary` | Binary dilation and erosion |
| `image_mean_filter` | Mean/average filtering |
| `image_variance_filter` | Local variance filtering |
| `image_gaussian_blur` | Gaussian blur convolution |
| `image_sharpen` | Unsharp masking |
| `image_sobel_3d` | 3D Sobel edge detection |
| `image_convolve_separable` | Separable 3D convolution |
| `image_stencil` | Custom stencil/convolution |
| `image_anisotropic_diffusion` | Perona-Malik edge-preserving diffusion |
| `image_edge_preserve_smooth` | Edge-preserving smoothing |
| `image_watershed` | Marker-based watershed segmentation |
| `image_watershed_simple` | Simple watershed segmentation |
| `image_slic` | SLIC superpixel segmentation |
| `image_connected_threshold` | Connected threshold region growing |
| `image_connected_components_3d` | 3D connected component labeling |
| `image_connected_components_count` | Connected component counting |
| `image_connected_filter` | Connected component filtering (largest/by size) |
| `image_flood_fill` | Flood fill on ImageData |
| `image_island_remove` | Remove small connected components from binary image |
| `image_canny_approx` | Approximate Canny edge detection |
| `image_difference_of_gaussians` | Difference of Gaussians blob detection |
| `image_distance_transform_chamfer` | Borgefors chamfer distance transform |
| `image_distance_to_surface` | Distance from voxels to PolyData surface |
| `image_chamfer_distance` | 3-4-5 Borgefors chamfer distance |
| `image_moments` | Spatial moments (center of mass, variance) |
| `image_cumulative` | Cumulative sums and integral image |
| `image_integral_invariant` | Integral invariant computations |
| `image_colorize` | Scalar-to-RGB colormapping and lookup table |
| `image_compose` | Vector field composition/decomposition |
| `image_channel_math` | Channel-wise math operations |
| `image_magnitude` | Vector magnitude computation |
| `image_phase_congruency` | Phase congruency feature detection |
| `image_feature_extraction` | Harris corner and feature detection |
| `image_histogram_match` | Histogram specification (match to reference) |
| `image_histogram_compute` | ImageData scalar histogram computation |
| `image_gabor` | Gabor filter for oriented texture analysis |
| `image_template_match` | NCC-based template matching |
| `image_ssim` | Structural Similarity Index (SSIM) |
| `image_psnr` | PSNR and MAE image quality metrics |
| `image_correlation` | Pearson correlation and covariance |
| `image_cosine_similarity` | Cosine similarity, Euclidean, L1 distance |
| `image_knn_classify` | K-nearest-neighbor scalar classification |
| `image_label_stats` | Per-label region statistics |
| `image_label_boundary` | Label region boundaries and inter-label contact |
| `image_region_props` | Per-label region props (area, centroid, bbox, mean) |
| `image_statistics_map` | Local statistics map computation |
| `image_percentile` | Percentile queries and percentile-based clipping |
| `image_percentile_filter` | Generalized percentile filter |
| `image_max_pool` | Max/min/avg pooling operations |
| `image_pyramid` | Gaussian and Laplacian image pyramids |
| `image_resize_nearest` | Nearest-neighbor image resizing |
| `image_slice_extract` | Axis-aligned slice extraction |
| `image_extract_slice` | 2D slice extraction from 3D ImageData |
| `image_extract_channel` | Extract single channel from multi-component array |
| `image_extract_component` | Extract component from multi-component array |
| `image_stack` | Stack 2D slices to 3D volume |
| `image_profile` | 1D profile extraction (row/column/diagonal) |
| `image_max_projection` | Maximum intensity projection |
| `image_min_projection` | Minimum intensity projection |
| `image_sum_projection` | Sum intensity projection |
| `image_gradient_magnitude` | Gradient magnitude via central differences |
| `image_gradient_vector` | 3-component gradient vector |
| `image_gradient_direction` | Gradient direction and orientation |
| `image_laplacian` (filters) | Discrete Laplacian of scalar field |
| `image_divergence` | Divergence of vector fields |
| `image_curl` | Curl magnitude of vector fields |
| `image_rle_compress` | Run-length encoding and value analysis |
| `image_quantize` | Scalar quantization and dithered quantization |
| `image_invert` | Invert ImageData values |
| `image_clamp` | Clamp values to range |
| `image_normalize_range` | Normalize to [0,1] range |
| `image_remap_range` | Linear scalar range remapping |
| `image_mask_by_scalar` | Scalar-range masking |
| `image_mask_apply` | Mask application and creation |
| `image_resample_simple` | Simple nearest-neighbor resampling |
| `image_power_law` | Gamma correction and sigmoid contrast |
| `image_math` | Element-wise ImageData math operations |
| `image_abs` | Absolute value, sqrt, log, exp on ImageData |
| `image_logical` | Logical operations (AND/OR/XOR/NOT) on binary ImageData |
| `image_local_binary_pattern` | Local Binary Pattern texture descriptor |
| `image_entropy` | Local Shannon entropy for texture analysis |
| `image_register_translation` | 2D image translation registration via cross-correlation |
| `image_thin` | Zhang-Suen morphological thinning |
| `image_hessian` | Hessian determinant for blob/saddle detection |

#### Data Conversion & Other (additional)
| `poly_data_bounds` | PolyData bounding box queries |
| `poly_data_normals_flip` | Flip polygon normals |
| `poly_data_transform_filter` | Direct geometric transformations on PolyData |
| `poly_data_to_table` | Convert PolyData to Table |
| `angle_between` | Dihedral angles between adjacent triangles |
| `point_set_registration` | Translation-only registration + error |
| `point_cloud_downsample` | Voxel and random point cloud downsampling |
| `face_varying` | Convert cell data to per-face-vertex point data |
| `image_dilate_erode` | Morphological dilation and erosion on ImageData |
| `cell_data_to_point_data_avg` | Average cell data values to shared points |
| `mesh_dual_graph` | Face adjacency dual graph as lines between centroids |
| `mesh_sharp_edges` | Sharp edge detection by dihedral angle |

### Implemented — VTK C++ Filter Categories with Coverage

#### Flow Visualization (VTK/Filters/FlowPaths — 21 filters)
- [x] `particle_tracer` — Time-dependent particle advection through vector fields (steady + temporal RK4)
- [x] `particle_path_filter` — RK4 particle pathlines with speed/time/particle ID arrays
- [x] `streak_line` — Streakline computation for unsteady flows (RK2, multi-seed)
- [x] `temporal_stream_tracer` — Temporally interpolated streamlines between consecutive fields
- [x] `lagrangian_particle_tracker` — Euler-step tracking with age, distance, speed accumulation
- [x] `bsp_tree` — BSP tree for fast spatial cell location with nearest-cell queries
- [x] `interpolated_velocity_field` — VelocityField trait, ImageData + composite multi-block fields, RK4 integration
- [x] `vector_field_topology` — Critical points (source/sink/saddle/center), vorticity computation
- [x] `even_spaced_streamlines_2d` — Evenly-spaced streamline seeding (Jobard-Lefer, RK2)

#### Statistics (VTK/Filters/Statistics — 25 filters)
- [x] `descriptive_statistics` — Mean, variance, skewness, kurtosis, percentiles, correlation matrix
- [x] `correlative_statistics` — Covariance/Pearson/Spearman correlation matrices
- [x] `multi_correlative_statistics` — Multivariate stats, Mahalanobis distance, outlier detection
- [x] `principal_component_analysis` — PCA with power iteration, projection, explained variance
- [x] `kmeans_statistics` — K-means clustering (k-means++ init, points + tables)
- [x] `contingency_statistics` — Contingency tables with chi-squared and Cramér's V
- [x] `order_statistics` — Quantiles, ranks (with ties), IQR outlier detection
- [x] `compute_quartiles` — Quartiles, IQR, fences, five-number summary
- [x] `high_dimensional_histogram` — N-D histogram with bin centers and 2D table output
- [x] `auto_correlative_statistics` — Autocorrelation + partial autocorrelation (Durbin-Levinson)

#### Parallel & Distributed (VTK/Filters/Parallel — 62 filters)
- [x] `parallel_for` — Rayon parallel map/transform/reduce over points and voxels
- [x] `ghost_cells` — Ghost cell generation for domain decomposition + mesh partitioning
- [x] `process_id_scalars` — Process/partition ID labeling + tag-and-merge
- [x] `redistribute_data_set` — Spatial-sort partitioning of points and cells with balance metric
- [x] `d3_decomposition` — K-d tree recursive bisection for points and cells
- [x] `aggregate_dataset` — Merge PolyData/Tables, vertical concat, cross-table statistics
- [x] `piece_request` — Streaming/chunked data: compute pieces, extract sub-regions, stream processing

#### Adaptive Mesh Refinement (VTK/Filters/AMR — 7 filters)
- [x] `amr_resample_filter` — Resample ImageData onto HTG coarse cells with averaging
- [x] `amr_slice_filter` — Axis-aligned slice through HTG as PolyData quads
- [x] `amr_connectivity` — Flood-fill connected regions on HTG coarse grid with threshold
- [x] `amr_to_multi_block` — HTG to MultiBlock + MultiBlock to uniform ImageData

#### HyperTree Grid (VTK/Filters/HyperTree — 29 filters)
- [x] `hyper_tree_grid_contour` — Isovalue contour on HTG coarse cell interfaces
- [x] `hyper_tree_grid_slice` — Arbitrary-plane slice through HTG with polygon intersection
- [x] `hyper_tree_grid_to_unstructured` — Convert HTG coarse cells to quads (2D) or hexahedra (3D)
- [x] `hyper_tree_grid_threshold` — Scalar threshold with mask array on HTG cell data
- [x] `hyper_tree_grid_gradient` — Central-difference gradient (X/Y/Z + magnitude) on HTG
- [x] `hyper_tree_grid_cell_centers` — Cell centers + geometry + to_image_data conversion
- [x] `hyper_tree_grid_depth_limiter` — Resolution limiting, depth stats, uniform grid conversion
- [x] `hyper_tree_grid_geometry` — Extract coarse-grid surface as PolyData quads

#### Point Cloud Processing (VTK/Filters/Points — 42 filters, partially covered)
- [x] `point_interpolator` — Kernel-based point interpolation (Gaussian, Shepard, IDW, SPH cubic spline)
- [x] `sph_interpolator` — SPH interpolation (cubic spline, Wendland C2, quintic) + density
- [x] `point_occupancy` — Voxel occupancy grid from points (count + binary)
- [x] `point_smoothing` — K-NN geometric position + scalar smoothing
- [x] `signed_distance_field_par` — Parallel SDF from oriented point clouds (rayon)
- [x] `fit_implicit_function` — Fit quadric surface + plane to point cloud (least-squares)
- [x] `extract_enclosed_points` — Ray-cast inside/outside test against closed surface
- [x] `extract_hierarchical_bins` — Morton-code octree binning with extract and count
- [x] `radii_to_scalars` — Distance from center/centroid, normalized radii
- [x] `voronoi_gaussian_kernel` — Voronoi (IDW), Gaussian, and Epanechnikov kernel interpolation + KDE

#### Modeling (VTK/Filters/Modeling — 36 filters, partially covered)
- [x] `collision_detection` — Collision detection between triangle meshes (AABB + SAT)
- [x] `select_polygon` — Select/extract cells and points inside 2D polygon
- [x] `cookie_cutter` — Cut surface by polygon or circle outline
- [x] `project_sphere` — Radial projection onto sphere with spherical coordinates
- [x] `earth` source — Earth globe geometry with latitude/longitude data + equirectangular UV
- [x] `sector` source — Sector/pie-slice geometry (filled + annular)
- [x] `outline_filter` — Wireframe bounding box for PolyData/ImageData/RectilinearGrid/UnstructuredGrid

#### Extraction (VTK/Filters/Extraction — 30 filters, partially covered)
- [x] `extract_time_steps` — Extract by index/time/range, subsample every Nth, snap to nearest
- [x] `extract_block` — Extract block by index/name, merge/flatten multi-block
- [x] `extract_data_over_time` — Extract scalars + position at closest point across time steps
- [x] `extract_rectilinear_grid` — Sub-grid extraction by index or coordinate ranges
- [x] `extract_unstructured_grid` — Extract by cell type, predicate, or index list
- [x] `extract_geometry` — Extract surface geometry from any dataset type to PolyData

#### Hybrid (VTK/Filters/Hybrid — 26 filters)
- [x] `weighted_transform` — Multi-transform blending + skeletal per-vertex bone weights
- [x] `depth_sort_poly_data` — Painter's algorithm depth sorting (in depth_peeling module)
- [x] `image_to_points` — ImageData to point cloud with optional scalar threshold
- [x] `poly_data_silhouette` — View-dependent silhouette edges (parallel + perspective)
- [x] `temporal_snap_to_time_step` — Snap to nearest available time step (in extract_time_steps)
- [x] `force_directed_layout` — Fruchterman-Reingold graph layout with degree data
- [x] `render_large_image` — Tile computation, RGBA compositing, and frustum adjustment for tiled rendering

#### Geometry (VTK/Filters/Geometry — 23 filters, partially covered)
- [x] `linear_to_quadratic` — Upgrade linear triangles to quadratic (6-node) with midpoint interpolation
- [x] `mark_boundary` — Mark boundary points/cells/edges with data arrays
- [x] `structured_grid_partitioner` — Partition along longest dimension with balance metric
- [x] `unstructured_grid_geometry` — Extract boundary/all/typed faces from UnstructuredGrid

#### Texture (VTK/Filters/Texture — 8 filters)
- [x] `texture_map_to_cylinder` — Cylindrical texture projection (manual + auto axis)
- [x] `texture_map_to_sphere_full` — Spherical texture mapping with seam correction (vertex duplication)
- [x] `transform_texture_coords` — Translate/scale/rotate/flip/tile/wrap/clamp UV coords
- [x] `threshold_texture_coords` — Scalar-to-UV mapping (linear, range, binary threshold)

#### Meshing (VTK/Filters/Meshing — 13 filters)
- [x] `quadric_decimation` — Garland-Heckbert with boundary preservation and topology checks
- [x] `constrained_smooth` — Boundary-fixed, displacement-limited, and normal-only smoothing
- [x] `incremental_octree` — Dynamic insert octree with nearest-neighbor and radius queries

#### Generic Data Handling (VTK/Filters/Generic — 9 filters)
- [x] `generic_clip` — Clip any dataset type (dispatches to PolyData conversion + clip_by_plane)
- [x] `generic_contour` — Contour any dataset type (dispatches to PolyData + contour)
- [x] `generic_cutter` — Cut any dataset type (dispatches to PolyData + slice_by_plane)
- [x] `generic_glyph` — Glyph any dataset type (dispatches to PolyData + glyph)

### Not Yet Implemented — VTK C++ Filter Categories

#### Temporal (VTK/Filters/Temporal — 4 filters)
- [ ] `CriticalTime` — time at which field exceeds threshold
- [ ] `DataObjectMeshCache` — cache mesh topology across time steps
- [ ] `ForceStaticMesh` — bypass mesh updates for static-topology data
- [ ] `TemporalSmoothing` — temporal averaging of data arrays

#### Reduction (VTK/Filters/Reduction — 6 filters)
- [ ] `ToImplicitArrayFilter` — compress explicit arrays to implicit/analytic representation
- [ ] `ToAffineArrayStrategy` / `ToConstantArrayStrategy` — specialized compression strategies
- [ ] `ToImplicitRamerDouglasPeuckerStrategy` — polyline simplification as implicit array

#### ReebGraph (VTK/Filters/ReebGraph — 5 filters)
- [ ] `UnstructuredGridToReebGraphFilter` — Reeb graph from scalar field topology
- [ ] `ReebGraphSimplificationFilter` — simplify Reeb graph by persistence
- [ ] `ReebGraphSurfaceSkeletonFilter` / `ReebGraphVolumeSkeletonFilter` — surface/volume skeleton extraction
- [ ] `ReebGraphToJoinSplitTreeFilter` — join/split tree decomposition

#### Tensor (VTK/Filters/Tensor — 2 filters)
- [ ] `TensorPrincipalInvariants` — principal tensor invariants (I1, I2, I3, J2, J3)
- [ ] `YieldCriteria` — von Mises, Tresca, and other yield criteria from stress tensors

#### Verdict Mesh Quality (VTK/Filters/Verdict — 5 filters)
- [ ] `MeshQuality` — comprehensive Verdict-library mesh quality metrics (50+ metrics)
- [ ] `BoundaryMeshQuality` — boundary-specific quality metrics
- [ ] `CellQuality` / `CellSizeFilter` — Verdict-based cell quality and size
- [ ] `MatrixMathFilter` — matrix eigenvalues, determinant, inverse operations

#### Selection (VTK/Filters/Selection — 3 filters)
- [ ] `CellDistanceSelector` — select cells within distance of seed cells
- [ ] `KdTreeSelector` — spatial selection via k-d tree queries
- [ ] `LinearSelector` — select cells along a line/ray

#### Topology (VTK/Filters/Topology — 1 filter)
- [ ] `FiberSurface` — extract fiber surfaces from bivariate fields

#### CellGrid / Discontinuous Galerkin (VTK/Filters/CellGrid — 46 filters)
- [ ] Entire DG/CellGrid subsystem — high-order DG cells, operators, evaluation, transcription

#### VTK Core Filters Not Yet Covered
- [ ] `HedgeHog` — oriented line glyphs from vector fields
- [ ] `FlyingEdges2D` — fast 2D contour extraction
- [ ] `SurfaceNets2D` — smooth 2D contour extraction
- [ ] `ClipClosedSurface` — clip mesh and cap with closed cross-section faces
- [ ] `BoxClipDataSet` — clip by oriented box volume
- [ ] `ClipVolume` — clip volumetric data by implicit function
- [ ] `PolyDataEdgeConnectivity` — edge-connectivity-based region extraction
- [ ] `StaticCleanPolyData` / `StaticCleanUnstructuredGrid` — topology-preserving cleaning
- [ ] `ResampleWithDataSet` — resample using dataset geometry as probe
- [ ] `ContourTriangulator` — triangulate contour polygons from slice output

#### VTK General Filters Not Yet Covered
- [ ] `ClipClosedSurface` — closed-surface clipping with cap generation
- [ ] `OBBTree` — oriented bounding box tree for spatial queries
- [ ] `OBBDicer` — OBB-based spatial decomposition
- [ ] `BooleanOperationPolyDataFilter` / `LoopBooleanPolyDataFilter` — VTK-native boolean operations
- [ ] `TableFFT` — FFT on table data
- [ ] `TableToStructuredGrid` — convert table columns to structured grid
- [ ] `TemporalPathLineFilter` / `TemporalStatistics` — temporal data processing
- [ ] `SphericalHarmonics` — spherical harmonic computation
- [ ] `UncertaintyTubeFilter` — tube rendering for uncertain trajectories
- [ ] `YoungsMaterialInterface` — multi-material interface reconstruction
- [ ] `WarpLens` — lens distortion correction
- [ ] `ContourTriangulator` — contour polygon triangulation

#### VTK Modeling Filters Not Yet Covered
- [ ] `ImprintFilter` — imprint one mesh onto another
- [ ] `ContourLoopExtraction` — extract closed contour loops
- [ ] `FitToHeightMapFilter` — fit mesh to height map
- [ ] `VolumeOfRevolutionFilter` — volumetric revolution solid
- [ ] `TrimmedExtrusionFilter` — extrusion trimmed by surface

#### VTK Points Filters Not Yet Covered
- [ ] `EuclideanClusterExtraction` — DBSCAN-like point cloud clustering
- [ ] `PCACurvatureEstimation` — PCA-based curvature from point clouds
- [ ] `DensifyPointCloudFilter` — densify sparse point clouds
- [ ] `PoissonDiskSampler` — blue-noise point sampling
- [ ] `VoxelGrid` — voxel-based point cloud downsampling
- [ ] `ConnectedPointsFilter` — connected components for point clouds

#### VTK Imaging Filters Not Yet Covered (~142 C++ classes, partially covered)
- [ ] `vtkImageFFT` / `vtkImageRFFT` — true FFT/inverse FFT (current impl is DFT-based)
- [ ] `vtkImageIdealHighPass` / `vtkImageIdealLowPass` — frequency-domain ideal filters
- [ ] `vtkImageButterworthHighPass` / `vtkImageButterworthLowPass` — Butterworth frequency filters
- [ ] `vtkImageCheckerboard` — checkerboard comparison of two images
- [ ] `vtkImageStencil` / `vtkImageToImageStencil` — stencil-based region processing
- [ ] `vtkImageMapToColors` — map scalar image through color lookup table
- [ ] `vtkImageMapToWindowLevelColors` — window/level → RGBA mapping

---

## I/O Formats (23 / ~70+ in VTK)

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
| glTF Binary `.glb` | vtk-io-gltf | yes | yes | glTF 2.0 binary (positions, normals, triangles) |
| VTK XML `.vtm` | vtk-io-xml | yes | yes | MultiBlock index file with block loading |
| EnSight Gold | vtk-io-ensight | yes | yes | ASCII case + geo + scalar/vector files |
| XDMF `.xdmf` | vtk-io-xdmf | | yes | XML with inline data for PolyData + ImageData |
| OFF `.off` | vtk-io-off | yes | yes | Object File Format ASCII + COFF vertex colors |
| DXF `.dxf` | vtk-io-dxf | yes | yes | AutoCAD DXF (3DFACE + LINE entities) |
| GeoJSON `.geojson` | vtk-io-geojson | yes | yes | Point, LineString, Polygon features |
| CSV/TSV `.csv/.tsv` | vtk-io-csv | yes | yes | Flexible delimiter, quoted fields, point import |
| BYU `.byu` | vtk-io-byu | yes | yes | BYU Movie format, 1-based polygon connectivity |
| LAS `.las` | vtk-io-las | yes | yes | LAS 1.2, intensity + classification, scale/offset |
| Facet `.facet` | vtk-io-facet | yes | yes | Simple vertex + triangle connectivity format |
| SEG-Y `.sgy` | vtk-io-segy | yes | | SEG-Y rev 1 seismic traces as ImageData or points |
| Tecplot `.dat` | vtk-io-tecplot | yes | yes | Tecplot ASCII FEPOINT zones with variables |
| FITS `.fits` | vtk-io-fits | yes | yes | 2D/3D image, multiple BITPIX, BSCALE/BZERO |

### Not Yet Implemented — By Priority

#### High Priority (widely used scientific/engineering formats)
- [ ] Exodus II `.exo` — HDF5-based FEM format (requires HDF5 C library)
- [ ] CGNS `.cgns` — CFD General Notation System (requires CGNS/HDF5 C library)
- [ ] NetCDF `.nc` — Network Common Data Format (requires NetCDF C library)

#### Medium Priority (domain-specific formats)
- [ ] Alembic `.abc` — Alembic interchange format (requires Alembic C++ library)
- [ ] OpenVDB `.vdb` — Sparse volumetric data (requires OpenVDB C++ library)
- [ ] USD `.usd/.usda/.usdc` — Universal Scene Description (requires USD C++ library)
- [ ] AMR formats — BoxLib, Chombo, SAMRAI AMR data
- [ ] CityGML `.gml` — 3D city model data
- [ ] DICOM `.dcm` — Medical imaging (DCMTK C++ library)
- [ ] MINC `.mnc` — Medical imaging (NetCDF-based)

#### Lower Priority (require heavy C/C++ dependencies)
- [ ] GDAL — Geospatial raster/vector via GDAL C library
- [ ] PDAL — Point cloud via PDAL C library
- [ ] OCCT/STEP/IGES — CAD formats via OpenCASCADE C++ library
- [ ] FFMPEG — Video I/O via FFMPEG C library
- [ ] MySQL/PostgreSQL/ODBC — Database I/O

---

## Rendering

### Implemented

- [x] `Camera` — position, focal point, FOV, orbit, dolly, pan, unproject, look_at, standard views (front/top/right/iso), right/up vectors
- [x] `Scene` / `Actor` — fluent builders, position/scale/visibility, `reset_camera()`, `with_axes()`, `with_silhouette()`
- [x] `Renderer` trait — backend-agnostic rendering abstraction
- [x] `ColorMap` — 15 presets: jet, viridis, plasma, inferno, magma, cividis, turbo, cool_to_warm, grayscale, black_body, blue_red, rd_yl_bu, rainbow_desaturated, parula, spectral + `by_name()` lookup
- [x] `Coloring::Solid` / `Coloring::ScalarMap` — per-actor coloring modes
- [x] wgpu backend — Blinn-Phong shading with per-actor materials, multi-light support, scalar color mapping
- [x] `Representation` modes — Surface, Wireframe, Points rendering via separate GPU pipelines
- [x] Mouse interaction — orbit (left drag), pan (middle drag), zoom (scroll)
- [x] Depth buffer

- [x] `Light` — directional, point, spot, and ambient light types with color/intensity, wired into GPU shader
- [x] `Scene.lights` — multi-light support (up to 8), headlight default
- [x] `Material` — per-actor surface properties wired into shader: ambient/diffuse/specular, flat shading (dpdx/dpdy), backface culling

- [x] Edge overlay (surface + edges) — `material.edge_visibility` + `edge_color` rendered as line overlay
- [x] Transparency / alpha blending — per-actor `opacity`, opaque-first then translucent with alpha blending
- [x] Texture mapping (2D textures on surfaces) — `Coloring::TextureMap` with `Texture`, CPU-side UV sampling
- [x] Volume rendering — CPU ray casting + GPU ray marching (3D texture, transfer function LUT, proxy cube geometry)
- [x] Axes / orientation widget — `AxesWidget` with camera-rotated XYZ arrows, labels, 2D overlay
- [x] Scalar bar (color legend) — `ScalarBar` widget with color bands, tick marks, 2D overlay pipeline
- [x] Text / label rendering (2D overlay) — built-in 5x7 bitmap font, used by scalar bar and axes widget
- [x] Screenshot / offscreen rendering — `render_to_image()` returns RGBA pixel buffer
- [x] Picking — CPU ray-cast picking with `pick()` returning actor/cell/point/position, `Camera::unproject()`
- [x] Anti-aliasing (MSAA) — 4x multisampling with resolve target
- [x] PBR materials — Cook-Torrance BRDF with metallic/roughness, `Material::pbr_metal()` / `pbr_dielectric()`
- [x] Silhouette / outline rendering — view-dependent silhouette edge extraction + line rendering
- [x] LOD (level of detail) — `LodSet` with distance-based mesh selection, integrated into renderer
- [x] Instanced rendering for glyphs — `InstancedGlyphs` with position/scale/color per instance, CPU flatten
- [x] Animation — `CameraAnimation` with `Track<T>` keyframes, easing functions (linear, quad, cubic)
- [x] Distance fog — linear/exponential/exponential² fog with configurable color, near/far, density
- [x] GPU color-ID picking — `GpuPicker` renders actor/cell IDs to offscreen buffer for pixel-perfect selection
- [x] Measurement tools — `measure()` for surface area, edge lengths, `point_distance()`, `angle_at_vertex()`, `triangle_area()`
- [x] JSON scene export — `scene_to_json()` / `scene_to_json_string()` for scene configuration
- [x] Shadow mapping — `ShadowConfig` with resolution, bias, softness, light VP matrix computation
- [x] Skybox — `Skybox::Solid`/`Gradient`/`ThreeStop` + presets (studio, sky, white, black)
- [x] Bloom — `BloomConfig` with threshold, intensity, radius, Gaussian kernel generation
- [x] Annotations — `Label3D`, `DistanceRuler`, `AngleProtractor` for 3D measurement overlays
- [x] Stereo rendering — `StereoConfig` with side-by-side, anaglyph, top/bottom modes + eye separation
- [x] Screenshot — `save_ppm()`, `save_bmp()`, `save_tga()` for saving rendered images (no deps)
- [x] Temporal data — `TemporalDataSet` with time steps, interpolation, bracketing
- [x] Streamline seeding — `seed_line()`, `seed_plane()`, `seed_sphere()`, `seed_circle()`
- [x] Transfer function editor — `TransferFunctionEditor` with add/remove/move control points
- [x] Viewport system — `Viewport` with split-screen presets (half, quad grid), pixel conversion, hit testing
- [x] Camera turntable + zoom — `CameraAnimation::turntable()` / `zoom()` for presentations
- [x] Extract by cell type — `extract_cells_by_type()`, `cell_type_counts()` for UnstructuredGrid
- [x] Table statistics — `describe_column()`, `correlation()`, `linear_regression()` for data analysis
- [x] GPU clip planes — up to 6 clip planes per scene for section views, `ClipPlane::x()`/`y()`/`z()`
- [x] Scene serialization — `save_scene_config()` / `load_scene_config()` for camera, lights, clip planes

### Not Yet Implemented — Rendering

#### High Priority
- [ ] Clip plane rendering — render clipped geometry with capped cross-section faces
- [ ] Multi-viewport rendering — render multiple views in split-screen layout (Viewport exists, not wired to GPU)
- [x] Depth peeling — CPU depth peeling layers + depth-sorted mesh for painter's algorithm
- [ ] Shadow rendering — shadow map pass wired into GPU pipeline (ShadowConfig exists, not wired)
- [ ] Bloom post-process — bloom pass wired into GPU pipeline (BloomConfig exists, not wired)
- [ ] Stereo rendering — stereo pass wired into GPU pipeline (StereoConfig exists, not wired)
- [ ] Point sprites — GPU point rendering with configurable size/shape
- [ ] Line width — GPU line width control

#### Medium Priority
- [ ] Ray tracing — hardware/software ray tracing for reflections and global illumination
- [ ] Environment mapping — cube map / IBL environment lighting
- [ ] Screen-space ambient occlusion (SSAO) — post-process ambient occlusion
- [ ] Glyph instancing on GPU — instanced rendering of glyph meshes (currently CPU flatten)
- [ ] Texture atlas — efficient multi-texture rendering
- [ ] TrueType font rendering — FreeType or equivalent for crisp text at any size
- [ ] Color bar improvements — gradient rendering, custom tick formatting
- [ ] Axes cube widget — interactive axes cube (like ParaView)
- [ ] 2D rendering context — 2D drawing API for charts, plots, annotations

#### Lower Priority
- [ ] VR/XR support — OpenXR-based head-mounted display rendering
- [ ] Remote rendering — server-side rendering with image streaming
- [ ] Parallel rendering — distributed rendering with compositing
- [ ] Path tracing — Monte Carlo path tracing for photorealistic rendering
- [ ] Subdivision surface rendering — GPU tessellation of subdivision surfaces
- [ ] Impostor rendering — billboards for distant objects
- [ ] Depth-of-field — post-process depth-of-field blur

---

## Infrastructure

### Implemented

- [x] Workspace with 48 crates + `vtk` convenience crate (`use vtk::prelude::*`)
- [x] 10 examples (triangle, shapes, isosurface, scalar_colors, showcase, pipeline_demo, volume, mesh_info, headless_render, bench_filters)
- [x] ~9100 unit tests (incl. cross-format I/O roundtrip integration tests, doctests, proptests, Send+Sync assertions)
- [x] Clippy-clean (`-D warnings`)

- [x] Pipeline system — `Pipeline` with chained filters, lazy evaluation, caching, and invalidation
- [x] Parallel filter execution (rayon) — `compute_normals_par`, `elevation_par`, `smooth_par`, `marching_cubes_par`
- [x] WASM / web target support — all non-GPU crates compile for wasm32-unknown-unknown
- [x] Python bindings (PyO3) — `vtk-python` crate with PolyData, sources, filters, I/O
- [x] Benchmarks — `cargo run --release --example bench_filters` (sphere, normals, marching cubes, etc.)
- [x] Documentation (rustdoc with examples) — doc examples on DataArray, PolyData, ImageData, Camera, ColorMap
- [x] Property / fuzz testing — proptest for DataArray, CellArray, Points, PolyData roundtrips
- [x] Streaming / chunked processing — mmap_data utilities, chunked read/write, memory estimation
- [x] Parallel pipeline — rayon-based parallel branches, map, merge, chunked processing

### Not Yet Implemented — Infrastructure

#### High Priority
- [ ] CI/CD — automated testing, linting, and benchmarking
- [ ] Published documentation — hosted rustdoc with examples

#### Medium Priority
- [ ] MPI support — distributed memory parallelism via MPI bindings
- [ ] GPU compute filters — wgpu compute shaders for data-parallel filters
- [ ] Memory-mapped data — mmap large files for zero-copy I/O
- [ ] Data streaming — progressive loading and rendering of large datasets
- [ ] Plugin system — dynamic loading of filter/reader/writer plugins

#### Lower Priority
- [ ] In-situ visualization — co-processing with simulation codes
- [ ] Web viewer — WASM-based interactive 3D viewer in browser
- [ ] Jupyter integration — inline rendering in Jupyter notebooks via Python bindings
- [ ] Language bindings — Java, C, JavaScript/TypeScript bindings
