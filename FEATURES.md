# FEATURES.md

Feature tracking for vtk-rs — a pure Rust reimplementation of VTK.

Last updated: 2026-03-26 | Tests: 1960 | Clippy: clean

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
- [x] `CellLocator` — BVH-based cell locator for closest-cell and radius queries
- [ ] Higher-order / Lagrange / Bezier cell support
- [ ] `Molecule` data type

---

## Geometry Sources (33 / ~40 in VTK)

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

### Not Yet Implemented

---

## Processing Filters (432 / ~300+ in VTK)

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

### High-Priority — Not Yet Implemented

### Lower-Priority — Not Yet Implemented

- [x] `temporal_interpolator` — Interpolate between time steps

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
- [x] `ColorMap` — jet, viridis, cool-to-warm, grayscale, plasma, inferno, turbo, black-body, blue-red
- [x] `Coloring::Solid` / `Coloring::ScalarMap` — per-actor coloring modes
- [x] wgpu backend — Phong shading, smooth normals, scalar color mapping
- [x] Mouse interaction — orbit (left drag), zoom (scroll)
- [x] Depth buffer

- [x] `Light` — directional, point, spot, and ambient light types with color/intensity
- [x] `Scene.lights` — multi-light support with headlight default
- [x] `Material` — surface material properties (ambient/diffuse/specular/shininess, edge, flat shading)

### Not Yet Implemented

- [ ] Edge overlay (surface + edges)
- [ ] Transparency / alpha blending
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
- [x] 1588 unit tests
- [x] Clippy-clean (`-D warnings`)

### Not Yet Implemented

- [ ] Pipeline system (lazy evaluation, caching, automatic updates)
- [ ] Parallel filter execution (rayon integration)
- [ ] WASM / web target support
- [ ] Python bindings (PyO3)
- [ ] Benchmarks
- [ ] Documentation (rustdoc with examples)
- [ ] Property / fuzz testing
