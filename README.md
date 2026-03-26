# vtk-rs

A pure Rust reimplementation of [VTK](https://vtk.org/) (The Visualization Toolkit) — the widely-used C++ library for 3D graphics, image processing, and scientific visualization.

This is **not** an FFI binding to the C++ library. It is a ground-up Rust implementation of VTK's core concepts, designed to be idiomatic, safe, and fast.

> **Status:** Early development. Core data structures, I/O (6 formats), geometry sources (6), mesh filters (18), scalar color mapping, and wgpu rendering are functional. The API is unstable.

## Quick Start

```rust
use vtk_data::PolyData;
use vtk_filters::sources;
use vtk_io_legacy::LegacyWriter;

// Generate a sphere
let sphere = sources::sphere(&sources::sphere::SphereParams {
    radius: 1.0,
    ..Default::default()
});

// Write to VTK legacy format
let writer = LegacyWriter::ascii();
writer.write_poly_data("sphere.vtk".as_ref(), &sphere).unwrap();
```

## Crates

### Core

| Crate | Description |
|-------|-------------|
| **vtk-types** | `Scalar` trait, `ScalarType` / `CellType` enums, `VtkError`, `BoundingBox` |
| **vtk-data** | `DataArray<T>`, `AnyDataArray`, `CellArray`, `Points`, `FieldData`, `DataSetAttributes`, `PolyData`, `ImageData`, `UnstructuredGrid` |

### Filters

| Crate | Description |
|-------|-------------|
| **vtk-filters** | 17 geometry sources + 65 processing filters (see [Filters](#filters) below) |

### I/O

| Crate | Formats | Read | Write |
|-------|---------|:----:|:-----:|
| **vtk-io-legacy** | VTK legacy `.vtk` (PolyData + ImageData + UnstructuredGrid), ASCII & binary | yes | yes |
| **vtk-io-xml** | VTK XML `.vtp` / `.vtu` / `.vti` / `.vtr` / `.vts`, ASCII | yes | yes |
| **vtk-io-stl** | STL `.stl`, ASCII & binary | yes | yes |
| **vtk-io-obj** | Wavefront `.obj` | yes | yes |
| **vtk-io-ply** | Stanford PLY `.ply`, ASCII & binary | yes | yes |

### Rendering

| Crate | Description |
|-------|-------------|
| **vtk-render** | Backend-agnostic: `Camera`, `Scene`, `Actor`, `Renderer` trait, `ColorMap` (jet, viridis, cool-to-warm, grayscale) |
| **vtk-render-wgpu** | [wgpu](https://wgpu.rs/) backend with Phong shading, smooth normals, scalar color mapping, mouse orbit/zoom |

## Building

```bash
cargo build                              # build all crates
cargo test --workspace                   # run all tests (296 tests)
cargo clippy --workspace -- -D warnings  # lint
```

## Examples

```bash
cargo run --example triangle      # create a triangle, write .vtk, open render window
cargo run --example shapes        # render sphere, cube, cone, cylinder, arrow
cargo run --example isosurface    # marching cubes on a gyroid, write .stl, render
cargo run --example scalar_colors # elevation + jet/viridis/cool-to-warm colormaps
```

All examples open an interactive 3D window (requires a GPU). Left-click drag to orbit, scroll to zoom.

## Data Structures

The core data model mirrors VTK's design, translated to idiomatic Rust:

| Type | Description |
|------|-------------|
| `DataArray<T>` | Contiguous array of N-component tuples, generic over scalar type |
| `AnyDataArray` | Type-erased enum over all `DataArray<T>` variants (f32, f64, i8..u64) |
| `CellArray` | Cell topology via offsets + connectivity arrays (matches VTK's `vtkCellArray`) |
| `PolyData` | Polygonal mesh: 4 cell arrays (vertices, lines, polygons, triangle strips) + point/cell data |
| `ImageData` | Regular grid with implicit coordinates from extent, spacing, and origin |
| `UnstructuredGrid` | Arbitrary mixed-cell mesh with explicit connectivity and per-cell types |

Traits `DataObject` and `DataSet` replace VTK's class hierarchy.

## Filters

**Sources:** `sphere`, `cube`, `cone`, `cylinder`, `plane`, `arrow`, `disk`, `line`, `point_source`, `regular_polygon`, `arc`, `superquadric`, `platonic_solid`, `frustum`, `parametric_function`, `bounding_box_source`, `axes`

**Processing filters (65):**

| Filter | Description |
|--------|-------------|
| `normals` | Compute smooth vertex normals (Newell's method + averaging) |
| `triangulate` | Convert quads/polygons/strips to triangles (fan triangulation) |
| `append` | Merge multiple PolyData with point index renumbering |
| `clean` | Merge duplicate points (spatial hash), remove degenerate cells |
| `transform` | Apply 4x4 matrix to points and normals |
| `marching_cubes` | Extract isosurface from scalar field on ImageData |
| `clip` | Clip mesh by plane, splitting crossing triangles |
| `slice` | Cut mesh by plane, producing intersection line segments |
| `contour` | Extract contour lines at scalar isovalues (2D analogue of marching cubes) |
| `elevation` | Compute scalar measuring projection along an axis |
| `threshold` | Extract cells by scalar value range |
| `decimate` | Quadric error metric mesh simplification |
| `smooth` | Laplacian smoothing with boundary preservation |
| `subdivide` | Loop subdivision (each triangle becomes 4) |
| `warp` | Displace vertices by scalar (along normals) or vector field |
| `connectivity` | Extract connected components (union-find) |
| `extract_surface` | Boundary surface of UnstructuredGrid (tetra, hex, wedge, pyramid) |
| `feature_edges` | Extract boundary, feature, manifold, and non-manifold edges |
| `reflect` | Mirror mesh across a coordinate plane, with optional input copy |
| `shrink` | Shrink cells toward their centroids |
| `orient` | Consistent polygon winding orientation (BFS edge traversal) |
| `cell_centers` | Generate vertex points at cell centroids |
| `extract_edges` | Extract all unique edges as line segments |
| `cell_data_to_point_data` | Convert cell attributes to point attributes by averaging |
| `point_data_to_cell_data` | Convert point attributes to cell attributes by averaging |
| `generate_ids` | Generate sequential PointIds and/or CellIds arrays |
| `mask_points` | Subsample points (every Nth or random) |
| `curvatures` | Discrete Gaussian and mean curvature computation |
| `gradient` | Least-squares gradient of a scalar field |
| `mass_properties` | Surface area, volume, and centroid of closed meshes |
| `densify` | Subdivide polygon edges exceeding a max length |
| `hull` | 3D convex hull (incremental algorithm) |
| `texture_map` | Generate texture coords (plane projection or spherical mapping) |
| `glyph` | Place scaled copies of a mesh at each input point |
| `tube` | Generate tubes with caps around line cells |
| `delaunay_2d` | 2D Delaunay triangulation (Bowyer-Watson) |
| `spline` | Catmull-Rom spline interpolation along polylines |
| `clip_data_set` | Clip UnstructuredGrid by plane |
| `windowed_sinc_smooth` | Low-pass windowed sinc mesh smoothing |
| `probe` | Interpolate source data at probe points |
| `stream_tracer` | RK4 streamline integration through vector fields |
| `voxel_modeller` | Convert PolyData to binary voxel ImageData |
| `sample_function` | Evaluate scalar function on ImageData grid |
| `integrate_attributes` | Integrate point data over surface area |
| `distance_poly_data` | Minimum distance from target points to source surface |
| `implicit_modeller` | Distance field from PolyData on ImageData grid |
| `tensor_glyph` | Place tensor-transformed glyphs at input points |
| `resample` | Resample source PolyData onto target ImageData grid |
| `calculator` | Expression-based scalar/vector field computation |
| `select_enclosed_points` | Ray-casting inside/outside classification |
| `extract_cells` | Extract cells by index or predicate |
| `icp` | Iterative Closest Point rigid registration |
| `extract_points` | Extract points by index or scalar range |
| `signed_distance` | Signed distance field from closed surface |
| `cell_size` | Compute area/length of polygon/line cells |
| `extrude` | Linear extrusion of 2D geometry along a direction |
| `cell_quality` | Aspect ratio, min/max angle, area quality metrics |
| `reverse_sense` | Flip polygon winding and normals |
| `strip` | Convert triangles to triangle strips |
| `interpolate` | Inverse distance weighting interpolation |
| `fill_holes` | Close open boundary loops with fan triangulation |
| `center_of_mass` | Point and area-weighted center of mass |
| `project_points` | Project points onto nearest surface location |
| `subdivide_midpoint` | Midpoint subdivision (triangles/quads) |
| `random_attributes` | Generate random scalar/vector data |
| `image_to_poly_data` | Convert ImageData surface to PolyData |

All filters are plain functions — no pipeline system required:

```rust
use vtk_filters::{normals, triangulate, transform, append, clean};

let sphere = vtk_filters::sources::sphere(&Default::default());

// Compute smooth vertex normals
let with_normals = normals::compute_normals(&sphere);

// Convert quads to triangles
let triangulated = triangulate::triangulate(&with_normals);

// Apply a transformation
let moved = transform::transform(&triangulated, &transform::translation(5.0, 0.0, 0.0));

// Merge multiple meshes
let merged = append::append(&[&triangulated, &moved]);

// Remove duplicate points
let cleaned = clean::clean(&merged, &clean::CleanParams::default());
```

Isosurface extraction from a scalar field:

```rust
use vtk_data::{ImageData, DataSet};
use vtk_filters::marching_cubes;

let image = ImageData::with_dimensions(50, 50, 50);
let scalars: Vec<f64> = (0..image.num_points())
    .map(|i| {
        let p = image.point(i);
        (p[0] * p[0] + p[1] * p[1] + p[2] * p[2]).sqrt()
    })
    .collect();
let isosurface = marching_cubes::marching_cubes(&image, &scalars, 25.0);
```

Scalar visualization with color maps:

```rust
use vtk_filters::{sources, elevation};
use vtk_render::{Actor, ColorMap, Scene};

let sphere = sources::sphere(&Default::default());
let colored = elevation::elevation_z(&sphere);
let mut scene = Scene::new();
scene.add_actor(
    Actor::new(colored).with_scalar_coloring(ColorMap::viridis(), None)
);
```

## Design Principles

- **Ownership over refcounting** — Rust ownership replaces VTK's `vtkObjectBase` reference counting. No `Arc` by default.
- **Enum-based type erasure** — `AnyDataArray` uses an enum over a closed set of scalar types instead of trait objects, enabling exhaustive matching and zero vtable overhead.
- **Traits over inheritance** — `DataObject` and `DataSet` traits replace VTK's deep class hierarchies.
- **Filters as functions** — No pipeline infrastructure needed. Compose filters by calling functions.
- **wgpu rendering** — Cross-platform WebGPU backend with Phong shading, smooth normals, and scalar-to-color mapping.
- **Multiple I/O formats** — VTK legacy (ASCII/binary), VTK XML, STL, OBJ, PLY.

## License

Same as the original VTK
