# vtk-rs

A pure Rust reimplementation of [VTK](https://vtk.org/) (The Visualization Toolkit). Not an FFI binding — a ground-up Rust implementation of VTK's core concepts.

**~300K lines of Rust | 5300+ source files | 9300+ tests | 22 I/O formats | 4000+ filters | wgpu rendering**

Note that this package is yet to undergo a fair amount of testing!

## Quick Start

```rust
use vtk::prelude::*;

// Generate a sphere, compute normals, render
let mesh = vtk::filters::core::quick::sphere_with_elevation();
let scene = Scene::from_poly_data(mesh);
```

```toml
[dependencies]
vtk = "0.1"
```

## Build Times

| What you need | Crates | Clean build |
|---------------|--------|-------------|
| Data structures only | `vtk-data` | **2s** |
| + STL/OBJ I/O | + `vtk-io-stl`, `vtk-io-obj` | **6s** |
| + Core filters (sources, normals, clip, pipeline) | + `vtk-filters` | **12s** |
| + Rendering (wgpu) | + `vtk-render-wgpu` | **20s** |
| Everything (`vtk` crate) | all | **29s** |

The heavy filter crates (`vtk-filters-image` 3000+ modules, `vtk-filters-mesh` 800+ modules) are feature-gated in `vtk-filters` so they don't compile unless you need them.

## Module Structure

```
vtk::types         Scalar, ScalarType, CellType, VtkError, BoundingBox, math, color
vtk::data          PolyData, ImageData, UnstructuredGrid, DataArray, CellArray, KdTree, ...
vtk::filters       Sources (62), processing (4000+), pipeline, convert, topology
vtk::io            22 formats: VTK, STL, OBJ, PLY, XML, glTF, OFF, DXF, GeoJSON, CSV, ...
vtk::render        Camera, Scene, Actor, Material, Light, ColorMap (15 presets), ...
vtk::render_wgpu   wgpu GPU backend: MSAA, PBR, shadows, SSAO, bloom, volume rendering
vtk::parallel      MPI-ready: spatial decomposition, ghost exchange, collective ops
```

## Data Structures

| Type | Description |
|------|-------------|
| `PolyData` | Polygonal mesh (triangles, quads, lines, vertices) |
| `ImageData` | Regular grid with implicit coordinates |
| `UnstructuredGrid` | Mixed-cell mesh with explicit connectivity |
| `RectilinearGrid` | Axis-aligned grid with non-uniform spacing |
| `StructuredGrid` | Curvilinear grid with explicit points |
| `DataArray<T>` | N-component typed array for point/cell data |
| `AnyDataArray` | Type-erased enum over all `DataArray<T>` variants |
| `Table` | Columnar data for analysis |
| `CellGrid` | Discontinuous Galerkin high-order cells |
| `MultiBlockDataSet` | Composite dataset |

## Geometry Sources (62)

sphere, cube, cone, cylinder, plane, arrow, disk, line, torus, helix, ellipsoid, capsule, geodesic_sphere, icosphere, superquadric, platonic_solid, frustum, spring, grid, text_3d, wavelet, mobius, klein_bottle, trefoil_knot, boy_surface, seashell, terrain, gear, star, and 30+ more.

## I/O Formats (22)

| Format | Extension | Read | Write |
|--------|-----------|:----:|:-----:|
| VTK Legacy | `.vtk` | yes | yes |
| VTK XML | `.vtp/.vtu/.vti/.vtr/.vts/.vtm` | yes | yes |
| STL | `.stl` | yes | yes |
| OBJ | `.obj` | yes | yes |
| PLY | `.ply` | yes | yes |
| glTF | `.glb` | yes | yes |
| OFF | `.off` | yes | yes |
| DXF | `.dxf` | yes | yes |
| GeoJSON | `.geojson` | yes | yes |
| CSV/TSV | `.csv/.tsv` | yes | yes |
| EnSight | `.case` | yes | yes |
| FITS | `.fits` | yes | yes |
| LAS | `.las` | yes | yes |
| SEG-Y | `.sgy` | yes | |
| Tecplot | `.dat` | yes | yes |
| BYU | `.byu` | yes | yes |
| Facet | `.facet` | yes | yes |
| XDMF | `.xdmf` | | yes |
| DICOM | `.dcm` | yes | |
| CityGML | `.gml` | yes | |
| Video | `.mp4` etc | | yes* |

*Video requires `ffmpeg` feature + system ffmpeg libs. Also available: HDF5-based formats (Exodus, CGNS, NetCDF, MINC, AMR) via `vtk-io-hdf5` with system `libhdf5-dev`, and GDAL (GeoTIFF, Shapefile) via `vtk-io-gdal` with system `libgdal-dev`.

## Rendering

wgpu-based GPU rendering with:

- Blinn-Phong and Cook-Torrance PBR shading
- Shadow mapping with 3x3 PCF
- Screen-space ambient occlusion (SSAO)
- Bloom post-processing
- Depth-of-field
- GPU volume rendering (ray marching)
- Stereo rendering (side-by-side, anaglyph, top/bottom)
- Multi-viewport split-screen
- Silhouette edges, wireframe, point rendering
- 15 color map presets (jet, viridis, plasma, inferno, etc.)
- Scalar bar, axes widget, axes cube, annotations
- GPU color-ID picking
- LOD, instanced glyphs, clip planes (6 max)
- Fog (linear, exponential)
- Offscreen rendering, PPM/BMP/TGA screenshot export
- CPU ray tracer and Monte Carlo path tracer
- TrueType font rendering (feature-gated `truetype`)

## Examples

```bash
cargo run --example triangle        # basic PolyData + render window
cargo run --example shapes          # sphere, cube, cone, cylinder, arrow
cargo run --example isosurface      # marching cubes on gyroid
cargo run --example scalar_colors   # elevation + colormap visualization
cargo run --example showcase        # PBR, transparency, axes, scalar bar
cargo run --example pipeline_demo   # filter pipeline + multi-format I/O
cargo run --example volume          # GPU volume rendering
```

## Design Principles

- **Ownership over refcounting** — no `Arc` by default; Rust ownership replaces VTK's reference counting
- **Enum-based type erasure** — `AnyDataArray` is an enum, not `Box<dyn Trait>`
- **Traits over inheritance** — `DataObject` and `DataSet` traits replace class hierarchies
- **Filters as functions** — plain `fn(&PolyData) -> PolyData`, composable without a pipeline
- **Pipeline optional** — `Pipeline` struct available for lazy evaluation + caching when needed
- **Feature-gated heavy deps** — image/mesh filters, ffmpeg, hdf5, gdal, fontdue are all optional

## Feature Flags

```toml
# Core filters only (12s build):
vtk-filters = "0.1"

# Add specific filter groups:
vtk-filters = { version = "0.1", features = ["image", "mesh"] }

# Everything:
vtk-filters = { version = "0.1", features = ["all"] }

# Optional system-library features:
vtk-io-hdf5 = { version = "0.1", features = ["exodus", "cgns", "netcdf"] }
vtk-io-gdal = { version = "0.1", features = ["gdal"] }
vtk-io-video = { version = "0.1", features = ["ffmpeg"] }
vtk-render = { version = "0.1", features = ["truetype"] }
```
