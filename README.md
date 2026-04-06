# vtk-pure-rs

A pure Rust reimplementation of [VTK](https://vtk.org/) (The Visualization Toolkit). Not an FFI binding — a ground-up Rust implementation of VTK's core concepts.

**~300K lines of Rust | 5184 features | 575+ tests | 22 I/O formats | wgpu rendering**

Most functions tested to give the same output as C++ version. But see TODO.md for some that are not; example files needed to do full testing

A few features have not been included; contact if interested in these

## Performance vs VTK C++ 9.6

Benchmarked 141 operations against VTK C++ 9.6. **On average 17% faster than C++** (0.83x ratio).

| Category | Count | % |
|---|---|---|
| **Faster than C++** | **109** | **77%** |
| Within 2x | 23 | 16% |
| 2-3x slower | 5 | 4% |
| >3x slower | 4 | 3% |

Biggest wins: signed_distance (32x), poly_data_distance (36x), boolean_union (6x), normals (10x), reflect (11x), shrink_large (10x), connectivity_large (7x), surface_nets (7x).

Uses `Arc<Vec<T>>` copy-on-write storage for zero-copy clone (matching VTK's `ShallowCopy` semantics), `target-cpu=native`, and rayon parallelism for flying edges.

## Test Coverage

- **575 tests** (434 validation + 141 performance)
- **366/396 (92%)** non-extra features tested against VTK C++ reference output
- **30 features** remain untested (exotic I/O, GPU, HyperTreeGrid, Reeb graph internals) — see [TODO.md](TODO.md)

## Quick Start

```rust
use vtk_pure_rs::data::*;
use vtk_pure_rs::filters::core::sources::sphere::{sphere, SphereParams};

let mesh = sphere(&SphereParams::default());
println!("points: {}, cells: {}", mesh.points.len(), mesh.polys.num_cells());
```

```toml
[dependencies]
vtk-pure-rs = "0.2"
```

## Build Times

| What you need | Features | Clean release build |
|---------------|----------|---------------------|
| Core filters + common I/O | default | **16s** |
| + all non-extra filter groups + all I/O | + `filters-smooth`, `filters-transform`, etc. + `io-all` | **25s** |
| + extra sources + image + mesh (everything minus GPU) | all non-GPU features | **1m 50s** |

The heavy filter modules (`filters-image` 3000+ modules, `filters-mesh` 800+ modules) are feature-gated so they don't compile unless you need them.

## Feature Counts

| Category | Non-extra | Extra | Total |
|----------|-----------|-------|-------|
| Sources | 64 | 350 | 414 |
| Core Filters | 26 | — | 26 |
| Geometry | 71 | — | 71 |
| Extract | 20 | — | 20 |
| Filter Data | 31 | — | 31 |
| Points | 24 | — | 24 |
| Grid | 19 | — | 19 |
| Transform | 17 | — | 17 |
| Cell | 16 | — | 16 |
| Clip | 11 | — | 11 |
| Statistics | 11 | — | 11 |
| Flow | 11 | — | 11 |
| Distance | 11 | — | 11 |
| Smooth | 10 | — | 10 |
| Subdivide | 8 | — | 8 |
| Normals | 8 | — | 8 |
| Texture | 7 | — | 7 |
| GPU | 5 | — | 5 |
| Boolean | 4 | — | 4 |
| Image | — | 3015 | 3015 |
| Mesh | — | 816 | 816 |
| Core Image | — | 187 | 187 |
| Core Mesh | — | 420 | 420 |
| I/O Formats | 22 | — | 22 |
| **Total** | **396** | **4788** | **5184** |

See [FEATURES.md](FEATURES.md) for the full annotated feature list.

## Module Structure

```
vtk_pure_rs::types      Scalar, ScalarType, CellType, VtkError, BoundingBox, math, color
vtk_pure_rs::data       PolyData, ImageData, UnstructuredGrid, DataArray, CellArray, KdTree, ...
vtk_pure_rs::filters    Sources (64+350), processing (4000+), pipeline, convert, topology
vtk_pure_rs::io         22 formats: VTK, STL, OBJ, PLY, XML, glTF, OFF, DXF, GeoJSON, CSV, ...
vtk_pure_rs::render     Camera, Scene, Actor, Material, Light, ColorMap (15 presets), ...
vtk_pure_rs::render_wgpu  wgpu GPU backend: MSAA, PBR, shadows, SSAO, bloom, volume rendering
```

## Data Structures

| Type | Description |
|------|-------------|
| `PolyData` | Polygonal mesh (triangles, quads, lines, vertices) |
| `ImageData` | Regular grid with implicit coordinates |
| `UnstructuredGrid` | Mixed-cell mesh with explicit connectivity |
| `RectilinearGrid` | Axis-aligned grid with non-uniform spacing |
| `StructuredGrid` | Curvilinear grid with explicit points |
| `DataArray<T>` | N-component typed array with `Arc<Vec<T>>` copy-on-write storage |
| `AnyDataArray` | Type-erased enum over all `DataArray<T>` variants |
| `Table` | Columnar data for analysis |
| `CellGrid` | Discontinuous Galerkin high-order cells |
| `MultiBlockDataSet` | Composite dataset |

## Geometry Sources (64 base + 350 extra)

sphere, cube, cone, cylinder, plane, arrow, disk, line, torus, helix, ellipsoid, capsule, geodesic_sphere, icosphere, superquadric, platonic_solid, frustum, spring, grid, text_3d, wavelet, mobius, klein_bottle, trefoil_knot, boy_surface, seashell, terrain, gear, star, and 30+ more.

Extra sources (behind `sources-extra` feature): airplane, amphora, castle_tower, dna_helix, lighthouse, rocket, space_station, wind_turbine, and 340+ more architectural/scientific/artistic models.

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

*Video requires `ffmpeg` feature. HDF5-based formats (Exodus, CGNS, NetCDF, MINC, AMR) available via `io-hdf5` feature with system `libhdf5-dev`. GDAL (GeoTIFF, Shapefile) via `io-gdal` feature with system `libgdal-dev`.

## Rendering

wgpu-based GPU rendering (behind `render-wgpu` feature):

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

- **Copy-on-write storage** — `Arc<Vec<T>>` in DataArray/CellArray gives zero-copy clone with automatic CoW on mutation, matching VTK's `ShallowCopy` semantics
- **Enum-based type erasure** — `AnyDataArray` is an enum, not `Box<dyn Trait>`
- **Traits over inheritance** — `DataObject` and `DataSet` traits replace class hierarchies
- **Filters as functions** — plain `fn(&PolyData) -> PolyData`, composable without a pipeline
- **Pipeline optional** — `Pipeline` struct available for lazy evaluation + caching when needed
- **Feature-gated heavy deps** — image/mesh filters, ffmpeg, hdf5, gdal, fontdue are all optional
- **Native CPU targeting** — `.cargo/config.toml` sets `target-cpu=native` for optimal SIMD

## Feature Flags

```toml
[dependencies]
vtk-pure-rs = "0.2"                                          # core filters only
vtk-pure-rs = { version = "0.2", features = ["filters-all"] } # all filters (excl. GPU)
vtk-pure-rs = { version = "0.2", features = ["full"] }        # everything

# Individual filter groups:
vtk-pure-rs = { version = "0.2", features = [
    "filters-smooth",      # smoothing filters
    "filters-transform",   # transform/warp/mirror/extrude
    "filters-subdivide",   # subdivision surfaces
    "filters-cell",        # cell operations
    "filters-boolean",     # boolean mesh operations
    "filters-distance",    # distance/collision/hausdorff
    "filters-image",       # 3000+ image processing filters [heavy]
    "filters-mesh",        # 800+ mesh processing filters [heavy]
    "sources-extra",       # 350 extra geometry sources
    "io-all",              # all 22 I/O formats
    "render-wgpu",         # wgpu GPU rendering
    "parallel",            # MPI-ready parallel ops
    "truetype",            # TrueType font rendering
] }
```
