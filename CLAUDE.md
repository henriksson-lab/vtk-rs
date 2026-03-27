# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

vtk-rs is a pure Rust reimplementation of VTK (The Visualization Toolkit). The `VTK/` directory contains the C++ 9.6.0 source as reference — this is **not** an FFI bindings project.

## Build / Test / Lint

```bash
cargo build                              # build all crates
cargo test --workspace --exclude vtk-python  # run all tests (~2350+)
cargo clippy --workspace -- -D warnings  # lint
cargo run --example triangle             # render a triangle (requires GPU)
cargo run --example shapes               # render sphere/cube/cone/cylinder/arrow
cargo run --example isosurface           # marching cubes on a gyroid field, writes gyroid.stl
cargo run --example scalar_colors        # elevation-colored sphere/cone/plane with colormaps
cargo run --example showcase             # PBR, transparency, axes, scalar bar, keyboard controls
cargo run --example pipeline_demo        # filter pipeline + multi-format I/O + topology analysis
cargo run --release --example bench_filters  # performance benchmarks
```

## Workspace Structure

```
crates/
  vtk-types/       # Scalar trait, ScalarType, CellType (46 types), VtkError, BoundingBox, ImplicitFunction, higher_order
  vtk-data/        # DataArray, CellArray, Points, FieldData, DataSetAttributes, PolyData, ImageData, UnstructuredGrid, RectilinearGrid, StructuredGrid, MultiBlockDataSet, Table, KdTree, OctreePointLocator, CellLocator, Selection, Graph, Tree, ExplicitStructuredGrid, HyperTreeGrid, Molecule
  vtk-filters/     # 575 filters + 40 sources + pipeline + convert + topology + merge + selection_extract + io_utils
  vtk-io-legacy/   # VTK legacy format (.vtk) ASCII/binary reader and writer
  vtk-io-stl/      # STL format ASCII/binary reader and writer
  vtk-io-obj/      # Wavefront OBJ format reader and writer
  vtk-io-xml/      # VTK XML format reader/writer (.vtp, .vtu, .vti, .vtr, .vts, .vtm) — ASCII + binary
  vtk-io-ply/      # Stanford PLY format ASCII and binary reader and writer
  vtk-io-gltf/     # glTF 2.0 binary (.glb) reader and writer
  vtk-io-ensight/  # EnSight Gold ASCII reader/writer + LS-DYNA keyword reader
  vtk-io-xdmf/     # XDMF writer (inline XML data for PolyData + ImageData)
  vtk-render/      # Backend-agnostic: Camera, Scene, Actor, Renderer, ColorMap (15 presets), Material (Phong + PBR), Light, ScalarBar, AxesWidget, Picker, LOD, InstancedGlyphs, VolumeActor, TransferFunction, Animation, ClipPlane, SilhouetteConfig, Texture, Fog, Measurement
  vtk-render-wgpu/ # wgpu: MSAA, Phong + PBR shader, edge overlay, wireframe/points, transparency, offscreen, silhouette, overlay (scalar bar + axes + bitmap font), clip planes, fog, per-actor model matrix, GPU volume rendering, GPU color-ID picking
  vtk-python/      # PyO3 Python bindings (cdylib) — PolyData, sources, filters, I/O
examples/
  triangle.rs      # Basic PolyData + render window
  shapes.rs        # Multiple geometric primitives
  isosurface.rs    # Marching cubes isosurface extraction
  scalar_colors.rs # Elevation filter + colormap visualization
  showcase.rs      # PBR, transparency, axes, scalar bar, keyboard interaction
  pipeline_demo.rs # Filter pipeline + multi-format I/O + topology analysis
  volume.rs        # GPU volume rendering with transfer function
  bench_filters.rs # Performance benchmarks (normals, marching cubes, etc.)
VTK/               # C++ 9.6.0 reference source (read-only)
```

**Dependency graph:** `vtk-render-wgpu → vtk-render → vtk-data → vtk-types`, `vtk-filters → vtk-data + vtk-io-*`, `vtk-io-{legacy,stl,obj,xml,ply,gltf,ensight,xdmf} → vtk-data → vtk-types`

### Key vtk-filters modules

**Sources (40):** sphere, cube, cone, cylinder, plane, arrow, disk, line, point_source, regular_polygon, arc, superquadric, platonic_solid, frustum, parametric, bounding_box_source, axes, torus, helix, ellipsoid, spring, capsule, geodesic_sphere, grid, text_3d, wavelet, circle, mobius, star, noise_field, ring, klein_bottle, trefoil_knot, cross, boy_surface, spiral, icosphere, mobius_strip, gear, grid_2d

**Infrastructure:** pipeline (lazy evaluation + caching), convert (dataset conversions), topology (manifold/euler/boundary analysis), merge (combine meshes), selection_extract (apply selections), io_utils (auto-format read/write)

**Filters (570+):** normals, triangulate, append, clean, transform, marching_cubes, clip, elevation, threshold, decimate, smooth, subdivide, warp, connectivity, extract_surface, and 555+ more covering mesh processing, image processing, and data manipulation.

## Key Design Decisions

- **Ownership over refcounting** — no `Arc` by default; Rust ownership replaces VTK's `vtkObjectBase` reference counting.
- **Enum-based type erasure** — `AnyDataArray` is an enum over `DataArray<f32>`, `DataArray<f64>`, etc. (closed set of scalar types). Prefer this over `Box<dyn Trait>`.
- **Traits over inheritance** — `DataObject` and `DataSet` traits replace VTK's class hierarchy. `DataSet` has default methods: `center()`, `diagonal()`, `is_empty()`.
- **Pipeline system** — `Pipeline` struct with lazy evaluation, caching, and invalidation. Builder API with `with_normals()`, `with_decimate()`, etc.
- **Scalar visualization** — `Actor` supports `Coloring::ScalarMap` with `ColorMap` (15 presets + `by_name()` lookup). Active scalars from point data are mapped through the color map.
- **PBR rendering** — Blinn-Phong and Cook-Torrance PBR (metallic/roughness) selectable per-actor via `material.pbr` flag.
- **wgpu rendering** — 4x MSAA, per-actor model matrix (position/scale), 6 clip planes, flat shading (dpdx/dpdy), backface culling, silhouette edges.
- **VTK legacy format** uses version 4.2 (legacy cell format: `npts id0 id1 ...`). Binary is big-endian.
- **Prelude modules** — `use vtk_data::prelude::*` and `use vtk_render::prelude::*` for quick starts.
- **WASM compatible** — all non-GPU crates compile for wasm32-unknown-unknown.

## VTK C++ Reference Architecture

The C++ source at `VTK/` is organized as a modular toolkit. Key reference files:

- `VTK/Common/DataModel/vtkPolyData.h` — 4-CellArray structure (verts/lines/polys/strips)
- `VTK/Common/DataModel/vtkCellArray.h` — offsets+connectivity design
- `VTK/IO/Legacy/vtkDataWriter.cxx` — legacy format output details
- `VTK/IO/Legacy/vtkDataReader.cxx` — legacy format parsing
- `VTK/Filters/Sources/vtkSphereSource.h` — geometry source pattern
