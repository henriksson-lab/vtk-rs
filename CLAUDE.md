# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

vtk-rs is a pure Rust reimplementation of VTK (The Visualization Toolkit). The `VTK/` directory contains the C++ 9.6.0 source as reference — this is **not** an FFI bindings project.

## Build / Test / Lint

```bash
cargo build                              # build all crates
cargo test --workspace --exclude vtk-python  # run all tests (~2350+)
cargo test -p vtk-data                   # test a single crate
cargo test -p vtk-filters test_normals   # run a single test by name
cargo clippy --workspace -- -D warnings  # lint
cargo run --example showcase             # PBR, transparency, axes, scalar bar, keyboard controls
cargo run --release --example bench_filters  # performance benchmarks
```

## Workspace Structure

The workspace has ~48 crates organized in layers:

**Core:** `vtk-types` → `vtk-data` — foundational types and data model
**Filters:** `vtk-filters` (main) + 18 split crates (`vtk-filters-image`, `vtk-filters-mesh`, `vtk-filters-extract`, `vtk-filters-transform`, `vtk-filters-subdivide`, `vtk-filters-clip`, `vtk-filters-smooth`, `vtk-filters-cell`, `vtk-filters-points`, `vtk-filters-statistics`, `vtk-filters-texture`, `vtk-filters-flow`, `vtk-filters-boolean`, `vtk-filters-grid`, `vtk-filters-data`, `vtk-filters-distance`, `vtk-filters-normals`, `vtk-filters-geometry`) + `vtk-filters-gpu` (wgpu compute)
**I/O:** `vtk-io-{legacy,stl,obj,xml,ply,off,dxf,geojson,csv,byu,las,facet,segy,tecplot,fits,gltf,ensight,xdmf}`
**Rendering:** `vtk-render` (backend-agnostic) → `vtk-render-wgpu` (GPU backend)
**Bindings:** `vtk-python` (PyO3)
**Umbrella:** `vtk` — re-exports `vtk-types`, `vtk-data`, `vtk-filters`, `vtk-render` with `use vtk::prelude::*`

**Dependency graph:** `vtk-render-wgpu → vtk-render → vtk-data → vtk-types`, `vtk-filters → vtk-data + vtk-io-*`, all `vtk-io-*` → `vtk-data → vtk-types`

### Key vtk-filters modules

**Sources (~42):** sphere, cube, cone, cylinder, plane, arrow, disk, line, point_source, regular_polygon, arc, superquadric, platonic_solid, frustum, parametric, bounding_box_source, axes, torus, helix, ellipsoid, spring, capsule, geodesic_sphere, grid, text_3d, wavelet, circle, mobius, star, noise_field, ring, klein_bottle, trefoil_knot, cross, boy_surface, spiral, icosphere, mobius_strip, gear, grid_2d, earth, sector

**Infrastructure:** pipeline (lazy evaluation + caching), convert (dataset conversions), topology (manifold/euler/boundary analysis), merge (combine meshes), selection_extract (apply selections), io_utils (auto-format read/write)

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
