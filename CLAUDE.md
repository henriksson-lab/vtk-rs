# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

vtk-rs is a pure Rust reimplementation of VTK (The Visualization Toolkit). The `VTK/` directory contains the C++ 9.6.0 source as reference — this is **not** an FFI bindings project.

## Build / Test / Lint

```bash
cargo build                          # build all crates
cargo test --workspace               # run all tests
cargo clippy --workspace -- -D warnings  # lint
cargo run --example triangle         # render a triangle (requires GPU)
cargo run --example shapes           # render sphere/cube/cone/cylinder/arrow
cargo run --example isosurface       # marching cubes on a gyroid field, writes gyroid.stl
cargo run --example scalar_colors    # elevation-colored sphere/cone/plane with colormaps
```

## Workspace Structure

```
crates/
  vtk-types/       # Scalar trait, ScalarType enum, CellType, VtkError, BoundingBox
  vtk-data/        # DataArray, CellArray, Points, FieldData, DataSetAttributes, PolyData, ImageData, UnstructuredGrid
  vtk-filters/     # Sources and filters — see below
  vtk-io-legacy/   # VTK legacy format (.vtk) ASCII/binary reader and writer
  vtk-io-stl/      # STL format ASCII/binary reader and writer
  vtk-io-obj/      # Wavefront OBJ format reader and writer
  vtk-io-xml/      # VTK XML format reader and writer (.vtp, .vtu, .vti)
  vtk-io-ply/      # Stanford PLY format ASCII and binary reader and writer
  vtk-render/      # Backend-agnostic: Camera, Scene, Actor, Renderer trait, ColorMap
  vtk-render-wgpu/ # wgpu implementation of the Renderer trait
examples/
  triangle.rs      # Creates PolyData, writes .vtk file, opens render window
  shapes.rs        # Renders sphere, cube, cone, cylinder, arrow side by side
  isosurface.rs    # Marching cubes isosurface extraction, writes gyroid.stl
  scalar_colors.rs # Elevation filter + colormap-based scalar visualization
VTK/               # C++ 9.6.0 reference source (read-only)
```

**Dependency graph:** `vtk-render-wgpu → vtk-render → vtk-data → vtk-types`, `vtk-filters → vtk-data → vtk-types`, `vtk-io-{legacy,stl,obj,xml,ply} → vtk-data → vtk-types`

### vtk-filters modules

**Sources:** sphere, cube, cone, cylinder, plane, arrow, disk, line, point_source, regular_polygon, arc
**Filters:** normals, triangulate, append, clean, transform, marching_cubes, clip, elevation, threshold, decimate, smooth, subdivide, warp, connectivity, extract_surface, slice, glyph, tube, feature_edges, reflect, shrink, contour, cell_centers, extract_edges, orient, attribute_convert, generate_ids, mask_points, curvatures, gradient, mass_properties, densify, hull, texture_map

`vtk-io-legacy` also supports ImageData (STRUCTURED_POINTS) and UnstructuredGrid read/write via `write_image_data` / `read_image_data` and `write_unstructured_grid` / `read_unstructured_grid`.

## Key Design Decisions

- **Ownership over refcounting** — no `Arc` by default; Rust ownership replaces VTK's `vtkObjectBase` reference counting.
- **Enum-based type erasure** — `AnyDataArray` is an enum over `DataArray<f32>`, `DataArray<f64>`, etc. (closed set of scalar types). Prefer this over `Box<dyn Trait>`.
- **Traits over inheritance** — `DataObject` and `DataSet` traits replace VTK's class hierarchy.
- **No pipeline system yet** — filters are plain functions. A pipeline crate can be added later.
- **Scalar visualization** — `Actor` supports `Coloring::ScalarMap` with `ColorMap` (jet, viridis, cool_to_warm, grayscale). Active scalars from point data are mapped through the color map.
- **wgpu rendering** — cross-platform WebGPU backend; rendering abstractions are backend-agnostic in `vtk-render`.
- **VTK legacy format** uses version 4.2 (legacy cell format: `npts id0 id1 ...`). Binary is big-endian.

## VTK C++ Reference Architecture

The C++ source at `VTK/` is organized as a modular toolkit. Key reference files:

- `VTK/Common/DataModel/vtkPolyData.h` — 4-CellArray structure (verts/lines/polys/strips)
- `VTK/Common/DataModel/vtkCellArray.h` — offsets+connectivity design
- `VTK/IO/Legacy/vtkDataWriter.cxx` — legacy format output details
- `VTK/IO/Legacy/vtkDataReader.cxx` — legacy format parsing
- `VTK/Filters/Sources/vtkSphereSource.h` — geometry source pattern
