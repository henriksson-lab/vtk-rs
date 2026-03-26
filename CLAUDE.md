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
  vtk-data/        # DataArray, CellArray, Points, FieldData, DataSetAttributes, PolyData, ImageData, UnstructuredGrid, RectilinearGrid, StructuredGrid, MultiBlockDataSet, Table, KdTree, Selection
  vtk-filters/     # Sources and filters — see below
  vtk-io-legacy/   # VTK legacy format (.vtk) ASCII/binary reader and writer
  vtk-io-stl/      # STL format ASCII/binary reader and writer
  vtk-io-obj/      # Wavefront OBJ format reader and writer
  vtk-io-xml/      # VTK XML format reader and writer (.vtp, .vtu, .vti, .vtr, .vts)
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

**Sources:** sphere, cube, cone, cylinder, plane, arrow, disk, line, point_source, regular_polygon, arc, superquadric, platonic_solid, frustum, parametric_function, bounding_box_source, axes, torus, helix, ellipsoid, spring, capsule, geodesic_sphere, grid, text_3d, wavelet, circle, mobius
**Filters:** normals, triangulate, append, clean, transform, marching_cubes, clip, elevation, threshold, decimate, smooth, subdivide, warp, connectivity, extract_surface, slice, glyph, tube, feature_edges, reflect, shrink, contour, cell_centers, extract_edges, orient, attribute_convert, generate_ids, mask_points, curvatures, gradient, mass_properties, densify, hull, texture_map, delaunay_2d, spline, clip_data_set, windowed_sinc_smooth, probe, stream_tracer, voxel_modeller, sample_function, integrate_attributes, distance_poly_data, implicit_modeller, tensor_glyph, resample, calculator, select_enclosed_points, extract_cells, icp, extract_points, signed_distance, cell_size, extrude, cell_quality, reverse_sense, strip, interpolate, fill_holes, center_of_mass, project_points, subdivide_midpoint, random_attributes, image_to_poly_data, flying_edges, hausdorff, color_transfer, extrude_normals, cell_normals, merge_points, bounding_box_filter, data_array_math, edge_lengths, extract_region, warp_implicit, delaunay_3d, boolean, ribbon, rotation_extrude, geometry_filter, banded_contour, ruled_surface, multi_threshold, data_set_triangulate, point_sampler, outline, histogram, attribute_smooth, dicer, intersection_poly_data, quadric_clustering, adaptive_subdivide, poly_data_to_image_data, normal_estimation, array_rename, vertex_glue, iso_volume, surface_nets, voronoi_2d, median_smooth, extract_component, pass_arrays, winding_number, decimate_boundary, scalar_range, point_merge, rectilinear_to_poly_data, structured_to_poly_data, offset_surface, table_to_poly_data, image_threshold, image_gaussian_smooth, image_gradient, image_dilate_erode, programmable_filter, field_data_to_attribute, image_resample, laplacian, aggregate, subdivide_butterfly, poly_data_distance, image_flip, cell_data_to_point_data_avg, compute_area, image_math, deep_copy, extract_surface_by_normal, select_by_scalar, image_crop, mesh_quality, resample_to_image, poly_data_bounds, clip_by_scalar, image_distance_transform, close_holes, image_statistics, sample_along_line, collapse_edges, image_pad, unstructured_to_poly_data, point_cloud_density, remesh, image_sobel, poly_line_to_strip, centroid_filter, image_median, triangulate_strips, image_connected_components, temporal_interpolator, separate_cells, image_normalize, poly_data_summary, image_logical, dual_mesh, spherical_coordinates, image_histogram_equalize, extract_largest, poly_data_normals_flip, image_stencil, coplanar_faces, split_by_array, image_variance, point_to_vertex, image_downsample, angle_between, sample_on_sphere, poly_data_transform_filter

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
