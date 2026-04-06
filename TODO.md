# TODO: Features Missing Tests

30 non-extra features remain untested (out of 396 total = 92% covered).

## Why these are untested

| Reason | Count | Features |
|--------|-------|----------|
| Exotic I/O (no test files) | 12 | dicom, ensight, fits, gdal, geojson, gltf, hdf5, segy, tecplot, video, xdmf, xml |
| Reeb graph / Graph type | 4 | force_directed_layout, reeb_join_split, reeb_simplify, reeb_skeleton |
| HyperTreeGrid type needed | 4 | amr_to_multi_block, hyper_tree_grid_depth_limiter, hyper_tree_grid_filters, hyper_tree_grid_slice |
| Internal infrastructure | 4 | generic_filters, io_utils, mmap_data, plugin |
| GPU (wgpu runtime) | 1 | context |
| Complex setup / private API | 5 | cell_data_to_point_data_avg, interpolated_velocity_field, order_statistics, seed_strategy, streak_line |

## Detailed List

### Exotic I/O Formats

**Reason:** Need sample test files (DICOM volumes, EnSight cases, FITS images, etc.)

- `I/O Formats/dicom`
- `I/O Formats/ensight`
- `I/O Formats/fits`
- `I/O Formats/gdal`
- `I/O Formats/geojson`
- `I/O Formats/gltf`
- `I/O Formats/hdf5`
- `I/O Formats/segy`
- `I/O Formats/tecplot`
- `I/O Formats/video`
- `I/O Formats/xdmf`
- `I/O Formats/xml`

### Reeb Graph / Graph Algorithms

**Reason:** Require `Graph` or `ReebGraph` type construction not exposed via simple API

- `Geometry/force_directed_layout`
- `Geometry/reeb_join_split`
- `Geometry/reeb_simplify`
- `Geometry/reeb_skeleton`

### HyperTreeGrid Filters

**Reason:** Require `HyperTreeGrid` construction which has complex AMR setup

- `Grid/amr_to_multi_block`
- `Grid/hyper_tree_grid_depth_limiter`
- `Grid/hyper_tree_grid_filters`
- `Grid/hyper_tree_grid_slice`

### Internal Infrastructure

**Reason:** Plugin system, memory-mapped I/O, internal dispatch — not user-facing algorithms

- `Core Filters/generic_filters`
- `Core Filters/io_utils`
- `Core Filters/mmap_data`
- `Core Filters/plugin`

### GPU Filters

**Reason:** Require wgpu device initialization at runtime

- `GPU/context`

### Complex Setup / Private API

**Reason:** Require private struct construction, multi-timestep data, or specialized inputs

- `Cell/cell_data_to_point_data_avg`
- `Flow/interpolated_velocity_field`
- `Statistics/order_statistics`
- `Geometry/seed_strategy`
- `Flow/streak_line`

## How to test these

- **Exotic I/O**: Add sample files to `tests/data/` and roundtrip test
- **HyperTreeGrid**: Implement `HyperTreeGrid::from_uniform()` constructor for testing
- **Reeb graph**: Expose `Graph` type or add `from_mesh()` constructor
- **GPU**: Feature-gate tests behind `#[cfg(feature = "filters-gpu")]` with wgpu init
- **Infrastructure**: These are internal and don't need user-facing tests
