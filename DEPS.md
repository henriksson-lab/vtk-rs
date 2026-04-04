# DEPS.md

Dependencies needed to implement remaining features.

Last updated: 2026-04-04 | 41 unchecked features remaining

---

## I/O Formats (15 items)

| Format | Rust Crate | System Package | Notes |
|--------|-----------|----------------|-------|
| Exodus II `.exo` | `hdf5` | `libhdf5-dev` | HDF5 C library via FFI |
| CGNS `.cgns` | `hdf5` | `libhdf5-dev` | CGNS is layered on HDF5 |
| NetCDF `.nc` | `netcdf` | `libnetcdf-dev` | C library FFI |
| Alembic `.abc` | None stable | Alembic C++ SDK | Would need custom C++ FFI |
| OpenVDB `.vdb` | None stable | OpenVDB C++ SDK | Would need custom C++ FFI |
| USD `.usd/.usda/.usdc` | None stable | USD C++ SDK | Massive C++ dependency |
| AMR formats | `hdf5` | `libhdf5-dev` | BoxLib/Chombo use HDF5 |
| CityGML `.gml` | `quick-xml` | **None** | Pure Rust XML parsing |
| DICOM `.dcm` | `dicom-object` | **None** | Pure Rust DICOM library |
| MINC `.mnc` | `netcdf` | `libnetcdf-dev` | NetCDF-based format |
| GDAL | `gdal` | `libgdal-dev` | Large C library |
| PDAL | None stable | PDAL C++ SDK | No good Rust bindings |
| OCCT/STEP/IGES | None stable | OpenCASCADE C++ SDK | Huge C++ dependency |
| FFMPEG | `ffmpeg-next` | `libavcodec-dev libavformat-dev` | C library FFI |
| MySQL/PostgreSQL/ODBC | `sqlx` | **None** | Pure Rust async SQL |

### Implementable now (no C deps)

```toml
quick-xml = "0.37"       # CityGML reader
dicom-object = "0.7"     # DICOM medical imaging
sqlx = { version = "0.8", features = ["runtime-tokio", "postgres", "mysql"] }
```

### Requires system libraries

```toml
hdf5 = "0.8"             # Exodus, CGNS, AMR (needs libhdf5-dev)
netcdf = "0.9"           # NetCDF, MINC (needs libnetcdf-dev)
gdal = "0.17"            # Geospatial (needs libgdal-dev)
ffmpeg-next = "7"        # Video I/O (needs libavcodec-dev)
```

```bash
# Ubuntu/Debian
apt install libhdf5-dev libnetcdf-dev libgdal-dev libavcodec-dev libavformat-dev

# Fedora
dnf install hdf5-devel netcdf-devel gdal-devel ffmpeg-devel
```

---

## Rendering (8 items)

| Feature | Rust Crate | System Package | Notes |
|---------|-----------|----------------|-------|
| Ray tracing | **None** | **None** | Pure Rust math (CPU) or wgpu compute (GPU) |
| Texture atlas | **None** | **None** | Rectangle packing algorithm |
| TrueType fonts | `fontdue` | **None** | Pure Rust font rasterizer |
| 2D rendering context | **None** | **None** | Vertex generation for shapes |
| VR/XR | `openxr` | OpenXR runtime | System runtime required |
| Remote rendering | `tokio` + `image` | **None** | WebSocket + PNG streaming |
| Parallel rendering | `rayon` (already dep) | **None** | CPU compositing |
| Path tracing | **None** | **None** | Pure Rust math |
| Subdivision surfaces | `wgpu` (already dep) | **None** | WGSL compute shader |
| Impostor rendering | **None** | **None** | Billboard quad generation |

### Recommended additions

```toml
fontdue = "0.9"           # TrueType font rasterization (pure Rust)
image = "0.25"            # PNG encoding for remote rendering
tokio = { version = "1", features = ["net", "rt-multi-thread"] }
```

---

## Infrastructure (10 items)

| Feature | Rust Crate | System Package | Notes |
|---------|-----------|----------------|-------|
| CI/CD | N/A | GitHub Actions | Config files, not code |
| Published docs | N/A | N/A | `cargo doc` + hosting |
| MPI | `mpi` | `libopenmpi-dev` | C library FFI |
| GPU compute filters | `wgpu` (already dep) | **None** | WGSL compute shaders |
| Memory-mapped data | `memmap2` | **None** | Pure Rust |
| Data streaming | **None** | **None** | Iterator/channel based |
| Plugin system | `libloading` | **None** | Pure Rust dynamic loading |
| In-situ | Depends on sim framework | Varies | Usually C/Fortran FFI |
| Web viewer | `wasm-bindgen` + `web-sys` | **None** | Compile to WASM |
| Jupyter | `pyo3` (already dep) | **None** | Python widget |
| Language bindings | `cbindgen`/`jni`/`napi` | **None** | Per-language |

### Recommended additions

```toml
memmap2 = "0.9"           # Memory-mapped file I/O
libloading = "0.8"        # Dynamic plugin loading
```

---

## Niche Algorithms (3 items)

| Feature | Deps | Notes |
|---------|------|-------|
| FiberSurface | **None** | Bivariate field topology — pure math |
| DG/CellGrid subsystem (×2) | **None** | High-order basis functions — pure math, ~5K lines |

---

## Priority for next implementation round

### Tier 1: Pure Rust, high impact, easy
1. CityGML reader (`quick-xml`)
2. DICOM reader (`dicom-object`)
3. TrueType fonts (`fontdue`)
4. Memory-mapped I/O (`memmap2`)
5. Impostor rendering (no deps)
6. Texture atlas (no deps)

### Tier 2: Pure Rust, medium effort
7. 2D rendering context (no deps)
8. GPU compute filters (wgpu compute shaders)
9. Remote rendering (`tokio` + image streaming)
10. FiberSurface (pure math)

### Tier 3: Requires system libraries
11. HDF5 I/O: Exodus + CGNS + AMR (`hdf5` crate)
12. NetCDF + MINC (`netcdf` crate)
13. Video I/O (`ffmpeg-next`)
14. Geospatial I/O (`gdal`)

### Tier 4: Major projects
15. DG/CellGrid subsystem (~5K lines)
16. VR/XR (`openxr`)
17. Path tracing (CPU ray tracer)
18. MPI distributed computing
