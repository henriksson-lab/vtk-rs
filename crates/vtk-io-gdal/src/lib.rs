//! GDAL-based geospatial raster and vector I/O for vtk-rs.
//!
//! Feature-gated behind `gdal` — requires system GDAL library:
//!
//! ```bash
//! apt install libgdal-dev    # Ubuntu/Debian
//! dnf install gdal-devel     # Fedora
//! brew install gdal           # macOS
//! ```
//!
//! # Features
//!
//! - **Raster**: read GeoTIFF, DEM, satellite imagery → `ImageData`
//! - **Vector**: read Shapefile, GeoPackage, KML → `PolyData`
//!
//! # Example
//!
//! ```toml
//! vtk-io-gdal = { path = "crates/vtk-io-gdal", features = ["gdal"] }
//! ```

// Always available: geospatial types
pub mod types;

#[cfg(feature = "gdal")]
pub mod raster;
#[cfg(feature = "gdal")]
pub mod vector;
