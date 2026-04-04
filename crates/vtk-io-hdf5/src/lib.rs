//! HDF5-based I/O formats for vtk-rs.
//!
//! All formats are feature-gated and require system HDF5/NetCDF libraries:
//!
//! ```bash
//! # Ubuntu/Debian
//! apt install libhdf5-dev libnetcdf-dev
//!
//! # Fedora
//! dnf install hdf5-devel netcdf-devel
//! ```
//!
//! # Features
//!
//! | Feature   | Format            | System dep       |
//! |-----------|-------------------|------------------|
//! | `exodus`  | Exodus II `.exo`  | `libhdf5-dev`    |
//! | `cgns`    | CGNS `.cgns`      | `libhdf5-dev`    |
//! | `amr`     | BoxLib/Chombo AMR | `libhdf5-dev`    |
//! | `netcdf`  | NetCDF `.nc`      | `libnetcdf-dev`  |
//! | `minc`    | MINC `.mnc`       | `libnetcdf-dev`  |
//! | `all`     | All of the above  | both             |
//!
//! # Example
//!
//! ```toml
//! vtk-io-hdf5 = { path = "crates/vtk-io-hdf5", features = ["exodus"] }
//! ```

#[cfg(feature = "exodus")]
pub mod exodus;
#[cfg(feature = "cgns")]
pub mod cgns;
#[cfg(feature = "amr")]
pub mod amr;
#[cfg(feature = "netcdf")]
pub mod netcdf_io;
#[cfg(feature = "minc")]
pub mod minc;

// Always-available types (no HDF5 needed)
pub mod types;
