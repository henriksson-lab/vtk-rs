//! MPI-based distributed parallel computing for vtk-rs.
//!
//! Feature-gated behind `mpi` — requires system MPI library:
//!
//! ```bash
//! # Ubuntu/Debian
//! apt install libopenmpi-dev
//!
//! # Fedora
//! dnf install openmpi-devel
//! ```
//!
//! # Usage
//!
//! ```toml
//! vtk-parallel = { path = "crates/vtk-parallel", features = ["mpi"] }
//! ```
//!
//! # Design
//!
//! Provides distributed data decomposition, ghost cell exchange, and
//! collective operations on PolyData/ImageData across MPI ranks.

// Always available: data decomposition types and serial stubs
pub mod decomposition;
pub mod ghost_exchange;
pub mod collective;

// MPI backend (feature-gated)
#[cfg(feature = "mpi")]
pub mod mpi_backend;
