//! Object File Format (OFF) reader and writer for vtk-rs.
//!
//! OFF is a simple mesh format storing vertices and polygonal faces.
//! Supports ASCII OFF files with optional vertex colors (COFF).
//!
//! # Examples
//!
//! ```
//! use vtk_data::PolyData;
//! use vtk_io_off::{OffWriter, OffReader};
//!
//! let mesh = PolyData::from_triangles(
//!     vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]],
//!     vec![[0, 1, 2]],
//! );
//!
//! // Write
//! let mut buf = Vec::new();
//! OffWriter::new(&mut buf).write(&mesh).unwrap();
//!
//! // Read
//! let loaded = OffReader::new(&buf[..]).read().unwrap();
//! assert_eq!(loaded.points.len(), 3);
//! ```

mod writer;
mod reader;

pub use writer::OffWriter;
pub use reader::{OffReader, read_off_file, write_off_file};
