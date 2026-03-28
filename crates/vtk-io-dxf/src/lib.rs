//! AutoCAD DXF format reader and writer for vtk-rs.
//!
//! Supports reading/writing 3DFACE entities and LINE entities from ASCII DXF files.
//! This is a subset of the full DXF specification focused on mesh geometry.

mod writer;
mod reader;

pub use writer::DxfWriter;
pub use reader::{DxfReader, read_dxf_file, write_dxf_file};
