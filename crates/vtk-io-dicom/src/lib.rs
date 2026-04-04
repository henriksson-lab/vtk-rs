//! Basic DICOM image reader for vtk-rs (no external dependencies).
//!
//! Reads uncompressed DICOM files (.dcm) with explicit VR little-endian
//! transfer syntax. Extracts pixel data as ImageData.
//!
//! Supports: 8-bit and 16-bit grayscale, single-frame images.
//! For full DICOM support (compressed, multi-frame), use the `dicom-object` crate.

mod reader;

pub use reader::{read_dicom, DicomInfo};
