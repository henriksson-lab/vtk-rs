//! Video I/O for vtk-rs.
//!
//! Write rendered frame sequences as video files. Supports two backends:
//!
//! - **Built-in PPM sequence writer** (always available, no deps)
//! - **FFmpeg encoder** (feature-gated behind `ffmpeg`, requires system ffmpeg libs)
//!
//! # Usage
//!
//! ```no_run
//! use vtk_io_video::{FrameSequence, write_ppm_sequence};
//!
//! let mut seq = FrameSequence::new(640, 480);
//! // Add RGBA frames from render_to_image()
//! seq.add_frame(vec![0u8; 640 * 480 * 4]);
//! seq.add_frame(vec![0u8; 640 * 480 * 4]);
//!
//! // Write as numbered PPM files (always available)
//! write_ppm_sequence(&seq, std::path::Path::new("/tmp/frames")).unwrap();
//! ```

mod frame;
mod ppm_writer;
#[cfg(feature = "ffmpeg")]
mod ffmpeg_writer;

pub use frame::FrameSequence;
pub use ppm_writer::write_ppm_sequence;
#[cfg(feature = "ffmpeg")]
pub use ffmpeg_writer::write_video;
