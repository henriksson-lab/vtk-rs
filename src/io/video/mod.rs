#![allow(unexpected_cfgs)]
mod frame;
mod ppm_writer;
#[cfg(feature = "ffmpeg")]
mod ffmpeg_writer;
pub use frame::FrameSequence;
pub use ppm_writer::write_ppm_sequence;
#[cfg(feature = "ffmpeg")]
pub use ffmpeg_writer::write_video;
