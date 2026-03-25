use thiserror::Error;

#[derive(Debug, Error)]
pub enum VtkError {
    #[error("I/O error: {0}")]
    Io(#[from] std::io::Error),

    #[error("parse error: {0}")]
    Parse(String),

    #[error("invalid data: {0}")]
    InvalidData(String),

    #[error("unsupported feature: {0}")]
    Unsupported(String),

    #[error("index out of bounds: {index} (size: {size})")]
    IndexOutOfBounds { index: usize, size: usize },
}
