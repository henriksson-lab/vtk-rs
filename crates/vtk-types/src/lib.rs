mod scalar;
mod error;
mod bounds;
mod cell_type;

pub use scalar::{Scalar, ScalarType};
pub use error::VtkError;
pub use bounds::BoundingBox;
pub use cell_type::CellType;
