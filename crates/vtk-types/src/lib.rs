//! Foundational types for the vtk-rs visualization toolkit.
//!
//! Provides scalar traits, cell types, error types, bounding boxes,
//! implicit functions, and higher-order cell evaluation.

mod scalar;
mod error;
mod bounds;
mod cell_type;
pub mod implicit;
pub mod higher_order;
pub mod color;
pub mod math;
pub mod progress;
pub mod sparse;
pub mod ode;
pub mod timer;

pub use scalar::{Scalar, ScalarType};
pub use error::VtkError;
pub use bounds::BoundingBox;
pub use cell_type::CellType;
pub use implicit::{ImplicitFunction, ImplicitPlane, ImplicitSphere, ImplicitBox};
