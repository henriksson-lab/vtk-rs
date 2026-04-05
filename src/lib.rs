//! vtk-rs — a pure Rust visualization toolkit.

// Always available
pub mod types;
pub mod data;
pub mod filters;

// I/O (feature-gated)
#[cfg(feature = "io-common")]
pub mod io;

// Rendering (feature-gated)
#[cfg(feature = "render")]
pub mod render;
#[cfg(feature = "render-wgpu")]
pub mod render_wgpu;

// Parallel (feature-gated)
#[cfg(feature = "parallel")]
pub mod parallel;

// Convenience re-exports (always available)
pub use types::{Scalar, ScalarType, VtkError, BoundingBox, CellType};
pub use data::{DataArray, AnyDataArray, CellArray, Points, FieldData, DataSetAttributes,
    DataObject, DataSet, PolyData, ImageData, UnstructuredGrid, RectilinearGrid,
    StructuredGrid, MultiBlockDataSet, Table, KdTree, Selection};

pub mod prelude {
    pub use crate::data::prelude::*;
    pub use crate::{VtkError, BoundingBox, CellType};
}
