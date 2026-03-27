//! Core data structures for the vtk-rs visualization toolkit.
//!
//! This crate provides the fundamental data model types analogous to VTK's
//! data model: meshes, grids, arrays, and spatial data structures.
//!
//! # Quick Start
//!
//! ```
//! use vtk_data::prelude::*;
//!
//! // Create a triangle mesh
//! let mesh = PolyData::from_triangles(
//!     vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]],
//!     vec![[0, 1, 2]],
//! );
//! assert_eq!(mesh.points.len(), 3);
//!
//! // Create a scalar field on a grid
//! let grid = ImageData::from_function(
//!     [10, 10, 10], [0.1, 0.1, 0.1], [0.0, 0.0, 0.0],
//!     "density", |x, y, z| (x*x + y*y + z*z).sqrt(),
//! );
//!
//! // Create a data table
//! let table = Table::new()
//!     .with_column(AnyDataArray::F64(DataArray::from_vec("x", vec![1.0, 2.0], 1)));
//! ```
//!
//! # Key Types
//!
//! - [`PolyData`] — polygonal mesh (triangles, quads, lines, vertices)
//! - [`ImageData`] — regular grid with implicit coordinates
//! - [`UnstructuredGrid`] — mixed-cell mesh with explicit connectivity
//! - [`RectilinearGrid`] — axis-aligned grid with non-uniform spacing
//! - [`StructuredGrid`] — curvilinear grid with explicit points
//! - [`DataArray`] — typed tuple array for point/cell data
//! - [`Table`] — columnar data for analysis

pub mod prelude;
mod data_array;
mod cell_array;
mod points;
mod field_data;
mod attributes;
mod traits;
mod poly_data;
mod image_data;
mod unstructured_grid;
mod rectilinear_grid;
mod structured_grid;
mod multi_block;
mod table;
pub mod kd_tree;
pub mod octree;
mod selection;
pub mod cell_locator;
pub mod graph;
mod explicit_structured_grid;
mod hyper_tree_grid;
pub mod temporal;
pub mod table_stats;
pub mod poly_data_builder;
pub mod spatial_hash;
mod molecule;

pub use data_array::{DataArray, AnyDataArray, ArrayStatistics, DataArrayTupleIter};
pub use cell_array::CellArray;
pub use points::{Points, PointsIter};
pub use field_data::FieldData;
pub use attributes::DataSetAttributes;
pub use traits::{DataObject, DataSet};
pub use poly_data::PolyData;
pub use image_data::ImageData;
pub use unstructured_grid::UnstructuredGrid;
pub use rectilinear_grid::RectilinearGrid;
pub use structured_grid::StructuredGrid;
pub use multi_block::{MultiBlockDataSet, Block};
pub use table::Table;
pub use kd_tree::KdTree;
pub use selection::{Selection, SelectionNode, SelectionContentType, SelectionFieldType};
pub use octree::OctreePointLocator;
pub use cell_locator::CellLocator;
pub use graph::{Graph, Tree};
pub use explicit_structured_grid::ExplicitStructuredGrid;
pub use hyper_tree_grid::HyperTreeGrid;
pub use molecule::Molecule;
pub use poly_data_builder::PolyDataBuilder;
pub use temporal::TemporalDataSet;
