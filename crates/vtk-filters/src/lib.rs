//! Geometry sources, processing filters, and utilities for vtk-rs.
//!
//! # Quick Start
//!
//! ```
//! use vtk_data::PolyData;
//! use vtk_filters::sources::sphere::{sphere, SphereParams};
//! use vtk_filters::normals::compute_normals;
//! use vtk_filters::elevation::elevation_z;
//! use vtk_filters::pipeline::Pipeline;
//!
//! // Generate a sphere and process it
//! let src = sphere(&SphereParams::default());
//! let mut pipe = Pipeline::new(src)
//!     .with_normals()
//!     .with_elevation_z();
//! let result = pipe.output();
//! assert!(result.points.len() > 0);
//! assert!(result.point_data().scalars().is_some());
//! ```

pub mod sources;

// Sub-crate re-exports
pub use vtk_filters_image as image;
pub use vtk_filters_mesh as mesh;
pub use vtk_filters_extract as extract;
pub use vtk_filters_transform as transform_crate;
pub use vtk_filters_subdivide as subdivide_crate;
pub use vtk_filters_clip as clip_crate;
pub use vtk_filters_smooth as smooth_crate;
pub use vtk_filters_cell as cell;
pub use vtk_filters_points as points;
pub use vtk_filters_statistics as statistics;
pub use vtk_filters_texture as texture;
pub use vtk_filters_flow as flow;
pub use vtk_filters_boolean as boolean_crate;
pub use vtk_filters_grid as grid;
pub use vtk_filters_data as data;
pub use vtk_filters_distance as distance;
pub use vtk_filters_normals as normals_crate;
pub use vtk_filters_geometry as geometry;

// Modules kept in vtk-filters (have cross-references or are infrastructure)
pub mod pipeline;
pub mod convert;
pub mod selection_extract;
pub mod topology;
pub mod io_utils;
pub mod merge;
pub mod append;
pub mod quick;
pub mod generic_filters;
pub mod extract_geometry;
pub mod parallel_pipeline;
pub mod mmap_data;
pub mod extract_largest;
pub mod mass_properties;
pub mod flying_edges;
pub mod decimate;
pub mod process_id_scalars;
pub mod procedural_displacement;
pub mod extract_block;
pub mod contour;
pub mod aggregate_dataset;
pub mod marching_cubes;
pub mod slice;
pub mod piece_request;
pub mod plugin;

// Re-export commonly used modules from sub-crates for backward compatibility
pub use vtk_filters_normals::normals;
pub use vtk_filters_geometry::triangulate;
pub use vtk_filters_geometry::clean;
pub use vtk_filters_transform::transform;
pub use vtk_filters_clip::clip;
pub use vtk_filters_geometry::elevation;
pub use vtk_filters_clip::threshold;
pub use vtk_filters_smooth::smooth;
pub use vtk_filters_subdivide::subdivide;
pub use vtk_filters_transform::warp;
pub use vtk_filters_geometry::connectivity;
pub use vtk_filters_extract::extract_surface;
pub use vtk_filters_geometry::glyph;
pub use vtk_filters_geometry::tube;
pub use vtk_filters_geometry::feature_edges;
pub use vtk_filters_transform::reflect;
pub use vtk_filters_cell::shrink;
pub use vtk_filters_extract::extract_edges;
pub use vtk_filters_normals::orient;
pub use vtk_filters_extract::extract_component;
