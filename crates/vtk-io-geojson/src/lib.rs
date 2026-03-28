//! GeoJSON format reader and writer for vtk-rs.
//!
//! Supports Point, LineString, Polygon, MultiPoint, MultiLineString,
//! and MultiPolygon geometry types. Properties are stored as point/cell data.

mod writer;
mod reader;

pub use writer::write_geojson;
pub use reader::read_geojson;
