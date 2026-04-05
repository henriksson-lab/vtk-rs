//! vtk-rs — a pure Rust visualization toolkit.
//!
//! Single-crate API that re-exports the entire vtk-rs workspace.
//!
//! # Quick Start
//! ```
//! use vtk::prelude::*;
//! let mesh = vtk::quick::sphere_with_elevation();
//! let scene = Scene::from_poly_data(mesh);
//! assert_eq!(scene.num_actors(), 1);
//! ```

pub use vtk_types::{Scalar, ScalarType, VtkError, BoundingBox, CellType};
pub use vtk_types::{ImplicitFunction, ImplicitPlane, ImplicitSphere, ImplicitBox};

pub mod types { pub use vtk_types::*; }

pub use vtk_data::{
    DataArray, AnyDataArray, CellArray, Points, FieldData, DataSetAttributes,
    DataObject, DataSet, PolyData, ImageData, UnstructuredGrid, RectilinearGrid,
    StructuredGrid, MultiBlockDataSet, Table, KdTree, Selection,
};
pub mod data { pub use vtk_data::*; }

pub use vtk_filters::quick;
pub use vtk_filters::pipeline;
pub mod filters { pub use vtk_filters::*; }

mod io_utils;
pub mod io {
    pub use crate::io_utils::*;
    pub mod legacy { pub use vtk_io_legacy::*; }
    pub mod stl { pub use vtk_io_stl::*; }
    pub mod obj { pub use vtk_io_obj::*; }
    pub mod xml { pub use vtk_io_xml::*; }
    pub mod ply { pub use vtk_io_ply::*; }
    pub mod gltf { pub use vtk_io_gltf::*; }
    pub mod off { pub use vtk_io_off::*; }
}

pub use vtk_render::{
    Camera, Scene, Actor, Coloring, Representation, Material, Light, LightType,
    ColorMap, Renderer, Fog, FogMode, ShadowConfig, BloomConfig, StereoConfig,
    StereoMode, Texture, TextureAtlas, ClipPlane, Viewport, SsaoConfig,
};
pub mod render { pub use vtk_render::*; }

pub mod math { pub use vtk_types::math::*; }
pub mod color { pub use vtk_types::color::*; }

pub mod prelude {
    pub use vtk_data::prelude::*;
    pub use vtk_render::prelude::*;
    pub use crate::{VtkError, BoundingBox, CellType};
}

#[cfg(test)]
mod tests {
    use super::prelude::*;

    #[test]
    fn quick_start() {
        let mesh = crate::quick::sphere_with_elevation();
        let scene = Scene::from_poly_data(mesh);
        assert_eq!(scene.num_actors(), 1);
    }

    #[test]
    fn full_pipeline() {
        use crate::pipeline::Pipeline;
        let mesh = crate::quick::sphere();
        let mut pipe = Pipeline::new(mesh).with_normals().with_elevation_z();
        let result = pipe.output();
        assert!(result.point_data().normals().is_some());
    }
}
