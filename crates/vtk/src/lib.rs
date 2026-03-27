//! Convenience crate that re-exports the entire vtk-rs toolkit.
//!
//! # Quick Start
//!
//! ```
//! use vtk::prelude::*;
//!
//! // Generate a sphere with elevation coloring
//! let mesh = vtk::quick::sphere_with_elevation();
//! assert!(mesh.points.len() > 0);
//!
//! // Create a scene
//! let scene = Scene::from_poly_data(mesh);
//! assert_eq!(scene.num_actors(), 1);
//! ```

/// Re-export core types
pub use vtk_types;
pub use vtk_data;
pub use vtk_filters;
pub use vtk_render;

/// Quick source generation functions (no parameter structs needed).
pub use vtk_filters::quick;

/// I/O utilities: auto-detect format by extension.
///
/// ```no_run
/// vtk::io::write_poly_data(std::path::Path::new("mesh.stl"), &vtk::quick::sphere()).unwrap();
/// ```
pub mod io {
    pub use vtk_filters::io_utils::*;
}

/// Math utilities: vectors, interpolation, noise, coordinate transforms.
pub mod math {
    pub use vtk_types::math::*;
}

/// Color utilities: RGB/HSV/hex conversion, blending.
pub mod color {
    pub use vtk_types::color::*;
}

/// Prelude: import everything you need with `use vtk::prelude::*`.
pub mod prelude {
    pub use vtk_data::prelude::*;
    pub use vtk_render::prelude::*;
}

#[cfg(test)]
mod tests {
    use super::prelude::*;

    #[test]
    fn quick_start() {
        let mesh = crate::quick::sphere_with_elevation();
        let scene = Scene::from_poly_data(mesh);
        assert_eq!(scene.num_actors(), 1);
        assert!(scene.camera.distance() > 0.0);
    }

    #[test]
    fn full_pipeline() {
        use vtk_filters::pipeline::Pipeline;

        let mesh = crate::quick::sphere();
        let mut pipe = Pipeline::new(mesh)
            .with_normals()
            .with_elevation_z();
        let result = pipe.output();
        assert!(result.point_data().normals().is_some());
        assert!(result.point_data().scalars().is_some());
    }
}
