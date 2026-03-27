use vtk_types::VtkError;

use crate::Scene;

/// Trait for rendering backends.
pub trait Renderer {
    /// Render the scene.
    fn render(&mut self, scene: &Scene) -> Result<(), VtkError>;

    /// Handle window resize.
    fn resize(&mut self, width: u32, height: u32);

    /// Render the scene to an RGBA image buffer (offscreen).
    ///
    /// Returns raw RGBA pixel data (4 bytes per pixel, row-major).
    /// Default implementation returns `Unsupported`.
    fn render_to_image(
        &mut self,
        _scene: &Scene,
        _width: u32,
        _height: u32,
    ) -> Result<Vec<u8>, VtkError> {
        Err(VtkError::Unsupported(
            "offscreen rendering not supported by this backend".into(),
        ))
    }
}
