use vtk_types::VtkError;

use crate::Scene;

/// Trait for rendering backends.
pub trait Renderer {
    /// Render the scene.
    fn render(&mut self, scene: &Scene) -> Result<(), VtkError>;

    /// Handle window resize.
    fn resize(&mut self, width: u32, height: u32);
}
