mod camera;
mod scene;
mod renderer;
pub mod color_map;

pub use camera::Camera;
pub use scene::{Actor, Coloring, Representation, Scene};
pub use renderer::Renderer;
pub use color_map::ColorMap;
