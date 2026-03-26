mod camera;
mod scene;
mod renderer;
pub mod color_map;
mod light;
mod material;

pub use camera::Camera;
pub use scene::{Actor, Coloring, Representation, Scene};
pub use renderer::Renderer;
pub use color_map::ColorMap;
pub use light::{Light, LightType};
pub use material::Material;
