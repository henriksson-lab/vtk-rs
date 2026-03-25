use vtk_data::PolyData;

use crate::Camera;

/// How to render the surface.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub enum Representation {
    #[default]
    Surface,
    Wireframe,
    Points,
}

/// A renderable entity in the scene.
#[derive(Debug, Clone)]
pub struct Actor {
    pub data: PolyData,
    pub color: [f32; 3],
    pub opacity: f32,
    pub representation: Representation,
}

impl Actor {
    pub fn new(data: PolyData) -> Self {
        Self {
            data,
            color: [1.0, 1.0, 1.0],
            opacity: 1.0,
            representation: Representation::Surface,
        }
    }

    pub fn with_color(mut self, r: f32, g: f32, b: f32) -> Self {
        self.color = [r, g, b];
        self
    }
}

/// A 3D scene containing actors, a camera, and background color.
#[derive(Debug, Clone)]
pub struct Scene {
    pub actors: Vec<Actor>,
    pub camera: Camera,
    pub background: [f32; 4],
}

impl Default for Scene {
    fn default() -> Self {
        Self {
            actors: Vec::new(),
            camera: Camera::default(),
            background: [0.1, 0.1, 0.1, 1.0],
        }
    }
}

impl Scene {
    pub fn new() -> Self {
        Self::default()
    }

    pub fn add_actor(&mut self, actor: Actor) {
        self.actors.push(actor);
    }
}
