use vtk_data::PolyData;

use crate::{Camera, ColorMap, Light, Material};

/// How to render the surface.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub enum Representation {
    #[default]
    Surface,
    Wireframe,
    Points,
}

/// How to determine vertex colors.
#[derive(Debug, Clone)]
pub enum Coloring {
    /// Uniform color for the entire actor.
    Solid([f32; 3]),
    /// Map active scalars from point data through a color map.
    ScalarMap {
        color_map: ColorMap,
        /// Scalar range [min, max]. If None, auto-computed from data.
        range: Option<[f64; 2]>,
    },
}

impl Default for Coloring {
    fn default() -> Self {
        Coloring::Solid([1.0, 1.0, 1.0])
    }
}

/// A renderable entity in the scene.
#[derive(Debug, Clone)]
pub struct Actor {
    pub data: PolyData,
    pub coloring: Coloring,
    pub opacity: f32,
    pub representation: Representation,
    pub material: Material,
}

impl Actor {
    pub fn new(data: PolyData) -> Self {
        Self {
            data,
            coloring: Coloring::default(),
            opacity: 1.0,
            representation: Representation::Surface,
            material: Material::default(),
        }
    }

    pub fn with_color(mut self, r: f32, g: f32, b: f32) -> Self {
        self.coloring = Coloring::Solid([r, g, b]);
        self
    }

    pub fn with_scalar_coloring(mut self, color_map: ColorMap, range: Option<[f64; 2]>) -> Self {
        self.coloring = Coloring::ScalarMap { color_map, range };
        self
    }
}

/// A 3D scene containing actors, lights, a camera, and background color.
#[derive(Debug, Clone)]
pub struct Scene {
    pub actors: Vec<Actor>,
    pub lights: Vec<Light>,
    pub camera: Camera,
    pub background: [f32; 4],
}

impl Default for Scene {
    fn default() -> Self {
        Self {
            actors: Vec::new(),
            lights: vec![Light::headlight()],
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

    pub fn add_light(&mut self, light: Light) {
        self.lights.push(light);
    }

    /// Remove all lights and start fresh.
    pub fn clear_lights(&mut self) {
        self.lights.clear();
    }
}
