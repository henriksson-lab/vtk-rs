//! Backend-agnostic rendering abstractions for vtk-rs.
//!
//! Provides [`Scene`], [`Actor`], [`Camera`], [`Material`], [`Light`],
//! and [`ColorMap`] types that describe what to render. The actual GPU
//! rendering is implemented by the `vtk-render-wgpu` crate.
//!
//! # Quick Start
//!
//! ```
//! use vtk_data::PolyData;
//! use vtk_render::{Actor, Scene, Material};
//!
//! let mesh = PolyData::from_triangles(
//!     vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]],
//!     vec![[0, 1, 2]],
//! );
//!
//! let scene = Scene::new()
//!     .with_actor(Actor::new(mesh).with_color(0.8, 0.2, 0.2))
//!     .with_background(0.1, 0.1, 0.1)
//!     .with_axes();
//! ```

pub mod prelude;
mod camera;
mod scene;
mod renderer;
pub mod color_map;
mod light;
mod material;
pub mod picker;
pub mod scalar_bar;
pub mod axes_widget;
pub mod silhouette;
pub mod fog;
pub mod texture;
pub mod lod;
pub mod instancing;
pub mod volume;
pub mod animation;
pub mod clip_plane;
pub mod scene_io;
pub mod measurement;
pub mod scene_json;
pub mod shadow;
pub mod skybox;
pub mod bloom;
pub mod annotation;
pub mod stereo;
pub mod screenshot;
pub mod transfer_editor;
pub mod viewport;
pub mod ssao;
pub mod dof;
pub mod environment;
pub mod axes_cube;

pub use camera::Camera;
pub use scene::{Actor, Coloring, Representation, Scene};
pub use renderer::Renderer;
pub use color_map::ColorMap;
pub use light::{Light, LightType};
pub use material::Material;
pub use picker::{pick, PickResult};
pub use scalar_bar::{ScalarBar, ScalarBarOrientation};
pub use axes_widget::AxesWidget;
pub use silhouette::SilhouetteConfig;
pub use fog::{Fog, FogMode};
pub use shadow::ShadowConfig;
pub use skybox::Skybox;
pub use bloom::BloomConfig;
pub use annotation::{Annotations, Label3D, DistanceRuler, AngleProtractor};
pub use stereo::{StereoConfig, StereoMode};
pub use texture::Texture;
pub use lod::{LodSet, LodLevel};
pub use instancing::{InstancedGlyphs, GlyphInstance};
pub use volume::{VolumeActor, TransferFunction};
pub use animation::{CameraAnimation, Track, Keyframe, Easing};
pub use clip_plane::ClipPlane;
pub use viewport::Viewport;
pub use ssao::SsaoConfig;
pub use dof::DofConfig;
pub use environment::EnvironmentMap;
pub use axes_cube::AxesCube;
