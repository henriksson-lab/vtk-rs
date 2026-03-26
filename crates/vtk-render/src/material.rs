/// Surface material properties for Phong/PBR shading.
///
/// Analogous to VTK's `vtkProperty`. Controls how light interacts
/// with the surface.
#[derive(Debug, Clone)]
pub struct Material {
    /// Ambient light contribution. Default: 0.1
    pub ambient: f64,
    /// Diffuse light contribution. Default: 0.7
    pub diffuse: f64,
    /// Specular highlight contribution. Default: 0.3
    pub specular: f64,
    /// Specular power (shininess). Higher = tighter highlights. Default: 32.0
    pub specular_power: f64,
    /// Specular highlight color. Default: white [1, 1, 1]
    pub specular_color: [f32; 3],
    /// Edge color for wireframe/edge overlay. Default: black [0, 0, 0]
    pub edge_color: [f32; 3],
    /// Edge visibility (for edge overlay mode). Default: false
    pub edge_visibility: bool,
    /// Line width for wireframe/edge rendering. Default: 1.0
    pub line_width: f32,
    /// Point size for point rendering. Default: 3.0
    pub point_size: f32,
    /// Whether to use flat shading (per-face normals). Default: false
    pub flat_shading: bool,
    /// Backface culling. Default: false
    pub backface_culling: bool,
}

impl Default for Material {
    fn default() -> Self {
        Self {
            ambient: 0.1,
            diffuse: 0.7,
            specular: 0.3,
            specular_power: 32.0,
            specular_color: [1.0, 1.0, 1.0],
            edge_color: [0.0, 0.0, 0.0],
            edge_visibility: false,
            line_width: 1.0,
            point_size: 3.0,
            flat_shading: false,
            backface_culling: false,
        }
    }
}

impl Material {
    /// Create a default material.
    pub fn new() -> Self {
        Self::default()
    }

    /// Create a matte (no specular) material.
    pub fn matte() -> Self {
        Self {
            specular: 0.0,
            ..Default::default()
        }
    }

    /// Create a shiny/metallic material.
    pub fn shiny() -> Self {
        Self {
            specular: 0.8,
            specular_power: 128.0,
            diffuse: 0.5,
            ..Default::default()
        }
    }

    /// Create a flat-shaded material (per-face coloring).
    pub fn flat() -> Self {
        Self {
            flat_shading: true,
            ..Default::default()
        }
    }

    /// Create a material with edge overlay.
    pub fn with_edges(mut self) -> Self {
        self.edge_visibility = true;
        self
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn default_material() {
        let m = Material::default();
        assert_eq!(m.ambient, 0.1);
        assert_eq!(m.diffuse, 0.7);
        assert!(!m.flat_shading);
    }

    #[test]
    fn matte_no_specular() {
        let m = Material::matte();
        assert_eq!(m.specular, 0.0);
    }

    #[test]
    fn shiny_high_power() {
        let m = Material::shiny();
        assert_eq!(m.specular_power, 128.0);
    }

    #[test]
    fn with_edges_builder() {
        let m = Material::flat().with_edges();
        assert!(m.flat_shading);
        assert!(m.edge_visibility);
    }
}
