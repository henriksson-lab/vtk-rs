use vtk_data::PolyData;

/// An instance of a glyph placed at a position with optional scale and color.
#[derive(Debug, Clone)]
pub struct GlyphInstance {
    /// World-space position.
    pub position: [f32; 3],
    /// Uniform scale factor.
    pub scale: f32,
    /// Override color (RGB). If None, uses the glyph's original coloring.
    pub color: Option<[f32; 3]>,
}

/// A set of glyph instances sharing the same template mesh.
///
/// For efficient rendering of many copies of the same geometry at different
/// positions (e.g., arrow glyphs at vector field points, sphere glyphs at
/// point locations).
#[derive(Debug, Clone)]
pub struct InstancedGlyphs {
    /// Template mesh to instance.
    pub template: PolyData,
    /// Per-instance transforms.
    pub instances: Vec<GlyphInstance>,
}

impl InstancedGlyphs {
    /// Create a new instanced glyph set.
    pub fn new(template: PolyData) -> Self {
        Self {
            template,
            instances: Vec::new(),
        }
    }

    /// Add an instance at the given position with unit scale.
    pub fn add(&mut self, position: [f32; 3]) {
        self.instances.push(GlyphInstance {
            position,
            scale: 1.0,
            color: None,
        });
    }

    /// Add an instance with position, scale, and color.
    pub fn add_with(&mut self, position: [f32; 3], scale: f32, color: [f32; 3]) {
        self.instances.push(GlyphInstance {
            position,
            scale,
            color: Some(color),
        });
    }

    /// Flatten all instances into a single PolyData by copying and transforming
    /// the template mesh for each instance.
    ///
    /// This is a CPU-side approach suitable for moderate instance counts.
    /// For very large counts, GPU instancing should be used.
    pub fn flatten(&self) -> PolyData {
        let tpl_npts = self.template.points.len();
        if tpl_npts == 0 || self.instances.is_empty() {
            return PolyData::new();
        }

        let mut result = PolyData::new();

        for inst in &self.instances {
            let base = result.points.len() as i64;

            // Copy and transform points
            for i in 0..tpl_npts {
                let p = self.template.points.get(i);
                result.points.push([
                    p[0] * inst.scale as f64 + inst.position[0] as f64,
                    p[1] * inst.scale as f64 + inst.position[1] as f64,
                    p[2] * inst.scale as f64 + inst.position[2] as f64,
                ]);
            }

            // Copy cells with offset indices
            for cell in self.template.polys.iter() {
                let offset_cell: Vec<i64> = cell.iter().map(|&id| id + base).collect();
                result.polys.push_cell(&offset_cell);
            }
        }

        result
    }

    /// Number of instances.
    pub fn len(&self) -> usize {
        self.instances.len()
    }

    /// Whether the instance set is empty.
    pub fn is_empty(&self) -> bool {
        self.instances.is_empty()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn flatten_instances() {
        let template = PolyData::from_triangles(
            vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]],
            vec![[0, 1, 2]],
        );

        let mut glyphs = InstancedGlyphs::new(template);
        glyphs.add([0.0, 0.0, 0.0]);
        glyphs.add([5.0, 0.0, 0.0]);
        glyphs.add([0.0, 5.0, 0.0]);

        let result = glyphs.flatten();
        assert_eq!(result.points.len(), 9); // 3 instances * 3 points
        assert_eq!(result.polys.num_cells(), 3);

        // Second instance should be offset by (5,0,0)
        let p = result.points.get(3);
        assert!((p[0] - 5.0).abs() < 1e-10);
    }

    #[test]
    fn scaled_instances() {
        let template = PolyData::from_triangles(
            vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]],
            vec![[0, 1, 2]],
        );

        let mut glyphs = InstancedGlyphs::new(template);
        glyphs.add_with([0.0, 0.0, 0.0], 2.0, [1.0, 0.0, 0.0]);

        let result = glyphs.flatten();
        let p = result.points.get(1);
        assert!((p[0] - 2.0).abs() < 1e-10); // scaled by 2
    }

    #[test]
    fn empty_instances() {
        let template = PolyData::from_triangles(
            vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]],
            vec![[0, 1, 2]],
        );
        let glyphs = InstancedGlyphs::new(template);
        assert!(glyphs.is_empty());
        let result = glyphs.flatten();
        assert_eq!(result.points.len(), 0);
    }
}
