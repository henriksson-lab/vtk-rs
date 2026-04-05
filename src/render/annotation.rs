/// 3D annotation overlays for measurement and labeling.

/// A positioned text label in the scene.
#[derive(Debug, Clone)]
pub struct Label3D {
    /// World-space position where the label is anchored.
    pub position: [f64; 3],
    /// Text to display.
    pub text: String,
    /// Text color (RGB).
    pub color: [f32; 3],
    /// Text scale (NDC height). Default: 0.015
    pub scale: f32,
}

impl Label3D {
    pub fn new(position: [f64; 3], text: impl Into<String>) -> Self {
        Self {
            position,
            text: text.into(),
            color: [1.0, 1.0, 1.0],
            scale: 0.015,
        }
    }

    pub fn with_color(mut self, r: f32, g: f32, b: f32) -> Self {
        self.color = [r, g, b];
        self
    }
}

/// A distance measurement ruler between two 3D points.
#[derive(Debug, Clone)]
pub struct DistanceRuler {
    /// Start point in world coordinates.
    pub start: [f64; 3],
    /// End point in world coordinates.
    pub end: [f64; 3],
    /// Line color (RGB).
    pub color: [f32; 3],
    /// Whether to show the distance value as text.
    pub show_label: bool,
}

impl DistanceRuler {
    pub fn new(start: [f64; 3], end: [f64; 3]) -> Self {
        Self {
            start,
            end,
            color: [1.0, 1.0, 0.0],
            show_label: true,
        }
    }

    /// Compute the distance between start and end.
    pub fn distance(&self) -> f64 {
        let dx = self.end[0] - self.start[0];
        let dy = self.end[1] - self.start[1];
        let dz = self.end[2] - self.start[2];
        (dx * dx + dy * dy + dz * dz).sqrt()
    }

    /// Midpoint between start and end.
    pub fn midpoint(&self) -> [f64; 3] {
        [
            (self.start[0] + self.end[0]) / 2.0,
            (self.start[1] + self.end[1]) / 2.0,
            (self.start[2] + self.end[2]) / 2.0,
        ]
    }
}

/// An angle measurement protractor at three 3D points.
#[derive(Debug, Clone)]
pub struct AngleProtractor {
    /// First arm endpoint.
    pub point_a: [f64; 3],
    /// Vertex of the angle.
    pub vertex: [f64; 3],
    /// Second arm endpoint.
    pub point_b: [f64; 3],
    /// Line color (RGB).
    pub color: [f32; 3],
    /// Whether to show the angle value as text.
    pub show_label: bool,
}

impl AngleProtractor {
    pub fn new(a: [f64; 3], vertex: [f64; 3], b: [f64; 3]) -> Self {
        Self {
            point_a: a,
            vertex,
            point_b: b,
            color: [0.0, 1.0, 1.0],
            show_label: true,
        }
    }

    /// Compute the angle in degrees.
    pub fn angle_degrees(&self) -> f64 {
        crate::render::measurement::angle_at_vertex(self.point_a, self.vertex, self.point_b)
    }
}

/// Collection of all annotation overlays in a scene.
#[derive(Debug, Clone, Default)]
pub struct Annotations {
    pub labels: Vec<Label3D>,
    pub rulers: Vec<DistanceRuler>,
    pub protractors: Vec<AngleProtractor>,
}

impl Annotations {
    pub fn new() -> Self {
        Self::default()
    }

    pub fn add_label(&mut self, label: Label3D) {
        self.labels.push(label);
    }

    pub fn add_ruler(&mut self, ruler: DistanceRuler) {
        self.rulers.push(ruler);
    }

    pub fn add_protractor(&mut self, protractor: AngleProtractor) {
        self.protractors.push(protractor);
    }

    pub fn is_empty(&self) -> bool {
        self.labels.is_empty() && self.rulers.is_empty() && self.protractors.is_empty()
    }

    pub fn clear(&mut self) {
        self.labels.clear();
        self.rulers.clear();
        self.protractors.clear();
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn label() {
        let l = Label3D::new([1.0, 2.0, 3.0], "Hello").with_color(1.0, 0.0, 0.0);
        assert_eq!(l.text, "Hello");
        assert_eq!(l.color, [1.0, 0.0, 0.0]);
    }

    #[test]
    fn ruler_distance() {
        let r = DistanceRuler::new([0.0, 0.0, 0.0], [3.0, 4.0, 0.0]);
        assert!((r.distance() - 5.0).abs() < 1e-10);
        let mid = r.midpoint();
        assert!((mid[0] - 1.5).abs() < 1e-10);
    }

    #[test]
    fn protractor_angle() {
        let p = AngleProtractor::new([1.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 1.0, 0.0]);
        assert!((p.angle_degrees() - 90.0).abs() < 1e-6);
    }

    #[test]
    fn annotations_collection() {
        let mut a = Annotations::new();
        assert!(a.is_empty());
        a.add_label(Label3D::new([0.0; 3], "test"));
        a.add_ruler(DistanceRuler::new([0.0; 3], [1.0, 0.0, 0.0]));
        assert!(!a.is_empty());
        a.clear();
        assert!(a.is_empty());
    }
}
