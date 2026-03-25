/// A color map that maps scalar values in [0, 1] to RGB colors.
#[derive(Debug, Clone)]
pub struct ColorMap {
    /// Control points: (parameter, [r, g, b]) where parameter is in [0, 1].
    points: Vec<(f64, [f32; 3])>,
}

impl ColorMap {
    /// Create a color map from control points. Points must be sorted by parameter.
    pub fn new(points: Vec<(f64, [f32; 3])>) -> Self {
        assert!(points.len() >= 2, "need at least 2 control points");
        Self { points }
    }

    /// Map a value in [0, 1] to an RGB color by linear interpolation.
    pub fn map(&self, t: f64) -> [f32; 3] {
        let t = t.clamp(0.0, 1.0);

        if t <= self.points[0].0 {
            return self.points[0].1;
        }
        if t >= self.points[self.points.len() - 1].0 {
            return self.points[self.points.len() - 1].1;
        }

        for i in 0..self.points.len() - 1 {
            let (t0, c0) = self.points[i];
            let (t1, c1) = self.points[i + 1];
            if t >= t0 && t <= t1 {
                let frac = if (t1 - t0).abs() > 1e-10 {
                    ((t - t0) / (t1 - t0)) as f32
                } else {
                    0.0
                };
                return [
                    c0[0] + frac * (c1[0] - c0[0]),
                    c0[1] + frac * (c1[1] - c0[1]),
                    c0[2] + frac * (c1[2] - c0[2]),
                ];
            }
        }

        self.points[self.points.len() - 1].1
    }

    /// Map a raw scalar value using the given range [min, max].
    pub fn map_value(&self, value: f64, min: f64, max: f64) -> [f32; 3] {
        let t = if (max - min).abs() > 1e-10 {
            (value - min) / (max - min)
        } else {
            0.5
        };
        self.map(t)
    }

    // --- Built-in color maps ---

    /// Blue → Cyan → Green → Yellow → Red (rainbow / jet).
    pub fn jet() -> Self {
        Self::new(vec![
            (0.0, [0.0, 0.0, 1.0]),
            (0.25, [0.0, 1.0, 1.0]),
            (0.5, [0.0, 1.0, 0.0]),
            (0.75, [1.0, 1.0, 0.0]),
            (1.0, [1.0, 0.0, 0.0]),
        ])
    }

    /// Cool to warm (blue → white → red), perceptually balanced.
    pub fn cool_to_warm() -> Self {
        Self::new(vec![
            (0.0, [0.231, 0.298, 0.753]),
            (0.5, [0.865, 0.865, 0.865]),
            (1.0, [0.706, 0.016, 0.150]),
        ])
    }

    /// Viridis (dark purple → blue → green → yellow).
    pub fn viridis() -> Self {
        Self::new(vec![
            (0.0, [0.267, 0.004, 0.329]),
            (0.25, [0.282, 0.140, 0.458]),
            (0.5, [0.127, 0.566, 0.551]),
            (0.75, [0.544, 0.774, 0.247]),
            (1.0, [0.993, 0.906, 0.144]),
        ])
    }

    /// Grayscale (black → white).
    pub fn grayscale() -> Self {
        Self::new(vec![
            (0.0, [0.0, 0.0, 0.0]),
            (1.0, [1.0, 1.0, 1.0]),
        ])
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn jet_endpoints() {
        let cm = ColorMap::jet();
        let c0 = cm.map(0.0);
        assert_eq!(c0, [0.0, 0.0, 1.0]); // blue
        let c1 = cm.map(1.0);
        assert_eq!(c1, [1.0, 0.0, 0.0]); // red
    }

    #[test]
    fn interpolation_midpoint() {
        let cm = ColorMap::grayscale();
        let c = cm.map(0.5);
        assert!((c[0] - 0.5).abs() < 1e-5);
    }

    #[test]
    fn map_value_range() {
        let cm = ColorMap::grayscale();
        let c = cm.map_value(50.0, 0.0, 100.0);
        assert!((c[0] - 0.5).abs() < 1e-5);
    }

    #[test]
    fn clamp_out_of_range() {
        let cm = ColorMap::jet();
        let below = cm.map(-1.0);
        let above = cm.map(2.0);
        assert_eq!(below, [0.0, 0.0, 1.0]);
        assert_eq!(above, [1.0, 0.0, 0.0]);
    }
}
