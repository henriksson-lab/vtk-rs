/// Axis-aligned bounding box in 3D space.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct BoundingBox {
    pub x_min: f64,
    pub x_max: f64,
    pub y_min: f64,
    pub y_max: f64,
    pub z_min: f64,
    pub z_max: f64,
}

impl BoundingBox {
    /// An empty (inverted) bounding box that will expand on the first point added.
    pub fn empty() -> Self {
        Self {
            x_min: f64::INFINITY,
            x_max: f64::NEG_INFINITY,
            y_min: f64::INFINITY,
            y_max: f64::NEG_INFINITY,
            z_min: f64::INFINITY,
            z_max: f64::NEG_INFINITY,
        }
    }

    /// Returns true if no points have been added.
    pub fn is_empty(&self) -> bool {
        self.x_min > self.x_max
    }

    /// Expand the bounding box to include the given point.
    pub fn expand(&mut self, point: [f64; 3]) {
        self.x_min = self.x_min.min(point[0]);
        self.x_max = self.x_max.max(point[0]);
        self.y_min = self.y_min.min(point[1]);
        self.y_max = self.y_max.max(point[1]);
        self.z_min = self.z_min.min(point[2]);
        self.z_max = self.z_max.max(point[2]);
    }

    /// Center of the bounding box.
    pub fn center(&self) -> [f64; 3] {
        [
            (self.x_min + self.x_max) * 0.5,
            (self.y_min + self.y_max) * 0.5,
            (self.z_min + self.z_max) * 0.5,
        ]
    }

    /// Diagonal length of the bounding box.
    pub fn diagonal_length(&self) -> f64 {
        let dx = self.x_max - self.x_min;
        let dy = self.y_max - self.y_min;
        let dz = self.z_max - self.z_min;
        (dx * dx + dy * dy + dz * dz).sqrt()
    }
}

impl Default for BoundingBox {
    fn default() -> Self {
        Self::empty()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn empty_bounds() {
        let bb = BoundingBox::empty();
        assert!(bb.is_empty());
    }

    #[test]
    fn expand_and_center() {
        let mut bb = BoundingBox::empty();
        bb.expand([0.0, 0.0, 0.0]);
        bb.expand([2.0, 4.0, 6.0]);
        assert!(!bb.is_empty());
        assert_eq!(bb.center(), [1.0, 2.0, 3.0]);
    }
}
