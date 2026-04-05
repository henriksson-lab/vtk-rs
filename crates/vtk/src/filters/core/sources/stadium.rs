//! Stadium (discorectangle / oblong) geometry source.

use crate::data::{CellArray, Points, PolyData};

/// Generate a stadium (rounded rectangle) in the XY plane.
///
/// A stadium is a rectangle with semicircular caps on both ends.
pub fn stadium(length: f64, radius: f64, resolution: usize) -> PolyData {
    let n = resolution.max(4);
    let half_n = n / 2;
    let half_len = length / 2.0;

    let mut points = Points::<f64>::new();
    let mut polys = CellArray::new();

    // Center point for fan triangulation
    points.push([0.0, 0.0, 0.0]); // index 0

    // Right semicircle
    for i in 0..=half_n {
        let angle = -std::f64::consts::FRAC_PI_2 + std::f64::consts::PI * i as f64 / half_n as f64;
        points.push([half_len + radius * angle.cos(), radius * angle.sin(), 0.0]);
    }

    // Left semicircle
    for i in 0..=half_n {
        let angle = std::f64::consts::FRAC_PI_2 + std::f64::consts::PI * i as f64 / half_n as f64;
        points.push([-half_len + radius * angle.cos(), radius * angle.sin(), 0.0]);
    }

    // Fan triangulation from center
    let n_outline = points.len() - 1;
    for i in 0..n_outline {
        let a = (i + 1) as i64;
        let b = if i + 1 < n_outline { (i + 2) as i64 } else { 1 };
        polys.push_cell(&[0, a, b]);
    }

    let mut mesh = PolyData::new();
    mesh.points = points;
    mesh.polys = polys;
    mesh
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn basic() {
        let s = stadium(2.0, 0.5, 16);
        assert!(s.points.len() > 10);
        assert!(s.polys.num_cells() > 5);
    }

    #[test]
    fn zero_length() {
        // Becomes a circle
        let s = stadium(0.0, 1.0, 16);
        assert!(s.polys.num_cells() > 0);
    }
}
