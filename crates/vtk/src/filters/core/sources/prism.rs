//! Triangular and hexagonal prism geometry sources.

use crate::data::{CellArray, Points, PolyData};

/// Generate a triangular prism (3 rectangular sides + 2 triangular caps).
pub fn triangular_prism(radius: f64, height: f64, resolution: usize) -> PolyData {
    let n = resolution.max(3);
    let mut points = Points::<f64>::new();
    let mut polys = CellArray::new();

    // Bottom and top vertices
    for layer in 0..2 {
        let z = if layer == 0 { 0.0 } else { height };
        for i in 0..n {
            let angle = 2.0 * std::f64::consts::PI * i as f64 / n as f64;
            points.push([radius * angle.cos(), radius * angle.sin(), z]);
        }
    }

    // Bottom cap (fan)
    for i in 1..n - 1 {
        polys.push_cell(&[0, (i + 1) as i64, i as i64]);
    }
    // Top cap
    for i in 1..n - 1 {
        polys.push_cell(&[n as i64, (n + i) as i64, (n + i + 1) as i64]);
    }
    // Side quads
    for i in 0..n {
        let i0 = i as i64;
        let i1 = ((i + 1) % n) as i64;
        let i2 = i1 + n as i64;
        let i3 = i0 + n as i64;
        polys.push_cell(&[i0, i1, i2]);
        polys.push_cell(&[i0, i2, i3]);
    }

    let mut mesh = PolyData::new();
    mesh.points = points;
    mesh.polys = polys;
    mesh
}

/// Generate a hexagonal prism.
pub fn hexagonal_prism(radius: f64, height: f64) -> PolyData {
    triangular_prism(radius, height, 6)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn tri_prism() {
        let p = triangular_prism(1.0, 2.0, 3);
        assert_eq!(p.points.len(), 6); // 3 bottom + 3 top
        assert!(p.polys.num_cells() > 0);
    }

    #[test]
    fn hex_prism() {
        let p = hexagonal_prism(1.0, 2.0);
        assert_eq!(p.points.len(), 12); // 6 bottom + 6 top
    }
}
