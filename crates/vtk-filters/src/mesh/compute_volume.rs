//! Compute volume of a closed triangulated mesh.

use vtk_data::PolyData;

/// Compute signed volume of a closed triangle mesh using divergence theorem.
pub fn signed_volume(mesh: &PolyData) -> f64 {
    let mut vol = 0.0;
    for cell in mesh.polys.iter() {
        if cell.len() < 3 { continue; }
        let a = mesh.points.get(cell[0] as usize);
        for i in 1..cell.len() - 1 {
            let b = mesh.points.get(cell[i] as usize);
            let c = mesh.points.get(cell[i + 1] as usize);
            vol += a[0] * (b[1] * c[2] - b[2] * c[1])
                 + a[1] * (b[2] * c[0] - b[0] * c[2])
                 + a[2] * (b[0] * c[1] - b[1] * c[0]);
        }
    }
    vol / 6.0
}

/// Compute unsigned volume (absolute value).
pub fn volume(mesh: &PolyData) -> f64 {
    signed_volume(mesh).abs()
}

/// Check if mesh has consistent winding (positive volume).
pub fn has_consistent_winding(mesh: &PolyData) -> bool {
    signed_volume(mesh) > 0.0
}

/// Compute centroid of a closed triangle mesh.
pub fn volume_centroid(mesh: &PolyData) -> [f64; 3] {
    let mut cx = 0.0;
    let mut cy = 0.0;
    let mut cz = 0.0;
    let mut vol = 0.0;
    for cell in mesh.polys.iter() {
        if cell.len() < 3 { continue; }
        let a = mesh.points.get(cell[0] as usize);
        for i in 1..cell.len() - 1 {
            let b = mesh.points.get(cell[i] as usize);
            let c = mesh.points.get(cell[i + 1] as usize);
            let v = a[0] * (b[1] * c[2] - b[2] * c[1])
                  + a[1] * (b[2] * c[0] - b[0] * c[2])
                  + a[2] * (b[0] * c[1] - b[1] * c[0]);
            vol += v;
            cx += (a[0] + b[0] + c[0]) * v;
            cy += (a[1] + b[1] + c[1]) * v;
            cz += (a[2] + b[2] + c[2]) * v;
        }
    }
    if vol.abs() < 1e-30 { return [0.0, 0.0, 0.0]; }
    [cx / (4.0 * vol), cy / (4.0 * vol), cz / (4.0 * vol)]
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::sources::cube::{cube, CubeParams};
    use crate::triangulate::triangulate;
    #[test]
    fn test_cube_volume() {
        let c = cube(&CubeParams::default());
        let tri = triangulate(&c);
        let v = volume(&tri);
        assert!((v - 1.0).abs() < 0.1, "volume = {v}");
    }
    #[test]
    fn test_centroid() {
        let c = cube(&CubeParams::default());
        let tri = triangulate(&c);
        let ctr = volume_centroid(&tri);
        assert!(ctr[0].abs() < 0.1);
        assert!(ctr[1].abs() < 0.1);
        assert!(ctr[2].abs() < 0.1);
    }
}
