use std::f64::consts::PI;

use crate::data::{CellArray, DataArray, Points, PolyData};

/// Parameters for generating Boy's surface (an immersion of the real projective plane in 3D).
pub struct BoySurfaceParams {
    pub center: [f64; 3],
    pub radius: f64,
    pub resolution: usize,
}

impl Default for BoySurfaceParams {
    fn default() -> Self {
        Self {
            center: [0.0, 0.0, 0.0],
            radius: 1.0,
            resolution: 32,
        }
    }
}

/// Generate Boy's surface using Apéry's parametrization.
///
/// The parametrization maps (u, v) with u in [0, PI] and v in [0, PI] to 3D,
/// where the surface is an immersion of the real projective plane.
pub fn boy_surface(params: &BoySurfaceParams) -> PolyData {
    let n = params.resolution.max(3);
    let [cx, cy, cz] = params.center;
    let r = params.radius;

    let mut points = Points::new();
    let mut normals_data = DataArray::<f64>::new("Normals", 3);
    let mut polys = CellArray::new();

    // Apéry's parametrization of Boy's surface:
    // u in [0, PI], v in [0, PI]
    // We use the standard form with complex coordinates.
    let eps = 1e-10;

    for j in 0..=n {
        let v = PI * j as f64 / n as f64;
        for i in 0..=n {
            let u = PI * i as f64 / n as f64;

            let (x, y, z) = apery(u, v);

            points.push([cx + r * x, cy + r * y, cz + r * z]);

            // Compute normal via finite differences
            let du = eps;
            let dv = eps;
            let (x1, y1, z1) = apery(u + du, v);
            let (x2, y2, z2) = apery(u, v + dv);
            let tu = [(x1 - x) / du, (y1 - y) / du, (z1 - z) / du];
            let tv = [(x2 - x) / dv, (y2 - y) / dv, (z2 - z) / dv];
            let nx = tu[1] * tv[2] - tu[2] * tv[1];
            let ny = tu[2] * tv[0] - tu[0] * tv[2];
            let nz = tu[0] * tv[1] - tu[1] * tv[0];
            let len = (nx * nx + ny * ny + nz * nz).sqrt().max(1e-12);
            normals_data.push_tuple(&[nx / len, ny / len, nz / len]);
        }
    }

    // Create quads connecting the grid
    let cols = n + 1;
    for j in 0..n {
        for i in 0..n {
            let p0 = (j * cols + i) as i64;
            let p1 = (j * cols + i + 1) as i64;
            let p2 = ((j + 1) * cols + i + 1) as i64;
            let p3 = ((j + 1) * cols + i) as i64;
            polys.push_cell(&[p0, p1, p2, p3]);
        }
    }

    let mut pd = PolyData::new();
    pd.points = points;
    pd.polys = polys;
    pd.point_data_mut().add_array(normals_data.into());
    pd.point_data_mut().set_active_normals("Normals");
    pd
}

/// Apéry's parametrization of Boy's surface.
fn apery(u: f64, v: f64) -> (f64, f64, f64) {
    let cu = u.cos();
    let su = u.sin();
    let cv = v.cos();
    let _sv = v.sin();
    let s2v = (2.0 * v).sin();
    let _c2v = (2.0 * v).cos();
    let _s2u = (2.0 * u).sin();

    let sqrt2 = std::f64::consts::SQRT_2;

    // Apéry's parametrization
    let denom = 2.0 - sqrt2 * s2v * (3.0 * u).sin();

    let x = (sqrt2 * cv * cv * (2.0 * u).cos() + cu * s2v) / denom;
    let y = (sqrt2 * cv * cv * (2.0 * u).sin() - su * s2v) / denom;
    let z = 3.0 * cv * cv / denom;

    (x, y, z)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn default_boy_surface() {
        let pd = boy_surface(&BoySurfaceParams::default());
        let n = 32usize;
        assert_eq!(pd.points.len(), (n + 1) * (n + 1));
        assert_eq!(pd.polys.num_cells(), n * n);
        assert!(pd.point_data().normals().is_some());
    }

    #[test]
    fn minimal_boy_surface() {
        let pd = boy_surface(&BoySurfaceParams {
            resolution: 3,
            ..Default::default()
        });
        assert!(pd.points.len() > 0);
        assert!(pd.polys.num_cells() > 0);
    }

    #[test]
    fn custom_center_and_radius() {
        let pd = boy_surface(&BoySurfaceParams {
            center: [1.0, 2.0, 3.0],
            radius: 2.0,
            resolution: 4,
            ..Default::default()
        });
        assert!(pd.points.len() > 0);
    }
}
