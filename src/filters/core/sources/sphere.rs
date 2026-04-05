use std::f64::consts::PI;

use crate::data::{CellArray, DataArray, Points, PolyData};

/// Parameters for generating a UV sphere.
pub struct SphereParams {
    pub center: [f64; 3],
    pub radius: f64,
    pub theta_resolution: usize,
    pub phi_resolution: usize,
}

impl Default for SphereParams {
    fn default() -> Self {
        Self {
            center: [0.0, 0.0, 0.0],
            radius: 0.5,
            theta_resolution: 16,
            phi_resolution: 16,
        }
    }
}

/// Generate a UV sphere as PolyData with normals.
pub fn sphere(params: &SphereParams) -> PolyData {
    let n_theta = params.theta_resolution.max(3);
    let n_phi = params.phi_resolution.max(3);
    let [cx, cy, cz] = params.center;
    let r = params.radius;

    let mut points = Points::new();
    let mut normals = DataArray::<f64>::new("Normals", 3);
    let mut polys = CellArray::new();

    // North pole
    points.push([cx, cy, cz + r]);
    normals.push_tuple(&[0.0, 0.0, 1.0]);

    // Interior rings (from top to bottom, excluding poles)
    for j in 1..n_phi {
        let phi = PI * j as f64 / n_phi as f64;
        let sp = phi.sin();
        let cp = phi.cos();

        for i in 0..n_theta {
            let theta = 2.0 * PI * i as f64 / n_theta as f64;
            let st = theta.sin();
            let ct = theta.cos();

            let nx = sp * ct;
            let ny = sp * st;
            let nz = cp;

            points.push([cx + r * nx, cy + r * ny, cz + r * nz]);
            normals.push_tuple(&[nx, ny, nz]);
        }
    }

    // South pole
    points.push([cx, cy, cz - r]);
    normals.push_tuple(&[0.0, 0.0, -1.0]);

    let south_pole = points.len() as i64 - 1;

    // North cap triangles (pole to first ring)
    for i in 0..n_theta {
        let next = (i + 1) % n_theta;
        polys.push_cell(&[0, (1 + i) as i64, (1 + next) as i64]);
    }

    // Body quads (ring to ring)
    for j in 0..n_phi - 2 {
        let ring_start = 1 + j * n_theta;
        let next_ring_start = 1 + (j + 1) * n_theta;

        for i in 0..n_theta {
            let next = (i + 1) % n_theta;
            polys.push_cell(&[
                (ring_start + i) as i64,
                (next_ring_start + i) as i64,
                (next_ring_start + next) as i64,
                (ring_start + next) as i64,
            ]);
        }
    }

    // South cap triangles (last ring to south pole)
    let last_ring_start = 1 + (n_phi - 2) * n_theta;
    for i in 0..n_theta {
        let next = (i + 1) % n_theta;
        polys.push_cell(&[
            (last_ring_start + i) as i64,
            south_pole,
            (last_ring_start + next) as i64,
        ]);
    }

    let mut pd = PolyData::new();
    pd.points = points;
    pd.polys = polys;
    pd.point_data_mut().add_array(normals.into());
    pd.point_data_mut().set_active_normals("Normals");
    pd
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn default_sphere() {
        let pd = sphere(&SphereParams::default());
        // poles (2) + rings (14 * 16)
        assert_eq!(pd.points.len(), 2 + 15 * 16);
        assert!(pd.polys.num_cells() > 0);
        assert!(pd.point_data().normals().is_some());
    }

    #[test]
    fn minimal_sphere() {
        let pd = sphere(&SphereParams {
            theta_resolution: 3,
            phi_resolution: 3,
            ..Default::default()
        });
        assert!(pd.points.len() > 0);
        assert!(pd.polys.num_cells() > 0);
    }
}
