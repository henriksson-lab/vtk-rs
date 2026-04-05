use std::f64::consts::PI;

use crate::data::{CellArray, DataArray, Points, PolyData};

/// Parameters for generating a trefoil knot tube.
pub struct TrefoilKnotParams {
    pub center: [f64; 3],
    pub tube_radius: f64,
    /// Number of segments along the knot curve.
    pub resolution: usize,
    /// Number of segments around the tube cross-section.
    pub tube_resolution: usize,
}

impl Default for TrefoilKnotParams {
    fn default() -> Self {
        Self {
            center: [0.0, 0.0, 0.0],
            tube_radius: 0.1,
            resolution: 128,
            tube_resolution: 16,
        }
    }
}

/// Evaluate the trefoil knot curve at parameter t in [0, 2pi).
/// x = sin(t) + 2 sin(2t)
/// y = cos(t) - 2 cos(2t)
/// z = -sin(3t)
fn trefoil_point(t: f64) -> [f64; 3] {
    [
        t.sin() + 2.0 * (2.0 * t).sin(),
        t.cos() - 2.0 * (2.0 * t).cos(),
        -(3.0 * t).sin(),
    ]
}

fn normalize(v: [f64; 3]) -> [f64; 3] {
    let len = (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt().max(1e-12);
    [v[0] / len, v[1] / len, v[2] / len]
}

fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

fn sub(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

/// Generate a trefoil knot tube as PolyData.
pub fn trefoil_knot(params: &TrefoilKnotParams) -> PolyData {
    let n_curve = params.resolution.max(6);
    let n_tube = params.tube_resolution.max(3);
    let [cx, cy, cz] = params.center;
    let tr = params.tube_radius;

    let mut points = Points::new();
    let mut normals = DataArray::<f64>::new("Normals", 3);
    let mut polys = CellArray::new();

    // Compute Frenet frames along the curve and generate tube vertices
    for j in 0..n_curve {
        let t = 2.0 * PI * j as f64 / n_curve as f64;
        let eps = 1e-5;
        let p = trefoil_point(t);
        let p_next = trefoil_point(t + eps);
        let p_prev = trefoil_point(t - eps);

        // Tangent
        let tangent = normalize(sub(p_next, p_prev));

        // Approximate normal from second derivative
        let d2 = [
            p_next[0] - 2.0 * p[0] + p_prev[0],
            p_next[1] - 2.0 * p[1] + p_prev[1],
            p_next[2] - 2.0 * p[2] + p_prev[2],
        ];

        // Binormal = tangent x d2, then normal = binormal x tangent
        let binormal = normalize(cross(tangent, d2));
        let normal = cross(binormal, tangent);

        // Generate ring of tube vertices
        for i in 0..n_tube {
            let angle = 2.0 * PI * i as f64 / n_tube as f64;
            let ca = angle.cos();
            let sa = angle.sin();

            // Direction in the normal plane
            let dx = ca * normal[0] + sa * binormal[0];
            let dy = ca * normal[1] + sa * binormal[1];
            let dz = ca * normal[2] + sa * binormal[2];

            points.push([
                cx + p[0] + tr * dx,
                cy + p[1] + tr * dy,
                cz + p[2] + tr * dz,
            ]);
            normals.push_tuple(&[dx, dy, dz]);
        }
    }

    // Generate quads connecting adjacent rings
    for j in 0..n_curve {
        let j_next = (j + 1) % n_curve;
        for i in 0..n_tube {
            let i_next = (i + 1) % n_tube;
            let p00 = (j * n_tube + i) as i64;
            let p10 = (j * n_tube + i_next) as i64;
            let p01 = (j_next * n_tube + i) as i64;
            let p11 = (j_next * n_tube + i_next) as i64;
            polys.push_cell(&[p00, p10, p11, p01]);
        }
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
    fn default_trefoil_knot() {
        let pd = trefoil_knot(&TrefoilKnotParams::default());
        assert_eq!(pd.points.len(), 128 * 16);
        assert_eq!(pd.polys.num_cells(), 128 * 16);
        assert!(pd.point_data().normals().is_some());
    }

    #[test]
    fn minimal_trefoil_knot() {
        let pd = trefoil_knot(&TrefoilKnotParams {
            resolution: 6,
            tube_resolution: 3,
            ..Default::default()
        });
        assert_eq!(pd.points.len(), 18);
        assert_eq!(pd.polys.num_cells(), 18);
    }
}
