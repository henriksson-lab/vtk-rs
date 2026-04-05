use std::f64::consts::PI;

use crate::data::{CellArray, DataArray, Points, PolyData};

/// Parameters for generating a superquadric (superellipsoid).
pub struct SuperquadricParams {
    /// Phi roundness exponent. Default: 1.0 (sphere-like)
    pub phi_roundness: f64,
    /// Theta roundness exponent. Default: 1.0 (sphere-like)
    pub theta_roundness: f64,
    /// Scale factors [x, y, z]. Default: [1, 1, 1]
    pub scale: [f64; 3],
    /// Center. Default: [0, 0, 0]
    pub center: [f64; 3],
    /// Resolution in theta direction. Default: 16
    pub theta_resolution: usize,
    /// Resolution in phi direction. Default: 16
    pub phi_resolution: usize,
}

impl Default for SuperquadricParams {
    fn default() -> Self {
        Self {
            phi_roundness: 1.0,
            theta_roundness: 1.0,
            scale: [1.0, 1.0, 1.0],
            center: [0.0, 0.0, 0.0],
            theta_resolution: 16,
            phi_resolution: 16,
        }
    }
}

/// Generate a superquadric surface.
///
/// The parametric equations use signed power functions to produce a variety
/// of shapes: spheres, cubes (rounded), cylinders, octrahedra, etc.
pub fn superquadric(params: &SuperquadricParams) -> PolyData {
    let n_theta = params.theta_resolution.max(3);
    let n_phi = params.phi_resolution.max(3);

    let mut points = Points::new();
    let mut normals = DataArray::<f64>::new("Normals", 3);
    let mut polys = CellArray::new();

    let e1 = params.phi_roundness.max(0.01);
    let e2 = params.theta_roundness.max(0.01);

    // Generate vertices
    for j in 0..=n_phi {
        let phi = -PI / 2.0 + PI * j as f64 / n_phi as f64;
        let cp = phi.cos();
        let sp = phi.sin();

        for i in 0..=n_theta {
            let theta = -PI + 2.0 * PI * i as f64 / n_theta as f64;
            let ct = theta.cos();
            let st = theta.sin();

            let x = params.scale[0] * signed_pow(cp, e1) * signed_pow(ct, e2) + params.center[0];
            let y = params.scale[1] * signed_pow(cp, e1) * signed_pow(st, e2) + params.center[1];
            let z = params.scale[2] * signed_pow(sp, e1) + params.center[2];

            // Normal (approximation: derivative of parametric form)
            let nx = signed_pow(cp, 2.0 - e1) * signed_pow(ct, 2.0 - e2) / params.scale[0];
            let ny = signed_pow(cp, 2.0 - e1) * signed_pow(st, 2.0 - e2) / params.scale[1];
            let nz = signed_pow(sp, 2.0 - e1) / params.scale[2];
            let len = (nx * nx + ny * ny + nz * nz).sqrt();
            let n = if len > 1e-10 {
                [nx / len, ny / len, nz / len]
            } else {
                [0.0, 0.0, 1.0]
            };

            points.push([x, y, z]);
            normals.push_tuple(&n);
        }
    }

    // Generate quads
    let row = n_theta + 1;
    for j in 0..n_phi {
        for i in 0..n_theta {
            let p0 = (j * row + i) as i64;
            let p1 = p0 + 1;
            let p2 = p0 + row as i64 + 1;
            let p3 = p0 + row as i64;
            polys.push_cell(&[p0, p1, p2, p3]);
        }
    }

    let mut pd = PolyData::new();
    pd.points = points;
    pd.polys = polys;
    pd.point_data_mut().add_array(normals.into());
    pd.point_data_mut().set_active_normals("Normals");
    pd
}

fn signed_pow(base: f64, exp: f64) -> f64 {
    base.signum() * base.abs().powf(exp)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn default_superquadric() {
        let pd = superquadric(&SuperquadricParams::default());
        // (16+1) * (16+1) = 289 points, 16*16 = 256 quads
        assert_eq!(pd.points.len(), 289);
        assert_eq!(pd.polys.num_cells(), 256);
    }

    #[test]
    fn boxy_superquadric() {
        let pd = superquadric(&SuperquadricParams {
            phi_roundness: 0.1,
            theta_roundness: 0.1,
            ..Default::default()
        });
        assert!(pd.points.len() > 0);
        assert!(pd.polys.num_cells() > 0);
    }
}
