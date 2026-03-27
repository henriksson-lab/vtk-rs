use std::f64::consts::PI;

use vtk_data::{CellArray, DataArray, Points, PolyData};

/// Parameters for generating a Möbius strip surface with normals.
pub struct MobiusStripParams {
    /// Radius of the center circle. Default: 1.0
    pub radius: f64,
    /// Half-width of the strip. Default: 0.3
    pub width: f64,
    /// Number of segments around the loop. Default: 64
    pub resolution: usize,
    /// Number of subdivisions across the strip width. Default: 4
    pub width_resolution: usize,
    /// Center of the strip. Default: [0, 0, 0]
    pub center: [f64; 3],
}

impl Default for MobiusStripParams {
    fn default() -> Self {
        Self {
            radius: 1.0,
            width: 0.3,
            resolution: 64,
            width_resolution: 4,
            center: [0.0, 0.0, 0.0],
        }
    }
}

/// Generate a Möbius strip as a triangulated PolyData surface with normals.
///
/// This generates the classic parametric Möbius band with a half-twist,
/// triangulated for clean rendering. Unlike the simpler `mobius` source,
/// this version supports configurable width resolution and generates normals.
pub fn mobius_strip(params: &MobiusStripParams) -> PolyData {
    let n_u = params.resolution.max(8);
    let n_v = params.width_resolution.max(1);
    let r = params.radius;
    let w = params.width;
    let [cx, cy, cz] = params.center;

    let mut points = Points::new();
    let mut normals = DataArray::<f64>::new("Normals", 3);
    let mut polys = CellArray::new();

    // Generate points on the parametric surface.
    // u goes around the loop [0, 2*PI], v goes across the width [-1, 1].
    // x = (R + w*v*cos(u/2)) * cos(u)
    // y = (R + w*v*cos(u/2)) * sin(u)
    // z = w * v * sin(u/2)
    for i in 0..=n_u {
        let u = 2.0 * PI * i as f64 / n_u as f64;
        let cu = u.cos();
        let su = u.sin();
        let ch = (u / 2.0).cos();
        let sh = (u / 2.0).sin();

        for j in 0..=n_v {
            let v = -1.0 + 2.0 * j as f64 / n_v as f64;

            let x = cx + (r + w * v * ch) * cu;
            let y = cy + (r + w * v * ch) * su;
            let z = cz + w * v * sh;
            points.push([x, y, z]);

            // Approximate normal via cross product of partial derivatives.
            // du: partial derivative with respect to u
            let du_x = -(r + w * v * ch) * su + w * v * (-0.5 * sh) * cu;
            let du_y = (r + w * v * ch) * cu + w * v * (-0.5 * sh) * su;
            let du_z = w * v * 0.5 * ch;

            // dv: partial derivative with respect to v
            let dv_x = w * ch * cu;
            let dv_y = w * ch * su;
            let dv_z = w * sh;

            // normal = du x dv
            let nx = du_y * dv_z - du_z * dv_y;
            let ny = du_z * dv_x - du_x * dv_z;
            let nz = du_x * dv_y - du_y * dv_x;
            let len = (nx * nx + ny * ny + nz * nz).sqrt().max(1e-12);
            normals.push_tuple(&[nx / len, ny / len, nz / len]);
        }
    }

    // Triangulate: two triangles per quad.
    let row = n_v + 1;
    for i in 0..n_u {
        for j in 0..n_v {
            let a = (i * row + j) as i64;
            let b = (i * row + j + 1) as i64;
            let c = ((i + 1) * row + j + 1) as i64;
            let d = ((i + 1) * row + j) as i64;
            polys.push_cell(&[a, b, c]);
            polys.push_cell(&[a, c, d]);
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
    fn default_mobius_strip() {
        let pd = mobius_strip(&MobiusStripParams::default());
        // (64+1) * (4+1) = 325 points
        assert_eq!(pd.points.len(), 325);
        // 64 * 4 * 2 = 512 triangles
        assert_eq!(pd.polys.num_cells(), 512);
        assert!(pd.point_data().normals().is_some());
    }

    #[test]
    fn minimal_mobius_strip() {
        let pd = mobius_strip(&MobiusStripParams {
            resolution: 8,
            width_resolution: 1,
            ..Default::default()
        });
        // (8+1) * (1+1) = 18 points
        assert_eq!(pd.points.len(), 18);
        // 8 * 1 * 2 = 16 triangles
        assert_eq!(pd.polys.num_cells(), 16);
    }

    #[test]
    fn custom_center() {
        let pd = mobius_strip(&MobiusStripParams {
            center: [1.0, 2.0, 3.0],
            resolution: 8,
            width_resolution: 1,
            ..Default::default()
        });
        assert!(pd.points.len() > 0);
    }
}
