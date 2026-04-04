//! Generate tubes around polylines where radius varies by an uncertainty scalar.
//!
//! Each line cell in the input PolyData is expanded into a tube mesh where the
//! radius at each vertex is `base_radius * uncertainty_value`.

use std::f64::consts::PI;

use vtk_data::{CellArray, Points, PolyData};

/// Generate uncertainty tubes around polylines.
///
/// The input must have line cells and a point data array named `"Uncertainty"`
/// with scalar values that modulate the tube radius.
///
/// # Arguments
/// * `input` - PolyData with line cells
/// * `base_radius` - Base tube radius (multiplied by uncertainty at each vertex)
/// * `sides` - Number of facets around the tube circumference (minimum 3)
pub fn uncertainty_tube(input: &PolyData, base_radius: f64, sides: usize) -> PolyData {
    let sides = sides.max(3);

    let mut out_points = Points::<f64>::new();
    let mut out_polys = CellArray::new();

    // Get uncertainty array
    let uncertainty_arr = input.point_data().get_array("Uncertainty");

    for cell in input.lines.iter() {
        if cell.len() < 2 {
            continue;
        }

        let mut rings: Vec<Vec<usize>> = Vec::new();

        for seg_idx in 0..cell.len() {
            let vid = cell[seg_idx] as usize;
            let p = input.points.get(vid);

            // Get uncertainty value for this vertex
            let unc = if let Some(arr) = uncertainty_arr {
                let mut buf = [0.0f64];
                arr.tuple_as_f64(vid, &mut buf);
                buf[0]
            } else {
                1.0
            };
            let radius = base_radius * unc;

            // Compute tangent direction
            let tangent = if seg_idx == 0 {
                let pn = input.points.get(cell[1] as usize);
                normalize([pn[0] - p[0], pn[1] - p[1], pn[2] - p[2]])
            } else if seg_idx == cell.len() - 1 {
                let pp = input.points.get(cell[seg_idx - 1] as usize);
                normalize([p[0] - pp[0], p[1] - pp[1], p[2] - pp[2]])
            } else {
                let pp = input.points.get(cell[seg_idx - 1] as usize);
                let pn = input.points.get(cell[seg_idx + 1] as usize);
                normalize([pn[0] - pp[0], pn[1] - pp[1], pn[2] - pp[2]])
            };

            let (u, v) = perpendicular_frame(tangent);

            let mut ring = Vec::with_capacity(sides);
            for s in 0..sides {
                let angle = 2.0 * PI * s as f64 / sides as f64;
                let ct = angle.cos();
                let st = angle.sin();

                let normal = [
                    ct * u[0] + st * v[0],
                    ct * u[1] + st * v[1],
                    ct * u[2] + st * v[2],
                ];

                let idx = out_points.len();
                out_points.push([
                    p[0] + radius * normal[0],
                    p[1] + radius * normal[1],
                    p[2] + radius * normal[2],
                ]);
                ring.push(idx);
            }
            rings.push(ring);
        }

        // Connect adjacent rings with quads
        for r in 0..rings.len() - 1 {
            for s in 0..sides {
                let sn = (s + 1) % sides;
                out_polys.push_cell(&[
                    rings[r][s] as i64,
                    rings[r][sn] as i64,
                    rings[r + 1][sn] as i64,
                    rings[r + 1][s] as i64,
                ]);
            }
        }
    }

    let mut output = PolyData::new();
    output.points = out_points;
    output.polys = out_polys;
    output
}

fn normalize(v: [f64; 3]) -> [f64; 3] {
    let len = (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt();
    if len < 1e-12 {
        [1.0, 0.0, 0.0]
    } else {
        [v[0] / len, v[1] / len, v[2] / len]
    }
}

fn perpendicular_frame(tangent: [f64; 3]) -> ([f64; 3], [f64; 3]) {
    // Choose a vector not parallel to tangent
    let seed = if tangent[0].abs() < 0.9 {
        [1.0, 0.0, 0.0]
    } else {
        [0.0, 1.0, 0.0]
    };

    // u = normalize(seed cross tangent)
    let u = normalize(cross(seed, tangent));
    // v = tangent cross u
    let v = cross(tangent, u);
    (u, v)
}

fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

#[cfg(test)]
mod tests {
    use super::*;
    use vtk_data::{AnyDataArray, DataArray};

    #[test]
    fn uncertainty_tube_basic() {
        // Create a simple line along the x-axis with varying uncertainty
        let mut pd = PolyData::new();
        pd.points.push([0.0, 0.0, 0.0]);
        pd.points.push([1.0, 0.0, 0.0]);
        pd.points.push([2.0, 0.0, 0.0]);
        pd.lines.push_cell(&[0, 1, 2]);

        let unc = DataArray::from_vec("Uncertainty", vec![1.0, 2.0, 0.5], 1);
        pd.point_data_mut().add_array(AnyDataArray::F64(unc));

        let result = uncertainty_tube(&pd, 0.1, 6);

        // 3 rings of 6 points = 18 points
        assert_eq!(result.points.len(), 18);
        // 2 segments * 6 quads = 12 quads
        assert_eq!(result.polys.num_cells(), 12);

        // Check that the second ring has larger radius than first
        // Ring 0 center at (0,0,0), radius 0.1*1.0 = 0.1
        // Ring 1 center at (1,0,0), radius 0.1*2.0 = 0.2
        let p0 = result.points.get(0); // first ring point
        let r0 = ((p0[1] * p0[1]) + (p0[2] * p0[2])).sqrt();
        let p6 = result.points.get(6); // second ring point
        let r6 = (((p6[0] - 1.0) * (p6[0] - 1.0)) + (p6[1] * p6[1]) + (p6[2] * p6[2])).sqrt();
        // r6 should be roughly double r0
        assert!((r6 / r0 - 2.0).abs() < 0.1);
    }
}
