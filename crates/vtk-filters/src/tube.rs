use std::f64::consts::PI;

use vtk_data::{CellArray, DataArray, Points, PolyData};

/// Generate a tube around each line in a PolyData.
///
/// Each line segment is expanded into a tube with the given `radius` and
/// `sides` (number of facets around the circumference).
pub fn tube(input: &PolyData, radius: f64, sides: usize) -> PolyData {
    let sides = sides.max(3);

    let mut out_points = Points::<f64>::new();
    let mut out_normals = DataArray::<f64>::new("Normals", 3);
    let mut out_polys = CellArray::new();

    for cell in input.lines.iter() {
        if cell.len() < 2 {
            continue;
        }

        // Generate tube along this polyline
        let mut rings: Vec<Vec<usize>> = Vec::new();

        for seg_idx in 0..cell.len() {
            let p = input.points.get(cell[seg_idx] as usize);

            // Compute tangent direction
            let tangent = if seg_idx == 0 {
                let pn = input.points.get(cell[1] as usize);
                normalize_3([pn[0] - p[0], pn[1] - p[1], pn[2] - p[2]])
            } else if seg_idx == cell.len() - 1 {
                let pp = input.points.get(cell[seg_idx - 1] as usize);
                normalize_3([p[0] - pp[0], p[1] - pp[1], p[2] - pp[2]])
            } else {
                let pp = input.points.get(cell[seg_idx - 1] as usize);
                let pn = input.points.get(cell[seg_idx + 1] as usize);
                normalize_3([pn[0] - pp[0], pn[1] - pp[1], pn[2] - pp[2]])
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
                out_normals.push_tuple(&normal);
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

        // Cap start
        let base = out_points.len() as i64;
        let first_center = input.points.get(cell[0] as usize);
        out_points.push(first_center);
        let tangent0 = {
            let p0 = input.points.get(cell[0] as usize);
            let p1 = input.points.get(cell[1] as usize);
            normalize_3([p0[0] - p1[0], p0[1] - p1[1], p0[2] - p1[2]])
        };
        out_normals.push_tuple(&tangent0);
        for s in 0..sides {
            let sn = (s + 1) % sides;
            out_polys.push_cell(&[base, rings[0][sn] as i64, rings[0][s] as i64]);
        }

        // Cap end
        let base = out_points.len() as i64;
        let last_center = input.points.get(cell[cell.len() - 1] as usize);
        out_points.push(last_center);
        let tangent_end = {
            let pp = input.points.get(cell[cell.len() - 2] as usize);
            let pl = input.points.get(cell[cell.len() - 1] as usize);
            normalize_3([pl[0] - pp[0], pl[1] - pp[1], pl[2] - pp[2]])
        };
        out_normals.push_tuple(&tangent_end);
        let last_ring = &rings[rings.len() - 1];
        for s in 0..sides {
            let sn = (s + 1) % sides;
            out_polys.push_cell(&[base, last_ring[s] as i64, last_ring[sn] as i64]);
        }
    }

    let mut pd = PolyData::new();
    pd.points = out_points;
    pd.polys = out_polys;
    pd.point_data_mut().add_array(out_normals.into());
    pd.point_data_mut().set_active_normals("Normals");
    pd
}

fn perpendicular_frame(dir: [f64; 3]) -> ([f64; 3], [f64; 3]) {
    let seed = if dir[0].abs() < 0.9 {
        [1.0, 0.0, 0.0]
    } else {
        [0.0, 1.0, 0.0]
    };
    let u = normalize_3(cross_3(seed, dir));
    let v = cross_3(dir, u);
    (u, v)
}

fn cross_3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

fn normalize_3(v: [f64; 3]) -> [f64; 3] {
    let len = (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt();
    if len > 1e-10 {
        [v[0] / len, v[1] / len, v[2] / len]
    } else {
        [0.0, 0.0, 1.0]
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn tube_around_single_line() {
        let mut pd = PolyData::new();
        pd.points.push([0.0, 0.0, 0.0]);
        pd.points.push([0.0, 0.0, 1.0]);
        pd.lines.push_cell(&[0, 1]);

        let result = tube(&pd, 0.1, 8);
        // 2 rings * 8 sides + 2 cap centers = 18 points
        assert_eq!(result.points.len(), 18);
        // 8 side quads + 8 start cap tris + 8 end cap tris = 24
        assert_eq!(result.polys.num_cells(), 24);
    }

    #[test]
    fn tube_around_polyline() {
        let mut pd = PolyData::new();
        pd.points.push([0.0, 0.0, 0.0]);
        pd.points.push([1.0, 0.0, 0.0]);
        pd.points.push([1.0, 1.0, 0.0]);
        pd.lines.push_cell(&[0, 1, 2]);

        let result = tube(&pd, 0.05, 6);
        // 3 rings * 6 sides + 2 caps = 20 points
        assert_eq!(result.points.len(), 20);
        // 2*6 side quads + 6+6 cap tris = 24
        assert_eq!(result.polys.num_cells(), 24);
    }
}
