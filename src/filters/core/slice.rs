use crate::data::{CellArray, Points, PolyData};

/// Slice a PolyData mesh with a plane, producing intersection line segments.
///
/// Returns a PolyData containing line cells where the mesh intersects the plane.
/// The plane is defined by a point on the plane and the plane normal.
pub fn slice_by_plane(
    input: &PolyData,
    origin: [f64; 3],
    normal: [f64; 3],
) -> PolyData {
    let nc = input.polys.num_cells();
    if nc == 0 {
        return PolyData::new();
    }

    // Pre-compute signed distances for all points (avoids redundant per-cell recomputation)
    let np = input.points.len();
    let mut dists = Vec::with_capacity(np);
    for i in 0..np {
        let p = input.points.get(i);
        dists.push(
            (p[0] - origin[0]) * normal[0]
                + (p[1] - origin[1]) * normal[1]
                + (p[2] - origin[2]) * normal[2],
        );
    }

    // Pre-sized flat buffers for output
    let mut pts_flat: Vec<f64> = Vec::with_capacity(nc * 6); // ~2 points per intersected cell
    let mut line_conn: Vec<i64> = Vec::with_capacity(nc * 2);
    let mut line_off: Vec<i64> = Vec::with_capacity(nc + 1);
    line_off.push(0);

    let offsets = input.polys.offsets();
    let conn = input.polys.connectivity();

    for ci in 0..nc {
        let start = offsets[ci] as usize;
        let end = offsets[ci + 1] as usize;
        let cell = &conn[start..end];
        let n = cell.len();
        if n < 3 {
            continue;
        }

        // Find edge crossings for this cell
        let mut c0 = [0.0f64; 3];
        let mut c1 = [0.0f64; 3];
        let mut num_crossings = 0u32;

        for i in 0..n {
            let j = (i + 1) % n;
            let ai = cell[i] as usize;
            let aj = cell[j] as usize;
            let di = unsafe { *dists.get_unchecked(ai) };
            let dj = unsafe { *dists.get_unchecked(aj) };

            if (di >= 0.0) != (dj >= 0.0) {
                let t = di / (di - dj);
                let pi = input.points.get(ai);
                let pj = input.points.get(aj);
                let pt = [
                    pi[0] + t * (pj[0] - pi[0]),
                    pi[1] + t * (pj[1] - pi[1]),
                    pi[2] + t * (pj[2] - pi[2]),
                ];
                if num_crossings == 0 {
                    c0 = pt;
                    num_crossings = 1;
                } else if num_crossings == 1 {
                    c1 = pt;
                    num_crossings = 2;
                }
                // For triangles we always get exactly 2 crossings
            } else if di.abs() < 1e-10 && dj.abs() >= 1e-10 {
                let pt = input.points.get(ai);
                if num_crossings == 0 {
                    c0 = pt;
                    num_crossings = 1;
                } else if num_crossings == 1 {
                    c1 = pt;
                    num_crossings = 2;
                }
            }
        }

        if num_crossings == 2 {
            let idx = (pts_flat.len() / 3) as i64;
            pts_flat.push(c0[0]); pts_flat.push(c0[1]); pts_flat.push(c0[2]);
            pts_flat.push(c1[0]); pts_flat.push(c1[1]); pts_flat.push(c1[2]);
            line_conn.push(idx);
            line_conn.push(idx + 1);
            line_off.push(line_conn.len() as i64);
        }
    }

    let mut pd = PolyData::new();
    pd.points = Points::from_flat_vec(pts_flat);
    pd.lines = CellArray::from_raw(line_off, line_conn);
    pd
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn slice_triangle_through_middle() {
        let pd = PolyData::from_triangles(
            vec![[-1.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]],
            vec![[0, 1, 2]],
        );

        // Slice with plane at x=0, normal +X
        let result = slice_by_plane(&pd, [0.0, 0.0, 0.0], [1.0, 0.0, 0.0]);
        assert_eq!(result.lines.num_cells(), 1);
        assert_eq!(result.points.len(), 2);
    }

    #[test]
    fn slice_misses_triangle() {
        let pd = PolyData::from_triangles(
            vec![[1.0, 0.0, 0.0], [2.0, 0.0, 0.0], [1.5, 0.0, 1.0]],
            vec![[0, 1, 2]],
        );

        // Slice with plane at x=0 — triangle is entirely on positive side
        let result = slice_by_plane(&pd, [0.0, 0.0, 0.0], [1.0, 0.0, 0.0]);
        assert_eq!(result.lines.num_cells(), 0);
    }

    #[test]
    fn slice_multiple_triangles() {
        let pd = PolyData::from_triangles(
            vec![
                [-1.0, -1.0, 0.0], [1.0, -1.0, 0.0], [0.0, -1.0, 1.0],
                [-1.0, 1.0, 0.0], [1.0, 1.0, 0.0], [0.0, 1.0, 1.0],
            ],
            vec![[0, 1, 2], [3, 4, 5]],
        );

        let result = slice_by_plane(&pd, [0.0, 0.0, 0.0], [1.0, 0.0, 0.0]);
        assert_eq!(result.lines.num_cells(), 2);
    }
}
