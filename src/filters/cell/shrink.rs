use crate::data::{CellArray, Points, PolyData};

/// Shrink each cell toward its centroid.
///
/// Each polygon's vertices are moved toward the polygon's center by
/// `factor` (0.0 = all vertices collapse to center, 1.0 = no change).
/// Each cell gets its own copy of its vertices, so the output has
/// more points than the input (no shared vertices).
pub fn shrink(input: &PolyData, factor: f64) -> PolyData {
    let factor = factor.clamp(0.0, 1.0);

    let offsets = input.polys.offsets();
    let conn = input.polys.connectivity();
    let nc = input.polys.num_cells();
    let pts = input.points.as_flat_slice();

    // Pre-allocated flat output: each cell vertex becomes a unique point.
    // Uses raw offsets/connectivity for zero-overhead cell traversal.
    let total_out_pts = conn.len();
    let mut out_flat = Vec::with_capacity(total_out_pts * 3);
    let mut out_off = Vec::with_capacity(nc + 1);
    let mut out_conn = Vec::with_capacity(total_out_pts);
    out_off.push(0i64);

    let mut pt_idx: i64 = 0;

    for ci in 0..nc {
        let start = offsets[ci] as usize;
        let end = offsets[ci + 1] as usize;
        let n = (end - start) as f64;
        if n < 1.0 { continue; }

        // Compute centroid
        let mut cx = 0.0f64;
        let mut cy = 0.0f64;
        let mut cz = 0.0f64;
        for idx in start..end {
            let b = conn[idx] as usize * 3;
            cx += pts[b];
            cy += pts[b + 1];
            cz += pts[b + 2];
        }
        cx /= n;
        cy /= n;
        cz /= n;

        // Create new vertices
        for idx in start..end {
            let b = conn[idx] as usize * 3;
            out_flat.push(cx + factor * (pts[b] - cx));
            out_flat.push(cy + factor * (pts[b + 1] - cy));
            out_flat.push(cz + factor * (pts[b + 2] - cz));
            out_conn.push(pt_idx);
            pt_idx += 1;
        }
        out_off.push(pt_idx);
    }

    let mut pd = PolyData::new();
    pd.points = Points::from_flat_vec(out_flat);
    pd.polys = CellArray::from_raw(out_off, out_conn);
    pd
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn shrink_factor_one_preserves_geometry() {
        let pd = PolyData::from_triangles(
            vec![[0.0, 0.0, 0.0], [3.0, 0.0, 0.0], [0.0, 3.0, 0.0]],
            vec![[0, 1, 2]],
        );
        let result = shrink(&pd, 1.0);
        assert_eq!(result.points.len(), 3);
        assert_eq!(result.polys.num_cells(), 1);
        let p0 = result.points.get(0);
        assert!((p0[0]).abs() < 1e-10);
    }

    #[test]
    fn shrink_factor_zero_collapses_to_centroid() {
        let pd = PolyData::from_triangles(
            vec![[0.0, 0.0, 0.0], [3.0, 0.0, 0.0], [0.0, 3.0, 0.0]],
            vec![[0, 1, 2]],
        );
        let result = shrink(&pd, 0.0);
        assert_eq!(result.points.len(), 3);
        // All points should be at centroid (1.0, 1.0, 0.0)
        for i in 0..3 {
            let p = result.points.get(i);
            assert!((p[0] - 1.0).abs() < 1e-10);
            assert!((p[1] - 1.0).abs() < 1e-10);
            assert!((p[2]).abs() < 1e-10);
        }
    }

    #[test]
    fn shrink_duplicates_shared_points() {
        let pd = PolyData::from_triangles(
            vec![
                [0.0, 0.0, 0.0],
                [1.0, 0.0, 0.0],
                [0.5, 1.0, 0.0],
                [1.5, 1.0, 0.0],
            ],
            vec![[0, 1, 2], [1, 3, 2]],
        );
        let result = shrink(&pd, 0.5);
        // Each triangle gets its own 3 points
        assert_eq!(result.points.len(), 6);
        assert_eq!(result.polys.num_cells(), 2);
    }
}
