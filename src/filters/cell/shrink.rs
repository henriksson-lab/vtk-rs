use crate::data::{CellArray, Points, PolyData};

/// Shrink each cell toward its centroid.
///
/// Each polygon's vertices are moved toward the polygon's center by
/// `factor` (0.0 = all vertices collapse to center, 1.0 = no change).
/// Each cell gets its own copy of its vertices, so the output has
/// more points than the input (no shared vertices).
pub fn shrink(input: &PolyData, factor: f64) -> PolyData {
    let factor = factor.clamp(0.0, 1.0);

    let mut out_points = Points::<f64>::new();
    let mut out_polys = CellArray::new();

    for cell in input.polys.iter() {
        if cell.is_empty() {
            continue;
        }

        // Compute centroid
        let mut cx = 0.0;
        let mut cy = 0.0;
        let mut cz = 0.0;
        for &id in cell {
            let p = input.points.get(id as usize);
            cx += p[0];
            cy += p[1];
            cz += p[2];
        }
        let n = cell.len() as f64;
        cx /= n;
        cy /= n;
        cz /= n;

        // Create new vertices, moved toward centroid
        let base = out_points.len() as i64;
        let mut ids = Vec::with_capacity(cell.len());
        for &id in cell {
            let p = input.points.get(id as usize);
            out_points.push([
                cx + factor * (p[0] - cx),
                cy + factor * (p[1] - cy),
                cz + factor * (p[2] - cz),
            ]);
            ids.push(base + ids.len() as i64);
        }
        out_polys.push_cell(&ids);
    }

    let mut pd = PolyData::new();
    pd.points = out_points;
    pd.polys = out_polys;
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
