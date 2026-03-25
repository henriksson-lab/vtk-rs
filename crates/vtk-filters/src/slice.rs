use vtk_data::{CellArray, Points, PolyData};

/// Slice a PolyData mesh with a plane, producing intersection line segments.
///
/// Returns a PolyData containing line cells where the mesh intersects the plane.
/// The plane is defined by a point on the plane and the plane normal.
pub fn slice_by_plane(
    input: &PolyData,
    origin: [f64; 3],
    normal: [f64; 3],
) -> PolyData {
    let mut points = Points::<f64>::new();
    let mut lines = CellArray::new();

    for cell in input.polys.iter() {
        if cell.len() < 3 {
            continue;
        }

        // Compute signed distance from each vertex to the plane
        let dists: Vec<f64> = cell
            .iter()
            .map(|&id| {
                let p = input.points.get(id as usize);
                (p[0] - origin[0]) * normal[0]
                    + (p[1] - origin[1]) * normal[1]
                    + (p[2] - origin[2]) * normal[2]
            })
            .collect();

        // Find edge crossings
        let mut crossings = Vec::new();
        let n = cell.len();
        for i in 0..n {
            let j = (i + 1) % n;
            let di = dists[i];
            let dj = dists[j];

            if (di >= 0.0) != (dj >= 0.0) {
                // Edge crosses the plane
                let t = di / (di - dj);
                let pi = input.points.get(cell[i] as usize);
                let pj = input.points.get(cell[j] as usize);
                let intersection = [
                    pi[0] + t * (pj[0] - pi[0]),
                    pi[1] + t * (pj[1] - pi[1]),
                    pi[2] + t * (pj[2] - pi[2]),
                ];
                let idx = points.len() as i64;
                points.push(intersection);
                crossings.push(idx);
            } else if di.abs() < 1e-10 && dj.abs() >= 1e-10 {
                // Vertex i is exactly on the plane
                let idx = points.len() as i64;
                points.push(input.points.get(cell[i] as usize));
                crossings.push(idx);
            }
        }

        // For a planar cut through a polygon, we expect exactly 2 crossings
        if crossings.len() == 2 {
            lines.push_cell(&[crossings[0], crossings[1]]);
        }
    }

    let mut pd = PolyData::new();
    pd.points = points;
    pd.lines = lines;
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
