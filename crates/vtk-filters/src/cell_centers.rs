use vtk_data::{CellArray, Points, PolyData};

/// Generate points at the centroid of each polygon cell.
///
/// The output is a PolyData with one vertex cell per input polygon.
/// Useful for placing labels or glyphs at cell centers.
pub fn cell_centers(input: &PolyData) -> PolyData {
    let mut out_points = Points::<f64>::new();
    let mut out_verts = CellArray::new();

    for cell in input.polys.iter() {
        if cell.is_empty() {
            continue;
        }
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
        let idx = out_points.len() as i64;
        out_points.push([cx / n, cy / n, cz / n]);
        out_verts.push_cell(&[idx]);
    }

    // Also handle line cells
    for cell in input.lines.iter() {
        if cell.is_empty() {
            continue;
        }
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
        let idx = out_points.len() as i64;
        out_points.push([cx / n, cy / n, cz / n]);
        out_verts.push_cell(&[idx]);
    }

    let mut pd = PolyData::new();
    pd.points = out_points;
    pd.verts = out_verts;
    pd
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn centers_of_triangles() {
        let pd = PolyData::from_triangles(
            vec![
                [0.0, 0.0, 0.0],
                [3.0, 0.0, 0.0],
                [0.0, 3.0, 0.0],
                [3.0, 3.0, 0.0],
            ],
            vec![[0, 1, 2], [1, 3, 2]],
        );
        let result = cell_centers(&pd);
        assert_eq!(result.points.len(), 2);
        assert_eq!(result.verts.num_cells(), 2);

        let c0 = result.points.get(0);
        assert!((c0[0] - 1.0).abs() < 1e-10);
        assert!((c0[1] - 1.0).abs() < 1e-10);
    }

    #[test]
    fn centers_of_lines() {
        let mut pd = PolyData::new();
        pd.points.push([0.0, 0.0, 0.0]);
        pd.points.push([2.0, 0.0, 0.0]);
        pd.lines.push_cell(&[0, 1]);

        let result = cell_centers(&pd);
        assert_eq!(result.points.len(), 1);
        let c = result.points.get(0);
        assert!((c[0] - 1.0).abs() < 1e-10);
    }
}
