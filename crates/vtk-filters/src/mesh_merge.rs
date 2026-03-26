use vtk_data::{CellArray, PolyData};

/// Merge two PolyData meshes into one by concatenating points and cells.
///
/// Cell point indices from the second mesh are adjusted by the number of points
/// in the first mesh. All cell types (verts, lines, polys, strips) are merged.
pub fn merge_poly_data(a: &PolyData, b: &PolyData) -> PolyData {
    let mut output = PolyData::new();

    // Copy all points from A
    for i in 0..a.points.len() {
        output.points.push(a.points.get(i));
    }
    // Copy all points from B
    for i in 0..b.points.len() {
        output.points.push(b.points.get(i));
    }

    let offset: i64 = a.points.len() as i64;

    // Copy cells from A directly
    copy_cells(&a.verts, &mut output.verts, 0);
    copy_cells(&a.lines, &mut output.lines, 0);
    copy_cells(&a.polys, &mut output.polys, 0);
    copy_cells(&a.strips, &mut output.strips, 0);

    // Copy cells from B with offset
    copy_cells(&b.verts, &mut output.verts, offset);
    copy_cells(&b.lines, &mut output.lines, offset);
    copy_cells(&b.polys, &mut output.polys, offset);
    copy_cells(&b.strips, &mut output.strips, offset);

    output
}

fn copy_cells(src: &CellArray, dst: &mut CellArray, offset: i64) {
    for cell in src.iter() {
        let remapped: Vec<i64> = cell.iter().map(|&id| id + offset).collect();
        dst.push_cell(&remapped);
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn merge_two_triangles() {
        let a = PolyData::from_triangles(
            vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]],
            vec![[0, 1, 2]],
        );
        let b = PolyData::from_triangles(
            vec![[2.0, 0.0, 0.0], [3.0, 0.0, 0.0], [2.0, 1.0, 0.0]],
            vec![[0, 1, 2]],
        );
        let result = merge_poly_data(&a, &b);
        assert_eq!(result.points.len(), 6);
        assert_eq!(result.polys.num_cells(), 2);
        // Second triangle's indices should be offset by 3
        let cell1 = result.polys.cell(1);
        assert_eq!(cell1, &[3, 4, 5]);
    }

    #[test]
    fn merge_with_empty() {
        let a = PolyData::from_triangles(
            vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]],
            vec![[0, 1, 2]],
        );
        let b = PolyData::new();
        let result = merge_poly_data(&a, &b);
        assert_eq!(result.points.len(), 3);
        assert_eq!(result.polys.num_cells(), 1);
    }

    #[test]
    fn merge_preserves_all_points() {
        let mut a = PolyData::new();
        a.points.push([0.0, 0.0, 0.0]);
        a.points.push([1.0, 0.0, 0.0]);

        let mut b = PolyData::new();
        b.points.push([5.0, 5.0, 5.0]);

        let result = merge_poly_data(&a, &b);
        assert_eq!(result.points.len(), 3);
        let p = result.points.get(2);
        assert!((p[0] - 5.0).abs() < 1e-10);
        assert!((p[1] - 5.0).abs() < 1e-10);
        assert!((p[2] - 5.0).abs() < 1e-10);
    }
}
