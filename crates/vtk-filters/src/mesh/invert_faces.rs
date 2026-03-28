use vtk_data::{CellArray, PolyData};

/// Reverse the winding order of all polygon faces, effectively flipping
/// the face normals.
///
/// Unlike `reverse_sense`, this does not flip point normal arrays --
/// it only reverses the vertex order in each polygon cell.
pub fn invert_faces(input: &PolyData) -> PolyData {
    let mut out_polys = CellArray::new();

    for cell in input.polys.iter() {
        let reversed: Vec<i64> = cell.iter().rev().copied().collect();
        out_polys.push_cell(&reversed);
    }

    let mut pd = input.clone();
    pd.polys = out_polys;
    pd
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn invert_single_triangle() {
        let pd = PolyData::from_triangles(
            vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]],
            vec![[0, 1, 2]],
        );
        let result = invert_faces(&pd);
        assert_eq!(result.polys.cell(0), &[2, 1, 0]);
    }

    #[test]
    fn invert_preserves_point_count() {
        let pd = PolyData::from_triangles(
            vec![
                [0.0, 0.0, 0.0],
                [1.0, 0.0, 0.0],
                [0.5, 1.0, 0.0],
                [2.0, 0.0, 0.0],
            ],
            vec![[0, 1, 2], [1, 3, 2]],
        );
        let result = invert_faces(&pd);
        assert_eq!(result.points.len(), 4);
        assert_eq!(result.polys.num_cells(), 2);
    }

    #[test]
    fn double_invert_is_identity() {
        let pd = PolyData::from_triangles(
            vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]],
            vec![[0, 1, 2]],
        );
        let result = invert_faces(&invert_faces(&pd));
        assert_eq!(result.polys.cell(0), &[0, 1, 2]);
    }
}
