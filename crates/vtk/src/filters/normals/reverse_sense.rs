use crate::data::{CellArray, PolyData};

/// Reverse the winding order of all polygon cells and flip normals.
///
/// This effectively flips the "inside" and "outside" of a surface.
pub fn reverse_sense(input: &PolyData) -> PolyData {
    let mut out_polys = CellArray::new();

    for cell in input.polys.iter() {
        let reversed: Vec<i64> = cell.iter().rev().copied().collect();
        out_polys.push_cell(&reversed);
    }

    let mut pd = input.clone();
    pd.polys = out_polys;

    // Flip normals if present
    let flipped_normals = pd.point_data().normals().map(|normals| {
        let nc = normals.num_components();
        let nt = normals.num_tuples();
        let name = normals.name().to_string();
        let mut flipped = Vec::with_capacity(nt * nc);
        let mut buf = vec![0.0f64; nc];
        for i in 0..nt {
            normals.tuple_as_f64(i, &mut buf);
            for v in &buf {
                flipped.push(-v);
            }
        }
        (name, flipped, nc)
    });

    if let Some((name, flipped, nc)) = flipped_normals {
        let arr = crate::data::AnyDataArray::F64(
            crate::data::DataArray::from_vec(&name, flipped, nc),
        );
        pd.point_data_mut().add_array(arr);
        pd.point_data_mut().set_active_normals(&name);
    }

    pd
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn reverse_triangle() {
        let pd = PolyData::from_triangles(
            vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.5, 1.0, 0.0]],
            vec![[0, 1, 2]],
        );
        let result = reverse_sense(&pd);
        assert_eq!(result.polys.cell(0), &[2, 1, 0]);
    }

    #[test]
    fn reverse_multiple() {
        let pd = PolyData::from_triangles(
            vec![
                [0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.5, 1.0, 0.0],
                [2.0, 0.0, 0.0],
            ],
            vec![[0, 1, 2], [1, 3, 2]],
        );
        let result = reverse_sense(&pd);
        assert_eq!(result.polys.cell(0), &[2, 1, 0]);
        assert_eq!(result.polys.cell(1), &[2, 3, 1]);
    }

    #[test]
    fn preserves_points() {
        let pd = PolyData::from_triangles(
            vec![[1.0, 2.0, 3.0], [4.0, 5.0, 6.0], [7.0, 8.0, 9.0]],
            vec![[0, 1, 2]],
        );
        let result = reverse_sense(&pd);
        assert_eq!(result.points.len(), 3);
        assert_eq!(result.points.get(0), [1.0, 2.0, 3.0]);
    }
}
