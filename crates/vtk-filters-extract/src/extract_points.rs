use vtk_data::{CellArray, Points, PolyData};

/// Extract points by index from a PolyData, producing vertex cells.
pub fn extract_points(input: &PolyData, point_indices: &[usize]) -> PolyData {
    let mut out_points = Points::<f64>::new();
    let mut out_verts = CellArray::new();

    for (new_idx, &pi) in point_indices.iter().enumerate() {
        if pi < input.points.len() {
            out_points.push(input.points.get(pi));
            out_verts.push_cell(&[new_idx as i64]);
        }
    }

    let mut pd = PolyData::new();
    pd.points = out_points;
    pd.verts = out_verts;
    pd
}

/// Extract points where a scalar condition is met.
pub fn extract_points_by_scalar(
    input: &PolyData,
    scalar_name: &str,
    min_val: f64,
    max_val: f64,
) -> PolyData {
    let arr = match input.point_data().get_array(scalar_name) {
        Some(a) => a,
        None => return PolyData::new(),
    };

    let mut indices = Vec::new();
    let mut buf = [0.0f64];
    for i in 0..arr.num_tuples().min(input.points.len()) {
        arr.tuple_as_f64(i, &mut buf);
        if buf[0] >= min_val && buf[0] <= max_val {
            indices.push(i);
        }
    }

    extract_points(input, &indices)
}

#[cfg(test)]
mod tests {
    use super::*;
    use vtk_data::DataArray;

    #[test]
    fn extract_by_index() {
        let mut pd = PolyData::new();
        pd.points.push([0.0, 0.0, 0.0]);
        pd.points.push([1.0, 0.0, 0.0]);
        pd.points.push([2.0, 0.0, 0.0]);
        pd.points.push([3.0, 0.0, 0.0]);

        let result = extract_points(&pd, &[1, 3]);
        assert_eq!(result.points.len(), 2);
        assert_eq!(result.verts.num_cells(), 2);
        let p = result.points.get(0);
        assert!((p[0] - 1.0).abs() < 1e-10);
    }

    #[test]
    fn extract_by_scalar() {
        let mut pd = PolyData::new();
        for i in 0..5 {
            pd.points.push([i as f64, 0.0, 0.0]);
        }
        let scalars = DataArray::from_vec("temp", vec![10.0, 20.0, 30.0, 40.0, 50.0], 1);
        pd.point_data_mut().add_array(scalars.into());

        let result = extract_points_by_scalar(&pd, "temp", 25.0, 45.0);
        assert_eq!(result.points.len(), 2); // points with temp 30, 40
    }
}
