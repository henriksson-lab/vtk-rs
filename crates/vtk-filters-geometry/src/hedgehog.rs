//! HedgeHog filter: oriented line glyphs from vector fields.
//!
//! For each point with a vector data array, creates a line segment
//! from the point in the direction of the vector, scaled by a factor.

use vtk_data::{AnyDataArray, CellArray, Points, PolyData};

/// Generate hedgehog (oriented line) glyphs from a vector field.
///
/// Each point produces a line from the point position to
/// `position + vector * scale_factor`.
pub fn hedgehog(input: &PolyData, vector_name: &str, scale_factor: f64) -> PolyData {
    let vectors = match input.point_data().get_array(vector_name) {
        Some(arr) => arr,
        None => return PolyData::new(),
    };

    let n = input.points.len();
    let nc = vectors.num_components();
    if nc < 3 {
        return PolyData::new();
    }

    let mut points = Points::<f64>::new();
    let mut lines = CellArray::new();

    for i in 0..n {
        let p = input.points.get(i);
        let mut v = [0.0f64; 3];
        vectors.tuple_as_f64(i, &mut v);

        let tip = [
            p[0] + v[0] * scale_factor,
            p[1] + v[1] * scale_factor,
            p[2] + v[2] * scale_factor,
        ];

        let base_idx = points.len();
        points.push(p);
        points.push(tip);
        lines.push_cell(&[base_idx as i64, (base_idx + 1) as i64]);
    }

    let mut result = PolyData::new();
    result.points = points;
    result.lines = lines;
    result
}

#[cfg(test)]
mod tests {
    use super::*;
    use vtk_data::DataArray;

    #[test]
    fn hedgehog_basic() {
        let mut pd = PolyData::from_triangles(
            vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]],
            vec![[0, 1, 2]],
        );
        pd.point_data_mut().add_array(AnyDataArray::F64(
            DataArray::from_vec("vectors", vec![
                1.0, 0.0, 0.0,
                0.0, 1.0, 0.0,
                0.0, 0.0, 1.0,
            ], 3),
        ));
        let result = hedgehog(&pd, "vectors", 2.0);
        assert_eq!(result.points.len(), 6); // 3 points × 2 (base + tip)
        assert_eq!(result.lines.num_cells(), 3);
        // Check first tip
        let tip = result.points.get(1);
        assert!((tip[0] - 2.0).abs() < 1e-10);
    }

    #[test]
    fn hedgehog_missing_array() {
        let pd = PolyData::from_triangles(
            vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]],
            vec![[0, 1, 2]],
        );
        let result = hedgehog(&pd, "nonexistent", 1.0);
        assert_eq!(result.points.len(), 0);
    }
}
