use vtk_data::{Points, PolyData};

/// Scale a mesh uniformly around a center point.
///
/// Each point is moved relative to `center` by the given `factor`.
pub fn scale_uniform(input: &PolyData, factor: f64, center: [f64; 3]) -> PolyData {
    scale_nonuniform(input, [factor, factor, factor], center)
}

/// Scale a mesh non-uniformly around a center point.
///
/// Each point is moved relative to `center` by the given per-axis `factors`.
pub fn scale_nonuniform(input: &PolyData, factors: [f64; 3], center: [f64; 3]) -> PolyData {
    let mut new_points = Points::new();
    for i in 0..input.points.len() {
        let p = input.points.get(i);
        let x: f64 = center[0] + (p[0] - center[0]) * factors[0];
        let y: f64 = center[1] + (p[1] - center[1]) * factors[1];
        let z: f64 = center[2] + (p[2] - center[2]) * factors[2];
        new_points.push([x, y, z]);
    }

    let mut pd = input.clone();
    pd.points = new_points;
    pd
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn uniform_scale_from_origin() {
        let pd = PolyData::from_triangles(
            vec![[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
            vec![[0, 1, 2]],
        );
        let result = scale_uniform(&pd, 2.0, [0.0, 0.0, 0.0]);
        let p0 = result.points.get(0);
        assert!((p0[0] - 2.0).abs() < 1e-10);
        assert!((p0[1] - 0.0).abs() < 1e-10);
        let p1 = result.points.get(1);
        assert!((p1[1] - 2.0).abs() < 1e-10);
    }

    #[test]
    fn uniform_scale_from_center() {
        let pd = PolyData::from_triangles(
            vec![[2.0, 0.0, 0.0], [0.0, 2.0, 0.0], [0.0, 0.0, 0.0]],
            vec![[0, 1, 2]],
        );
        let center = [1.0, 1.0, 0.0];
        let result = scale_uniform(&pd, 3.0, center);
        // Point [2,0,0] -> center + 3*(p - center) = [1,1,0] + 3*[1,-1,0] = [4,-2,0]
        let p0 = result.points.get(0);
        assert!((p0[0] - 4.0).abs() < 1e-10);
        assert!((p0[1] - (-2.0)).abs() < 1e-10);
    }

    #[test]
    fn nonuniform_scale() {
        let pd = PolyData::from_triangles(
            vec![[1.0, 1.0, 1.0], [2.0, 2.0, 2.0], [3.0, 3.0, 3.0]],
            vec![[0, 1, 2]],
        );
        let result = scale_nonuniform(&pd, [2.0, 1.0, 0.5], [0.0, 0.0, 0.0]);
        let p0 = result.points.get(0);
        assert!((p0[0] - 2.0).abs() < 1e-10); // x*2
        assert!((p0[1] - 1.0).abs() < 1e-10); // y*1
        assert!((p0[2] - 0.5).abs() < 1e-10); // z*0.5
    }
}
