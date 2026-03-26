use vtk_data::{AnyDataArray, DataArray, PolyData};

/// Quality metric to compute.
#[derive(Debug, Clone, Copy)]
pub enum QualityMetric {
    /// Aspect ratio (longest edge / shortest edge).
    AspectRatio,
    /// Minimum interior angle (degrees).
    MinAngle,
    /// Maximum interior angle (degrees).
    MaxAngle,
    /// Area of the polygon.
    Area,
}

/// Compute a quality metric for each polygon cell.
///
/// Adds a "Quality" scalar array to cell data.
pub fn cell_quality(input: &PolyData, metric: QualityMetric) -> PolyData {
    let mut values = Vec::new();

    for cell in input.polys.iter() {
        let n = cell.len();
        if n < 3 {
            values.push(0.0);
            continue;
        }

        let pts: Vec<[f64; 3]> = cell.iter().map(|&id| input.points.get(id as usize)).collect();

        let val = match metric {
            QualityMetric::AspectRatio => {
                let edges = edge_lengths(&pts);
                let min_e = edges.iter().cloned().fold(f64::MAX, f64::min);
                let max_e = edges.iter().cloned().fold(0.0f64, f64::max);
                if min_e > 1e-20 { max_e / min_e } else { f64::MAX }
            }
            QualityMetric::MinAngle => {
                let angles = interior_angles(&pts);
                angles.iter().cloned().fold(f64::MAX, f64::min).to_degrees()
            }
            QualityMetric::MaxAngle => {
                let angles = interior_angles(&pts);
                angles.iter().cloned().fold(0.0f64, f64::max).to_degrees()
            }
            QualityMetric::Area => {
                polygon_area(&pts)
            }
        };
        values.push(val);
    }

    let mut pd = input.clone();
    pd.cell_data_mut().add_array(AnyDataArray::F64(
        DataArray::from_vec("Quality", values, 1),
    ));
    pd
}

fn edge_lengths(pts: &[[f64; 3]]) -> Vec<f64> {
    let n = pts.len();
    (0..n).map(|i| {
        let j = (i + 1) % n;
        let d = [pts[j][0] - pts[i][0], pts[j][1] - pts[i][1], pts[j][2] - pts[i][2]];
        (d[0] * d[0] + d[1] * d[1] + d[2] * d[2]).sqrt()
    }).collect()
}

fn interior_angles(pts: &[[f64; 3]]) -> Vec<f64> {
    let n = pts.len();
    (0..n).map(|i| {
        let prev = if i == 0 { n - 1 } else { i - 1 };
        let next = (i + 1) % n;
        let a = [pts[prev][0] - pts[i][0], pts[prev][1] - pts[i][1], pts[prev][2] - pts[i][2]];
        let b = [pts[next][0] - pts[i][0], pts[next][1] - pts[i][1], pts[next][2] - pts[i][2]];
        let la = (a[0] * a[0] + a[1] * a[1] + a[2] * a[2]).sqrt();
        let lb = (b[0] * b[0] + b[1] * b[1] + b[2] * b[2]).sqrt();
        if la < 1e-20 || lb < 1e-20 { return 0.0; }
        let cos_angle = ((a[0] * b[0] + a[1] * b[1] + a[2] * b[2]) / (la * lb)).clamp(-1.0, 1.0);
        cos_angle.acos()
    }).collect()
}

fn polygon_area(pts: &[[f64; 3]]) -> f64 {
    if pts.len() < 3 { return 0.0; }
    let mut total = 0.0;
    for i in 1..pts.len() - 1 {
        let e1 = [pts[i][0] - pts[0][0], pts[i][1] - pts[0][1], pts[i][2] - pts[0][2]];
        let e2 = [pts[i+1][0] - pts[0][0], pts[i+1][1] - pts[0][1], pts[i+1][2] - pts[0][2]];
        let c = [e1[1]*e2[2]-e1[2]*e2[1], e1[2]*e2[0]-e1[0]*e2[2], e1[0]*e2[1]-e1[1]*e2[0]];
        total += 0.5 * (c[0]*c[0] + c[1]*c[1] + c[2]*c[2]).sqrt();
    }
    total
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn equilateral_aspect_ratio() {
        let pd = PolyData::from_triangles(
            vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.5, (3.0f64).sqrt() / 2.0, 0.0]],
            vec![[0, 1, 2]],
        );
        let result = cell_quality(&pd, QualityMetric::AspectRatio);
        let arr = result.cell_data().get_array("Quality").unwrap();
        let mut val = [0.0f64];
        arr.tuple_as_f64(0, &mut val);
        assert!((val[0] - 1.0).abs() < 1e-10); // equilateral: all edges equal
    }

    #[test]
    fn right_triangle_angles() {
        let pd = PolyData::from_triangles(
            vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]],
            vec![[0, 1, 2]],
        );
        let result = cell_quality(&pd, QualityMetric::MinAngle);
        let arr = result.cell_data().get_array("Quality").unwrap();
        let mut val = [0.0f64];
        arr.tuple_as_f64(0, &mut val);
        assert!((val[0] - 45.0).abs() < 1e-6);
    }

    #[test]
    fn triangle_area_metric() {
        let pd = PolyData::from_triangles(
            vec![[0.0, 0.0, 0.0], [2.0, 0.0, 0.0], [0.0, 3.0, 0.0]],
            vec![[0, 1, 2]],
        );
        let result = cell_quality(&pd, QualityMetric::Area);
        let arr = result.cell_data().get_array("Quality").unwrap();
        let mut val = [0.0f64];
        arr.tuple_as_f64(0, &mut val);
        assert!((val[0] - 3.0).abs() < 1e-10);
    }
}
