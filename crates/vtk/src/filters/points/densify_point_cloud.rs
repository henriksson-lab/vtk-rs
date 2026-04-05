//! DensifyPointCloudFilter — add interpolated points between distant neighbors.

use crate::data::{CellArray, Points, PolyData};

/// For each pair of points within `radius` that are farther apart than
/// `min_distance`, add a midpoint. Output is a denser PolyData with
/// vertex cells for all points (original + interpolated).
pub fn densify_point_cloud(input: &PolyData, radius: f64, min_distance: f64) -> PolyData {
    let n = input.points.len();
    let mut new_pts = Points::<f64>::new();

    // Copy original points
    for i in 0..n {
        new_pts.push(input.points.get(i));
    }

    let pts: Vec<[f64; 3]> = (0..n).map(|i| input.points.get(i)).collect();
    let radius2 = radius * radius;
    let min_dist2 = min_distance * min_distance;

    // For each pair within radius, add midpoints if they are far enough apart
    for i in 0..n {
        for j in (i + 1)..n {
            let dx = pts[i][0] - pts[j][0];
            let dy = pts[i][1] - pts[j][1];
            let dz = pts[i][2] - pts[j][2];
            let d2 = dx * dx + dy * dy + dz * dz;

            if d2 <= radius2 && d2 >= min_dist2 {
                let mid = [
                    (pts[i][0] + pts[j][0]) * 0.5,
                    (pts[i][1] + pts[j][1]) * 0.5,
                    (pts[i][2] + pts[j][2]) * 0.5,
                ];
                new_pts.push(mid);
            }
        }
    }

    // Build vertex cells for all points
    let total = new_pts.len();
    let mut verts = CellArray::new();
    for i in 0..total {
        verts.push_cell(&[i as i64]);
    }

    let mut result = PolyData::new();
    result.points = new_pts;
    result.verts = verts;
    result
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn adds_midpoints() {
        let mut pd = PolyData::new();
        pd.points.push([0.0, 0.0, 0.0]);
        pd.points.push([2.0, 0.0, 0.0]);

        // radius=3, min_distance=1: the pair is within radius and above min_distance
        let result = densify_point_cloud(&pd, 3.0, 1.0);
        assert_eq!(result.points.len(), 3); // 2 original + 1 midpoint
        let mid = result.points.get(2);
        assert!((mid[0] - 1.0).abs() < 1e-10);
    }

    #[test]
    fn no_midpoint_for_close_pair() {
        let mut pd = PolyData::new();
        pd.points.push([0.0, 0.0, 0.0]);
        pd.points.push([0.1, 0.0, 0.0]);

        // min_distance=1: pair is too close
        let result = densify_point_cloud(&pd, 10.0, 1.0);
        assert_eq!(result.points.len(), 2);
    }

    #[test]
    fn vertex_cells_for_all() {
        let mut pd = PolyData::new();
        pd.points.push([0.0, 0.0, 0.0]);
        pd.points.push([2.0, 0.0, 0.0]);

        let result = densify_point_cloud(&pd, 3.0, 1.0);
        assert_eq!(result.verts.num_cells(), result.points.len());
    }
}
