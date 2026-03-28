use vtk_data::PolyData;

/// Result of computing distance statistics between two meshes.
#[derive(Debug, Clone, Copy)]
pub struct MeshDistanceResult {
    /// Maximum distance from any point in A to its closest point in B.
    pub max_a_to_b: f64,
    /// Maximum distance from any point in B to its closest point in A.
    pub max_b_to_a: f64,
    /// Mean distance from points in A to their closest points in B.
    pub mean_a_to_b: f64,
    /// Mean distance from points in B to their closest points in A.
    pub mean_b_to_a: f64,
}

/// Compute symmetric Hausdorff-like distance statistics between two meshes.
///
/// For each point in mesh A, finds the closest point in mesh B (brute force)
/// and vice versa. Returns max and mean distances in both directions.
pub fn mesh_distance_stats(a: &PolyData, b: &PolyData) -> MeshDistanceResult {
    let (max_ab, mean_ab) = directed_distance(a, b);
    let (max_ba, mean_ba) = directed_distance(b, a);
    MeshDistanceResult {
        max_a_to_b: max_ab,
        max_b_to_a: max_ba,
        mean_a_to_b: mean_ab,
        mean_b_to_a: mean_ba,
    }
}

/// Compute directed distance stats from A to B.
/// Returns (max_distance, mean_distance).
fn directed_distance(a: &PolyData, b: &PolyData) -> (f64, f64) {
    let na: usize = a.points.len();
    let nb: usize = b.points.len();
    if na == 0 || nb == 0 {
        return (0.0, 0.0);
    }

    let mut max_d: f64 = 0.0;
    let mut sum_d: f64 = 0.0;

    for i in 0..na {
        let pa = a.points.get(i);
        let mut min_d2: f64 = f64::MAX;
        for j in 0..nb {
            let pb = b.points.get(j);
            let dx: f64 = pa[0] - pb[0];
            let dy: f64 = pa[1] - pb[1];
            let dz: f64 = pa[2] - pb[2];
            let d2: f64 = dx * dx + dy * dy + dz * dz;
            if d2 < min_d2 {
                min_d2 = d2;
            }
        }
        let d: f64 = min_d2.sqrt();
        if d > max_d {
            max_d = d;
        }
        sum_d += d;
    }

    (max_d, sum_d / na as f64)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn identical_meshes() {
        let mut a = PolyData::new();
        a.points.push([0.0, 0.0, 0.0]);
        a.points.push([1.0, 0.0, 0.0]);
        a.points.push([0.0, 1.0, 0.0]);

        let result = mesh_distance_stats(&a, &a);
        assert!(result.max_a_to_b < 1e-10);
        assert!(result.max_b_to_a < 1e-10);
        assert!(result.mean_a_to_b < 1e-10);
        assert!(result.mean_b_to_a < 1e-10);
    }

    #[test]
    fn known_distance() {
        let mut a = PolyData::new();
        a.points.push([0.0, 0.0, 0.0]);

        let mut b = PolyData::new();
        b.points.push([3.0, 4.0, 0.0]);

        let result = mesh_distance_stats(&a, &b);
        assert!((result.max_a_to_b - 5.0).abs() < 1e-10);
        assert!((result.max_b_to_a - 5.0).abs() < 1e-10);
        assert!((result.mean_a_to_b - 5.0).abs() < 1e-10);
    }

    #[test]
    fn asymmetric_distances() {
        let mut a = PolyData::new();
        a.points.push([0.0, 0.0, 0.0]);

        let mut b = PolyData::new();
        b.points.push([0.0, 0.0, 0.0]);
        b.points.push([10.0, 0.0, 0.0]);

        let result = mesh_distance_stats(&a, &b);
        // A->B: point (0,0,0) closest to (0,0,0) = 0
        assert!(result.max_a_to_b < 1e-10);
        // B->A: point (10,0,0) closest to (0,0,0) = 10
        assert!((result.max_b_to_a - 10.0).abs() < 1e-10);
        // mean B->A: (0 + 10) / 2 = 5
        assert!((result.mean_b_to_a - 5.0).abs() < 1e-10);
    }
}
