use vtk_data::PolyData;

/// Compute the Hausdorff distance between two point sets.
///
/// The Hausdorff distance is `max(d(A→B), d(B→A))` where
/// `d(A→B) = max over a in A of (min over b in B of dist(a, b))`.
///
/// Returns `(hausdorff_distance, mean_distance_a_to_b, mean_distance_b_to_a)`.
pub fn hausdorff_distance(a: &PolyData, b: &PolyData) -> (f64, f64, f64) {
    let (max_ab, mean_ab) = directed_hausdorff(a, b);
    let (max_ba, mean_ba) = directed_hausdorff(b, a);
    (max_ab.max(max_ba), mean_ab, mean_ba)
}

/// Compute the directed Hausdorff distance from A to B.
/// Returns (max_distance, mean_distance).
fn directed_hausdorff(a: &PolyData, b: &PolyData) -> (f64, f64) {
    let na = a.points.len();
    let nb = b.points.len();
    if na == 0 || nb == 0 {
        return (0.0, 0.0);
    }

    let mut max_d = 0.0f64;
    let mut sum_d = 0.0f64;

    for i in 0..na {
        let pa = a.points.get(i);
        let mut min_d2 = f64::MAX;
        for j in 0..nb {
            let pb = b.points.get(j);
            let d2 = (pa[0] - pb[0]) * (pa[0] - pb[0])
                + (pa[1] - pb[1]) * (pa[1] - pb[1])
                + (pa[2] - pb[2]) * (pa[2] - pb[2]);
            if d2 < min_d2 {
                min_d2 = d2;
            }
        }
        let d = min_d2.sqrt();
        max_d = max_d.max(d);
        sum_d += d;
    }

    (max_d, sum_d / na as f64)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn identical_sets() {
        let mut pd = PolyData::new();
        pd.points.push([0.0, 0.0, 0.0]);
        pd.points.push([1.0, 0.0, 0.0]);
        let (h, _, _) = hausdorff_distance(&pd, &pd);
        assert!(h < 1e-10);
    }

    #[test]
    fn known_distance() {
        let mut a = PolyData::new();
        a.points.push([0.0, 0.0, 0.0]);

        let mut b = PolyData::new();
        b.points.push([3.0, 4.0, 0.0]);

        let (h, _, _) = hausdorff_distance(&a, &b);
        assert!((h - 5.0).abs() < 1e-10);
    }

    #[test]
    fn asymmetric_sets() {
        let mut a = PolyData::new();
        a.points.push([0.0, 0.0, 0.0]);
        a.points.push([1.0, 0.0, 0.0]);

        let mut b = PolyData::new();
        b.points.push([0.0, 0.0, 0.0]);
        b.points.push([1.0, 0.0, 0.0]);
        b.points.push([5.0, 0.0, 0.0]); // extra far point

        let (h, mean_ab, mean_ba) = hausdorff_distance(&a, &b);
        // Hausdorff should be driven by b→a: point (5,0,0) is 4 from nearest in a
        assert!((h - 4.0).abs() < 1e-10);
        assert!(mean_ab < mean_ba); // a→b is closer on average
    }
}
