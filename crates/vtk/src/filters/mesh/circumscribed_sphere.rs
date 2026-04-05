use crate::data::PolyData;

/// Compute the smallest circumscribed (bounding) sphere of a mesh using
/// Ritter's algorithm.
///
/// Returns `(radius, center)`.
pub fn bounding_sphere(input: &PolyData) -> (f64, [f64; 3]) {
    let n: usize = input.points.len();
    if n == 0 {
        return (0.0, [0.0, 0.0, 0.0]);
    }
    if n == 1 {
        return (0.0, input.points.get(0));
    }

    // Step 1: Find an initial diameter by picking the two most distant points
    // along a coordinate axis, then refining.
    let p0 = input.points.get(0);

    // Find point farthest from p0
    let mut far1_idx: usize = 0;
    let mut far1_dist: f64 = 0.0;
    for i in 0..n {
        let d: f64 = dist_sq(&p0, &input.points.get(i));
        if d > far1_dist {
            far1_dist = d;
            far1_idx = i;
        }
    }

    // Find point farthest from far1
    let far1 = input.points.get(far1_idx);
    let mut far2_idx: usize = 0;
    let mut far2_dist: f64 = 0.0;
    for i in 0..n {
        let d: f64 = dist_sq(&far1, &input.points.get(i));
        if d > far2_dist {
            far2_dist = d;
            far2_idx = i;
        }
    }

    // Initial sphere: centered between far1 and far2
    let a = input.points.get(far1_idx);
    let b = input.points.get(far2_idx);
    let mut cx: f64 = (a[0] + b[0]) * 0.5;
    let mut cy: f64 = (a[1] + b[1]) * 0.5;
    let mut cz: f64 = (a[2] + b[2]) * 0.5;
    let mut radius: f64 = (far2_dist).sqrt() * 0.5;

    // Step 2: Ritter's expansion — grow sphere to include all points
    for i in 0..n {
        let p = input.points.get(i);
        let dx: f64 = p[0] - cx;
        let dy: f64 = p[1] - cy;
        let dz: f64 = p[2] - cz;
        let d: f64 = (dx * dx + dy * dy + dz * dz).sqrt();

        if d > radius {
            let new_radius: f64 = (radius + d) * 0.5;
            let shift: f64 = new_radius - radius;
            let inv_d: f64 = 1.0 / d;
            cx += dx * inv_d * shift;
            cy += dy * inv_d * shift;
            cz += dz * inv_d * shift;
            radius = new_radius;
        }
    }

    // Step 3: Second pass to ensure all points are enclosed (refine)
    for i in 0..n {
        let p = input.points.get(i);
        let dx: f64 = p[0] - cx;
        let dy: f64 = p[1] - cy;
        let dz: f64 = p[2] - cz;
        let d: f64 = (dx * dx + dy * dy + dz * dz).sqrt();

        if d > radius {
            let new_radius: f64 = (radius + d) * 0.5;
            let shift: f64 = new_radius - radius;
            let inv_d: f64 = 1.0 / d;
            cx += dx * inv_d * shift;
            cy += dy * inv_d * shift;
            cz += dz * inv_d * shift;
            radius = new_radius;
        }
    }

    (radius, [cx, cy, cz])
}

fn dist_sq(a: &[f64; 3], b: &[f64; 3]) -> f64 {
    let dx: f64 = a[0] - b[0];
    let dy: f64 = a[1] - b[1];
    let dz: f64 = a[2] - b[2];
    dx * dx + dy * dy + dz * dz
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn single_point() {
        let mut pd = PolyData::new();
        pd.points.push([3.0, 4.0, 5.0]);
        let (r, c) = bounding_sphere(&pd);
        assert_eq!(r, 0.0);
        assert_eq!(c, [3.0, 4.0, 5.0]);
    }

    #[test]
    fn two_points() {
        let mut pd = PolyData::new();
        pd.points.push([0.0, 0.0, 0.0]);
        pd.points.push([6.0, 0.0, 0.0]);
        let (r, c) = bounding_sphere(&pd);
        assert!((r - 3.0).abs() < 1e-10, "r={}", r);
        assert!((c[0] - 3.0).abs() < 1e-10);
        assert!(c[1].abs() < 1e-10);
        assert!(c[2].abs() < 1e-10);
    }

    #[test]
    fn cube_vertices() {
        let mut pd = PolyData::new();
        for &x in &[-1.0f64, 1.0] {
            for &y in &[-1.0f64, 1.0] {
                for &z in &[-1.0f64, 1.0] {
                    pd.points.push([x, y, z]);
                }
            }
        }
        let (r, c) = bounding_sphere(&pd);
        // Circumscribed sphere of unit cube centered at origin: r = sqrt(3)
        let expected: f64 = 3.0f64.sqrt();
        assert!((r - expected).abs() < 0.1, "r={}, expected={}", r, expected);
        // Center should be near origin
        assert!(c[0].abs() < 0.1);
        assert!(c[1].abs() < 0.1);
        assert!(c[2].abs() < 0.1);
    }
}
