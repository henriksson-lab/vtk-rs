use crate::data::{CellArray, Points, PolyData};

/// Parameters for generating a random point cloud.
pub struct PointSourceParams {
    /// Number of points to generate. Default: 100
    pub number_of_points: usize,
    /// Center of the point cloud. Default: [0, 0, 0]
    pub center: [f64; 3],
    /// Radius of the bounding sphere. Default: 0.5
    pub radius: f64,
    /// Random seed. Default: 42
    pub seed: u64,
}

impl Default for PointSourceParams {
    fn default() -> Self {
        Self {
            number_of_points: 100,
            center: [0.0, 0.0, 0.0],
            radius: 0.5,
            seed: 42,
        }
    }
}

/// Generate a random point cloud within a sphere.
///
/// Points are uniformly distributed within the sphere using rejection sampling
/// with a deterministic pseudo-random generator.
pub fn point_source(params: &PointSourceParams) -> PolyData {
    let mut points = Points::new();
    let mut verts = CellArray::new();
    let mut state = params.seed;

    let mut count = 0;
    while count < params.number_of_points {
        // Generate 3 random values in [-1, 1]
        let x = next_random(&mut state) * 2.0 - 1.0;
        let y = next_random(&mut state) * 2.0 - 1.0;
        let z = next_random(&mut state) * 2.0 - 1.0;

        // Rejection sampling: keep if inside unit sphere
        if x * x + y * y + z * z <= 1.0 {
            points.push([
                params.center[0] + x * params.radius,
                params.center[1] + y * params.radius,
                params.center[2] + z * params.radius,
            ]);
            verts.push_cell(&[count as i64]);
            count += 1;
        }
    }

    let mut pd = PolyData::new();
    pd.points = points;
    pd.verts = verts;
    pd
}

fn next_random(state: &mut u64) -> f64 {
    // xorshift64
    *state ^= *state << 13;
    *state ^= *state >> 7;
    *state ^= *state << 17;
    (*state & 0xFFFFFFFF) as f64 / 0xFFFFFFFF_u64 as f64
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn default_point_source() {
        let pd = point_source(&PointSourceParams::default());
        assert_eq!(pd.points.len(), 100);
        assert_eq!(pd.verts.num_cells(), 100);
    }

    #[test]
    fn points_within_radius() {
        let pd = point_source(&PointSourceParams {
            number_of_points: 50,
            radius: 1.0,
            center: [0.0, 0.0, 0.0],
            seed: 123,
        });
        for i in 0..pd.points.len() {
            let p = pd.points.get(i);
            let dist = (p[0] * p[0] + p[1] * p[1] + p[2] * p[2]).sqrt();
            assert!(dist <= 1.0 + 1e-10, "point {} at distance {}", i, dist);
        }
    }

    #[test]
    fn centered_point_source() {
        let pd = point_source(&PointSourceParams {
            number_of_points: 10,
            center: [5.0, 5.0, 5.0],
            radius: 0.1,
            seed: 99,
        });
        for i in 0..pd.points.len() {
            let p = pd.points.get(i);
            assert!((p[0] - 5.0).abs() <= 0.1 + 1e-10);
            assert!((p[1] - 5.0).abs() <= 0.1 + 1e-10);
            assert!((p[2] - 5.0).abs() <= 0.1 + 1e-10);
        }
    }
}
