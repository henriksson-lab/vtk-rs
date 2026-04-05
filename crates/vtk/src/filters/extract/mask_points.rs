use crate::data::{CellArray, Points, PolyData};

/// Subsample points from a PolyData by keeping every Nth point.
///
/// The output contains only vertex cells (one per kept point).
/// Useful for reducing point clouds before glyphing.
pub fn mask_points(input: &PolyData, every_nth: usize) -> PolyData {
    let every_nth = every_nth.max(1);
    let mut out_points = Points::<f64>::new();
    let mut out_verts = CellArray::new();

    let mut idx = 0;
    let mut i = 0;
    while i < input.points.len() {
        out_points.push(input.points.get(i));
        out_verts.push_cell(&[idx]);
        idx += 1;
        i += every_nth;
    }

    let mut pd = PolyData::new();
    pd.points = out_points;
    pd.verts = out_verts;
    pd
}

/// Subsample points randomly using a deterministic seed.
///
/// Keeps approximately `ratio` fraction of points (0.0–1.0).
pub fn mask_points_random(input: &PolyData, ratio: f64, seed: u64) -> PolyData {
    let ratio = ratio.clamp(0.0, 1.0);
    let mut out_points = Points::<f64>::new();
    let mut out_verts = CellArray::new();

    // Simple deterministic hash-based selection
    let mut state = seed;
    let mut idx = 0i64;
    for i in 0..input.points.len() {
        // xorshift64
        state ^= state << 13;
        state ^= state >> 7;
        state ^= state << 17;
        let r = (state & 0xFFFFFFFF) as f64 / 0xFFFFFFFF_u64 as f64;
        if r < ratio {
            out_points.push(input.points.get(i));
            out_verts.push_cell(&[idx]);
            idx += 1;
        }
    }

    let mut pd = PolyData::new();
    pd.points = out_points;
    pd.verts = out_verts;
    pd
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn every_second_point() {
        let mut pd = PolyData::new();
        for i in 0..10 {
            pd.points.push([i as f64, 0.0, 0.0]);
        }
        let result = mask_points(&pd, 2);
        assert_eq!(result.points.len(), 5); // 0, 2, 4, 6, 8
        assert_eq!(result.verts.num_cells(), 5);
        let p = result.points.get(2);
        assert!((p[0] - 4.0).abs() < 1e-10);
    }

    #[test]
    fn every_point() {
        let mut pd = PolyData::new();
        for i in 0..5 {
            pd.points.push([i as f64, 0.0, 0.0]);
        }
        let result = mask_points(&pd, 1);
        assert_eq!(result.points.len(), 5);
    }

    #[test]
    fn random_mask() {
        let mut pd = PolyData::new();
        for i in 0..100 {
            pd.points.push([i as f64, 0.0, 0.0]);
        }
        let result = mask_points_random(&pd, 0.5, 42);
        // Should keep roughly half
        assert!(result.points.len() > 20);
        assert!(result.points.len() < 80);
    }
}
