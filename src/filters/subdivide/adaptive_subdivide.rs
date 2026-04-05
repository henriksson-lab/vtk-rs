use crate::data::{CellArray, Points, PolyData};
use std::collections::HashMap;

/// Adaptively subdivide triangles that exceed a maximum edge length.
///
/// Iteratively splits triangles whose longest edge exceeds `max_edge_length`
/// by inserting midpoints on the longest edge. Repeats up to `max_passes` times.
pub fn adaptive_subdivide(input: &PolyData, max_edge_length: f64, max_passes: usize) -> PolyData {
    let mut points = input.points.clone();
    let mut current_tris: Vec<[i64; 3]> = Vec::new();

    // Collect initial triangles (fan-triangulate)
    for cell in input.polys.iter() {
        if cell.len() < 3 {
            continue;
        }
        for i in 1..cell.len() - 1 {
            current_tris.push([cell[0], cell[i], cell[i + 1]]);
        }
    }

    let max_len2 = max_edge_length * max_edge_length;

    for _ in 0..max_passes {
        let mut new_tris: Vec<[i64; 3]> = Vec::new();
        let mut midpoint_cache: HashMap<(i64, i64), i64> = HashMap::new();
        let mut any_split = false;

        for &[a, b, c] in &current_tris {
            let pa = points.get(a as usize);
            let pb = points.get(b as usize);
            let pc = points.get(c as usize);

            let d_ab = dist2(pa, pb);
            let d_bc = dist2(pb, pc);
            let d_ca = dist2(pc, pa);

            let max_d = d_ab.max(d_bc).max(d_ca);

            if max_d <= max_len2 {
                new_tris.push([a, b, c]);
                continue;
            }

            any_split = true;

            // Split the longest edge
            if d_ab >= d_bc && d_ab >= d_ca {
                let m = get_midpoint(&mut points, &mut midpoint_cache, a, b);
                new_tris.push([a, m, c]);
                new_tris.push([m, b, c]);
            } else if d_bc >= d_ca {
                let m = get_midpoint(&mut points, &mut midpoint_cache, b, c);
                new_tris.push([a, b, m]);
                new_tris.push([a, m, c]);
            } else {
                let m = get_midpoint(&mut points, &mut midpoint_cache, c, a);
                new_tris.push([a, b, m]);
                new_tris.push([m, b, c]);
            }
        }

        current_tris = new_tris;
        if !any_split {
            break;
        }
    }

    let mut polys = CellArray::new();
    for &[a, b, c] in &current_tris {
        polys.push_cell(&[a, b, c]);
    }

    let mut pd = PolyData::new();
    pd.points = points;
    pd.polys = polys;
    pd
}

fn get_midpoint(
    points: &mut Points<f64>,
    cache: &mut HashMap<(i64, i64), i64>,
    a: i64,
    b: i64,
) -> i64 {
    let key = if a < b { (a, b) } else { (b, a) };
    if let Some(&mid) = cache.get(&key) {
        return mid;
    }
    let pa = points.get(a as usize);
    let pb = points.get(b as usize);
    let mid_pt = [
        (pa[0] + pb[0]) * 0.5,
        (pa[1] + pb[1]) * 0.5,
        (pa[2] + pb[2]) * 0.5,
    ];
    let idx = points.len() as i64;
    points.push(mid_pt);
    cache.insert(key, idx);
    idx
}

fn dist2(a: [f64; 3], b: [f64; 3]) -> f64 {
    let dx = a[0] - b[0];
    let dy = a[1] - b[1];
    let dz = a[2] - b[2];
    dx * dx + dy * dy + dz * dz
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn subdivide_large_triangle() {
        let mut pd = PolyData::new();
        pd.points.push([0.0, 0.0, 0.0]);
        pd.points.push([10.0, 0.0, 0.0]);
        pd.points.push([5.0, 10.0, 0.0]);
        pd.polys.push_cell(&[0, 1, 2]);

        let result = adaptive_subdivide(&pd, 3.0, 10);
        assert!(result.polys.num_cells() > 1);
        assert!(result.points.len() > 3);
    }

    #[test]
    fn no_subdivide_small_triangle() {
        let mut pd = PolyData::new();
        pd.points.push([0.0, 0.0, 0.0]);
        pd.points.push([0.1, 0.0, 0.0]);
        pd.points.push([0.05, 0.1, 0.0]);
        pd.polys.push_cell(&[0, 1, 2]);

        let result = adaptive_subdivide(&pd, 1.0, 10);
        assert_eq!(result.polys.num_cells(), 1);
        assert_eq!(result.points.len(), 3);
    }

    #[test]
    fn max_passes_limits() {
        let mut pd = PolyData::new();
        pd.points.push([0.0, 0.0, 0.0]);
        pd.points.push([100.0, 0.0, 0.0]);
        pd.points.push([50.0, 100.0, 0.0]);
        pd.polys.push_cell(&[0, 1, 2]);

        let r1 = adaptive_subdivide(&pd, 1.0, 1);
        let r5 = adaptive_subdivide(&pd, 1.0, 5);
        assert!(r5.polys.num_cells() > r1.polys.num_cells());
    }

    #[test]
    fn empty_input() {
        let pd = PolyData::new();
        let result = adaptive_subdivide(&pd, 1.0, 10);
        assert_eq!(result.polys.num_cells(), 0);
    }

    #[test]
    fn shared_edges_use_same_midpoint() {
        let mut pd = PolyData::new();
        pd.points.push([0.0, 0.0, 0.0]);
        pd.points.push([10.0, 0.0, 0.0]);
        pd.points.push([10.0, 10.0, 0.0]);
        pd.points.push([0.0, 10.0, 0.0]);
        pd.polys.push_cell(&[0, 1, 2]);
        pd.polys.push_cell(&[0, 2, 3]);

        let result = adaptive_subdivide(&pd, 5.0, 1);
        // Shared edge 0-2 should create only one midpoint
        // No duplicate points for shared edges
        let n_pts = result.points.len();
        // With 4 original + shared midpoints, should be reasonable
        assert!(n_pts <= 10);
    }
}
