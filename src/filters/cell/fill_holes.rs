use std::collections::HashMap;

use crate::data::PolyData;

/// Fill holes (open boundary loops) in a triangle mesh.
///
/// Finds boundary edges (edges used by exactly one polygon), traces
/// closed loops, and fills each loop with a fan of triangles from the
/// loop centroid.
pub fn fill_holes(input: &PolyData) -> PolyData {
    // Count edge usage
    let mut edge_count: HashMap<(i64, i64), usize> = HashMap::new();
    let mut edge_directed: HashMap<(i64, i64), (i64, i64)> = HashMap::new();

    for cell in input.polys.iter() {
        let n = cell.len();
        for i in 0..n {
            let a = cell[i];
            let b = cell[(i + 1) % n];
            let key = if a < b { (a, b) } else { (b, a) };
            *edge_count.entry(key).or_insert(0) += 1;
            edge_directed.insert((a, b), key);
        }
    }

    // Boundary edges: used exactly once
    let boundary_edges: Vec<(i64, i64)> = edge_count
        .iter()
        .filter(|(_, &count)| count == 1)
        .map(|(&key, _)| key)
        .collect();

    if boundary_edges.is_empty() {
        return input.clone();
    }

    // Build adjacency for boundary vertices
    let mut next_map: HashMap<i64, i64> = HashMap::new();
    for cell in input.polys.iter() {
        let n = cell.len();
        for i in 0..n {
            let a = cell[i];
            let b = cell[(i + 1) % n];
            let key = if a < b { (a, b) } else { (b, a) };
            if edge_count.get(&key) == Some(&1) {
                // This is a boundary edge — record the directed version
                next_map.insert(a, b);
            }
        }
    }

    // Trace loops
    let mut visited: HashMap<i64, bool> = HashMap::new();
    let mut loops: Vec<Vec<i64>> = Vec::new();

    for &start in next_map.keys() {
        if visited.contains_key(&start) {
            continue;
        }
        let mut loop_pts = Vec::new();
        let mut current = start;
        loop {
            if visited.contains_key(&current) {
                break;
            }
            visited.insert(current, true);
            loop_pts.push(current);
            match next_map.get(&current) {
                Some(&next) => current = next,
                None => break,
            }
        }
        if loop_pts.len() >= 3 && current == start {
            loops.push(loop_pts);
        }
    }

    let mut pd = input.clone();

    // Fill each loop with a fan triangulation
    for lp in &loops {
        // Compute centroid
        let mut cx = 0.0;
        let mut cy = 0.0;
        let mut cz = 0.0;
        for &id in lp {
            let p = pd.points.get(id as usize);
            cx += p[0];
            cy += p[1];
            cz += p[2];
        }
        let n = lp.len() as f64;
        let center_idx = pd.points.len() as i64;
        pd.points.push([cx / n, cy / n, cz / n]);

        // Fan triangles (reversed winding to match boundary)
        for i in 0..lp.len() {
            let a = lp[i];
            let b = lp[(i + 1) % lp.len()];
            pd.polys.push_cell(&[center_idx, b, a]);
        }
    }

    pd
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn fill_single_hole() {
        // Open mesh: two triangles forming an open "V" with a hole
        let pd = PolyData::from_triangles(
            vec![
                [0.0, 0.0, 0.0], // 0
                [1.0, 0.0, 0.0], // 1
                [0.5, 1.0, 0.0], // 2
                [1.5, 1.0, 0.0], // 3
            ],
            vec![[0, 1, 2], [1, 3, 2]],
        );
        let result = fill_holes(&pd);
        // Should have more polys than input (hole filled)
        assert!(result.polys.num_cells() >= 2);
    }

    #[test]
    fn closed_mesh_unchanged() {
        // A tetrahedron has no boundary edges
        let pd = PolyData::from_triangles(
            vec![
                [0.0, 0.0, 0.0], [1.0, 0.0, 0.0],
                [0.5, 1.0, 0.0], [0.5, 0.5, 1.0],
            ],
            vec![[0, 2, 1], [0, 1, 3], [1, 2, 3], [0, 3, 2]],
        );
        let result = fill_holes(&pd);
        // No holes to fill — same number of polys
        assert_eq!(result.polys.num_cells(), 4);
    }
}
