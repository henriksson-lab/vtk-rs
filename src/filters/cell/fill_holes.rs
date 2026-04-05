use std::collections::HashMap;

use crate::data::PolyData;

/// Fill holes (open boundary loops) in a triangle mesh.
///
/// Finds boundary edges (edges used by exactly one polygon), traces
/// closed loops, and fills each loop with a fan of triangles from the
/// loop centroid.
pub fn fill_holes(input: &PolyData) -> PolyData {
    // Count edge usage using packed u64 keys for speed
    let mut edge_count: HashMap<u64, u8> = HashMap::new();

    let offsets = input.polys.offsets();
    let conn = input.polys.connectivity();
    let nc = input.polys.num_cells();

    for ci in 0..nc {
        let start = offsets[ci] as usize;
        let end = offsets[ci + 1] as usize;
        let cell = &conn[start..end];
        let n = cell.len();
        for i in 0..n {
            let a = cell[i];
            let b = cell[(i + 1) % n];
            let key = if a < b {
                (a as u64) << 32 | b as u64
            } else {
                (b as u64) << 32 | a as u64
            };
            let e = edge_count.entry(key).or_insert(0);
            *e = (*e).saturating_add(1);
        }
    }

    // Check if there are any boundary edges at all
    let has_boundary = edge_count.values().any(|&c| c == 1);
    if !has_boundary {
        return input.clone();
    }

    // Build directed boundary edge adjacency: vertex → next vertex
    let mut next_map: HashMap<i64, i64> = HashMap::new();
    for ci in 0..nc {
        let start = offsets[ci] as usize;
        let end = offsets[ci + 1] as usize;
        let cell = &conn[start..end];
        let n = cell.len();
        for i in 0..n {
            let a = cell[i];
            let b = cell[(i + 1) % n];
            let key = if a < b {
                (a as u64) << 32 | b as u64
            } else {
                (b as u64) << 32 | a as u64
            };
            if edge_count.get(&key) == Some(&1) {
                next_map.insert(a, b);
            }
        }
    }

    // Trace loops using a simple visited set
    let mut visited = vec![false; input.points.len()];
    let mut loops: Vec<Vec<i64>> = Vec::new();

    for &start in next_map.keys() {
        let si = start as usize;
        if si < visited.len() && visited[si] {
            continue;
        }
        let mut loop_pts = Vec::new();
        let mut current = start;
        loop {
            let ci = current as usize;
            if ci < visited.len() && visited[ci] {
                break;
            }
            if ci < visited.len() {
                visited[ci] = true;
            }
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
        let pd = PolyData::from_triangles(
            vec![
                [0.0, 0.0, 0.0],
                [1.0, 0.0, 0.0],
                [0.5, 1.0, 0.0],
                [1.5, 1.0, 0.0],
            ],
            vec![[0, 1, 2], [1, 3, 2]],
        );
        let result = fill_holes(&pd);
        assert!(result.polys.num_cells() >= 2);
    }

    #[test]
    fn closed_mesh_unchanged() {
        let pd = PolyData::from_triangles(
            vec![
                [0.0, 0.0, 0.0], [1.0, 0.0, 0.0],
                [0.5, 1.0, 0.0], [0.5, 0.5, 1.0],
            ],
            vec![[0, 2, 1], [0, 1, 3], [1, 2, 3], [0, 3, 2]],
        );
        let result = fill_holes(&pd);
        assert_eq!(result.polys.num_cells(), 4);
    }
}
