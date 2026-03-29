use std::collections::HashMap;

use vtk_data::{CellArray, PolyData};

/// Orient all polygons in a PolyData to have consistent winding order.
///
/// Uses a breadth-first traversal starting from the first polygon.
/// Neighboring polygons sharing an edge are flipped if their shared edge
/// has the same direction (indicating inconsistent winding). Only works
/// on manifold meshes (each edge shared by at most two polygons).
pub fn orient(input: &PolyData) -> PolyData {
    let n_cells = input.polys.num_cells();
    if n_cells == 0 {
        return input.clone();
    }

    // Build edge -> cell adjacency
    // Edge key: (min_id, max_id) -> list of (cell_idx, edge_direction_matches_key)
    let mut edge_cells: HashMap<(i64, i64), Vec<(usize, bool)>> = HashMap::new();

    let cells: Vec<Vec<i64>> = input.polys.iter().map(|c| c.to_vec()).collect();

    for (ci, cell) in cells.iter().enumerate() {
        let n = cell.len();
        for i in 0..n {
            let a = cell[i];
            let b = cell[(i + 1) % n];
            let key = if a < b { (a, b) } else { (b, a) };
            let forward = a < b;
            edge_cells.entry(key).or_default().push((ci, forward));
        }
    }

    // BFS to determine which cells need flipping
    let mut oriented: Vec<Option<bool>> = vec![None; n_cells]; // None=unvisited, Some(false)=keep, Some(true)=flip
    let mut queue = std::collections::VecDeque::new();

    // Process each connected component
    for start in 0..n_cells {
        if oriented[start].is_some() {
            continue;
        }
        oriented[start] = Some(false); // Keep first cell as-is
        queue.push_back(start);

        while let Some(ci) = queue.pop_front() {
            let cell = &cells[ci];
            let ci_flipped = oriented[ci].unwrap();
            let n = cell.len();

            for i in 0..n {
                let a = cell[i];
                let b = cell[(i + 1) % n];
                let key = if a < b { (a, b) } else { (b, a) };
                let ci_forward = a < b;
                // If ci is flipped, the actual direction is reversed
                let ci_actual_forward = ci_forward ^ ci_flipped;

                if let Some(neighbors) = edge_cells.get(&key) {
                    for &(ni, ni_forward) in neighbors {
                        if ni == ci || oriented[ni].is_some() {
                            continue;
                        }
                        // For consistent winding, neighboring cells should traverse
                        // the shared edge in opposite directions
                        let should_flip = ci_actual_forward == ni_forward;
                        oriented[ni] = Some(should_flip);
                        queue.push_back(ni);
                    }
                }
            }
        }
    }

    // Build output
    let mut out_polys = CellArray::new();
    for (ci, cell) in cells.iter().enumerate() {
        if oriented[ci] == Some(true) {
            let reversed: Vec<i64> = cell.iter().rev().copied().collect();
            out_polys.push_cell(&reversed);
        } else {
            out_polys.push_cell(cell);
        }
    }

    let mut pd = input.clone();
    pd.polys = out_polys;
    pd
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn already_consistent() {
        let pd = PolyData::from_triangles(
            vec![
                [0.0, 0.0, 0.0],
                [1.0, 0.0, 0.0],
                [0.5, 1.0, 0.0],
                [1.5, 1.0, 0.0],
            ],
            // Both CCW when viewed from +Z
            vec![[0, 1, 2], [1, 3, 2]],
        );
        let result = orient(&pd);
        assert_eq!(result.polys.num_cells(), 2);
        // Check that shared edge (1,2) is traversed in opposite directions
        let c0 = result.polys.cell(0);
        let c1 = result.polys.cell(1);
        // In c0: 1->2, in c1: should be 2->1 (opposite)
        let edge_forward_in_c0 = has_directed_edge(c0, 1, 2);
        let edge_forward_in_c1 = has_directed_edge(c1, 1, 2);
        assert_ne!(edge_forward_in_c0, edge_forward_in_c1);
    }

    #[test]
    fn fixes_inconsistent_winding() {
        let pd = PolyData::from_triangles(
            vec![
                [0.0, 0.0, 0.0],
                [1.0, 0.0, 0.0],
                [0.5, 1.0, 0.0],
                [1.5, 1.0, 0.0],
            ],
            // Second triangle has same winding direction for shared edge (inconsistent)
            vec![[0, 1, 2], [1, 2, 3]],
        );
        let result = orient(&pd);
        assert_eq!(result.polys.num_cells(), 2);
        // After orientation, shared edge should go in opposite directions
        let c0 = result.polys.cell(0);
        let c1 = result.polys.cell(1);
        let edge_forward_in_c0 = has_directed_edge(c0, 1, 2);
        let edge_forward_in_c1 = has_directed_edge(c1, 1, 2);
        assert_ne!(edge_forward_in_c0, edge_forward_in_c1);
    }

    fn has_directed_edge(cell: &[i64], a: i64, b: i64) -> bool {
        let n = cell.len();
        for i in 0..n {
            if cell[i] == a && cell[(i + 1) % n] == b {
                return true;
            }
        }
        false
    }
}
