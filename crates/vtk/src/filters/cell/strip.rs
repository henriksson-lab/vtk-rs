use crate::data::{CellArray, PolyData};

/// Convert triangle polygons to triangle strips.
///
/// Uses a simple greedy algorithm: picks an unvisited triangle, then
/// extends the strip by finding adjacent triangles sharing an edge.
pub fn to_triangle_strips(input: &PolyData) -> PolyData {
    let n_cells = input.polys.num_cells();
    if n_cells == 0 {
        return input.clone();
    }

    // Only works on triangles
    let cells: Vec<[i64; 3]> = input.polys.iter()
        .filter(|c| c.len() == 3)
        .map(|c| [c[0], c[1], c[2]])
        .collect();

    if cells.is_empty() {
        return input.clone();
    }

    // Build edge -> triangle adjacency
    use std::collections::HashMap;
    let mut edge_tris: HashMap<(i64, i64), Vec<usize>> = HashMap::new();
    for (ti, tri) in cells.iter().enumerate() {
        for e in 0..3 {
            let a = tri[e];
            let b = tri[(e + 1) % 3];
            let key = if a < b { (a, b) } else { (b, a) };
            edge_tris.entry(key).or_default().push(ti);
        }
    }

    let mut visited = vec![false; cells.len()];
    let mut out_strips = CellArray::new();

    for start in 0..cells.len() {
        if visited[start] {
            continue;
        }
        visited[start] = true;

        let mut strip = vec![cells[start][0], cells[start][1], cells[start][2]];

        // Try to extend the strip
        loop {
            let last_len = strip.len();
            let a = strip[last_len - 2];
            let b = strip[last_len - 1];
            let key = if a < b { (a, b) } else { (b, a) };

            let next = edge_tris.get(&key).and_then(|tris| {
                tris.iter().find(|&&ti| {
                    if visited[ti] { return false; }
                    let tri = &cells[ti];
                    // The third vertex (not a or b) is the new point
                    tri.iter().any(|&v| v != a && v != b)
                }).copied()
            });

            if let Some(ti) = next {
                visited[ti] = true;
                let tri = &cells[ti];
                let new_pt = tri.iter().find(|&&v| v != a && v != b).unwrap();
                strip.push(*new_pt);
            } else {
                break;
            }
        }

        out_strips.push_cell(&strip);
    }

    let mut pd = input.clone();
    pd.polys = CellArray::new();
    pd.strips = out_strips;
    pd
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn single_triangle_to_strip() {
        let pd = PolyData::from_triangles(
            vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.5, 1.0, 0.0]],
            vec![[0, 1, 2]],
        );
        let result = to_triangle_strips(&pd);
        assert_eq!(result.polys.num_cells(), 0);
        assert_eq!(result.strips.num_cells(), 1);
        assert_eq!(result.strips.cell(0).len(), 3);
    }

    #[test]
    fn two_adjacent_triangles() {
        let pd = PolyData::from_triangles(
            vec![
                [0.0, 0.0, 0.0], [1.0, 0.0, 0.0],
                [0.5, 1.0, 0.0], [1.5, 1.0, 0.0],
            ],
            vec![[0, 1, 2], [1, 3, 2]],
        );
        let result = to_triangle_strips(&pd);
        // Should form one strip of 4 vertices
        let total_strip_verts: usize = result.strips.iter().map(|s| s.len()).sum();
        assert!(total_strip_verts <= 5); // at most 4 in one strip or 3+3 in two
    }

    #[test]
    fn preserves_points() {
        let pd = PolyData::from_triangles(
            vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.5, 1.0, 0.0]],
            vec![[0, 1, 2]],
        );
        let result = to_triangle_strips(&pd);
        assert_eq!(result.points.len(), 3);
    }
}
