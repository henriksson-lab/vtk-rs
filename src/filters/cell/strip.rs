use std::collections::HashMap;
use crate::data::{CellArray, PolyData};

/// Convert triangle polygons to triangle strips.
///
/// Uses a greedy algorithm: picks an unvisited triangle, then
/// extends the strip by finding adjacent triangles sharing an edge.
pub fn to_triangle_strips(input: &PolyData) -> PolyData {
    let nc = input.polys.num_cells();
    if nc == 0 {
        return input.clone();
    }

    let offsets = input.polys.offsets();
    let conn = input.polys.connectivity();

    // Collect triangle indices and vertex ids using raw access
    let mut tri_verts: Vec<[i64; 3]> = Vec::with_capacity(nc);
    for ci in 0..nc {
        let start = offsets[ci] as usize;
        let end = offsets[ci + 1] as usize;
        if end - start == 3 {
            tri_verts.push([conn[start], conn[start + 1], conn[start + 2]]);
        }
    }

    let nt = tri_verts.len();
    if nt == 0 {
        return input.clone();
    }

    // Build per-triangle adjacency: for each triangle edge, store the neighboring triangle.
    // tri_adj[ti*3 + edge_idx] = neighbor tri index (u32::MAX if boundary).
    // Uses sorted edge list to pair triangles sharing each edge.
    let mut edge_list: Vec<(u64, u32, u8)> = Vec::with_capacity(nt * 3); // (key, tri, edge_idx)
    for (ti, tri) in tri_verts.iter().enumerate() {
        for e in 0u8..3 {
            let a = tri[e as usize];
            let b = tri[((e + 1) % 3) as usize];
            let key = if a < b { (a as u64) << 32 | b as u64 } else { (b as u64) << 32 | a as u64 };
            edge_list.push((key, ti as u32, e));
        }
    }
    edge_list.sort_unstable_by_key(|&(k, _, _)| k);

    let mut tri_adj = vec![u32::MAX; nt * 3]; // tri_adj[ti*3 + edge] = neighbor tri
    let n_edges = edge_list.len();
    let mut i = 0;
    while i < n_edges {
        let key = edge_list[i].0;
        if i + 1 < n_edges && edge_list[i + 1].0 == key {
            let (_, t0, e0) = edge_list[i];
            let (_, t1, e1) = edge_list[i + 1];
            tri_adj[t0 as usize * 3 + e0 as usize] = t1;
            tri_adj[t1 as usize * 3 + e1 as usize] = t0;
            i += 2;
        } else {
            i += 1;
        }
    }

    // Build edge -> triangle lookup for greedy extension (HashMap, but built from sorted data)
    let mut edge_tris: HashMap<u64, [u32; 2]> = HashMap::with_capacity(nt * 3 / 2);
    i = 0;
    while i < n_edges {
        let key = edge_list[i].0;
        let t0 = edge_list[i].1;
        let mut t1 = u32::MAX;
        if i + 1 < n_edges && edge_list[i + 1].0 == key {
            t1 = edge_list[i + 1].1;
            i += 1;
        }
        edge_tris.insert(key, [t0, t1]);
        i += 1;
    }

    let mut visited = vec![false; nt];

    // Pre-size strip output
    let mut strip_off: Vec<i64> = Vec::with_capacity(nt / 2 + 1);
    let mut strip_conn: Vec<i64> = Vec::with_capacity(nt * 2);
    strip_off.push(0);

    for start in 0..nt {
        if visited[start] {
            continue;
        }
        visited[start] = true;

        let tri = &tri_verts[start];
        let _strip_start = strip_conn.len();
        strip_conn.push(tri[0]);
        strip_conn.push(tri[1]);
        strip_conn.push(tri[2]);

        // Extend strip greedily
        loop {
            let len = strip_conn.len();
            let a = strip_conn[len - 2];
            let b = strip_conn[len - 1];
            let key = if a < b { (a as u64) << 32 | b as u64 } else { (b as u64) << 32 | a as u64 };

            let next = edge_tris.get(&key).and_then(|pair| {
                for &ti in pair {
                    if ti == u32::MAX { break; }
                    let ti = ti as usize;
                    if visited[ti] { continue; }
                    let t = &tri_verts[ti];
                    // Find the third vertex (not a or b)
                    for &v in t {
                        if v != a && v != b {
                            return Some((ti, v));
                        }
                    }
                }
                None
            });

            if let Some((ti, new_pt)) = next {
                visited[ti] = true;
                strip_conn.push(new_pt);
            } else {
                break;
            }
        }

        strip_off.push(strip_conn.len() as i64);
    }

    let mut pd = PolyData::new();
    pd.points = input.points.clone();
    pd.strips = CellArray::from_raw(strip_off, strip_conn);
    pd.lines = input.lines.clone();
    pd.verts = input.verts.clone();
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
        let total_strip_verts: usize = result.strips.iter().map(|s| s.len()).sum();
        assert!(total_strip_verts <= 5);
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
