use std::collections::{HashMap, HashSet};

use crate::data::{CellArray, Points, PolyData};

/// Extract closed boundary loops from a mesh.
///
/// Boundary edges are those used by exactly one polygon. The filter traces
/// connected boundary edges into closed loops and returns a PolyData with
/// line cells forming those loops, along with the number of loops found.
pub fn extract_boundary_loops(input: &PolyData) -> (PolyData, usize) {
    // Count edge usage
    let mut edge_count: HashMap<(i64, i64), usize> = HashMap::new();

    for cell in input.polys.iter() {
        let n: usize = cell.len();
        for i in 0..n {
            let a: i64 = cell[i];
            let b: i64 = cell[(i + 1) % n];
            let key = if a < b { (a, b) } else { (b, a) };
            *edge_count.entry(key).or_insert(0) += 1;
        }
    }

    // Build adjacency for boundary edges only
    let mut adj: HashMap<i64, Vec<i64>> = HashMap::new();
    for (&(a, b), &count) in &edge_count {
        if count == 1 {
            adj.entry(a).or_default().push(b);
            adj.entry(b).or_default().push(a);
        }
    }

    let mut out_points = Points::<f64>::new();
    let mut out_lines = CellArray::new();
    let mut visited: HashSet<i64> = HashSet::new();
    let mut pt_map: HashMap<i64, i64> = HashMap::new();
    let mut num_loops: usize = 0;

    let mut starts: Vec<i64> = adj.keys().copied().collect();
    starts.sort();

    for start in starts {
        if visited.contains(&start) {
            continue;
        }

        let mut loop_ids: Vec<i64> = Vec::new();
        let mut cur: i64 = start;

        loop {
            if visited.contains(&cur) && !loop_ids.is_empty() {
                break;
            }
            visited.insert(cur);

            // Map point
            let mapped: i64 = *pt_map.entry(cur).or_insert_with(|| {
                let idx: i64 = out_points.len() as i64;
                out_points.push(input.points.get(cur as usize));
                idx
            });
            loop_ids.push(mapped);

            let nexts = match adj.get(&cur) {
                Some(v) => v,
                None => break,
            };
            match nexts.iter().find(|&&n| !visited.contains(&n)) {
                Some(&n) => cur = n,
                None => {
                    // Close loop back to start if possible
                    if nexts.contains(&start) && loop_ids.len() >= 3 {
                        loop_ids.push(loop_ids[0]);
                    }
                    break;
                }
            }
        }

        if loop_ids.len() >= 3 {
            out_lines.push_cell(&loop_ids);
            num_loops += 1;
        }
    }

    let mut pd = PolyData::new();
    pd.points = out_points;
    pd.lines = out_lines;
    (pd, num_loops)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn single_triangle_one_loop() {
        let mut pd = PolyData::new();
        pd.points.push([0.0, 0.0, 0.0]);
        pd.points.push([1.0, 0.0, 0.0]);
        pd.points.push([0.5, 1.0, 0.0]);
        pd.polys.push_cell(&[0, 1, 2]);

        let (result, count) = extract_boundary_loops(&pd);
        assert_eq!(count, 1);
        assert_eq!(result.lines.num_cells(), 1);
        // The loop should contain 3 unique boundary vertices
        assert_eq!(result.points.len(), 3);
    }

    #[test]
    fn closed_tetrahedron_no_loops() {
        let mut pd = PolyData::new();
        pd.points.push([0.0, 0.0, 0.0]);
        pd.points.push([1.0, 0.0, 0.0]);
        pd.points.push([0.5, 1.0, 0.0]);
        pd.points.push([0.5, 0.5, 1.0]);
        pd.polys.push_cell(&[0, 1, 2]);
        pd.polys.push_cell(&[0, 1, 3]);
        pd.polys.push_cell(&[1, 2, 3]);
        pd.polys.push_cell(&[0, 2, 3]);

        let (_, count) = extract_boundary_loops(&pd);
        assert_eq!(count, 0);
    }

    #[test]
    fn two_separate_triangles_two_loops() {
        let mut pd = PolyData::new();
        // First triangle
        pd.points.push([0.0, 0.0, 0.0]);
        pd.points.push([1.0, 0.0, 0.0]);
        pd.points.push([0.5, 1.0, 0.0]);
        // Second triangle (separate)
        pd.points.push([3.0, 0.0, 0.0]);
        pd.points.push([4.0, 0.0, 0.0]);
        pd.points.push([3.5, 1.0, 0.0]);
        pd.polys.push_cell(&[0, 1, 2]);
        pd.polys.push_cell(&[3, 4, 5]);

        let (_, count) = extract_boundary_loops(&pd);
        assert_eq!(count, 2);
    }
}
