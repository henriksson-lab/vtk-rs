use std::collections::HashSet;

use crate::data::{CellArray, PolyData};

/// Extract all unique edges from a PolyData as line segments.
///
/// Each polygon edge appears exactly once in the output, regardless of
/// how many cells share it.
pub fn extract_edges(input: &PolyData) -> PolyData {
    let mut seen: HashSet<(i64, i64)> = HashSet::new();
    let mut out_lines = CellArray::new();

    let process_cell = |cell: &[i64], seen: &mut HashSet<(i64, i64)>, lines: &mut CellArray| {
        let n = cell.len();
        for i in 0..n {
            let a = cell[i];
            let b = cell[(i + 1) % n];
            let key = if a < b { (a, b) } else { (b, a) };
            if seen.insert(key) {
                lines.push_cell(&[a, b]);
            }
        }
    };

    for cell in input.polys.iter() {
        process_cell(cell, &mut seen, &mut out_lines);
    }
    for cell in input.strips.iter() {
        process_cell(cell, &mut seen, &mut out_lines);
    }
    // Line cells: edges are consecutive pairs
    for cell in input.lines.iter() {
        for i in 0..cell.len().saturating_sub(1) {
            let a = cell[i];
            let b = cell[i + 1];
            let key = if a < b { (a, b) } else { (b, a) };
            if seen.insert(key) {
                out_lines.push_cell(&[a, b]);
            }
        }
    }

    let mut pd = PolyData::new();
    pd.points = input.points.clone();
    pd.lines = out_lines;
    pd
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn edges_of_single_triangle() {
        let pd = PolyData::from_triangles(
            vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.5, 1.0, 0.0]],
            vec![[0, 1, 2]],
        );
        let result = extract_edges(&pd);
        assert_eq!(result.lines.num_cells(), 3);
    }

    #[test]
    fn shared_edge_appears_once() {
        let pd = PolyData::from_triangles(
            vec![
                [0.0, 0.0, 0.0],
                [1.0, 0.0, 0.0],
                [0.5, 1.0, 0.0],
                [0.5, -1.0, 0.0],
            ],
            vec![[0, 1, 2], [0, 3, 1]],
        );
        let result = extract_edges(&pd);
        // 2 triangles share edge (0,1): 3 + 3 - 1 = 5 unique edges
        assert_eq!(result.lines.num_cells(), 5);
    }

    #[test]
    fn edges_of_quad() {
        let mut pd = PolyData::new();
        pd.points.push([0.0, 0.0, 0.0]);
        pd.points.push([1.0, 0.0, 0.0]);
        pd.points.push([1.0, 1.0, 0.0]);
        pd.points.push([0.0, 1.0, 0.0]);
        pd.polys.push_cell(&[0, 1, 2, 3]);
        let result = extract_edges(&pd);
        assert_eq!(result.lines.num_cells(), 4);
    }
}
