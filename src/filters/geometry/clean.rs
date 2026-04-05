use std::collections::HashMap;

use crate::data::{CellArray, Points, PolyData};

/// Parameters for cleaning PolyData.
pub struct CleanParams {
    /// Tolerance for merging nearby points. Points within this distance are merged.
    pub tolerance: f64,
    /// If true, merge duplicate/nearby points.
    pub merge_points: bool,
    /// If true, remove degenerate cells (lines with <2 points, polys with <3 points).
    pub remove_degenerate: bool,
}

impl Default for CleanParams {
    fn default() -> Self {
        Self {
            tolerance: 1e-6,
            merge_points: true,
            remove_degenerate: true,
        }
    }
}

/// Clean a PolyData by merging duplicate points and removing degenerate cells.
pub fn clean(input: &PolyData, params: &CleanParams) -> PolyData {
    let (new_points, point_map) = if params.merge_points {
        merge_points(&input.points, params.tolerance)
    } else {
        // Identity mapping
        let map: Vec<usize> = (0..input.points.len()).collect();
        let pts = input.points.clone();
        (pts, map)
    };

    let mut output = PolyData::new();
    output.points = new_points;

    // Remap cells
    output.verts = remap_cells(&input.verts, &point_map, params.remove_degenerate, 1);
    output.lines = remap_cells(&input.lines, &point_map, params.remove_degenerate, 2);
    output.polys = remap_cells(&input.polys, &point_map, params.remove_degenerate, 3);
    output.strips = remap_cells(&input.strips, &point_map, params.remove_degenerate, 3);

    output
}

/// Merge points that are within `tolerance` distance of each other.
/// Returns (new points, mapping from old index → new index).
fn merge_points(points: &Points<f64>, tolerance: f64) -> (Points<f64>, Vec<usize>) {
    let tol2 = tolerance * tolerance;
    let n = points.len();
    let mut point_map = vec![0usize; n];
    let mut new_points = Points::new();

    // Simple O(n*m) algorithm. For large datasets, a spatial hash would be better.
    // We use a hash map with quantized coordinates for O(n) average case.
    let inv_tol = if tolerance > 0.0 { 1.0 / tolerance } else { 1e12 };

    let mut grid: HashMap<(i64, i64, i64), Vec<usize>> = HashMap::new();

    for (i, pm) in point_map.iter_mut().enumerate() {
        let p = points.get(i);
        let key = (
            (p[0] * inv_tol).round() as i64,
            (p[1] * inv_tol).round() as i64,
            (p[2] * inv_tol).round() as i64,
        );

        let mut found = None;

        // Check neighboring cells (3x3x3 neighborhood for robustness)
        'search: for dx in -1..=1 {
            for dy in -1..=1 {
                for dz in -1..=1 {
                    let nkey = (key.0 + dx, key.1 + dy, key.2 + dz);
                    if let Some(indices) = grid.get(&nkey) {
                        for &idx in indices {
                            let q: [f64; 3] = new_points.get(idx);
                            let d2 = (p[0] - q[0]) * (p[0] - q[0])
                                + (p[1] - q[1]) * (p[1] - q[1])
                                + (p[2] - q[2]) * (p[2] - q[2]);
                            if d2 <= tol2 {
                                found = Some(idx);
                                break 'search;
                            }
                        }
                    }
                }
            }
        }

        if let Some(existing) = found {
            *pm = existing;
        } else {
            let new_idx = new_points.len();
            new_points.push(p);
            *pm = new_idx;
            grid.entry(key).or_default().push(new_idx);
        }
    }

    (new_points, point_map)
}

/// Remap cell point indices and optionally remove degenerate cells.
fn remap_cells(
    cells: &CellArray,
    point_map: &[usize],
    remove_degenerate: bool,
    min_unique: usize,
) -> CellArray {
    let mut out = CellArray::new();

    for cell in cells.iter() {
        let remapped: Vec<i64> = cell
            .iter()
            .map(|&id| point_map[id as usize] as i64)
            .collect();

        if remove_degenerate {
            // Remove consecutive duplicates
            let mut deduped: Vec<i64> = Vec::with_capacity(remapped.len());
            for &id in &remapped {
                if deduped.last() != Some(&id) {
                    deduped.push(id);
                }
            }
            // Also check if first == last for polygons
            if deduped.len() > 1 && deduped.first() == deduped.last() {
                deduped.pop();
            }

            if deduped.len() >= min_unique {
                out.push_cell(&deduped);
            }
        } else {
            out.push_cell(&remapped);
        }
    }

    out
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn merge_duplicate_points() {
        let mut pd = PolyData::new();
        // Two triangles with duplicate points at indices 3,4,5 = copies of 0,1,2
        pd.points = Points::from_vec(vec![
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 0.0], // duplicate of 0
            [1.0, 0.0, 0.0], // duplicate of 1
            [0.0, 1.0, 0.0], // duplicate of 2
        ]);
        pd.polys.push_cell(&[0, 1, 2]);
        pd.polys.push_cell(&[3, 4, 5]);

        let result = clean(&pd, &CleanParams::default());
        assert_eq!(result.points.len(), 3);
        assert_eq!(result.polys.num_cells(), 2);
        // Both triangles should reference the same 3 points
        assert_eq!(result.polys.cell(0), &[0, 1, 2]);
        assert_eq!(result.polys.cell(1), &[0, 1, 2]);
    }

    #[test]
    fn remove_degenerate_cell() {
        let mut pd = PolyData::new();
        pd.points = Points::from_vec(vec![
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
        ]);
        // Degenerate triangle: only 2 unique points after collapse
        pd.polys.push_cell(&[0, 1, 0]);

        let result = clean(&pd, &CleanParams::default());
        assert_eq!(result.polys.num_cells(), 0); // degenerate removed
    }

    #[test]
    fn no_merge_mode() {
        let mut pd = PolyData::new();
        pd.points = Points::from_vec(vec![
            [0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0], // duplicate
        ]);
        pd.polys.push_cell(&[0, 1, 0]); // technically valid if not merging

        let result = clean(
            &pd,
            &CleanParams {
                merge_points: false,
                remove_degenerate: false,
                ..Default::default()
            },
        );
        assert_eq!(result.points.len(), 2); // no merging
        assert_eq!(result.polys.num_cells(), 1); // no removal
    }
}
