use std::collections::HashMap;

use crate::data::{Points, UnstructuredGrid};

/// Clip an UnstructuredGrid by a plane defined by a point and normal.
///
/// Keeps cells that are entirely in the half-space where
/// `dot(p - origin, normal) >= 0`. Cells that are partially clipped
/// are removed entirely (no splitting).
pub fn clip_data_set(
    input: &UnstructuredGrid,
    origin: [f64; 3],
    normal: [f64; 3],
) -> UnstructuredGrid {
    let n_points = input.points.len();

    // Classify each point
    let dists: Vec<f64> = (0..n_points)
        .map(|i| {
            let p = input.points.get(i);
            (p[0] - origin[0]) * normal[0]
                + (p[1] - origin[1]) * normal[1]
                + (p[2] - origin[2]) * normal[2]
        })
        .collect();

    let mut point_map: HashMap<usize, usize> = HashMap::new();
    let mut out_points = Points::<f64>::new();
    let mut out = UnstructuredGrid::new();

    let n_cells = input.cells().num_cells();
    for ci in 0..n_cells {
        let pts = input.cell_points(ci);
        let ct = input.cell_type(ci);

        // Keep cell only if all vertices are inside
        let all_inside = pts.iter().all(|&id| dists[id as usize] >= 0.0);
        if !all_inside {
            continue;
        }

        // Remap point indices
        let remapped: Vec<i64> = pts
            .iter()
            .map(|&id| {
                let uid = id as usize;
                *point_map.entry(uid).or_insert_with(|| {
                    let idx = out_points.len();
                    out_points.push(input.points.get(uid));
                    idx
                }) as i64
            })
            .collect();

        out.push_cell(ct, &remapped);
    }

    out.points = out_points;
    out
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::types::CellType;

    #[test]
    fn clip_keeps_inside_cells() {
        let mut grid = UnstructuredGrid::new();
        // Two triangles: one at x>0, one at x<0
        grid.points.push([1.0, 0.0, 0.0]);
        grid.points.push([2.0, 0.0, 0.0]);
        grid.points.push([1.5, 1.0, 0.0]);
        grid.points.push([-1.0, 0.0, 0.0]);
        grid.points.push([-2.0, 0.0, 0.0]);
        grid.points.push([-1.5, 1.0, 0.0]);

        grid.push_cell(CellType::Triangle, &[0, 1, 2]);
        grid.push_cell(CellType::Triangle, &[3, 4, 5]);

        // Clip at x=0, keep x>0
        let result = clip_data_set(&grid, [0.0, 0.0, 0.0], [1.0, 0.0, 0.0]);
        assert_eq!(result.cells().num_cells(), 1);
        assert_eq!(result.points.len(), 3);
    }

    #[test]
    fn clip_removes_all() {
        let mut grid = UnstructuredGrid::new();
        grid.points.push([-1.0, 0.0, 0.0]);
        grid.points.push([-2.0, 0.0, 0.0]);
        grid.points.push([-1.5, 1.0, 0.0]);
        grid.push_cell(CellType::Triangle, &[0, 1, 2]);

        let result = clip_data_set(&grid, [0.0, 0.0, 0.0], [1.0, 0.0, 0.0]);
        assert_eq!(result.cells().num_cells(), 0);
    }

    #[test]
    fn clip_mixed_cells() {
        let mut grid = UnstructuredGrid::new();
        grid.points.push([0.0, 0.0, 0.0]);
        grid.points.push([1.0, 0.0, 0.0]);
        grid.points.push([0.5, 1.0, 0.0]);
        grid.points.push([0.5, 0.5, 1.0]);
        grid.push_cell(CellType::Tetra, &[0, 1, 2, 3]);

        // Clip at x=0.5, normal=[1,0,0] — vertex 0 is at boundary (dist=0.5), but
        // all vertices have x >= 0, so the tetra should be kept
        let result = clip_data_set(&grid, [-0.1, 0.0, 0.0], [1.0, 0.0, 0.0]);
        assert_eq!(result.cells().num_cells(), 1);
    }
}
