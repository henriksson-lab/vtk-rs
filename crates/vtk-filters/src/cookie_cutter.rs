//! Cut a mesh surface by a polygon outline (cookie cutter).
//!
//! Removes cells that fall outside the cutting polygon, producing
//! a mesh trimmed to the polygon boundary.

use vtk_data::{CellArray, Points, PolyData};

/// Cut a mesh by a 2D polygon outline, keeping only cells inside.
///
/// The polygon is defined in the XY plane. Cells whose centroids lie
/// inside the polygon are kept; others are discarded.
pub fn cookie_cut(mesh: &PolyData, polygon: &[[f64; 2]]) -> PolyData {
    if polygon.len() < 3 || mesh.polys.num_cells() == 0 {
        return mesh.clone();
    }

    let mut new_points = Points::<f64>::new();
    let mut new_polys = CellArray::new();
    let mut point_map: std::collections::HashMap<usize, usize> = std::collections::HashMap::new();

    for cell in mesh.polys.iter() {
        if cell.is_empty() { continue; }
        let mut cx = 0.0;
        let mut cy = 0.0;
        for &pid in cell {
            let p = mesh.points.get(pid as usize);
            cx += p[0];
            cy += p[1];
        }
        let n = cell.len() as f64;
        if !point_in_polygon(cx / n, cy / n, polygon) { continue; }

        let mut new_ids = Vec::with_capacity(cell.len());
        for &pid in cell {
            let old = pid as usize;
            let new_idx = *point_map.entry(old).or_insert_with(|| {
                let idx = new_points.len();
                new_points.push(mesh.points.get(old));
                idx
            });
            new_ids.push(new_idx as i64);
        }
        new_polys.push_cell(&new_ids);
    }

    let mut result = PolyData::new();
    result.points = new_points;
    result.polys = new_polys;
    result
}

/// Cut by a circular region (keep cells inside circle in XY plane).
pub fn cookie_cut_circle(mesh: &PolyData, center: [f64; 2], radius: f64) -> PolyData {
    if mesh.polys.num_cells() == 0 { return mesh.clone(); }

    let r2 = radius * radius;
    let mut new_points = Points::<f64>::new();
    let mut new_polys = CellArray::new();
    let mut point_map: std::collections::HashMap<usize, usize> = std::collections::HashMap::new();

    for cell in mesh.polys.iter() {
        if cell.is_empty() { continue; }
        let mut cx = 0.0;
        let mut cy = 0.0;
        for &pid in cell {
            let p = mesh.points.get(pid as usize);
            cx += p[0];
            cy += p[1];
        }
        let n = cell.len() as f64;
        cx /= n;
        cy /= n;
        let dx = cx - center[0];
        let dy = cy - center[1];
        if dx * dx + dy * dy > r2 { continue; }

        let mut new_ids = Vec::with_capacity(cell.len());
        for &pid in cell {
            let old = pid as usize;
            let new_idx = *point_map.entry(old).or_insert_with(|| {
                let idx = new_points.len();
                new_points.push(mesh.points.get(old));
                idx
            });
            new_ids.push(new_idx as i64);
        }
        new_polys.push_cell(&new_ids);
    }

    let mut result = PolyData::new();
    result.points = new_points;
    result.polys = new_polys;
    result
}

fn point_in_polygon(px: f64, py: f64, polygon: &[[f64; 2]]) -> bool {
    let n = polygon.len();
    let mut inside = false;
    let mut j = n - 1;
    for i in 0..n {
        let yi = polygon[i][1];
        let yj = polygon[j][1];
        if (yi > py) != (yj > py) {
            let x_int = polygon[i][0] + (py - yi) / (yj - yi) * (polygon[j][0] - polygon[i][0]);
            if px < x_int { inside = !inside; }
        }
        j = i;
    }
    inside
}

#[cfg(test)]
mod tests {
    use super::*;
    use vtk_data::PolyData;

    fn make_grid_mesh() -> PolyData {
        // 4x4 grid of triangles in XY plane
        let mut pts = Vec::new();
        for y in 0..5 {
            for x in 0..5 {
                pts.push([x as f64, y as f64, 0.0]);
            }
        }
        let mut tris = Vec::new();
        for y in 0..4 {
            for x in 0..4 {
                let bl = y * 5 + x;
                tris.push([bl, bl + 1, bl + 6]);
                tris.push([bl, bl + 6, bl + 5]);
            }
        }
        PolyData::from_triangles(pts, tris)
    }

    #[test]
    fn cut_by_polygon() {
        let mesh = make_grid_mesh();
        let polygon = [[0.5, 0.5], [2.5, 0.5], [2.5, 2.5], [0.5, 2.5]];
        let result = cookie_cut(&mesh, &polygon);
        assert!(result.polys.num_cells() > 0);
        assert!(result.polys.num_cells() < mesh.polys.num_cells());
    }

    #[test]
    fn cut_by_circle() {
        let mesh = make_grid_mesh();
        let result = cookie_cut_circle(&mesh, [2.0, 2.0], 1.5);
        assert!(result.polys.num_cells() > 0);
        assert!(result.polys.num_cells() < mesh.polys.num_cells());
    }

    #[test]
    fn polygon_outside_mesh() {
        let mesh = make_grid_mesh();
        let polygon = [[100.0, 100.0], [200.0, 100.0], [200.0, 200.0]];
        let result = cookie_cut(&mesh, &polygon);
        assert_eq!(result.polys.num_cells(), 0);
    }

    #[test]
    fn polygon_covers_all() {
        let mesh = make_grid_mesh();
        let polygon = [[-1.0, -1.0], [10.0, -1.0], [10.0, 10.0], [-1.0, 10.0]];
        let result = cookie_cut(&mesh, &polygon);
        assert_eq!(result.polys.num_cells(), mesh.polys.num_cells());
    }
}
