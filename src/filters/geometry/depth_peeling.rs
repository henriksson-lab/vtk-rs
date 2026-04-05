//! Depth-sorted polygon layers for order-independent transparency.
//!
//! Splits a PolyData into front-to-back layers by depth sorting relative
//! to a view direction. Each layer can be rendered in order for correct
//! transparency blending. This is a CPU-side preparation for the depth
//! peeling rendering technique.

use crate::data::{AnyDataArray, CellArray, DataArray, Points, PolyData};

/// Split a mesh into depth-peeled layers by sorting faces front-to-back
/// relative to a view direction.
///
/// Returns a vector of PolyData, one per layer, sorted front-to-back.
/// `num_layers` controls maximum number of peeling passes.
pub fn depth_peel_layers(
    mesh: &PolyData,
    view_direction: [f64; 3],
    num_layers: usize,
) -> Vec<PolyData> {
    if mesh.polys.num_cells() == 0 || num_layers == 0 {
        return Vec::new();
    }

    // Compute depth (dot product of centroid with view direction) per cell
    let mut cell_depths: Vec<(usize, f64)> = Vec::new();
    for (ci, cell) in mesh.polys.iter().enumerate() {
        if cell.is_empty() { continue; }
        let mut cx = 0.0;
        let mut cy = 0.0;
        let mut cz = 0.0;
        for &pid in cell {
            let p = mesh.points.get(pid as usize);
            cx += p[0];
            cy += p[1];
            cz += p[2];
        }
        let n = cell.len() as f64;
        cx /= n;
        cy /= n;
        cz /= n;
        let depth = cx * view_direction[0] + cy * view_direction[1] + cz * view_direction[2];
        cell_depths.push((ci, depth));
    }

    // Sort by depth (front = smallest dot product for typical -Z view)
    cell_depths.sort_by(|a, b| a.1.partial_cmp(&b.1).unwrap_or(std::cmp::Ordering::Equal));

    // Split into layers of roughly equal size
    let cells_per_layer = (cell_depths.len() + num_layers - 1) / num_layers;

    let mut layers = Vec::new();
    for chunk in cell_depths.chunks(cells_per_layer) {
        let layer = extract_cells_to_poly_data(mesh, chunk);
        if layer.polys.num_cells() > 0 {
            layers.push(layer);
        }
    }

    layers
}

/// Sort all faces of a mesh by depth for painter's algorithm rendering.
///
/// Returns a new PolyData with faces reordered front-to-back relative to
/// the view direction, with a "Depth" cell data array.
pub fn depth_sort_mesh(
    mesh: &PolyData,
    view_direction: [f64; 3],
) -> PolyData {
    if mesh.polys.num_cells() == 0 {
        return mesh.clone();
    }

    let mut cell_depths: Vec<(usize, f64)> = Vec::new();
    for (ci, cell) in mesh.polys.iter().enumerate() {
        if cell.is_empty() { continue; }
        let mut cx = 0.0;
        let mut cy = 0.0;
        let mut cz = 0.0;
        for &pid in cell {
            let p = mesh.points.get(pid as usize);
            cx += p[0];
            cy += p[1];
            cz += p[2];
        }
        let n = cell.len() as f64;
        let depth = (cx / n) * view_direction[0]
            + (cy / n) * view_direction[1]
            + (cz / n) * view_direction[2];
        cell_depths.push((ci, depth));
    }

    cell_depths.sort_by(|a, b| a.1.partial_cmp(&b.1).unwrap_or(std::cmp::Ordering::Equal));

    let mut result = extract_cells_to_poly_data(mesh, &cell_depths);

    // Add depth as cell data
    let depths: Vec<f64> = cell_depths.iter().map(|(_, d)| *d).collect();
    result.cell_data_mut().add_array(AnyDataArray::F64(
        DataArray::from_vec("Depth", depths, 1),
    ));

    result
}

fn extract_cells_to_poly_data(mesh: &PolyData, cells: &[(usize, f64)]) -> PolyData {
    // Map old point indices to new ones
    let mut point_map: std::collections::HashMap<usize, usize> = std::collections::HashMap::new();
    let mut new_points = Points::<f64>::new();
    let mut new_polys = CellArray::new();

    let all_cells: Vec<Vec<i64>> = mesh.polys.iter().map(|c| c.to_vec()).collect();

    for &(ci, _) in cells {
        if ci >= all_cells.len() { continue; }
        let cell = &all_cells[ci];
        let mut new_ids = Vec::with_capacity(cell.len());
        for &old_id in cell {
            let old_idx = old_id as usize;
            let new_idx = *point_map.entry(old_idx).or_insert_with(|| {
                let idx = new_points.len();
                new_points.push(mesh.points.get(old_idx));
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

#[cfg(test)]
mod tests {
    use super::*;

    fn make_two_planes() -> PolyData {
        // Two quads at different Z depths
        PolyData::from_triangles(
            vec![
                // Front triangle (z=0)
                [0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.5, 1.0, 0.0],
                // Back triangle (z=1)
                [0.0, 0.0, 1.0], [1.0, 0.0, 1.0], [0.5, 1.0, 1.0],
            ],
            vec![[0, 1, 2], [3, 4, 5]],
        )
    }

    #[test]
    fn depth_peel_two_layers() {
        let mesh = make_two_planes();
        let layers = depth_peel_layers(&mesh, [0.0, 0.0, 1.0], 2);
        assert_eq!(layers.len(), 2);
        assert_eq!(layers[0].polys.num_cells(), 1);
        assert_eq!(layers[1].polys.num_cells(), 1);
    }

    #[test]
    fn depth_sort() {
        let mesh = make_two_planes();
        let sorted = depth_sort_mesh(&mesh, [0.0, 0.0, 1.0]);
        assert_eq!(sorted.polys.num_cells(), 2);
        assert!(sorted.cell_data().get_array("Depth").is_some());

        // First cell should have smaller depth
        let depth_arr = sorted.cell_data().get_array("Depth").unwrap();
        let mut d0 = [0.0f64];
        let mut d1 = [0.0f64];
        depth_arr.tuple_as_f64(0, &mut d0);
        depth_arr.tuple_as_f64(1, &mut d1);
        assert!(d0[0] <= d1[0]);
    }

    #[test]
    fn single_layer() {
        let mesh = make_two_planes();
        let layers = depth_peel_layers(&mesh, [0.0, 0.0, 1.0], 1);
        assert_eq!(layers.len(), 1);
        assert_eq!(layers[0].polys.num_cells(), 2);
    }

    #[test]
    fn empty_mesh() {
        let mesh = PolyData::new();
        let layers = depth_peel_layers(&mesh, [0.0, 0.0, 1.0], 4);
        assert!(layers.is_empty());
    }
}
