use std::collections::HashMap;

use crate::data::{CellArray, Points, PolyData};

/// Extract silhouette edges from a triangle mesh as seen from a given viewpoint.
///
/// An edge is a silhouette edge if it is shared by exactly two triangles where
/// one is front-facing and the other is back-facing with respect to the viewpoint.
/// Boundary edges (used by only one face) that are front-facing are also included.
///
/// Returns a PolyData with line cells representing the silhouette edges.
pub fn extract_silhouette(input: &PolyData, viewpoint: [f64; 3]) -> PolyData {
    let num_faces: usize = input.polys.num_cells();
    if num_faces == 0 {
        return PolyData::default();
    }

    // Compute face normals and determine front/back facing.
    let mut face_facing: Vec<bool> = Vec::with_capacity(num_faces); // true = front-facing

    for cell in input.polys.iter() {
        if cell.len() < 3 {
            face_facing.push(false);
            continue;
        }
        let p0 = input.points.get(cell[0] as usize);
        let p1 = input.points.get(cell[1] as usize);
        let p2 = input.points.get(cell[2] as usize);

        // Face normal via cross product.
        let e1 = [p1[0] - p0[0], p1[1] - p0[1], p1[2] - p0[2]];
        let e2 = [p2[0] - p0[0], p2[1] - p0[1], p2[2] - p0[2]];
        let nx: f64 = e1[1] * e2[2] - e1[2] * e2[1];
        let ny: f64 = e1[2] * e2[0] - e1[0] * e2[2];
        let nz: f64 = e1[0] * e2[1] - e1[1] * e2[0];

        // View direction from face centroid to viewpoint.
        let cx: f64 = (p0[0] + p1[0] + p2[0]) / 3.0;
        let cy: f64 = (p0[1] + p1[1] + p2[1]) / 3.0;
        let cz: f64 = (p0[2] + p1[2] + p2[2]) / 3.0;
        let vx: f64 = viewpoint[0] - cx;
        let vy: f64 = viewpoint[1] - cy;
        let vz: f64 = viewpoint[2] - cz;

        let dot: f64 = nx * vx + ny * vy + nz * vz;
        face_facing.push(dot > 0.0);
    }

    // Build edge-to-face adjacency. Key: (min_vertex, max_vertex).
    let mut edge_faces: HashMap<(i64, i64), Vec<usize>> = HashMap::new();
    for (fi, cell) in input.polys.iter().enumerate() {
        let n: usize = cell.len();
        for k in 0..n {
            let a: i64 = cell[k];
            let b: i64 = cell[(k + 1) % n];
            let key = if a < b { (a, b) } else { (b, a) };
            edge_faces.entry(key).or_default().push(fi);
        }
    }

    // Collect silhouette edges.
    let mut out_points: Points<f64> = Points::new();
    let mut out_lines: CellArray = CellArray::new();
    let mut point_map: HashMap<i64, i64> = HashMap::new();

    let mut get_or_insert = |vid: i64, pts: &mut Points<f64>, map: &mut HashMap<i64, i64>| -> i64 {
        if let Some(&id) = map.get(&vid) {
            return id;
        }
        let id: i64 = pts.len() as i64;
        pts.push(input.points.get(vid as usize));
        map.insert(vid, id);
        id
    };

    for (&(a, b), faces) in &edge_faces {
        let is_silhouette: bool = if faces.len() == 1 {
            // Boundary edge: include if the single face is front-facing.
            face_facing[faces[0]]
        } else if faces.len() == 2 {
            // Silhouette: one front-facing and one back-facing.
            face_facing[faces[0]] != face_facing[faces[1]]
        } else {
            false
        };

        if is_silhouette {
            let id_a: i64 = get_or_insert(a, &mut out_points, &mut point_map);
            let id_b: i64 = get_or_insert(b, &mut out_points, &mut point_map);
            out_lines.push_cell(&[id_a, id_b]);
        }
    }

    let mut pd = PolyData::new();
    pd.points = out_points;
    pd.lines = out_lines;
    pd
}

#[cfg(test)]
mod tests {
    use super::*;

    fn make_double_triangle() -> PolyData {
        // Two triangles sharing an edge (a "bowtie" seen from above).
        // Triangle 0: (0,0,0), (1,0,0), (0.5,1,0) — in XY plane, normal +Z
        // Triangle 1: (1,0,0), (0.5,1,0), (0.5,0.5,-1) — tilted away
        let mut points: Points<f64> = Points::new();
        points.push([0.0, 0.0, 0.0]);
        points.push([1.0, 0.0, 0.0]);
        points.push([0.5, 1.0, 0.0]);
        points.push([0.5, 0.5, -1.0]);
        let mut polys = CellArray::new();
        polys.push_cell(&[0, 1, 2]);
        polys.push_cell(&[1, 2, 3]);
        let mut pd = PolyData::new();
        pd.points = points;
        pd.polys = polys;
        pd
    }

    #[test]
    fn silhouette_from_above() {
        let mesh = make_double_triangle();
        // Viewpoint far above: first triangle front-facing, second back-facing.
        let result = extract_silhouette(&mesh, [0.5, 0.5, 10.0]);
        // Should have silhouette edges. The shared edge (1,2) should be silhouette
        // since triangle 0 faces up and triangle 1 faces away.
        assert!(result.lines.num_cells() > 0, "Expected silhouette edges");
    }

    #[test]
    fn single_triangle_boundary() {
        let mut points: Points<f64> = Points::new();
        points.push([0.0, 0.0, 0.0]);
        points.push([1.0, 0.0, 0.0]);
        points.push([0.0, 1.0, 0.0]);
        let mut polys = CellArray::new();
        polys.push_cell(&[0, 1, 2]);
        let mut mesh = PolyData::new();
        mesh.points = points;
        mesh.polys = polys;
        // Viewpoint above: all 3 boundary edges should be silhouette (front-facing).
        let result = extract_silhouette(&mesh, [0.0, 0.0, 10.0]);
        assert_eq!(result.lines.num_cells(), 3, "3 boundary edges");
    }

    #[test]
    fn empty_mesh() {
        let mesh = PolyData::default();
        let result = extract_silhouette(&mesh, [0.0, 0.0, 1.0]);
        assert_eq!(result.lines.num_cells(), 0);
    }
}
