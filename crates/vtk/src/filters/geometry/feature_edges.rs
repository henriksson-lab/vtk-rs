use std::collections::HashMap;

use crate::data::{CellArray, Points, PolyData};

/// Type of edge detected by the feature edges filter.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum EdgeType {
    /// Edge used by only one polygon (mesh boundary).
    Boundary,
    /// Edge where the dihedral angle exceeds the feature angle.
    Feature,
    /// Edge shared by exactly two polygons (smooth interior).
    Manifold,
    /// Edge shared by more than two polygons (non-manifold).
    NonManifold,
}

/// Parameters for feature edge extraction.
pub struct FeatureEdgesParams {
    /// Angle threshold in degrees. Edges with dihedral angle greater than
    /// this are classified as feature edges. Default: 30.0
    pub feature_angle: f64,
    /// Include boundary edges in output. Default: true
    pub boundary_edges: bool,
    /// Include feature edges in output. Default: true
    pub feature_edges: bool,
    /// Include manifold (interior smooth) edges. Default: false
    pub manifold_edges: bool,
    /// Include non-manifold edges. Default: true
    pub non_manifold_edges: bool,
}

impl Default for FeatureEdgesParams {
    fn default() -> Self {
        Self {
            feature_angle: 30.0,
            boundary_edges: true,
            feature_edges: true,
            manifold_edges: false,
            non_manifold_edges: true,
        }
    }
}

/// Extract feature, boundary, manifold, and non-manifold edges from a PolyData.
///
/// Returns a PolyData containing line cells for the selected edge types.
pub fn feature_edges(input: &PolyData, params: &FeatureEdgesParams) -> PolyData {
    // Build edge -> face list mapping
    // Key: sorted (min, max) point id pair
    let mut edge_faces: HashMap<(i64, i64), Vec<usize>> = HashMap::new();
    let mut face_normals: Vec<[f64; 3]> = Vec::new();

    for (face_idx, cell) in input.polys.iter().enumerate() {
        // Compute face normal
        let normal = polygon_normal(input, cell);
        face_normals.push(normal);

        let n = cell.len();
        for i in 0..n {
            let a = cell[i];
            let b = cell[(i + 1) % n];
            let key = if a < b { (a, b) } else { (b, a) };
            edge_faces.entry(key).or_default().push(face_idx);
        }
    }

    let cos_threshold = (params.feature_angle.to_radians()).cos();

    let mut point_map: HashMap<i64, usize> = HashMap::new();
    let mut out_points = Points::<f64>::new();
    let mut out_lines = CellArray::new();

    let map_point = |id: i64, pts: &Points<f64>, out_pts: &mut Points<f64>, pm: &mut HashMap<i64, usize>| -> i64 {
        *pm.entry(id).or_insert_with(|| {
            let idx = out_pts.len();
            out_pts.push(pts.get(id as usize));
            idx
        }) as i64
    };

    for (&(a, b), faces) in &edge_faces {
        let edge_type = if faces.len() == 1 {
            EdgeType::Boundary
        } else if faces.len() == 2 {
            let n1 = face_normals[faces[0]];
            let n2 = face_normals[faces[1]];
            let dot = n1[0] * n2[0] + n1[1] * n2[1] + n1[2] * n2[2];
            if dot < cos_threshold {
                EdgeType::Feature
            } else {
                EdgeType::Manifold
            }
        } else {
            EdgeType::NonManifold
        };

        let include = match edge_type {
            EdgeType::Boundary => params.boundary_edges,
            EdgeType::Feature => params.feature_edges,
            EdgeType::Manifold => params.manifold_edges,
            EdgeType::NonManifold => params.non_manifold_edges,
        };

        if include {
            let ma = map_point(a, &input.points, &mut out_points, &mut point_map);
            let mb = map_point(b, &input.points, &mut out_points, &mut point_map);
            out_lines.push_cell(&[ma, mb]);
        }
    }

    let mut pd = PolyData::new();
    pd.points = out_points;
    pd.lines = out_lines;
    pd
}

fn polygon_normal(input: &PolyData, cell: &[i64]) -> [f64; 3] {
    // Newell's method
    let mut nx = 0.0;
    let mut ny = 0.0;
    let mut nz = 0.0;
    let n = cell.len();
    for i in 0..n {
        let p = input.points.get(cell[i] as usize);
        let q = input.points.get(cell[(i + 1) % n] as usize);
        nx += (p[1] - q[1]) * (p[2] + q[2]);
        ny += (p[2] - q[2]) * (p[0] + q[0]);
        nz += (p[0] - q[0]) * (p[1] + q[1]);
    }
    let len = (nx * nx + ny * ny + nz * nz).sqrt();
    if len > 1e-10 {
        [nx / len, ny / len, nz / len]
    } else {
        [0.0, 0.0, 1.0]
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn boundary_edges_of_single_triangle() {
        let pd = PolyData::from_triangles(
            vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.5, 1.0, 0.0]],
            vec![[0, 1, 2]],
        );
        let result = feature_edges(&pd, &FeatureEdgesParams::default());
        // Single triangle: all 3 edges are boundary
        assert_eq!(result.lines.num_cells(), 3);
    }

    #[test]
    fn shared_edge_not_boundary() {
        let pd = PolyData::from_triangles(
            vec![
                [0.0, 0.0, 0.0],
                [1.0, 0.0, 0.0],
                [0.5, 1.0, 0.0],
                [0.5, -1.0, 0.0],
            ],
            vec![[0, 1, 2], [0, 3, 1]],
        );
        let params = FeatureEdgesParams {
            boundary_edges: true,
            feature_edges: false,
            manifold_edges: false,
            non_manifold_edges: false,
            ..Default::default()
        };
        let result = feature_edges(&pd, &params);
        // 2 triangles share edge (0,1), so 4 boundary edges remain
        assert_eq!(result.lines.num_cells(), 4);
    }

    #[test]
    fn feature_angle_detection() {
        // Two triangles meeting at 90 degrees
        let pd = PolyData::from_triangles(
            vec![
                [0.0, 0.0, 0.0],
                [1.0, 0.0, 0.0],
                [0.5, 1.0, 0.0],
                [0.5, 0.0, -1.0],
            ],
            vec![[0, 1, 2], [0, 1, 3]],
        );
        let params = FeatureEdgesParams {
            feature_angle: 45.0, // 90 deg > 45 deg threshold
            boundary_edges: false,
            feature_edges: true,
            manifold_edges: false,
            non_manifold_edges: false,
        };
        let result = feature_edges(&pd, &params);
        // The shared edge should be a feature edge
        assert_eq!(result.lines.num_cells(), 1);
    }
}
