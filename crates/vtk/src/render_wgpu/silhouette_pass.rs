use crate::data::PolyData;
use crate::render::SilhouetteConfig;

use crate::render_wgpu::mesh::Vertex;

/// Extract silhouette edges from a PolyData mesh viewed from a camera position.
///
/// A silhouette edge is one where one adjacent face is front-facing and the other
/// is back-facing relative to the view direction. Also includes boundary edges.
/// Returns line-list vertices and indices for rendering with the wireframe pipeline.
pub fn extract_silhouette_edges(
    pd: &PolyData,
    camera_pos: [f64; 3],
    config: &SilhouetteConfig,
) -> (Vec<Vertex>, Vec<u32>) {
    let color = config.color;
    let n_pts = pd.points.len();
    if n_pts == 0 {
        return (Vec::new(), Vec::new());
    }

    // Build edge-to-face adjacency
    let mut edge_faces: std::collections::HashMap<(usize, usize), Vec<usize>> =
        std::collections::HashMap::new();

    // Compute face normals
    let mut face_normals: Vec<[f64; 3]> = Vec::new();
    let mut face_centers: Vec<[f64; 3]> = Vec::new();

    for cell in pd.polys.iter() {
        let fi = face_normals.len();
        if cell.len() < 3 {
            face_normals.push([0.0, 0.0, 1.0]);
            face_centers.push([0.0, 0.0, 0.0]);
            continue;
        }

        // Face normal via cross product
        let p0 = pd.points.get(cell[0] as usize);
        let p1 = pd.points.get(cell[1] as usize);
        let p2 = pd.points.get(cell[2] as usize);
        let e1 = [p1[0] - p0[0], p1[1] - p0[1], p1[2] - p0[2]];
        let e2 = [p2[0] - p0[0], p2[1] - p0[1], p2[2] - p0[2]];
        let nx = e1[1] * e2[2] - e1[2] * e2[1];
        let ny = e1[2] * e2[0] - e1[0] * e2[2];
        let nz = e1[0] * e2[1] - e1[1] * e2[0];
        let len = (nx * nx + ny * ny + nz * nz).sqrt();
        if len > 1e-12 {
            face_normals.push([nx / len, ny / len, nz / len]);
        } else {
            face_normals.push([0.0, 0.0, 1.0]);
        }

        // Center
        let mut cx = 0.0;
        let mut cy = 0.0;
        let mut cz = 0.0;
        for &pid in cell {
            let p = pd.points.get(pid as usize);
            cx += p[0];
            cy += p[1];
            cz += p[2];
        }
        let nc = cell.len() as f64;
        face_centers.push([cx / nc, cy / nc, cz / nc]);

        // Register edges
        for i in 0..cell.len() {
            let a = cell[i] as usize;
            let b = cell[(i + 1) % cell.len()] as usize;
            let key = if a < b { (a, b) } else { (b, a) };
            edge_faces.entry(key).or_default().push(fi);
        }
    }

    // Determine front/back facing for each face
    let face_front: Vec<bool> = face_normals
        .iter()
        .zip(face_centers.iter())
        .map(|(normal, center)| {
            let view_dir = [
                camera_pos[0] - center[0],
                camera_pos[1] - center[1],
                camera_pos[2] - center[2],
            ];
            let dot = normal[0] * view_dir[0] + normal[1] * view_dir[1] + normal[2] * view_dir[2];
            dot > 0.0
        })
        .collect();

    // Extract silhouette edges
    let mut vertices: Vec<Vertex> = Vec::new();
    let mut indices: Vec<u32> = Vec::new();

    for (&(a, b), faces) in &edge_faces {
        let is_silhouette = if faces.len() == 1 {
            // Boundary edge — always silhouette
            face_front[faces[0]]
        } else if faces.len() == 2 {
            // Silhouette: one front, one back
            face_front[faces[0]] != face_front[faces[1]]
        } else {
            false
        };

        if is_silhouette {
            let base = vertices.len() as u32;
            let pa = pd.points.get(a);
            let pb = pd.points.get(b);
            vertices.push(Vertex {
                position: [pa[0] as f32, pa[1] as f32, pa[2] as f32],
                normal: [0.0, 0.0, 1.0],
                color,
            });
            vertices.push(Vertex {
                position: [pb[0] as f32, pb[1] as f32, pb[2] as f32],
                normal: [0.0, 0.0, 1.0],
                color,
            });
            indices.push(base);
            indices.push(base + 1);
        }
    }

    (vertices, indices)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn silhouette_of_triangle() {
        let pd = PolyData::from_triangles(
            vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]],
            vec![[0, 1, 2]],
        );
        let config = SilhouetteConfig {
            color: [0.0, 0.0, 0.0],
            enabled: true,
            ..Default::default()
        };
        let (verts, idxs) = extract_silhouette_edges(&pd, [0.5, 0.5, 5.0], &config);
        // Single triangle: all 3 edges are boundary edges, all silhouette
        assert_eq!(idxs.len(), 6); // 3 edges * 2 indices
        assert_eq!(verts.len(), 6); // 3 edges * 2 vertices
    }

    #[test]
    fn silhouette_of_two_coplanar_triangles() {
        // Two triangles sharing edge 1-2, same plane
        let pd = PolyData::from_triangles(
            vec![
                [0.0, 0.0, 0.0],
                [1.0, 0.0, 0.0],
                [0.5, 1.0, 0.0],
                [0.5, -1.0, 0.0],
            ],
            vec![[0, 1, 2], [0, 3, 1]],
        );
        let config = SilhouetteConfig {
            color: [0.0, 0.0, 0.0],
            enabled: true,
            ..Default::default()
        };
        let (_, idxs) = extract_silhouette_edges(&pd, [0.5, 0.0, 5.0], &config);
        // Shared edge 0-1 is NOT silhouette (both faces front-facing)
        // 4 boundary edges are silhouette
        assert_eq!(idxs.len(), 8); // 4 edges * 2
    }
}
