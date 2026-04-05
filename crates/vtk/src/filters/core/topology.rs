use std::collections::{HashMap, HashSet};
use crate::data::PolyData;

/// Topology analysis results for a PolyData mesh.
#[derive(Debug, Clone)]
pub struct TopologyInfo {
    /// Number of vertices (0-cells).
    pub num_points: usize,
    /// Number of edges.
    pub num_edges: usize,
    /// Number of faces (polygons).
    pub num_faces: usize,
    /// Number of boundary edges (shared by exactly one face).
    pub num_boundary_edges: usize,
    /// Number of non-manifold edges (shared by more than two faces).
    pub num_non_manifold_edges: usize,
    /// Euler characteristic (V - E + F).
    pub euler_characteristic: i64,
    /// Number of connected components.
    pub num_components: usize,
    /// Whether the mesh is manifold (every edge shared by 1 or 2 faces).
    pub is_manifold: bool,
    /// Whether the mesh is closed (no boundary edges).
    pub is_closed: bool,
    /// Whether all faces are triangles.
    pub is_triangle_mesh: bool,
    /// Genus (for closed manifold: (2 - euler) / 2).
    pub genus: Option<i64>,
}

/// Analyze the topology of a PolyData mesh.
pub fn analyze_topology(pd: &PolyData) -> TopologyInfo {
    let n_pts = pd.points.len();
    let n_faces = pd.polys.num_cells();

    // Build edge-to-face adjacency
    let mut edge_count: HashMap<(usize, usize), usize> = HashMap::new();
    let mut is_all_tris = true;

    for cell in pd.polys.iter() {
        if cell.len() != 3 {
            is_all_tris = false;
        }
        let len = cell.len();
        for i in 0..len {
            let a = cell[i] as usize;
            let b = cell[(i + 1) % len] as usize;
            let edge = if a < b { (a, b) } else { (b, a) };
            *edge_count.entry(edge).or_insert(0) += 1;
        }
    }

    let n_edges = edge_count.len();
    let n_boundary = edge_count.values().filter(|&&c| c == 1).count();
    let n_non_manifold = edge_count.values().filter(|&&c| c > 2).count();
    let is_manifold = n_non_manifold == 0;
    let is_closed = n_boundary == 0 && n_faces > 0;
    let euler = n_pts as i64 - n_edges as i64 + n_faces as i64;

    let genus = if is_closed && is_manifold {
        Some((2 - euler) / 2)
    } else {
        None
    };

    let num_components = count_components(pd);

    TopologyInfo {
        num_points: n_pts,
        num_edges: n_edges,
        num_faces: n_faces,
        num_boundary_edges: n_boundary,
        num_non_manifold_edges: n_non_manifold,
        euler_characteristic: euler,
        num_components,
        is_manifold,
        is_closed,
        is_triangle_mesh: is_all_tris,
        genus,
    }
}

/// Count connected components using union-find.
fn count_components(pd: &PolyData) -> usize {
    let n = pd.points.len();
    if n == 0 { return 0; }

    let mut parent: Vec<usize> = (0..n).collect();

    fn find(parent: &mut [usize], x: usize) -> usize {
        let mut x = x;
        while parent[x] != x {
            parent[x] = parent[parent[x]];
            x = parent[x];
        }
        x
    }

    fn union(parent: &mut [usize], a: usize, b: usize) {
        let ra = find(parent, a);
        let rb = find(parent, b);
        if ra != rb {
            parent[ra] = rb;
        }
    }

    for cell in pd.polys.iter() {
        for i in 1..cell.len() {
            union(&mut parent, cell[0] as usize, cell[i] as usize);
        }
    }

    // Count unique roots among used points
    let mut used_points: HashSet<usize> = HashSet::new();
    for cell in pd.polys.iter() {
        for &pid in cell {
            used_points.insert(pid as usize);
        }
    }

    let roots: HashSet<usize> = used_points.iter().map(|&p| find(&mut parent, p)).collect();
    roots.len().max(if used_points.is_empty() { 0 } else { 1 })
}

/// Find boundary edges (edges with exactly one adjacent face).
/// Returns pairs of point indices.
pub fn boundary_edges(pd: &PolyData) -> Vec<(usize, usize)> {
    let mut edge_count: HashMap<(usize, usize), usize> = HashMap::new();
    for cell in pd.polys.iter() {
        let len = cell.len();
        for i in 0..len {
            let a = cell[i] as usize;
            let b = cell[(i + 1) % len] as usize;
            let edge = if a < b { (a, b) } else { (b, a) };
            *edge_count.entry(edge).or_insert(0) += 1;
        }
    }
    edge_count.into_iter()
        .filter(|(_, count)| *count == 1)
        .map(|(edge, _)| edge)
        .collect()
}

/// Find boundary vertex indices (vertices on boundary edges).
pub fn boundary_vertices(pd: &PolyData) -> HashSet<usize> {
    let edges = boundary_edges(pd);
    let mut verts = HashSet::new();
    for (a, b) in edges {
        verts.insert(a);
        verts.insert(b);
    }
    verts
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn single_triangle() {
        let pd = PolyData::from_triangles(
            vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]],
            vec![[0, 1, 2]],
        );
        let info = analyze_topology(&pd);
        assert_eq!(info.num_points, 3);
        assert_eq!(info.num_edges, 3);
        assert_eq!(info.num_faces, 1);
        assert_eq!(info.num_boundary_edges, 3);
        assert!(info.is_manifold);
        assert!(!info.is_closed);
        assert!(info.is_triangle_mesh);
        assert_eq!(info.num_components, 1);
        assert_eq!(info.euler_characteristic, 1); // V-E+F = 3-3+1 = 1
    }

    #[test]
    fn two_components() {
        let pd = PolyData::from_triangles(
            vec![
                [0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0],
                [5.0, 0.0, 0.0], [6.0, 0.0, 0.0], [5.0, 1.0, 0.0],
            ],
            vec![[0, 1, 2], [3, 4, 5]],
        );
        let info = analyze_topology(&pd);
        assert_eq!(info.num_components, 2);
    }

    #[test]
    fn boundary_of_single_tri() {
        let pd = PolyData::from_triangles(
            vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]],
            vec![[0, 1, 2]],
        );
        let edges = boundary_edges(&pd);
        assert_eq!(edges.len(), 3);
        let verts = boundary_vertices(&pd);
        assert_eq!(verts.len(), 3);
    }

    #[test]
    fn shared_edge_not_boundary() {
        let pd = PolyData::from_triangles(
            vec![
                [0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.5, 1.0, 0.0], [0.5, -1.0, 0.0],
            ],
            vec![[0, 1, 2], [0, 3, 1]],
        );
        let edges = boundary_edges(&pd);
        // Edge 0-1 is shared by both triangles, so 4 boundary edges total
        assert_eq!(edges.len(), 4);
    }

    #[test]
    fn quad_mesh() {
        let pd = PolyData::from_quads(
            vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [1.0, 1.0, 0.0], [0.0, 1.0, 0.0]],
            vec![[0, 1, 2, 3]],
        );
        let info = analyze_topology(&pd);
        assert!(!info.is_triangle_mesh);
        assert!(info.is_manifold);
    }
}
