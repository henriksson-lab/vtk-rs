use std::collections::HashMap;

use crate::data::PolyData;

/// Result of computing mesh topology genus.
#[derive(Debug, Clone)]
pub struct TopologyResult {
    /// Number of vertices.
    pub vertices: usize,
    /// Number of unique edges.
    pub edges: usize,
    /// Number of faces (polygon cells).
    pub faces: usize,
    /// Euler characteristic: V - E + F.
    pub euler_characteristic: i64,
    /// Genus computed from Euler's formula for closed orientable surfaces:
    /// V - E + F = 2(1 - g), so g = 1 - (V - E + F) / 2.
    /// Only meaningful for closed manifold surfaces.
    pub genus: i64,
}

/// Compute the genus of a mesh using Euler's formula.
///
/// For a closed orientable surface, V - E + F = 2(1 - g), giving
/// g = 1 - (V - E + F) / 2.
///
/// Note: the result is only geometrically meaningful when the mesh is a
/// closed manifold surface (no boundary edges, no non-manifold edges).
pub fn compute_genus(input: &PolyData) -> TopologyResult {
    let v: usize = input.points.len();
    let f: usize = input.polys.num_cells();

    // Count unique edges
    let mut edge_set: HashMap<(i64, i64), bool> = HashMap::new();
    for cell in input.polys.iter() {
        let n: usize = cell.len();
        for i in 0..n {
            let a: i64 = cell[i];
            let b: i64 = cell[(i + 1) % n];
            let key = if a < b { (a, b) } else { (b, a) };
            edge_set.entry(key).or_insert(true);
        }
    }

    let e: usize = edge_set.len();
    let chi: i64 = v as i64 - e as i64 + f as i64;
    let g: i64 = 1 - chi / 2;

    TopologyResult {
        vertices: v,
        edges: e,
        faces: f,
        euler_characteristic: chi,
        genus: g,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn tetrahedron_genus_zero() {
        // Tetrahedron: V=4, E=6, F=4, chi=2, genus=0
        let mut pd = PolyData::new();
        pd.points.push([0.0, 0.0, 0.0]);
        pd.points.push([1.0, 0.0, 0.0]);
        pd.points.push([0.5, 1.0, 0.0]);
        pd.points.push([0.5, 0.5, 1.0]);
        pd.polys.push_cell(&[0, 1, 2]);
        pd.polys.push_cell(&[0, 1, 3]);
        pd.polys.push_cell(&[1, 2, 3]);
        pd.polys.push_cell(&[0, 2, 3]);

        let result = compute_genus(&pd);
        assert_eq!(result.vertices, 4);
        assert_eq!(result.edges, 6);
        assert_eq!(result.faces, 4);
        assert_eq!(result.euler_characteristic, 2);
        assert_eq!(result.genus, 0);
    }

    #[test]
    fn cube_genus_zero() {
        // Cube as 12 triangles: V=8, E=18, F=12, chi=2, genus=0
        let mut pd = PolyData::new();
        pd.points.push([0.0, 0.0, 0.0]); // 0
        pd.points.push([1.0, 0.0, 0.0]); // 1
        pd.points.push([1.0, 1.0, 0.0]); // 2
        pd.points.push([0.0, 1.0, 0.0]); // 3
        pd.points.push([0.0, 0.0, 1.0]); // 4
        pd.points.push([1.0, 0.0, 1.0]); // 5
        pd.points.push([1.0, 1.0, 1.0]); // 6
        pd.points.push([0.0, 1.0, 1.0]); // 7

        // 6 faces, each as 2 triangles = 12 triangles
        // bottom
        pd.polys.push_cell(&[0, 1, 2]);
        pd.polys.push_cell(&[0, 2, 3]);
        // top
        pd.polys.push_cell(&[4, 6, 5]);
        pd.polys.push_cell(&[4, 7, 6]);
        // front
        pd.polys.push_cell(&[0, 5, 1]);
        pd.polys.push_cell(&[0, 4, 5]);
        // back
        pd.polys.push_cell(&[2, 7, 3]);
        pd.polys.push_cell(&[2, 6, 7]);
        // left
        pd.polys.push_cell(&[0, 3, 7]);
        pd.polys.push_cell(&[0, 7, 4]);
        // right
        pd.polys.push_cell(&[1, 5, 6]);
        pd.polys.push_cell(&[1, 6, 2]);

        let result = compute_genus(&pd);
        assert_eq!(result.vertices, 8);
        assert_eq!(result.edges, 18);
        assert_eq!(result.faces, 12);
        assert_eq!(result.euler_characteristic, 2);
        assert_eq!(result.genus, 0);
    }

    #[test]
    fn single_triangle_open() {
        // Single triangle: V=3, E=3, F=1, chi=1
        // genus formula gives g = 1 - 1/2 = 0 (integer division)
        let mut pd = PolyData::new();
        pd.points.push([0.0, 0.0, 0.0]);
        pd.points.push([1.0, 0.0, 0.0]);
        pd.points.push([0.5, 1.0, 0.0]);
        pd.polys.push_cell(&[0, 1, 2]);

        let result = compute_genus(&pd);
        assert_eq!(result.vertices, 3);
        assert_eq!(result.edges, 3);
        assert_eq!(result.faces, 1);
        assert_eq!(result.euler_characteristic, 1);
    }
}
