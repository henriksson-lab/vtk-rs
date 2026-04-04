//! Oriented Bounding Box (OBB) tree for fast spatial queries on meshes.
//!
//! Builds a binary tree of oriented bounding boxes for ray intersection,
//! closest-point queries, and collision detection.

use crate::PolyData;

/// An oriented bounding box defined by a center, three axes, and half-extents.
#[derive(Debug, Clone)]
pub struct Obb {
    pub center: [f64; 3],
    pub axes: [[f64; 3]; 3],
    pub half_extents: [f64; 3],
}

impl Obb {
    /// Compute OBB from a set of points using PCA.
    pub fn from_points(points: &[[f64; 3]]) -> Self {
        if points.is_empty() {
            return Self {
                center: [0.0; 3],
                axes: [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
                half_extents: [0.0; 3],
            };
        }

        // Compute centroid
        let n = points.len() as f64;
        let mut center = [0.0; 3];
        for p in points {
            center[0] += p[0]; center[1] += p[1]; center[2] += p[2];
        }
        center[0] /= n; center[1] /= n; center[2] /= n;

        // Compute covariance matrix
        let mut cov = [[0.0f64; 3]; 3];
        for p in points {
            let d = [p[0] - center[0], p[1] - center[1], p[2] - center[2]];
            for i in 0..3 {
                for j in 0..3 {
                    cov[i][j] += d[i] * d[j];
                }
            }
        }
        for i in 0..3 { for j in 0..3 { cov[i][j] /= n; } }

        // Eigendecomposition via power iteration (3 iterations)
        let axes = eigen_axes_3x3(&cov);

        // Project points onto axes to get extents
        let mut mins = [f64::MAX; 3];
        let mut maxs = [f64::MIN; 3];
        for p in points {
            let d = [p[0] - center[0], p[1] - center[1], p[2] - center[2]];
            for (a, (mn, mx)) in axes.iter().zip(mins.iter_mut().zip(maxs.iter_mut())) {
                let proj = d[0] * a[0] + d[1] * a[1] + d[2] * a[2];
                *mn = mn.min(proj);
                *mx = mx.max(proj);
            }
        }

        let half_extents = [
            (maxs[0] - mins[0]) / 2.0,
            (maxs[1] - mins[1]) / 2.0,
            (maxs[2] - mins[2]) / 2.0,
        ];

        // Adjust center to OBB center
        let mid = [
            (maxs[0] + mins[0]) / 2.0,
            (maxs[1] + mins[1]) / 2.0,
            (maxs[2] + mins[2]) / 2.0,
        ];
        let obb_center = [
            center[0] + mid[0] * axes[0][0] + mid[1] * axes[1][0] + mid[2] * axes[2][0],
            center[1] + mid[0] * axes[0][1] + mid[1] * axes[1][1] + mid[2] * axes[2][1],
            center[2] + mid[0] * axes[0][2] + mid[1] * axes[1][2] + mid[2] * axes[2][2],
        ];

        Self { center: obb_center, axes, half_extents }
    }

    /// Test if a point is inside the OBB.
    pub fn contains(&self, point: [f64; 3]) -> bool {
        let d = [
            point[0] - self.center[0],
            point[1] - self.center[1],
            point[2] - self.center[2],
        ];
        for (i, axis) in self.axes.iter().enumerate() {
            let proj = (d[0] * axis[0] + d[1] * axis[1] + d[2] * axis[2]).abs();
            if proj > self.half_extents[i] {
                return false;
            }
        }
        true
    }

    /// Volume of the OBB.
    pub fn volume(&self) -> f64 {
        8.0 * self.half_extents[0] * self.half_extents[1] * self.half_extents[2]
    }
}

/// OBB tree node.
#[derive(Debug)]
enum ObbNode {
    Leaf {
        obb: Obb,
        cell_indices: Vec<usize>,
    },
    Internal {
        obb: Obb,
        left: Box<ObbNode>,
        right: Box<ObbNode>,
    },
}

/// Oriented Bounding Box tree for spatial queries.
#[derive(Debug)]
pub struct ObbTree {
    root: Option<ObbNode>,
}

impl ObbTree {
    /// Build an OBB tree from a PolyData mesh.
    pub fn build(poly_data: &PolyData, max_leaf_size: usize) -> Self {
        let num_cells = poly_data.polys.num_cells();
        if num_cells == 0 {
            return Self { root: None };
        }

        // Compute centroids for all cells
        let mut centroids = Vec::with_capacity(num_cells);
        for ci in 0..num_cells {
            let cell = poly_data.polys.cell(ci);
            let mut cx = 0.0;
            let mut cy = 0.0;
            let mut cz = 0.0;
            for &vid in cell {
                let p = poly_data.points.get(vid as usize);
                cx += p[0]; cy += p[1]; cz += p[2];
            }
            let n = cell.len() as f64;
            centroids.push([cx / n, cy / n, cz / n]);
        }

        let indices: Vec<usize> = (0..num_cells).collect();
        let root = Self::build_node(&centroids, indices, max_leaf_size);
        Self { root: Some(root) }
    }

    fn build_node(centroids: &[[f64; 3]], indices: Vec<usize>, max_leaf_size: usize) -> ObbNode {
        let points: Vec<[f64; 3]> = indices.iter().map(|&i| centroids[i]).collect();
        let obb = Obb::from_points(&points);

        if indices.len() <= max_leaf_size {
            return ObbNode::Leaf { obb, cell_indices: indices };
        }

        // Split along the longest OBB axis
        let split_axis = if obb.half_extents[0] >= obb.half_extents[1] && obb.half_extents[0] >= obb.half_extents[2] {
            0
        } else if obb.half_extents[1] >= obb.half_extents[2] {
            1
        } else {
            2
        };

        // Project centroids onto split axis and split at median
        let axis = obb.axes[split_axis];
        let mut projections: Vec<(usize, f64)> = indices
            .iter()
            .map(|&i| {
                let d = [
                    centroids[i][0] - obb.center[0],
                    centroids[i][1] - obb.center[1],
                    centroids[i][2] - obb.center[2],
                ];
                (i, d[0] * axis[0] + d[1] * axis[1] + d[2] * axis[2])
            })
            .collect();
        projections.sort_by(|a, b| a.1.partial_cmp(&b.1).unwrap_or(std::cmp::Ordering::Equal));

        let mid = projections.len() / 2;
        let left_indices: Vec<usize> = projections[..mid].iter().map(|(i, _)| *i).collect();
        let right_indices: Vec<usize> = projections[mid..].iter().map(|(i, _)| *i).collect();

        if left_indices.is_empty() || right_indices.is_empty() {
            return ObbNode::Leaf { obb, cell_indices: projections.iter().map(|(i, _)| *i).collect() };
        }

        let left = Box::new(Self::build_node(centroids, left_indices, max_leaf_size));
        let right = Box::new(Self::build_node(centroids, right_indices, max_leaf_size));

        ObbNode::Internal { obb, left, right }
    }

    /// Find all leaf cell indices whose OBB contains the given point.
    pub fn find_cells_containing(&self, point: [f64; 3]) -> Vec<usize> {
        let mut result = Vec::new();
        if let Some(ref root) = self.root {
            Self::query_node(root, point, &mut result);
        }
        result
    }

    fn query_node(node: &ObbNode, point: [f64; 3], result: &mut Vec<usize>) {
        match node {
            ObbNode::Leaf { obb, cell_indices } => {
                if obb.contains(point) {
                    result.extend_from_slice(cell_indices);
                }
            }
            ObbNode::Internal { obb, left, right } => {
                if obb.contains(point) {
                    Self::query_node(left, point, result);
                    Self::query_node(right, point, result);
                }
            }
        }
    }

    /// Count total leaf cells.
    pub fn num_cells(&self) -> usize {
        fn count(node: &ObbNode) -> usize {
            match node {
                ObbNode::Leaf { cell_indices, .. } => cell_indices.len(),
                ObbNode::Internal { left, right, .. } => count(left) + count(right),
            }
        }
        self.root.as_ref().map_or(0, count)
    }
}

/// Compute 3 orthogonal eigenvectors of a 3x3 symmetric matrix via power iteration.
fn eigen_axes_3x3(cov: &[[f64; 3]; 3]) -> [[f64; 3]; 3] {
    let v1 = power_iteration(cov, [1.0, 0.0, 0.0], 20);
    // Deflate
    let e1 = dot3(&mat_vec(cov, v1), v1);
    let mut cov2 = *cov;
    for i in 0..3 { for j in 0..3 { cov2[i][j] -= e1 * v1[i] * v1[j]; } }
    let v2_raw = power_iteration(&cov2, [0.0, 1.0, 0.0], 20);
    // Gram-Schmidt
    let proj = dot3(&v2_raw, v1);
    let mut v2 = [v2_raw[0] - proj * v1[0], v2_raw[1] - proj * v1[1], v2_raw[2] - proj * v1[2]];
    let len = (v2[0]*v2[0]+v2[1]*v2[1]+v2[2]*v2[2]).sqrt();
    if len > 1e-15 { v2[0] /= len; v2[1] /= len; v2[2] /= len; }
    else { v2 = [0.0, 1.0, 0.0]; }

    let v3 = cross3(v1, v2);
    [v1, v2, v3]
}

fn power_iteration(mat: &[[f64; 3]; 3], mut v: [f64; 3], iters: usize) -> [f64; 3] {
    for _ in 0..iters {
        v = mat_vec(mat, v);
        let len = (v[0]*v[0]+v[1]*v[1]+v[2]*v[2]).sqrt();
        if len < 1e-30 { return [1.0, 0.0, 0.0]; }
        v[0] /= len; v[1] /= len; v[2] /= len;
    }
    v
}

fn mat_vec(m: &[[f64; 3]; 3], v: [f64; 3]) -> [f64; 3] {
    [
        m[0][0]*v[0]+m[0][1]*v[1]+m[0][2]*v[2],
        m[1][0]*v[0]+m[1][1]*v[1]+m[1][2]*v[2],
        m[2][0]*v[0]+m[2][1]*v[1]+m[2][2]*v[2],
    ]
}

fn dot3(a: &[f64; 3], b: [f64; 3]) -> f64 { a[0]*b[0]+a[1]*b[1]+a[2]*b[2] }

fn cross3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2], a[0]*b[1]-a[1]*b[0]]
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn obb_from_axis_aligned_points() {
        let points = vec![[0.0, 0.0, 0.0], [2.0, 0.0, 0.0], [0.0, 1.0, 0.0], [2.0, 1.0, 1.0]];
        let obb = Obb::from_points(&points);
        assert!(obb.contains([1.0, 0.5, 0.25]));
        assert!(!obb.contains([5.0, 5.0, 5.0]));
        assert!(obb.volume() > 0.0);
    }

    #[test]
    fn obb_tree_build_and_query() {
        let pd = PolyData::from_triangles(
            vec![[0.0,0.0,0.0],[1.0,0.0,0.0],[0.0,1.0,0.5],[2.0,0.0,0.0],[2.0,1.0,0.5]],
            vec![[0,1,2],[1,3,4]],
        );
        let tree = ObbTree::build(&pd, 1);
        assert_eq!(tree.num_cells(), 2);
    }
}
