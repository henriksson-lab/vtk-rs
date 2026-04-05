use crate::data::{DataSetAttributes, FieldData};
use crate::types::BoundingBox;

/// A node in a hyper tree (octree/quadtree cell).
#[derive(Debug, Clone)]
struct HyperTreeNode {
    /// Whether this node is a leaf.
    is_leaf: bool,
    /// Index of the first child (children are stored contiguously).
    /// For a 3D tree: 8 children at indices first_child..first_child+8.
    first_child: usize,
    /// Global cell index for leaf nodes (used to look up cell data).
    global_id: usize,
}

/// A single hyper tree (recursive octree/quadtree).
#[derive(Debug, Clone)]
pub struct HyperTree {
    nodes: Vec<HyperTreeNode>,
    /// Branching factor: 4 for 2D (quadtree), 8 for 3D (octree).
    branch_factor: usize,
}

#[allow(dead_code)]
impl HyperTree {
    fn new(branch_factor: usize) -> Self {
        Self {
            nodes: vec![HyperTreeNode {
                is_leaf: true,
                first_child: 0,
                global_id: 0,
            }],
            branch_factor,
        }
    }

    /// Subdivide a leaf node, creating `branch_factor` children.
    /// Returns the index of the first child.
    fn subdivide(&mut self, node_idx: usize, next_global_id: &mut usize) -> usize {
        let first = self.nodes.len();
        self.nodes[node_idx].is_leaf = false;
        self.nodes[node_idx].first_child = first;

        for _ in 0..self.branch_factor {
            self.nodes.push(HyperTreeNode {
                is_leaf: true,
                first_child: 0,
                global_id: *next_global_id,
            });
            *next_global_id += 1;
        }
        first
    }

    /// Number of nodes in this tree.
    fn num_nodes(&self) -> usize {
        self.nodes.len()
    }

    /// Number of leaf nodes.
    fn num_leaves(&self) -> usize {
        self.nodes.iter().filter(|n| n.is_leaf).count()
    }

    /// Maximum depth of the tree.
    fn max_depth(&self) -> usize {
        self.depth_recursive(0)
    }

    fn depth_recursive(&self, node_idx: usize) -> usize {
        let node = &self.nodes[node_idx];
        if node.is_leaf {
            return 0;
        }
        let mut max_d = 0;
        for c in 0..self.branch_factor {
            let child_idx = node.first_child + c;
            if child_idx < self.nodes.len() {
                max_d = max_d.max(self.depth_recursive(child_idx));
            }
        }
        1 + max_d
    }
}

/// AMR-style hierarchical grid.
///
/// Analogous to VTK's `vtkHyperTreeGrid`. Stores a coarse grid where each
/// cell can be recursively subdivided into an octree (3D) or quadtree (2D).
///
/// This allows adaptive mesh refinement (AMR) with different resolution
/// in different parts of the domain.
#[derive(Debug, Clone)]
pub struct HyperTreeGrid {
    /// Number of dimensions (2 or 3).
    dimension: usize,
    /// Coarse grid dimensions [ni, nj, nk].
    grid_size: [usize; 3],
    /// Origin of the grid.
    origin: [f64; 3],
    /// Grid spacing for the coarse level.
    spacing: [f64; 3],
    /// One hyper tree per coarse cell.
    trees: Vec<Option<HyperTree>>,
    /// Global cell data (indexed by global_id across all trees).
    cell_data: DataSetAttributes,
    /// Field data.
    field_data: FieldData,
    /// Next available global cell ID.
    next_global_id: usize,
}

impl HyperTreeGrid {
    /// Create a new HyperTreeGrid.
    ///
    /// `grid_size` is the coarse grid dimensions [ni, nj, nk].
    /// For 2D grids, set nk=1.
    pub fn new(grid_size: [usize; 3], origin: [f64; 3], spacing: [f64; 3]) -> Self {
        let dimension = if grid_size[2] <= 1 { 2 } else { 3 };
        let n_cells = grid_size[0] * grid_size[1] * grid_size[2];
        Self {
            dimension,
            grid_size,
            origin,
            spacing,
            trees: vec![None; n_cells],
            cell_data: DataSetAttributes::new(),
            field_data: FieldData::new(),
            next_global_id: 0,
        }
    }

    /// Initialize a tree at coarse cell (i, j, k).
    /// Returns the global cell ID assigned to the root leaf.
    pub fn init_tree(&mut self, i: usize, j: usize, k: usize) -> usize {
        let idx = self.coarse_index(i, j, k);
        let bf = if self.dimension == 2 { 4 } else { 8 };
        let gid = self.next_global_id;
        self.next_global_id += 1;
        let mut tree = HyperTree::new(bf);
        tree.nodes[0].global_id = gid;
        self.trees[idx] = Some(tree);
        gid
    }

    /// Subdivide a leaf node in the tree at coarse cell (i, j, k).
    ///
    /// `node_index` is the index within the tree's node array (0 = root).
    /// Returns the index of the first child node.
    pub fn subdivide(&mut self, i: usize, j: usize, k: usize, node_index: usize) -> Option<usize> {
        let idx = self.coarse_index(i, j, k);
        let tree = self.trees[idx].as_mut()?;
        if !tree.nodes[node_index].is_leaf {
            return None;
        }
        Some(tree.subdivide(node_index, &mut self.next_global_id))
    }

    /// Coarse grid dimensions.
    pub fn grid_size(&self) -> [usize; 3] {
        self.grid_size
    }

    /// Number of dimensions (2 or 3).
    pub fn dimension(&self) -> usize {
        self.dimension
    }

    /// Number of coarse cells.
    pub fn num_coarse_cells(&self) -> usize {
        self.grid_size[0] * self.grid_size[1] * self.grid_size[2]
    }

    /// Total number of leaf cells across all trees (= actual data cells).
    pub fn num_cells(&self) -> usize {
        self.next_global_id
    }

    /// Maximum refinement depth across all trees.
    pub fn max_depth(&self) -> usize {
        self.trees.iter().filter_map(|t| t.as_ref()).map(|t| t.max_depth()).max().unwrap_or(0)
    }

    /// Number of initialized trees.
    pub fn num_trees(&self) -> usize {
        self.trees.iter().filter(|t| t.is_some()).count()
    }

    /// Compute bounding box of the coarse grid.
    pub fn bounds(&self) -> BoundingBox {
        BoundingBox {
            x_min: self.origin[0],
            x_max: self.origin[0] + self.grid_size[0] as f64 * self.spacing[0],
            y_min: self.origin[1],
            y_max: self.origin[1] + self.grid_size[1] as f64 * self.spacing[1],
            z_min: self.origin[2],
            z_max: self.origin[2] + self.grid_size[2] as f64 * self.spacing[2],
        }
    }

    pub fn cell_data(&self) -> &DataSetAttributes {
        &self.cell_data
    }

    pub fn cell_data_mut(&mut self) -> &mut DataSetAttributes {
        &mut self.cell_data
    }

    pub fn field_data(&self) -> &FieldData {
        &self.field_data
    }

    pub fn field_data_mut(&mut self) -> &mut FieldData {
        &mut self.field_data
    }

    fn coarse_index(&self, i: usize, j: usize, k: usize) -> usize {
        k * self.grid_size[0] * self.grid_size[1] + j * self.grid_size[0] + i
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn basic_2d_grid() {
        let htg = HyperTreeGrid::new([4, 4, 1], [0.0, 0.0, 0.0], [1.0, 1.0, 1.0]);
        assert_eq!(htg.dimension(), 2);
        assert_eq!(htg.num_coarse_cells(), 16);
        assert_eq!(htg.num_cells(), 0);
    }

    #[test]
    fn basic_3d_grid() {
        let htg = HyperTreeGrid::new([2, 2, 2], [0.0, 0.0, 0.0], [1.0, 1.0, 1.0]);
        assert_eq!(htg.dimension(), 3);
        assert_eq!(htg.num_coarse_cells(), 8);
    }

    #[test]
    fn init_and_subdivide() {
        let mut htg = HyperTreeGrid::new([2, 2, 1], [0.0, 0.0, 0.0], [1.0, 1.0, 1.0]);

        let gid = htg.init_tree(0, 0, 0);
        assert_eq!(gid, 0);
        assert_eq!(htg.num_cells(), 1);

        // Subdivide root: creates 4 children (2D quadtree)
        let first_child = htg.subdivide(0, 0, 0, 0).unwrap();
        assert_eq!(first_child, 1); // first child node index
        assert_eq!(htg.num_cells(), 5); // root no longer leaf, 4 new leaves

        assert_eq!(htg.max_depth(), 1);
    }

    #[test]
    fn multi_level_refinement() {
        let mut htg = HyperTreeGrid::new([1, 1, 1], [0.0, 0.0, 0.0], [1.0, 1.0, 1.0]);
        htg.init_tree(0, 0, 0);

        // Level 1
        htg.subdivide(0, 0, 0, 0).unwrap();
        assert_eq!(htg.max_depth(), 1);

        // Level 2: subdivide first child (node index 1)
        htg.subdivide(0, 0, 0, 1).unwrap();
        assert_eq!(htg.max_depth(), 2);
    }

    #[test]
    fn bounds() {
        let htg = HyperTreeGrid::new([4, 3, 2], [1.0, 2.0, 3.0], [0.5, 0.5, 0.5]);
        let bb = htg.bounds();
        assert_eq!(bb.x_min, 1.0);
        assert_eq!(bb.x_max, 3.0);
        assert_eq!(bb.y_min, 2.0);
        assert_eq!(bb.y_max, 3.5);
    }
}
