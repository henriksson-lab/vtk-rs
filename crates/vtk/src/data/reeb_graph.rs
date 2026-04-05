//! Reeb graph data structure for scalar field topology analysis.
//!
//! A Reeb graph captures the topological structure of a scalar field on a mesh
//! by tracking how level sets merge and split as the scalar value increases.

use crate::data::{CellArray, Points, PolyData};

/// Type of a critical point in the Reeb graph.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum NodeType {
    /// Local minimum: all neighbors have higher scalar values.
    Minimum,
    /// Local maximum: all neighbors have lower scalar values.
    Maximum,
    /// Saddle point: level set topology changes.
    Saddle,
    /// Regular point (not a critical point).
    Regular,
}

/// A node in the Reeb graph corresponding to a critical point.
#[derive(Debug, Clone)]
pub struct ReebNode {
    /// Index of the vertex in the original mesh.
    pub vertex_id: usize,
    /// Scalar value at this vertex.
    pub scalar_value: f64,
    /// Type of critical point.
    pub node_type: NodeType,
}

/// An arc in the Reeb graph connecting two critical points.
#[derive(Debug, Clone)]
pub struct ReebArc {
    /// Index of the source node in the Reeb graph.
    pub source: usize,
    /// Index of the target node in the Reeb graph.
    pub target: usize,
    /// Vertex indices of points along this arc (excluding endpoints).
    pub vertex_ids: Vec<usize>,
}

/// Reeb graph representation of scalar field topology.
///
/// Nodes represent critical points (minima, maxima, saddles) and arcs
/// represent monotone paths between them.
#[derive(Debug, Clone, Default)]
pub struct ReebGraph {
    /// Critical point nodes.
    pub nodes: Vec<ReebNode>,
    /// Arcs connecting critical points.
    pub arcs: Vec<ReebArc>,
}

impl ReebGraph {
    /// Create an empty Reeb graph.
    pub fn new() -> Self {
        Self::default()
    }

    /// Add a node to the graph and return its index.
    pub fn add_node(&mut self, vertex_id: usize, scalar_value: f64, node_type: NodeType) -> usize {
        let idx = self.nodes.len();
        self.nodes.push(ReebNode {
            vertex_id,
            scalar_value,
            node_type,
        });
        idx
    }

    /// Add an arc between two nodes and return its index.
    pub fn add_arc(&mut self, source: usize, target: usize, vertex_ids: Vec<usize>) -> usize {
        let idx = self.arcs.len();
        self.arcs.push(ReebArc {
            source,
            target,
            vertex_ids,
        });
        idx
    }

    /// Number of nodes.
    pub fn num_nodes(&self) -> usize {
        self.nodes.len()
    }

    /// Number of arcs.
    pub fn num_arcs(&self) -> usize {
        self.arcs.len()
    }

    /// Get node by index.
    pub fn node(&self, idx: usize) -> &ReebNode {
        &self.nodes[idx]
    }

    /// Get arc by index.
    pub fn arc(&self, idx: usize) -> &ReebArc {
        &self.arcs[idx]
    }

    /// Return indices of all nodes matching the given type.
    pub fn nodes_of_type(&self, node_type: NodeType) -> Vec<usize> {
        self.nodes
            .iter()
            .enumerate()
            .filter(|(_, n)| n.node_type == node_type)
            .map(|(i, _)| i)
            .collect()
    }

    /// Convert the Reeb graph to PolyData with arcs as line cells.
    ///
    /// Each node becomes a point, and each arc becomes a line cell
    /// connecting source to target through intermediate vertices.
    pub fn to_poly_data(&self) -> PolyData {
        let mut points = Points::<f64>::new();
        let mut lines = CellArray::new();

        // Add node points — node index i maps to point index i
        for node in &self.nodes {
            // Place nodes along x-axis at their scalar value, y=0
            points.push([node.scalar_value, 0.0, 0.0]);
        }

        // Add arcs as line cells
        for arc in &self.arcs {
            let mut cell = Vec::with_capacity(2 + arc.vertex_ids.len());
            cell.push(arc.source as i64);
            // Add intermediate points
            for &vid in &arc.vertex_ids {
                let pt_idx = points.len();
                // Place intermediate points at some interpolated position
                let _ = vid; // vertex_id from original mesh, use index as position
                let s0 = self.nodes[arc.source].scalar_value;
                let s1 = self.nodes[arc.target].scalar_value;
                let t = if arc.vertex_ids.len() > 0 {
                    (cell.len() as f64) / (arc.vertex_ids.len() as f64 + 1.0)
                } else {
                    0.5
                };
                points.push([s0 + t * (s1 - s0), 0.0, 0.0]);
                cell.push(pt_idx as i64);
            }
            cell.push(arc.target as i64);
            lines.push_cell(&cell);
        }

        let mut pd = PolyData::new();
        pd.points = points;
        pd.lines = lines;
        pd
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn basic_reeb_graph() {
        let mut rg = ReebGraph::new();
        let n0 = rg.add_node(0, 0.0, NodeType::Minimum);
        let n1 = rg.add_node(1, 1.0, NodeType::Maximum);
        rg.add_arc(n0, n1, vec![]);

        assert_eq!(rg.num_nodes(), 2);
        assert_eq!(rg.num_arcs(), 1);
        assert_eq!(rg.node(0).node_type, NodeType::Minimum);
        assert_eq!(rg.arc(0).source, 0);
        assert_eq!(rg.arc(0).target, 1);
    }

    #[test]
    fn nodes_of_type() {
        let mut rg = ReebGraph::new();
        rg.add_node(0, 0.0, NodeType::Minimum);
        rg.add_node(1, 0.5, NodeType::Saddle);
        rg.add_node(2, 1.0, NodeType::Maximum);
        rg.add_node(3, 0.2, NodeType::Minimum);

        let mins = rg.nodes_of_type(NodeType::Minimum);
        assert_eq!(mins.len(), 2);
        assert_eq!(mins, vec![0, 3]);

        let maxs = rg.nodes_of_type(NodeType::Maximum);
        assert_eq!(maxs.len(), 1);
    }

    #[test]
    fn to_poly_data() {
        let mut rg = ReebGraph::new();
        let n0 = rg.add_node(0, 0.0, NodeType::Minimum);
        let n1 = rg.add_node(1, 0.5, NodeType::Saddle);
        let n2 = rg.add_node(2, 1.0, NodeType::Maximum);
        rg.add_arc(n0, n1, vec![]);
        rg.add_arc(n1, n2, vec![]);

        let pd = rg.to_poly_data();
        assert_eq!(pd.points.len(), 3); // 3 nodes, no intermediate
        assert_eq!(pd.lines.num_cells(), 2);
    }
}
