/// A simple octree for spatial point queries.
///
/// Recursively subdivides space into octants. Supports nearest-neighbor
/// and radius search. Less optimal than KdTree for most cases but useful
/// for uniform point distributions.
#[derive(Debug, Clone)]
pub struct OctreePointLocator {
    root: Option<OctreeNode>,
    bounds: [f64; 6],
}

#[derive(Debug, Clone)]
enum OctreeNode {
    Leaf(Vec<(usize, [f64; 3])>),
    Branch(Box<[Option<OctreeNode>; 8]>),
}

impl OctreePointLocator {
    /// Build an octree from points with default settings.
    pub fn build(points: &[[f64; 3]]) -> Self {
        Self::build_with_params(points, 10, 16)
    }

    /// Build with custom parameters.
    pub fn build_with_params(
        points: &[[f64; 3]],
        max_points_per_leaf: usize,
        max_depth: usize,
    ) -> Self {
        if points.is_empty() {
            return Self { root: None, bounds: [0.0; 6] };
        }

        let mut bounds = [f64::MAX, f64::MIN, f64::MAX, f64::MIN, f64::MAX, f64::MIN];
        for p in points {
            bounds[0] = bounds[0].min(p[0]);
            bounds[1] = bounds[1].max(p[0]);
            bounds[2] = bounds[2].min(p[1]);
            bounds[3] = bounds[3].max(p[1]);
            bounds[4] = bounds[4].min(p[2]);
            bounds[5] = bounds[5].max(p[2]);
        }
        // Add small margin
        let eps = 1e-6;
        for i in (0..6).step_by(2) { bounds[i] -= eps; }
        for i in (1..6).step_by(2) { bounds[i] += eps; }

        let indexed: Vec<(usize, [f64; 3])> = points.iter().enumerate().map(|(i, &p)| (i, p)).collect();
        let root = build_node(indexed, &bounds, 0, max_depth, max_points_per_leaf);

        Self { root: Some(root), bounds }
    }

    /// Find nearest point. Returns (index, distance²).
    pub fn nearest(&self, query: [f64; 3]) -> Option<(usize, f64)> {
        let root = self.root.as_ref()?;
        let mut best_idx = 0;
        let mut best_d2 = f64::MAX;
        search_nearest(root, &self.bounds, query, &mut best_idx, &mut best_d2);
        Some((best_idx, best_d2))
    }

    /// Find all points within radius.
    pub fn find_within_radius(&self, query: [f64; 3], radius: f64) -> Vec<(usize, f64)> {
        let mut results = Vec::new();
        if let Some(root) = &self.root {
            let r2 = radius * radius;
            search_radius(root, &self.bounds, query, r2, &mut results);
        }
        results
    }
}

fn build_node(
    points: Vec<(usize, [f64; 3])>,
    bounds: &[f64; 6],
    depth: usize,
    max_depth: usize,
    max_per_leaf: usize,
) -> OctreeNode {
    if points.len() <= max_per_leaf || depth >= max_depth {
        return OctreeNode::Leaf(points);
    }

    let mid = [
        (bounds[0] + bounds[1]) * 0.5,
        (bounds[2] + bounds[3]) * 0.5,
        (bounds[4] + bounds[5]) * 0.5,
    ];

    let mut buckets: [Vec<(usize, [f64; 3])>; 8] = Default::default();
    for (idx, p) in points {
        let octant = ((p[0] >= mid[0]) as usize)
            | (((p[1] >= mid[1]) as usize) << 1)
            | (((p[2] >= mid[2]) as usize) << 2);
        buckets[octant].push((idx, p));
    }

    let children: [Option<OctreeNode>; 8] = std::array::from_fn(|i| {
        if buckets[i].is_empty() {
            None
        } else {
            let child_bounds = child_bounds(bounds, &mid, i);
            Some(build_node(std::mem::take(&mut buckets[i]), &child_bounds, depth + 1, max_depth, max_per_leaf))
        }
    });

    OctreeNode::Branch(Box::new(children))
}

fn child_bounds(bounds: &[f64; 6], mid: &[f64; 3], octant: usize) -> [f64; 6] {
    [
        if octant & 1 == 0 { bounds[0] } else { mid[0] },
        if octant & 1 == 0 { mid[0] } else { bounds[1] },
        if octant & 2 == 0 { bounds[2] } else { mid[1] },
        if octant & 2 == 0 { mid[1] } else { bounds[3] },
        if octant & 4 == 0 { bounds[4] } else { mid[2] },
        if octant & 4 == 0 { mid[2] } else { bounds[5] },
    ]
}

fn search_nearest(node: &OctreeNode, bounds: &[f64; 6], query: [f64; 3], best_idx: &mut usize, best_d2: &mut f64) {
    match node {
        OctreeNode::Leaf(pts) => {
            for &(idx, p) in pts {
                let d2 = dist2(query, p);
                if d2 < *best_d2 { *best_d2 = d2; *best_idx = idx; }
            }
        }
        OctreeNode::Branch(children) => {
            let mid = [(bounds[0]+bounds[1])*0.5, (bounds[2]+bounds[3])*0.5, (bounds[4]+bounds[5])*0.5];
            for (i, child) in children.iter().enumerate() {
                if let Some(c) = child {
                    let cb = child_bounds(bounds, &mid, i);
                    if min_dist2_to_box(query, &cb) < *best_d2 {
                        search_nearest(c, &cb, query, best_idx, best_d2);
                    }
                }
            }
        }
    }
}

fn search_radius(node: &OctreeNode, bounds: &[f64; 6], query: [f64; 3], r2: f64, results: &mut Vec<(usize, f64)>) {
    match node {
        OctreeNode::Leaf(pts) => {
            for &(idx, p) in pts {
                let d2 = dist2(query, p);
                if d2 <= r2 { results.push((idx, d2)); }
            }
        }
        OctreeNode::Branch(children) => {
            let mid = [(bounds[0]+bounds[1])*0.5, (bounds[2]+bounds[3])*0.5, (bounds[4]+bounds[5])*0.5];
            for (i, child) in children.iter().enumerate() {
                if let Some(c) = child {
                    let cb = child_bounds(bounds, &mid, i);
                    if min_dist2_to_box(query, &cb) <= r2 {
                        search_radius(c, &cb, query, r2, results);
                    }
                }
            }
        }
    }
}

fn min_dist2_to_box(p: [f64; 3], b: &[f64; 6]) -> f64 {
    let dx = (b[0] - p[0]).max(0.0).max(p[0] - b[1]);
    let dy = (b[2] - p[1]).max(0.0).max(p[1] - b[3]);
    let dz = (b[4] - p[2]).max(0.0).max(p[2] - b[5]);
    dx * dx + dy * dy + dz * dz
}

fn dist2(a: [f64; 3], b: [f64; 3]) -> f64 {
    (a[0]-b[0])*(a[0]-b[0]) + (a[1]-b[1])*(a[1]-b[1]) + (a[2]-b[2])*(a[2]-b[2])
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn nearest_neighbor() {
        let points = vec![[0.0,0.0,0.0],[1.0,0.0,0.0],[0.0,1.0,0.0],[5.0,5.0,5.0]];
        let tree = OctreePointLocator::build(&points);
        let (idx, _) = tree.nearest([0.1, 0.0, 0.0]).unwrap();
        assert_eq!(idx, 0);
    }

    #[test]
    fn radius_search() {
        let points = vec![[0.0,0.0,0.0],[0.5,0.0,0.0],[10.0,0.0,0.0]];
        let tree = OctreePointLocator::build(&points);
        let results = tree.find_within_radius([0.0, 0.0, 0.0], 1.0);
        assert_eq!(results.len(), 2);
    }

    #[test]
    fn empty() {
        let tree = OctreePointLocator::build(&[]);
        assert!(tree.nearest([0.0, 0.0, 0.0]).is_none());
    }
}
