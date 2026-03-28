use std::collections::BinaryHeap;
use std::cmp::Ordering;

use vtk_data::{CellArray, Points, PolyData};

/// Simplify a triangle mesh by collapsing edges using a quadric error metric.
///
/// `target_reduction` is the fraction of triangles to remove (0.0 = none, 1.0 = all).
/// The algorithm greedily collapses the lowest-error edge until the target is reached.
///
/// Input must be a triangle mesh (run `triangulate` first if needed).
pub fn decimate(input: &PolyData, target_reduction: f64) -> PolyData {
    let target_reduction = target_reduction.clamp(0.0, 0.99);
    let n_points = input.points.len();
    let n_cells = input.polys.num_cells();

    if n_cells == 0 || n_points < 4 {
        return input.clone();
    }

    let target_cells = ((1.0 - target_reduction) * n_cells as f64).ceil() as usize;
    let target_cells = target_cells.max(1);

    // Copy points
    let mut points: Vec<[f64; 3]> = (0..n_points).map(|i| input.points.get(i)).collect();

    // Copy triangles as mutable Vec<[usize; 3]>
    let mut triangles: Vec<Option<[usize; 3]>> = Vec::with_capacity(n_cells);
    for cell in input.polys.iter() {
        if cell.len() == 3 {
            triangles.push(Some([cell[0] as usize, cell[1] as usize, cell[2] as usize]));
        }
    }

    // Compute quadric matrices per vertex
    let mut quadrics = vec![[0.0f64; 10]; points.len()]; // symmetric 4x4 stored as 10 values
    for &[v0, v1, v2] in triangles.iter().flatten() {
        let q = triangle_quadric(points[v0], points[v1], points[v2]);
        add_quadric(&mut quadrics[v0], &q);
        add_quadric(&mut quadrics[v1], &q);
        add_quadric(&mut quadrics[v2], &q);
    }

    // Build edge set and priority queue
    let mut heap: BinaryHeap<EdgeCollapse> = BinaryHeap::new();
    let mut point_alive = vec![true; points.len()];
    // Union-find for collapsed points
    let mut parent: Vec<usize> = (0..points.len()).collect();

    // Collect edges from triangles
    let mut edges_seen = std::collections::HashSet::new();
    for &[v0, v1, v2] in triangles.iter().flatten() {
        for &(a, b) in &[(v0, v1), (v1, v2), (v2, v0)] {
            let edge = if a < b { (a, b) } else { (b, a) };
            if edges_seen.insert(edge) {
                let cost = edge_collapse_cost(&quadrics[edge.0], &quadrics[edge.1], points[edge.0], points[edge.1]);
                heap.push(EdgeCollapse { cost: -cost, v0: edge.0, v1: edge.1 });
            }
        }
    }

    let mut live_triangles = triangles.iter().filter(|t| t.is_some()).count();

    // Greedily collapse edges
    while live_triangles > target_cells {
        let Some(collapse) = heap.pop() else {
            break;
        };

        let a = find(&parent, collapse.v0);
        let b = find(&parent, collapse.v1);

        if a == b || !point_alive[a] || !point_alive[b] {
            continue;
        }

        // Merge b into a
        let new_pos = midpoint(points[a], points[b]);
        points[a] = new_pos;
        let mut merged_q = quadrics[a];
        add_quadric(&mut merged_q, &quadrics[b]);
        quadrics[a] = merged_q;

        parent[b] = a;
        point_alive[b] = false;

        // Update triangles: replace b with a, remove degenerate
        for tri in &mut triangles {
            if let Some(t) = tri {
                for v in t.iter_mut() {
                    if find(&parent, *v) == b || *v == b {
                        *v = a;
                    } else {
                        *v = find(&parent, *v);
                    }
                }
                // Degenerate check
                if t[0] == t[1] || t[1] == t[2] || t[0] == t[2] {
                    *tri = None;
                    live_triangles -= 1;
                }
            }
        }

        if live_triangles <= target_cells {
            break;
        }

        // Re-insert edges from remaining triangles touching a
        for &[v0, v1, v2] in triangles.iter().flatten() {
            if v0 == a || v1 == a || v2 == a {
                for &(x, y) in &[(v0, v1), (v1, v2), (v2, v0)] {
                    if x != y {
                        let cost = edge_collapse_cost(&quadrics[x], &quadrics[y], points[x], points[y]);
                        heap.push(EdgeCollapse { cost: -cost, v0: x, v1: y });
                    }
                }
            }
        }
    }

    // Build output: compact points and triangles
    let mut point_map = vec![usize::MAX; points.len()];
    let mut out_points = Points::<f64>::new();
    let mut out_polys = CellArray::new();

    for &[v0, v1, v2] in triangles.iter().flatten() {
        for &v in &[v0, v1, v2] {
            if point_map[v] == usize::MAX {
                point_map[v] = out_points.len();
                out_points.push(points[v]);
            }
        }
        out_polys.push_cell(&[
            point_map[v0] as i64,
            point_map[v1] as i64,
            point_map[v2] as i64,
        ]);
    }

    let mut output = PolyData::new();
    output.points = out_points;
    output.polys = out_polys;
    output
}

fn find(parent: &[usize], mut x: usize) -> usize {
    while parent[x] != x {
        x = parent[x];
    }
    x
}

#[derive(Debug)]
struct EdgeCollapse {
    cost: f64, // negated for max-heap → min behavior
    v0: usize,
    v1: usize,
}

impl PartialEq for EdgeCollapse {
    fn eq(&self, other: &Self) -> bool {
        self.cost == other.cost
    }
}
impl Eq for EdgeCollapse {}
impl PartialOrd for EdgeCollapse {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}
impl Ord for EdgeCollapse {
    fn cmp(&self, other: &Self) -> Ordering {
        self.cost.partial_cmp(&other.cost).unwrap_or(Ordering::Equal)
    }
}

/// Compute the plane equation coefficients for a triangle and return its quadric.
/// Quadric is stored as 10 values for the symmetric 4x4 matrix.
fn triangle_quadric(p0: [f64; 3], p1: [f64; 3], p2: [f64; 3]) -> [f64; 10] {
    let e1 = [p1[0] - p0[0], p1[1] - p0[1], p1[2] - p0[2]];
    let e2 = [p2[0] - p0[0], p2[1] - p0[1], p2[2] - p0[2]];
    let n = [
        e1[1] * e2[2] - e1[2] * e2[1],
        e1[2] * e2[0] - e1[0] * e2[2],
        e1[0] * e2[1] - e1[1] * e2[0],
    ];
    let len = (n[0] * n[0] + n[1] * n[1] + n[2] * n[2]).sqrt();
    if len < 1e-10 {
        return [0.0; 10];
    }
    let a = n[0] / len;
    let b = n[1] / len;
    let c = n[2] / len;
    let d = -(a * p0[0] + b * p0[1] + c * p0[2]);

    // Q = [a,b,c,d]^T * [a,b,c,d] (outer product), stored as symmetric:
    // [a*a, a*b, a*c, a*d, b*b, b*c, b*d, c*c, c*d, d*d]
    [a*a, a*b, a*c, a*d, b*b, b*c, b*d, c*c, c*d, d*d]
}

fn add_quadric(a: &mut [f64; 10], b: &[f64; 10]) {
    for i in 0..10 {
        a[i] += b[i];
    }
}

fn edge_collapse_cost(q0: &[f64; 10], q1: &[f64; 10], p0: [f64; 3], p1: [f64; 3]) -> f64 {
    // Combined quadric
    let mut q = [0.0f64; 10];
    for i in 0..10 {
        q[i] = q0[i] + q1[i];
    }
    // Evaluate at midpoint
    let mp = midpoint(p0, p1);
    evaluate_quadric(&q, mp)
}

fn evaluate_quadric(q: &[f64; 10], p: [f64; 3]) -> f64 {
    // Q stored as [a*a, a*b, a*c, a*d, b*b, b*c, b*d, c*c, c*d, d*d]
    // v = [x, y, z, 1]
    // error = v^T Q v
    let x = p[0];
    let y = p[1];
    let z = p[2];
    q[0]*x*x + 2.0*q[1]*x*y + 2.0*q[2]*x*z + 2.0*q[3]*x
    + q[4]*y*y + 2.0*q[5]*y*z + 2.0*q[6]*y
    + q[7]*z*z + 2.0*q[8]*z
    + q[9]
}

fn midpoint(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [(a[0]+b[0])*0.5, (a[1]+b[1])*0.5, (a[2]+b[2])*0.5]
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::sources;

    #[test]
    fn decimate_sphere() {
        let sphere = sources::sphere::sphere(&sources::sphere::SphereParams {
            theta_resolution: 16,
            phi_resolution: 16,
            ..Default::default()
        });
        // Triangulate first (sphere has quads)
        let tri = crate::triangulate::triangulate(&sphere);
        let original_cells = tri.polys.num_cells();

        let decimated = decimate(&tri, 0.5);

        assert!(
            decimated.polys.num_cells() < original_cells,
            "should have fewer cells: {} vs {}",
            decimated.polys.num_cells(),
            original_cells
        );
        assert!(decimated.polys.num_cells() > 0);
    }

    #[test]
    fn decimate_zero_reduction() {
        let pd = PolyData::from_triangles(
            vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]],
            vec![[0, 1, 2]],
        );
        let result = decimate(&pd, 0.0);
        assert_eq!(result.polys.num_cells(), 1);
    }
}
