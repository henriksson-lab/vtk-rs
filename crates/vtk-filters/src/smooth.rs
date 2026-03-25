use std::collections::HashSet;

use vtk_data::PolyData;

/// Laplacian smoothing of a PolyData mesh.
///
/// Each vertex is moved toward the average position of its neighbors.
/// `iterations` controls how many smoothing passes to apply.
/// `relaxation_factor` controls the step size (0.0 = no movement, 1.0 = move to average).
/// Boundary vertices (edges with only one adjacent polygon) are optionally fixed.
pub fn smooth(
    input: &PolyData,
    iterations: usize,
    relaxation_factor: f64,
    fix_boundary: bool,
) -> PolyData {
    let mut output = input.clone();
    let n = output.points.len();
    if n == 0 || iterations == 0 {
        return output;
    }

    // Build adjacency: for each point, collect its neighbor point indices
    let mut neighbors: Vec<HashSet<usize>> = vec![HashSet::new(); n];
    for cell in input.polys.iter() {
        let len = cell.len();
        for j in 0..len {
            let a = cell[j] as usize;
            let b = cell[(j + 1) % len] as usize;
            neighbors[a].insert(b);
            neighbors[b].insert(a);
        }
    }

    // Find boundary vertices if needed
    let boundary = if fix_boundary {
        find_boundary_vertices(input)
    } else {
        HashSet::new()
    };

    let factor = relaxation_factor.clamp(0.0, 1.0);

    for _ in 0..iterations {
        let mut new_positions = Vec::with_capacity(n);

        for (i, nbrs) in neighbors.iter().enumerate() {
            let p = output.points.get(i);
            if boundary.contains(&i) || nbrs.is_empty() {
                new_positions.push(p);
                continue;
            }

            // Average of neighbors
            let mut avg = [0.0f64; 3];
            let count = nbrs.len() as f64;
            for &nb in nbrs {
                let q = output.points.get(nb);
                avg[0] += q[0];
                avg[1] += q[1];
                avg[2] += q[2];
            }
            avg[0] /= count;
            avg[1] /= count;
            avg[2] /= count;

            // Blend
            new_positions.push([
                p[0] + factor * (avg[0] - p[0]),
                p[1] + factor * (avg[1] - p[1]),
                p[2] + factor * (avg[2] - p[2]),
            ]);
        }

        for (i, &pos) in new_positions.iter().enumerate() {
            output.points.set(i, pos);
        }
    }

    output
}

/// Find vertices on the mesh boundary (edges shared by only one polygon).
fn find_boundary_vertices(input: &PolyData) -> HashSet<usize> {
    use std::collections::HashMap;

    let mut edge_count: HashMap<(usize, usize), usize> = HashMap::new();
    for cell in input.polys.iter() {
        let len = cell.len();
        for j in 0..len {
            let a = cell[j] as usize;
            let b = cell[(j + 1) % len] as usize;
            let edge = if a < b { (a, b) } else { (b, a) };
            *edge_count.entry(edge).or_insert(0) += 1;
        }
    }

    let mut boundary = HashSet::new();
    for (&(a, b), &count) in &edge_count {
        if count == 1 {
            boundary.insert(a);
            boundary.insert(b);
        }
    }
    boundary
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn smooth_preserves_topology() {
        let pd = PolyData::from_triangles(
            vec![
                [0.0, 0.0, 0.0],
                [1.0, 0.0, 0.0],
                [0.5, 1.0, 0.0],
                [1.5, 1.0, 0.0],
            ],
            vec![[0, 1, 2], [1, 3, 2]],
        );

        let result = smooth(&pd, 5, 0.5, false);
        assert_eq!(result.points.len(), 4);
        assert_eq!(result.polys.num_cells(), 2);
    }

    #[test]
    fn smooth_moves_interior() {
        // Create a mesh with an interior vertex that can be smoothed
        let pd = PolyData::from_triangles(
            vec![
                [0.0, 0.0, 0.0],
                [2.0, 0.0, 0.0],
                [2.0, 2.0, 0.0],
                [0.0, 2.0, 0.0],
                [0.5, 0.5, 0.0], // off-center interior point
            ],
            vec![[0, 1, 4], [1, 2, 4], [2, 3, 4], [3, 0, 4]],
        );

        let before = pd.points.get(4);
        let result = smooth(&pd, 10, 0.5, false);
        let after = result.points.get(4);

        // Interior point should move toward center (1, 1, 0)
        assert!(
            (after[0] - 1.0).abs() < (before[0] - 1.0).abs(),
            "x should move toward center"
        );
    }

    #[test]
    fn zero_iterations_noop() {
        let pd = PolyData::from_triangles(
            vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]],
            vec![[0, 1, 2]],
        );
        let result = smooth(&pd, 0, 0.5, false);
        assert_eq!(result.points.get(0), pd.points.get(0));
    }
}
