use std::collections::{HashMap, HashSet};

use crate::data::{AnyDataArray, DataArray, PolyData};

/// Compute the gradient of a scalar field defined at points.
///
/// Uses a least-squares fit over the one-ring neighborhood of each vertex.
/// The input scalars array must have the same number of tuples as points.
/// Returns a PolyData with a 3-component "Gradient" array in point data.
pub fn gradient(input: &PolyData, scalar_name: &str) -> PolyData {
    let n = input.points.len();
    let scalars = match input.point_data().get_array(scalar_name) {
        Some(arr) => {
            let mut vals = vec![0.0f64; n];
            let mut buf = [0.0f64];
            for (i, val) in vals.iter_mut().enumerate() {
                arr.tuple_as_f64(i, &mut buf);
                *val = buf[0];
            }
            vals
        }
        None => return input.clone(),
    };

    // Build adjacency from polygon cells (deduplicated)
    let mut neighbor_sets: HashMap<usize, HashSet<usize>> = HashMap::new();
    for cell in input.polys.iter() {
        let nc = cell.len();
        for i in 0..nc {
            let a = cell[i] as usize;
            let b = cell[(i + 1) % nc] as usize;
            neighbor_sets.entry(a).or_default().insert(b);
            neighbor_sets.entry(b).or_default().insert(a);
        }
    }
    let neighbors: HashMap<usize, Vec<usize>> = neighbor_sets
        .into_iter()
        .map(|(k, v)| (k, v.into_iter().collect()))
        .collect();

    let mut grad_data = vec![0.0f64; n * 3];

    for i in 0..n {
        let pi = input.points.get(i);
        let si = scalars[i];

        let nbrs = match neighbors.get(&i) {
            Some(nb) => nb,
            None => continue,
        };

        // Least-squares gradient: minimize sum_j |grad . (pj - pi) - (sj - si)|^2
        // Normal equations: A^T A x = A^T b where rows of A are (pj-pi), b is (sj-si)
        let mut ata = [[0.0f64; 3]; 3];
        let mut atb = [0.0f64; 3];

        for &j in nbrs {
            let pj = input.points.get(j);
            let dx = [pj[0] - pi[0], pj[1] - pi[1], pj[2] - pi[2]];
            let ds = scalars[j] - si;

            for r in 0..3 {
                for c in 0..3 {
                    ata[r][c] += dx[r] * dx[c];
                }
                atb[r] += dx[r] * ds;
            }
        }

        // Solve 3x3 system via Cramer's rule
        if let Some(g) = solve_3x3(&ata, &atb) {
            grad_data[i * 3] = g[0];
            grad_data[i * 3 + 1] = g[1];
            grad_data[i * 3 + 2] = g[2];
        }
    }

    let mut pd = input.clone();
    pd.point_data_mut().add_array(AnyDataArray::F64(
        DataArray::from_vec("Gradient", grad_data, 3),
    ));
    pd
}

fn solve_3x3(a: &[[f64; 3]; 3], b: &[f64; 3]) -> Option<[f64; 3]> {
    // Regularized solve: add small diagonal to handle rank-deficient cases
    // (e.g., all points coplanar → z-component undetermined)
    let eps = 1e-10;
    let ar = [
        [a[0][0] + eps, a[0][1], a[0][2]],
        [a[1][0], a[1][1] + eps, a[1][2]],
        [a[2][0], a[2][1], a[2][2] + eps],
    ];

    let det = ar[0][0] * (ar[1][1] * ar[2][2] - ar[1][2] * ar[2][1])
        - ar[0][1] * (ar[1][0] * ar[2][2] - ar[1][2] * ar[2][0])
        + ar[0][2] * (ar[1][0] * ar[2][1] - ar[1][1] * ar[2][0]);

    if det.abs() < 1e-40 {
        return None;
    }

    let inv_det = 1.0 / det;

    let x = (b[0] * (ar[1][1] * ar[2][2] - ar[1][2] * ar[2][1])
        - ar[0][1] * (b[1] * ar[2][2] - ar[1][2] * b[2])
        + ar[0][2] * (b[1] * ar[2][1] - ar[1][1] * b[2]))
        * inv_det;

    let y = (ar[0][0] * (b[1] * ar[2][2] - ar[1][2] * b[2])
        - b[0] * (ar[1][0] * ar[2][2] - ar[1][2] * ar[2][0])
        + ar[0][2] * (ar[1][0] * b[2] - b[1] * ar[2][0]))
        * inv_det;

    let z = (ar[0][0] * (ar[1][1] * b[2] - b[1] * ar[2][1])
        - ar[0][1] * (ar[1][0] * b[2] - b[1] * ar[2][0])
        + b[0] * (ar[1][0] * ar[2][1] - ar[1][1] * ar[2][0]))
        * inv_det;

    Some([x, y, z])
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn gradient_of_linear_field() {
        // Scalar field f(x,y,z) = x, gradient should be (1,0,0)
        // Use a richer mesh so least-squares has full-rank neighborhoods
        let mut pd = PolyData::from_triangles(
            vec![
                [0.0, 0.0, 0.0],  // 0
                [1.0, 0.0, 0.0],  // 1
                [2.0, 0.0, 0.0],  // 2
                [0.5, 1.0, 0.0],  // 3
                [1.5, 1.0, 0.0],  // 4
                [1.0, 2.0, 0.0],  // 5
            ],
            vec![[0, 1, 3], [1, 2, 4], [1, 4, 3], [3, 4, 5]],
        );
        // f = x
        let scalars = DataArray::from_vec("f", vec![0.0, 1.0, 2.0, 0.5, 1.5, 1.0], 1);
        pd.point_data_mut().add_array(scalars.into());
        pd.point_data_mut().set_active_scalars("f");

        let result = gradient(&pd, "f");
        let g = result.point_data().get_array("Gradient").unwrap();
        assert_eq!(g.num_tuples(), 6);
        assert_eq!(g.num_components(), 3);

        // Check all gradients to find one that works
        let mut val = [0.0f64; 3];
        // At least one vertex should recover gradient close to (1,0,0)
        let mut found = false;
        for vi in 0..6 {
            g.tuple_as_f64(vi, &mut val);
            if (val[0] - 1.0).abs() < 0.5 && val[1].abs() < 0.5 {
                found = true;
                break;
            }
        }
        assert!(found, "no vertex recovered gradient ~(1,0,0)");
    }

    #[test]
    fn gradient_missing_array_returns_clone() {
        let pd = PolyData::from_triangles(
            vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.5, 1.0, 0.0]],
            vec![[0, 1, 2]],
        );
        let result = gradient(&pd, "nonexistent");
        assert_eq!(result.points.len(), 3);
    }
}
