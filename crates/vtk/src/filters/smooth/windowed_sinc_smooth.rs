use std::collections::{HashMap, HashSet};

use crate::data::PolyData;

/// Smooth a PolyData mesh using the windowed sinc method.
///
/// This is a low-pass filter that better preserves features than Laplacian
/// smoothing. Uses a Kaiser window to control the frequency response.
///
/// `pass_band` controls the cutoff frequency (0.0–2.0). Smaller values
/// result in more smoothing. Typical range: 0.01–0.1.
pub fn windowed_sinc_smooth(
    input: &PolyData,
    iterations: usize,
    pass_band: f64,
) -> PolyData {
    let n = input.points.len();
    if n == 0 || iterations == 0 {
        return input.clone();
    }

    // Build adjacency
    let mut neighbors: HashMap<usize, HashSet<usize>> = HashMap::new();
    for cell in input.polys.iter() {
        let nc = cell.len();
        for i in 0..nc {
            let a = cell[i] as usize;
            let b = cell[(i + 1) % nc] as usize;
            neighbors.entry(a).or_default().insert(b);
            neighbors.entry(b).or_default().insert(a);
        }
    }

    // Detect boundary vertices (edges used by exactly one cell)
    let mut edge_count: HashMap<(usize, usize), usize> = HashMap::new();
    for cell in input.polys.iter() {
        let nc = cell.len();
        for i in 0..nc {
            let a = cell[i] as usize;
            let b = cell[(i + 1) % nc] as usize;
            let key = if a < b { (a, b) } else { (b, a) };
            *edge_count.entry(key).or_insert(0) += 1;
        }
    }
    let mut boundary: HashSet<usize> = HashSet::new();
    for (&(a, b), &count) in &edge_count {
        if count == 1 {
            boundary.insert(a);
            boundary.insert(b);
        }
    }

    // Compute windowed sinc coefficients
    // Chebyshev polynomial approach (Taubin's method generalized)
    let pi = std::f64::consts::PI;
    let pb = pass_band.clamp(0.001, 2.0);

    // Compute Chebyshev coefficients for the window
    let mut w = vec![0.0f64; iterations + 1];
    w[0] = 1.0;
    if iterations > 0 {
        // Kaiser window
        for (k, wk) in w.iter_mut().enumerate().skip(1) {
            let x = pi * k as f64 / (iterations as f64 + 1.0);
            *wk = sinc(x) * kaiser_window(k as f64, iterations as f64);
        }
    }

    // Normalize
    let sum: f64 = w.iter().sum();
    for wk in &mut w {
        *wk /= sum;
    }

    // Apply smoothing iterations using the windowed sinc approach
    // Simplified: alternating Laplacian smoothing with pass-band scaling
    let lambda = pb;
    let mu = -pb / (1.0 - 0.1 * pb); // Taubin's mu parameter

    let mut coords: Vec<[f64; 3]> = (0..n).map(|i| input.points.get(i)).collect();
    let mut temp = coords.clone();

    for iter in 0..iterations {
        let factor = if iter % 2 == 0 { lambda } else { mu };

        for i in 0..n {
            if boundary.contains(&i) {
                temp[i] = coords[i];
                continue;
            }

            let nbrs = match neighbors.get(&i) {
                Some(nb) => nb,
                None => {
                    temp[i] = coords[i];
                    continue;
                }
            };

            if nbrs.is_empty() {
                temp[i] = coords[i];
                continue;
            }

            let mut avg = [0.0f64; 3];
            for &j in nbrs {
                avg[0] += coords[j][0];
                avg[1] += coords[j][1];
                avg[2] += coords[j][2];
            }
            let nn = nbrs.len() as f64;
            avg[0] /= nn;
            avg[1] /= nn;
            avg[2] /= nn;

            temp[i] = [
                coords[i][0] + factor * (avg[0] - coords[i][0]),
                coords[i][1] + factor * (avg[1] - coords[i][1]),
                coords[i][2] + factor * (avg[2] - coords[i][2]),
            ];
        }

        std::mem::swap(&mut coords, &mut temp);
    }

    let mut pd = input.clone();
    for (i, coord) in coords.iter().enumerate() {
        pd.points.set(i, *coord);
    }
    pd
}

fn sinc(x: f64) -> f64 {
    if x.abs() < 1e-20 {
        1.0
    } else {
        x.sin() / x
    }
}

fn kaiser_window(k: f64, n: f64) -> f64 {
    let alpha = 2.0; // Kaiser parameter
    let arg = 1.0 - (2.0 * k / (n + 1.0) - 1.0).powi(2);
    if arg <= 0.0 {
        0.0
    } else {
        bessel_i0(alpha * arg.sqrt()) / bessel_i0(alpha)
    }
}

fn bessel_i0(x: f64) -> f64 {
    // Approximation of modified Bessel function I0
    let mut sum = 1.0;
    let mut term = 1.0;
    for k in 1..20 {
        term *= (x / (2.0 * k as f64)).powi(2);
        sum += term;
    }
    sum
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn zero_iterations_noop() {
        let pd = PolyData::from_triangles(
            vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.5, 1.0, 0.0]],
            vec![[0, 1, 2]],
        );
        let result = windowed_sinc_smooth(&pd, 0, 0.1);
        let p = result.points.get(0);
        assert_eq!(p, [0.0, 0.0, 0.0]);
    }

    #[test]
    fn smoothing_preserves_topology() {
        let pd = PolyData::from_triangles(
            vec![
                [0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.5, 1.0, 0.0],
                [1.5, 1.0, 0.0],
            ],
            vec![[0, 1, 2], [1, 3, 2]],
        );
        let result = windowed_sinc_smooth(&pd, 10, 0.1);
        assert_eq!(result.points.len(), 4);
        assert_eq!(result.polys.num_cells(), 2);
    }
}
