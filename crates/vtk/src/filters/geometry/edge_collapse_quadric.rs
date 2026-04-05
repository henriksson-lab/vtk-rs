use crate::data::{CellArray, Points, PolyData};
use std::collections::{BinaryHeap, HashMap, HashSet};
use std::cmp::Ordering;

/// Simplify a triangle mesh using quadric error metrics (Garland-Heckbert).
///
/// Iteratively collapses edges with the lowest quadric error until
/// `target_ratio` of the original face count remains (0.0-1.0).
/// More accurate than `quadric_clustering` for preserving shape.
pub fn edge_collapse_quadric(input: &PolyData, target_ratio: f64) -> PolyData {
    let target = ((input.polys.num_cells() as f64 * target_ratio.clamp(0.01, 1.0)).ceil() as usize).max(1);

    let n = input.points.len();
    let mut pts: Vec<[f64; 3]> = (0..n).map(|i| input.points.get(i)).collect();

    let mut tris: Vec<[usize; 3]> = input.polys.iter()
        .filter(|c| c.len() >= 3)
        .map(|c| [c[0] as usize, c[1] as usize, c[2] as usize])
        .collect();

    // Compute per-vertex quadric matrices (symmetric 4x4 stored as 10 floats)
    let mut quadrics = vec![[0.0f64; 10]; n];
    for tri in &tris {
        let q = face_quadric(&pts[tri[0]], &pts[tri[1]], &pts[tri[2]]);
        for &v in tri { add_quadric(&mut quadrics[v], &q); }
    }

    // Build edge set
    let mut edges: HashSet<(usize, usize)> = HashSet::new();
    for tri in &tris {
        for k in 0..3 {
            let a = tri[k]; let b = tri[(k+1)%3];
            edges.insert(if a < b { (a,b) } else { (b,a) });
        }
    }

    // Merge map
    let mut remap: Vec<usize> = (0..n).collect();
    let find = |remap: &mut Vec<usize>, mut x: usize| -> usize {
        while remap[x] != x { remap[x] = remap[remap[x]]; x = remap[x]; }
        x
    };

    let mut current_faces = tris.len();

    // Greedy: sort edges by cost, collapse cheapest
    let mut edge_list: Vec<(f64, usize, usize)> = edges.iter().map(|&(a,b)| {
        let cost = edge_cost(&quadrics[a], &quadrics[b], &pts[a], &pts[b]);
        (cost, a, b)
    }).collect();
    edge_list.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap_or(Ordering::Equal));

    for &(_, a, b) in &edge_list {
        if current_faces <= target { break; }

        let ra = find(&mut remap, a);
        let rb = find(&mut remap, b);
        if ra == rb { continue; }

        // Collapse rb into ra
        let mid = [(pts[ra][0]+pts[rb][0])*0.5, (pts[ra][1]+pts[rb][1])*0.5, (pts[ra][2]+pts[rb][2])*0.5];
        pts[ra] = mid;
        let qb = quadrics[rb];
        add_quadric(&mut quadrics[ra], &qb);
        remap[rb] = ra;

        // Count removed faces
        let before = current_faces;
        tris.retain(|tri| {
            let v: Vec<usize> = tri.iter().map(|&v| find(&mut remap, v)).collect();
            v[0] != v[1] && v[1] != v[2] && v[0] != v[2]
        });
        current_faces = tris.len();
    }

    // Remap and compact
    let mut pt_map: HashMap<usize, i64> = HashMap::new();
    let mut out_points = Points::<f64>::new();
    let mut out_polys = CellArray::new();

    for tri in &tris {
        let mapped: Vec<i64> = tri.iter().map(|&v| {
            let rv = find(&mut remap, v);
            *pt_map.entry(rv).or_insert_with(|| {
                let idx = out_points.len() as i64;
                out_points.push(pts[rv]);
                idx
            })
        }).collect();
        out_polys.push_cell(&mapped);
    }

    let mut pd = PolyData::new();
    pd.points = out_points;
    pd.polys = out_polys;
    pd
}

fn face_quadric(v0: &[f64; 3], v1: &[f64; 3], v2: &[f64; 3]) -> [f64; 10] {
    let e1 = [v1[0]-v0[0], v1[1]-v0[1], v1[2]-v0[2]];
    let e2 = [v2[0]-v0[0], v2[1]-v0[1], v2[2]-v0[2]];
    let mut n = [e1[1]*e2[2]-e1[2]*e2[1], e1[2]*e2[0]-e1[0]*e2[2], e1[0]*e2[1]-e1[1]*e2[0]];
    let len = (n[0]*n[0]+n[1]*n[1]+n[2]*n[2]).sqrt();
    if len > 1e-15 { n[0] /= len; n[1] /= len; n[2] /= len; }
    let d = -(n[0]*v0[0] + n[1]*v0[1] + n[2]*v0[2]);
    // Q = [a b c d]^T * [a b c d] stored as upper triangle
    let (a, b, c) = (n[0], n[1], n[2]);
    [a*a, a*b, a*c, a*d, b*b, b*c, b*d, c*c, c*d, d*d]
}

fn add_quadric(dest: &mut [f64; 10], src: &[f64; 10]) {
    for i in 0..10 { dest[i] += src[i]; }
}

fn edge_cost(qa: &[f64; 10], qb: &[f64; 10], pa: &[f64; 3], pb: &[f64; 3]) -> f64 {
    let mid = [(pa[0]+pb[0])*0.5, (pa[1]+pb[1])*0.5, (pa[2]+pb[2])*0.5];
    let mut q = [0.0f64; 10];
    for i in 0..10 { q[i] = qa[i] + qb[i]; }
    eval_quadric(&q, &mid)
}

fn eval_quadric(q: &[f64; 10], v: &[f64; 3]) -> f64 {
    let (x, y, z) = (v[0], v[1], v[2]);
    q[0]*x*x + 2.0*q[1]*x*y + 2.0*q[2]*x*z + 2.0*q[3]*x
        + q[4]*y*y + 2.0*q[5]*y*z + 2.0*q[6]*y
        + q[7]*z*z + 2.0*q[8]*z + q[9]
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn reduces_face_count() {
        let mut pd = PolyData::new();
        for j in 0..5 {
            for i in 0..5 {
                pd.points.push([i as f64, j as f64, 0.0]);
            }
        }
        for j in 0..4 {
            for i in 0..4 {
                let a = (j*5+i) as i64;
                pd.polys.push_cell(&[a, a+1, a+6]);
                pd.polys.push_cell(&[a, a+6, a+5]);
            }
        }

        let result = edge_collapse_quadric(&pd, 0.5);
        assert!(result.polys.num_cells() < pd.polys.num_cells());
        assert!(result.polys.num_cells() > 0);
    }

    #[test]
    fn ratio_1_preserves() {
        let mut pd = PolyData::new();
        pd.points.push([0.0, 0.0, 0.0]);
        pd.points.push([1.0, 0.0, 0.0]);
        pd.points.push([0.5, 1.0, 0.0]);
        pd.polys.push_cell(&[0, 1, 2]);

        let result = edge_collapse_quadric(&pd, 1.0);
        assert_eq!(result.polys.num_cells(), 1);
    }

    #[test]
    fn empty_input() {
        let pd = PolyData::new();
        let result = edge_collapse_quadric(&pd, 0.5);
        assert_eq!(result.polys.num_cells(), 0);
    }
}
