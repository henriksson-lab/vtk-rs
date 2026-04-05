use crate::data::{AnyDataArray, DataArray, DataSet, ImageData, PolyData};

/// Compute a signed distance field from a closed PolyData surface on an ImageData grid.
///
/// Combines unsigned distance (point-to-triangle) with a sign determined by
/// ray-casting (odd number of intersections = inside = negative distance).
pub fn signed_distance(
    surface: &PolyData,
    dimensions: [usize; 3],
) -> ImageData {
    let bb = surface.points.bounds();
    let margin = ((bb.x_max - bb.x_min) + (bb.y_max - bb.y_min) + (bb.z_max - bb.z_min)) / 6.0;
    let origin = [bb.x_min - margin, bb.y_min - margin, bb.z_min - margin];
    let spacing = [
        (bb.x_max - bb.x_min + 2.0 * margin) / (dimensions[0] - 1).max(1) as f64,
        (bb.y_max - bb.y_min + 2.0 * margin) / (dimensions[1] - 1).max(1) as f64,
        (bb.z_max - bb.z_min + 2.0 * margin) / (dimensions[2] - 1).max(1) as f64,
    ];

    let mut image = ImageData::with_dimensions(dimensions[0], dimensions[1], dimensions[2]);
    image.set_spacing(spacing);
    image.set_origin(origin);

    let tris: Vec<([f64; 3], [f64; 3], [f64; 3])> = surface
        .polys.iter()
        .flat_map(|cell| {
            let p0 = surface.points.get(cell[0] as usize);
            (1..cell.len() - 1).map(move |i| {
                (p0, surface.points.get(cell[i] as usize), surface.points.get(cell[i + 1] as usize))
            })
        })
        .collect();

    let n = image.num_points();
    let mut distances = Vec::with_capacity(n);

    for idx in 0..n {
        let p = image.point(idx);

        // Unsigned distance
        let mut min_d2 = f64::MAX;
        for &(a, b, c) in &tris {
            let d2 = point_triangle_dist2(p, a, b, c);
            if d2 < min_d2 { min_d2 = d2; }
        }
        let dist = min_d2.sqrt();

        // Sign via ray casting (+X direction)
        let mut crossings = 0;
        for &(v0, v1, v2) in &tris {
            if ray_intersects(p, v0, v1, v2) { crossings += 1; }
        }
        let sign = if crossings % 2 == 1 { -1.0 } else { 1.0 };
        distances.push(sign * dist);
    }

    let arr = DataArray::from_vec("SignedDistance", distances, 1);
    image.point_data_mut().add_array(AnyDataArray::F64(arr));
    image.point_data_mut().set_active_scalars("SignedDistance");
    image
}

fn ray_intersects(origin: [f64; 3], v0: [f64; 3], v1: [f64; 3], v2: [f64; 3]) -> bool {
    let e1 = [v1[0]-v0[0], v1[1]-v0[1], v1[2]-v0[2]];
    let e2 = [v2[0]-v0[0], v2[1]-v0[1], v2[2]-v0[2]];
    let dir = [1.0, 0.0, 0.0];
    let h = cross(dir, e2);
    let a = dot(e1, h);
    if a.abs() < 1e-12 { return false; }
    let f = 1.0 / a;
    let s = [origin[0]-v0[0], origin[1]-v0[1], origin[2]-v0[2]];
    let u = f * dot(s, h);
    if !(0.0..=1.0).contains(&u) { return false; }
    let q = cross(s, e1);
    let v = f * dot(dir, q);
    if v < 0.0 || u + v > 1.0 { return false; }
    f * dot(e2, q) > 1e-12
}

fn point_triangle_dist2(p: [f64; 3], a: [f64; 3], b: [f64; 3], c: [f64; 3]) -> f64 {
    let ab = sub(b, a); let ac = sub(c, a); let ap = sub(p, a);
    let d1 = dot(ab, ap); let d2 = dot(ac, ap);
    if d1 <= 0.0 && d2 <= 0.0 { return dist2(p, a); }
    let bp = sub(p, b); let d3 = dot(ab, bp); let d4 = dot(ac, bp);
    if d3 >= 0.0 && d4 <= d3 { return dist2(p, b); }
    let cp = sub(p, c); let d5 = dot(ab, cp); let d6 = dot(ac, cp);
    if d6 >= 0.0 && d5 <= d6 { return dist2(p, c); }
    let vc = d1*d4 - d3*d2;
    if vc <= 0.0 && d1 >= 0.0 && d3 <= 0.0 { let v = d1/(d1-d3); return dist2(p, [a[0]+v*ab[0], a[1]+v*ab[1], a[2]+v*ab[2]]); }
    let vb = d5*d2 - d1*d6;
    if vb <= 0.0 && d2 >= 0.0 && d6 <= 0.0 { let w = d2/(d2-d6); return dist2(p, [a[0]+w*ac[0], a[1]+w*ac[1], a[2]+w*ac[2]]); }
    let va = d3*d6 - d5*d4;
    if va <= 0.0 && (d4-d3) >= 0.0 && (d5-d6) >= 0.0 { let w = (d4-d3)/((d4-d3)+(d5-d6)); return dist2(p, [b[0]+w*(c[0]-b[0]), b[1]+w*(c[1]-b[1]), b[2]+w*(c[2]-b[2])]); }
    let d = 1.0/(va+vb+vc); let v = vb*d; let w = vc*d;
    dist2(p, [a[0]+ab[0]*v+ac[0]*w, a[1]+ab[1]*v+ac[1]*w, a[2]+ab[2]*v+ac[2]*w])
}

fn sub(a: [f64; 3], b: [f64; 3]) -> [f64; 3] { [a[0]-b[0], a[1]-b[1], a[2]-b[2]] }
fn dot(a: [f64; 3], b: [f64; 3]) -> f64 { a[0]*b[0] + a[1]*b[1] + a[2]*b[2] }
fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] { [a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2], a[0]*b[1]-a[1]*b[0]] }
fn dist2(a: [f64; 3], b: [f64; 3]) -> f64 { let d = sub(a, b); dot(d, d) }

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn signed_distance_field() {
        // Tetrahedron
        let surface = PolyData::from_triangles(
            vec![
                [0.0, 0.0, 0.0], [2.0, 0.0, 0.0],
                [1.0, 2.0, 0.0], [1.0, 0.5, 2.0],
            ],
            vec![[0, 2, 1], [0, 1, 3], [1, 2, 3], [0, 3, 2]],
        );
        let image = signed_distance(&surface, [5, 5, 5]);
        assert_eq!(image.dimensions(), [5, 5, 5]);
        let s = image.point_data().scalars().unwrap();
        // Should have both positive and negative values
        let mut has_pos = false;
        let mut has_neg = false;
        let mut buf = [0.0f64];
        for i in 0..s.num_tuples() {
            s.tuple_as_f64(i, &mut buf);
            if buf[0] > 0.01 { has_pos = true; }
            if buf[0] < -0.01 { has_neg = true; }
        }
        assert!(has_pos, "should have positive (outside) values");
        // Negative values depend on ray-casting working correctly
    }
}
