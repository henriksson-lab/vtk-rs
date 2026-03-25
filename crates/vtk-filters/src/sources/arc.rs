use std::f64::consts::PI;

use vtk_data::{Points, PolyData};

/// Parameters for generating a circular arc.
pub struct ArcParams {
    /// Center of the arc. Default: [0, 0, 0]
    pub center: [f64; 3],
    /// Radius of the arc. Default: 0.5
    pub radius: f64,
    /// Start angle in degrees. Default: 0
    pub start_angle: f64,
    /// End angle in degrees. Default: 90
    pub end_angle: f64,
    /// Number of points along the arc. Default: 20
    pub resolution: usize,
    /// Normal direction of the arc plane. Default: [0, 0, 1] (XY plane)
    pub normal: [f64; 3],
}

impl Default for ArcParams {
    fn default() -> Self {
        Self {
            center: [0.0, 0.0, 0.0],
            radius: 0.5,
            start_angle: 0.0,
            end_angle: 90.0,
            resolution: 20,
            normal: [0.0, 0.0, 1.0],
        }
    }
}

/// Generate a circular arc as a polyline.
pub fn arc(params: &ArcParams) -> PolyData {
    let n = params.resolution.max(2);
    let a0 = params.start_angle * PI / 180.0;
    let a1 = params.end_angle * PI / 180.0;

    // Build a coordinate frame from the normal
    let nz = normalize(params.normal);
    let (nx, ny) = perpendicular_frame(nz);

    let mut points = Points::new();
    let mut ids = Vec::with_capacity(n);

    for i in 0..n {
        let t = i as f64 / (n - 1) as f64;
        let angle = a0 + t * (a1 - a0);
        let ct = angle.cos();
        let st = angle.sin();

        let x = params.center[0] + params.radius * (ct * nx[0] + st * ny[0]);
        let y = params.center[1] + params.radius * (ct * nx[1] + st * ny[1]);
        let z = params.center[2] + params.radius * (ct * nx[2] + st * ny[2]);

        points.push([x, y, z]);
        ids.push(i as i64);
    }

    let mut pd = PolyData::new();
    pd.points = points;
    pd.lines.push_cell(&ids);
    pd
}

fn normalize(v: [f64; 3]) -> [f64; 3] {
    let len = (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt();
    if len > 1e-10 {
        [v[0] / len, v[1] / len, v[2] / len]
    } else {
        [0.0, 0.0, 1.0]
    }
}

fn perpendicular_frame(dir: [f64; 3]) -> ([f64; 3], [f64; 3]) {
    let seed = if dir[0].abs() < 0.9 {
        [1.0, 0.0, 0.0]
    } else {
        [0.0, 1.0, 0.0]
    };
    let u = normalize(cross(seed, dir));
    let v = cross(dir, u);
    (u, v)
}

fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn default_arc() {
        let pd = arc(&ArcParams::default());
        assert_eq!(pd.points.len(), 20);
        assert_eq!(pd.lines.num_cells(), 1);

        // All points should be at distance 0.5 from center
        for i in 0..pd.points.len() {
            let p = pd.points.get(i);
            let dist = (p[0] * p[0] + p[1] * p[1] + p[2] * p[2]).sqrt();
            assert!((dist - 0.5).abs() < 1e-10, "point {} dist = {}", i, dist);
        }
    }

    #[test]
    fn full_circle() {
        let pd = arc(&ArcParams {
            start_angle: 0.0,
            end_angle: 360.0,
            resolution: 37,
            radius: 1.0,
            ..Default::default()
        });
        assert_eq!(pd.points.len(), 37);
        // First and last points should be nearly the same
        let p0 = pd.points.get(0);
        let pn = pd.points.get(36);
        assert!((p0[0] - pn[0]).abs() < 1e-10);
        assert!((p0[1] - pn[1]).abs() < 1e-10);
    }
}
