use vtk_data::{CellArray, Points, PolyData};

/// Parameters for generating a line source.
pub struct LineParams {
    /// Start point. Default: [-0.5, 0, 0]
    pub point1: [f64; 3],
    /// End point. Default: [0.5, 0, 0]
    pub point2: [f64; 3],
    /// Number of points along the line. Default: 2 (just endpoints)
    pub resolution: usize,
}

impl Default for LineParams {
    fn default() -> Self {
        Self {
            point1: [-0.5, 0.0, 0.0],
            point2: [0.5, 0.0, 0.0],
            resolution: 2,
        }
    }
}

/// Generate a line as PolyData with a single polyline cell.
pub fn line(params: &LineParams) -> PolyData {
    let n = params.resolution.max(2);
    let mut points = Points::new();
    let mut lines = CellArray::new();

    let mut ids = Vec::with_capacity(n);
    for i in 0..n {
        let t = i as f64 / (n - 1) as f64;
        points.push([
            params.point1[0] + t * (params.point2[0] - params.point1[0]),
            params.point1[1] + t * (params.point2[1] - params.point1[1]),
            params.point1[2] + t * (params.point2[2] - params.point1[2]),
        ]);
        ids.push(i as i64);
    }
    lines.push_cell(&ids);

    let mut pd = PolyData::new();
    pd.points = points;
    pd.lines = lines;
    pd
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn default_line() {
        let pd = line(&LineParams::default());
        assert_eq!(pd.points.len(), 2);
        assert_eq!(pd.lines.num_cells(), 1);
        let p0 = pd.points.get(0);
        let p1 = pd.points.get(1);
        assert!((p0[0] + 0.5).abs() < 1e-10);
        assert!((p1[0] - 0.5).abs() < 1e-10);
    }

    #[test]
    fn subdivided_line() {
        let pd = line(&LineParams {
            point1: [0.0, 0.0, 0.0],
            point2: [10.0, 0.0, 0.0],
            resolution: 11,
        });
        assert_eq!(pd.points.len(), 11);
        assert_eq!(pd.lines.num_cells(), 1);
        assert_eq!(pd.lines.cell(0).len(), 11);
        // Midpoint should be at x=5
        let mid = pd.points.get(5);
        assert!((mid[0] - 5.0).abs() < 1e-10);
    }
}
