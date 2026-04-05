use std::f64::consts::PI;

use crate::data::{DataArray, Points, PolyData};

/// Parameters for generating a regular polygon.
pub struct RegularPolygonParams {
    /// Number of sides. Default: 6 (hexagon)
    pub number_of_sides: usize,
    /// Radius of the circumscribed circle. Default: 0.5
    pub radius: f64,
    /// Center of the polygon. Default: [0, 0, 0]
    pub center: [f64; 3],
    /// If true, generate a filled polygon. If false, generate just the outline.
    /// Default: true
    pub generate_polygon: bool,
}

impl Default for RegularPolygonParams {
    fn default() -> Self {
        Self {
            number_of_sides: 6,
            radius: 0.5,
            center: [0.0, 0.0, 0.0],
            generate_polygon: true,
        }
    }
}

/// Generate a regular polygon in the XY plane.
pub fn regular_polygon(params: &RegularPolygonParams) -> PolyData {
    let n = params.number_of_sides.max(3);
    let cx = params.center[0];
    let cy = params.center[1];
    let cz = params.center[2];

    let mut points = Points::new();
    let mut normals = DataArray::<f64>::new("Normals", 3);

    for i in 0..n {
        let angle = 2.0 * PI * i as f64 / n as f64;
        points.push([
            cx + params.radius * angle.cos(),
            cy + params.radius * angle.sin(),
            cz,
        ]);
        normals.push_tuple(&[0.0, 0.0, 1.0]);
    }

    let mut pd = PolyData::new();
    pd.points = points;

    if params.generate_polygon {
        let ids: Vec<i64> = (0..n as i64).collect();
        pd.polys.push_cell(&ids);
    } else {
        // Outline: closed polyline
        let mut ids: Vec<i64> = (0..n as i64).collect();
        ids.push(0); // close the loop
        pd.lines.push_cell(&ids);
    }

    pd.point_data_mut().add_array(normals.into());
    pd.point_data_mut().set_active_normals("Normals");
    pd
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn default_hexagon() {
        let pd = regular_polygon(&RegularPolygonParams::default());
        assert_eq!(pd.points.len(), 6);
        assert_eq!(pd.polys.num_cells(), 1);
        assert_eq!(pd.polys.cell(0).len(), 6);
    }

    #[test]
    fn triangle() {
        let pd = regular_polygon(&RegularPolygonParams {
            number_of_sides: 3,
            radius: 1.0,
            ..Default::default()
        });
        assert_eq!(pd.points.len(), 3);
        // First point should be at (1, 0, 0)
        let p = pd.points.get(0);
        assert!((p[0] - 1.0).abs() < 1e-10);
        assert!(p[1].abs() < 1e-10);
    }

    #[test]
    fn outline_mode() {
        let pd = regular_polygon(&RegularPolygonParams {
            number_of_sides: 4,
            generate_polygon: false,
            ..Default::default()
        });
        assert_eq!(pd.polys.num_cells(), 0);
        assert_eq!(pd.lines.num_cells(), 1);
        // Closed polyline: 4 + 1 = 5 point references
        assert_eq!(pd.lines.cell(0).len(), 5);
    }
}
