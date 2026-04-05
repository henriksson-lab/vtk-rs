use std::f64::consts::PI;

use crate::data::{CellArray, DataArray, Points, PolyData};

/// Parameters for generating a spiral (helix projected onto a cone).
pub struct SpiralParams {
    pub center: [f64; 3],
    pub turns: f64,
    pub height: f64,
    pub base_radius: f64,
    pub top_radius: f64,
    pub resolution: usize,
}

impl Default for SpiralParams {
    fn default() -> Self {
        Self {
            center: [0.0, 0.0, 0.0],
            turns: 5.0,
            height: 1.0,
            base_radius: 0.5,
            top_radius: 0.0,
            resolution: 200,
        }
    }
}

/// Generate a spiral (helix on a cone) as PolyData with a polyline.
///
/// The spiral winds around the Z axis, starting at the base (z=0) and ending
/// at the top (z=height). The radius linearly interpolates from `base_radius`
/// at the bottom to `top_radius` at the top.
pub fn spiral(params: &SpiralParams) -> PolyData {
    let n = params.resolution.max(3);
    let [cx, cy, cz] = params.center;

    let mut points = Points::new();
    let mut normals_data = DataArray::<f64>::new("Normals", 3);
    let mut lines = CellArray::new();

    let total_angle = 2.0 * PI * params.turns;

    let mut line_ids: Vec<i64> = Vec::with_capacity(n + 1);

    for i in 0..=n {
        let t = i as f64 / n as f64;
        let angle = total_angle * t;
        let z = params.height * t;
        let radius = params.base_radius + (params.top_radius - params.base_radius) * t;

        let x = radius * angle.cos();
        let y = radius * angle.sin();

        points.push([cx + x, cy + y, cz + z]);

        // Normal points radially outward from the axis
        let len = (x * x + y * y).sqrt().max(1e-12);
        normals_data.push_tuple(&[x / len, y / len, 0.0]);

        line_ids.push(i as i64);
    }

    lines.push_cell(&line_ids);

    let mut pd = PolyData::new();
    pd.points = points;
    pd.lines = lines;
    pd.point_data_mut().add_array(normals_data.into());
    pd.point_data_mut().set_active_normals("Normals");
    pd
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn default_spiral() {
        let pd = spiral(&SpiralParams::default());
        assert_eq!(pd.points.len(), 201);
        assert_eq!(pd.lines.num_cells(), 1);
        assert!(pd.point_data().normals().is_some());
    }

    #[test]
    fn minimal_spiral() {
        let pd = spiral(&SpiralParams {
            resolution: 3,
            ..Default::default()
        });
        assert_eq!(pd.points.len(), 4);
        assert_eq!(pd.lines.num_cells(), 1);
    }

    #[test]
    fn cone_spiral() {
        let pd = spiral(&SpiralParams {
            turns: 3.0,
            height: 2.0,
            base_radius: 1.0,
            top_radius: 0.1,
            resolution: 100,
            ..Default::default()
        });
        assert_eq!(pd.points.len(), 101);
    }
}
