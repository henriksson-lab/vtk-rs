use vtk_data::PolyData;

use super::{cone, cylinder};
use crate::sources::cone::ConeParams;
use crate::sources::cylinder::CylinderParams;

/// Parameters for generating an arrow (cylinder shaft + cone tip).
pub struct ArrowParams {
    pub shaft_radius: f64,
    pub shaft_length: f64,
    pub tip_radius: f64,
    pub tip_length: f64,
    pub resolution: usize,
}

impl Default for ArrowParams {
    fn default() -> Self {
        Self {
            shaft_radius: 0.03,
            shaft_length: 0.7,
            tip_radius: 0.1,
            tip_length: 0.3,
            resolution: 16,
        }
    }
}

/// Generate an arrow pointing along +X as PolyData.
///
/// The total length is `shaft_length + tip_length`. The arrow is centered at the origin
/// with the base at `(-total/2, 0, 0)` and tip at `(+total/2, 0, 0)`.
pub fn arrow(params: &ArrowParams) -> PolyData {
    let total = params.shaft_length + params.tip_length;

    // Shaft: cylinder along Y axis, then we'll remap to X
    let shaft = cylinder(&CylinderParams {
        center: [0.0, 0.0, 0.0],
        height: params.shaft_length,
        radius: params.shaft_radius,
        resolution: params.resolution,
        capping: true,
    });

    // Tip: cone along X axis
    let tip = cone(&ConeParams {
        center: [0.0, 0.0, 0.0],
        height: params.tip_length,
        radius: params.tip_radius,
        direction: [1.0, 0.0, 0.0],
        resolution: params.resolution,
        capping: true,
    });

    // Transform shaft from Y-axis to X-axis and offset
    let shaft_offset_x = -total * 0.5 + params.shaft_length * 0.5;
    let mut shaft_points = Vec::new();
    for i in 0..shaft.points.len() {
        let [x, y, z] = shaft.points.get(i);
        // Rotate 90° around Z: (x,y,z) → (y,-x,z) maps Y→X
        // Actually: cylinder is along Y. To make it along X: swap x and y
        shaft_points.push([y + shaft_offset_x, x, z]);
    }

    // Transform tip and offset
    let tip_offset_x = -total * 0.5 + params.shaft_length + params.tip_length * 0.5;
    let mut tip_points = Vec::new();
    for i in 0..tip.points.len() {
        let [x, y, z] = tip.points.get(i);
        tip_points.push([x + tip_offset_x, y, z]);
    }

    // Merge into single PolyData
    let shaft_n = shaft_points.len();
    let mut merged = PolyData::new();

    for p in &shaft_points {
        merged.points.push(*p);
    }
    for p in &tip_points {
        merged.points.push(*p);
    }

    // Copy shaft polys
    for cell in shaft.polys.iter() {
        merged.polys.push_cell(cell);
    }

    // Copy tip polys with offset indices
    for cell in tip.polys.iter() {
        let offset_cell: Vec<i64> = cell.iter().map(|&id| id + shaft_n as i64).collect();
        merged.polys.push_cell(&offset_cell);
    }

    merged
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn default_arrow() {
        let pd = arrow(&ArrowParams::default());
        assert!(pd.points.len() > 0);
        assert!(pd.polys.num_cells() > 0);
    }
}
