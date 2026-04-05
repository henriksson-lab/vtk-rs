use crate::data::{CellArray, Points, PolyData};

/// Axis for axis-aligned plane slicing.
#[derive(Debug, Clone, Copy)]
pub enum SliceAxis {
    X,
    Y,
    Z,
}

/// Slice a mesh with an axis-aligned plane, returning line segments at the intersection.
///
/// For each polygon cell, computes signed distances from vertices to the plane,
/// finds edge crossings, and emits line segments connecting crossing points.
fn slice_by_axis(input: &PolyData, axis: SliceAxis, value: f64) -> PolyData {
    let mut points = Points::<f64>::new();
    let mut lines = CellArray::new();

    let axis_idx: usize = match axis {
        SliceAxis::X => 0,
        SliceAxis::Y => 1,
        SliceAxis::Z => 2,
    };

    for cell in input.polys.iter() {
        if cell.len() < 3 {
            continue;
        }

        let dists: Vec<f64> = cell
            .iter()
            .map(|&id| {
                let p = input.points.get(id as usize);
                p[axis_idx] - value
            })
            .collect();

        let mut crossings = Vec::new();
        let n: usize = cell.len();
        for i in 0..n {
            let j: usize = (i + 1) % n;
            let di: f64 = dists[i];
            let dj: f64 = dists[j];

            if (di >= 0.0) != (dj >= 0.0) {
                let t: f64 = di / (di - dj);
                let pi = input.points.get(cell[i] as usize);
                let pj = input.points.get(cell[j] as usize);
                let intersection = [
                    pi[0] + t * (pj[0] - pi[0]),
                    pi[1] + t * (pj[1] - pi[1]),
                    pi[2] + t * (pj[2] - pi[2]),
                ];
                let idx: i64 = points.len() as i64;
                points.push(intersection);
                crossings.push(idx);
            } else if di.abs() < 1e-10 && dj.abs() >= 1e-10 {
                let idx: i64 = points.len() as i64;
                points.push(input.points.get(cell[i] as usize));
                crossings.push(idx);
            }
        }

        if crossings.len() == 2 {
            lines.push_cell(&[crossings[0], crossings[1]]);
        }
    }

    let mut pd = PolyData::new();
    pd.points = points;
    pd.lines = lines;
    pd
}

/// Slice a mesh with the plane X = constant.
pub fn slice_by_plane_x(input: &PolyData, x: f64) -> PolyData {
    slice_by_axis(input, SliceAxis::X, x)
}

/// Slice a mesh with the plane Y = constant.
pub fn slice_by_plane_y(input: &PolyData, y: f64) -> PolyData {
    slice_by_axis(input, SliceAxis::Y, y)
}

/// Slice a mesh with the plane Z = constant.
pub fn slice_by_plane_z(input: &PolyData, z: f64) -> PolyData {
    slice_by_axis(input, SliceAxis::Z, z)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn slice_z_through_triangle() {
        // Triangle from z=0 to z=2
        let pd = PolyData::from_triangles(
            vec![[0.0, 0.0, 0.0], [2.0, 0.0, 0.0], [1.0, 0.0, 2.0]],
            vec![[0, 1, 2]],
        );
        let result = slice_by_plane_z(&pd, 1.0);
        assert_eq!(result.lines.num_cells(), 1, "expected one line segment");
        assert_eq!(result.points.len(), 2);

        // Both intersection points should have z == 1.0
        let p0 = result.points.get(0);
        let p1 = result.points.get(1);
        assert!((p0[2] - 1.0).abs() < 1e-10);
        assert!((p1[2] - 1.0).abs() < 1e-10);
    }

    #[test]
    fn slice_x_misses_triangle() {
        // Triangle entirely in 0 < x < 2
        let pd = PolyData::from_triangles(
            vec![[0.0, 0.0, 0.0], [2.0, 0.0, 0.0], [1.0, 1.0, 0.0]],
            vec![[0, 1, 2]],
        );
        let result = slice_by_plane_x(&pd, 5.0);
        assert_eq!(result.lines.num_cells(), 0);
    }

    #[test]
    fn slice_y_two_triangles() {
        let pd = PolyData::from_triangles(
            vec![
                [0.0, 0.0, 0.0],
                [1.0, 0.0, 0.0],
                [0.5, 2.0, 0.0],
                [2.0, 0.0, 0.0],
                [3.0, 0.0, 0.0],
                [2.5, 2.0, 0.0],
            ],
            vec![[0, 1, 2], [3, 4, 5]],
        );
        let result = slice_by_plane_y(&pd, 1.0);
        assert_eq!(result.lines.num_cells(), 2, "expected two line segments");
    }
}
