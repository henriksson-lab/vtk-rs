use crate::data::{CellArray, Points, PolyData};

/// Place a copy of a glyph mesh at each point of the input PolyData.
///
/// The glyph is translated to each input point. If `scale_by_scalar` is true
/// and active scalars exist, each glyph is uniformly scaled by the scalar value.
pub fn glyph(
    input: &PolyData,
    glyph_source: &PolyData,
    scale_factor: f64,
    scale_by_scalar: bool,
) -> PolyData {
    let n = input.points.len();
    let glyph_n_pts = glyph_source.points.len();

    if n == 0 || glyph_n_pts == 0 {
        return PolyData::new();
    }

    // Read scalar values if needed
    let scalars: Option<Vec<f64>> = if scale_by_scalar {
        input.point_data().scalars().map(|s| {
            let mut values = Vec::with_capacity(n);
            let mut buf = [0.0f64];
            for i in 0..n {
                s.tuple_as_f64(i, &mut buf);
                values.push(buf[0]);
            }
            values
        })
    } else {
        None
    };

    let mut out_points = Points::<f64>::new();
    let mut out_polys = CellArray::new();

    for i in 0..n {
        let center = input.points.get(i);
        let s = if let Some(ref vals) = scalars {
            vals[i] * scale_factor
        } else {
            scale_factor
        };

        let base = out_points.len() as i64;

        // Copy and transform glyph points
        for j in 0..glyph_n_pts {
            let gp = glyph_source.points.get(j);
            out_points.push([
                center[0] + gp[0] * s,
                center[1] + gp[1] * s,
                center[2] + gp[2] * s,
            ]);
        }

        // Copy glyph cells with offset
        for cell in glyph_source.polys.iter() {
            let remapped: Vec<i64> = cell.iter().map(|&id| id + base).collect();
            out_polys.push_cell(&remapped);
        }
    }

    let mut pd = PolyData::new();
    pd.points = out_points;
    pd.polys = out_polys;
    pd
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn glyph_at_two_points() {
        // Input: two points
        let mut input = PolyData::new();
        input.points.push([0.0, 0.0, 0.0]);
        input.points.push([5.0, 0.0, 0.0]);

        // Glyph: single triangle
        let glyph_src = PolyData::from_triangles(
            vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]],
            vec![[0, 1, 2]],
        );

        let result = glyph(&input, &glyph_src, 1.0, false);
        assert_eq!(result.points.len(), 6); // 3 per glyph * 2 points
        assert_eq!(result.polys.num_cells(), 2);

        // First glyph at origin
        assert_eq!(result.points.get(0), [0.0, 0.0, 0.0]);
        // Second glyph at (5, 0, 0)
        assert_eq!(result.points.get(3), [5.0, 0.0, 0.0]);
        assert_eq!(result.points.get(4), [6.0, 0.0, 0.0]);
    }

    #[test]
    fn glyph_with_scaling() {
        let mut input = PolyData::new();
        input.points.push([0.0, 0.0, 0.0]);

        let glyph_src = PolyData::from_triangles(
            vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]],
            vec![[0, 1, 2]],
        );

        let result = glyph(&input, &glyph_src, 2.0, false);
        // Glyph scaled by 2
        assert_eq!(result.points.get(1), [2.0, 0.0, 0.0]);
    }
}
