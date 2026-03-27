use vtk_data::{CellArray, DataArray, Points, PolyData};

/// Parameters for generating a flat 2D rectangular grid (triangulated).
pub struct Grid2dParams {
    /// Size in X direction. Default: 1.0
    pub x_size: f64,
    /// Size in Y direction. Default: 1.0
    pub y_size: f64,
    /// Number of divisions in X. Default: 10
    pub x_divisions: usize,
    /// Number of divisions in Y. Default: 10
    pub y_divisions: usize,
    /// Center of the grid. Default: [0, 0, 0]
    pub center: [f64; 3],
}

impl Default for Grid2dParams {
    fn default() -> Self {
        Self {
            x_size: 1.0,
            y_size: 1.0,
            x_divisions: 10,
            y_divisions: 10,
            center: [0.0, 0.0, 0.0],
        }
    }
}

/// Generate a flat 2D rectangular grid as triangulated PolyData with normals.
///
/// The grid lies in the XY plane at `center[2]`, spanning from
/// `center - size/2` to `center + size/2` in X and Y.
pub fn grid_2d(params: &Grid2dParams) -> PolyData {
    let nx = params.x_divisions.max(1);
    let ny = params.y_divisions.max(1);
    let [cx, cy, cz] = params.center;
    let x0 = cx - params.x_size / 2.0;
    let y0 = cy - params.y_size / 2.0;
    let dx = params.x_size / nx as f64;
    let dy = params.y_size / ny as f64;

    let mut points = Points::new();
    let mut normals = DataArray::<f64>::new("Normals", 3);
    let mut polys = CellArray::new();

    // Generate (nx+1) * (ny+1) grid points.
    for j in 0..=ny {
        for i in 0..=nx {
            points.push([x0 + i as f64 * dx, y0 + j as f64 * dy, cz]);
            normals.push_tuple(&[0.0, 0.0, 1.0]);
        }
    }

    // Triangulate each cell: two triangles per quad.
    let cols = nx + 1;
    for j in 0..ny {
        for i in 0..nx {
            let bl = (j * cols + i) as i64;
            let br = bl + 1;
            let tl = ((j + 1) * cols + i) as i64;
            let tr = tl + 1;

            polys.push_cell(&[bl, br, tr]);
            polys.push_cell(&[bl, tr, tl]);
        }
    }

    let mut pd = PolyData::new();
    pd.points = points;
    pd.polys = polys;
    pd.point_data_mut().add_array(normals.into());
    pd.point_data_mut().set_active_normals("Normals");
    pd
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn default_grid_2d() {
        let pd = grid_2d(&Grid2dParams::default());
        // (10+1) * (10+1) = 121 points
        assert_eq!(pd.points.len(), 121);
        // 10 * 10 * 2 = 200 triangles
        assert_eq!(pd.polys.num_cells(), 200);
        assert!(pd.point_data().normals().is_some());
    }

    #[test]
    fn minimal_grid_2d() {
        let pd = grid_2d(&Grid2dParams {
            x_divisions: 1,
            y_divisions: 1,
            ..Default::default()
        });
        assert_eq!(pd.points.len(), 4);
        assert_eq!(pd.polys.num_cells(), 2);
    }

    #[test]
    fn grid_is_flat() {
        let pd = grid_2d(&Grid2dParams {
            center: [0.0, 0.0, 5.0],
            ..Default::default()
        });
        for i in 0..pd.points.len() {
            assert!((pd.points.get(i)[2] - 5.0).abs() < 1e-12);
        }
    }

    #[test]
    fn grid_bounds() {
        let pd = grid_2d(&Grid2dParams {
            x_size: 2.0,
            y_size: 4.0,
            center: [1.0, 2.0, 0.0],
            ..Default::default()
        });
        let xs: Vec<f64> = (0..pd.points.len()).map(|i| pd.points.get(i)[0]).collect();
        let ys: Vec<f64> = (0..pd.points.len()).map(|i| pd.points.get(i)[1]).collect();
        let xmin = xs.iter().cloned().fold(f64::INFINITY, f64::min);
        let xmax = xs.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
        let ymin = ys.iter().cloned().fold(f64::INFINITY, f64::min);
        let ymax = ys.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
        assert!((xmin - 0.0).abs() < 1e-12);
        assert!((xmax - 2.0).abs() < 1e-12);
        assert!((ymin - 0.0).abs() < 1e-12);
        assert!((ymax - 4.0).abs() < 1e-12);
    }
}
