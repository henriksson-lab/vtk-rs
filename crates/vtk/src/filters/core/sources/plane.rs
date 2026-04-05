use crate::data::{CellArray, DataArray, Points, PolyData};

/// Parameters for generating a rectangular plane.
pub struct PlaneParams {
    pub origin: [f64; 3],
    pub point1: [f64; 3],
    pub point2: [f64; 3],
    pub x_resolution: usize,
    pub y_resolution: usize,
}

impl Default for PlaneParams {
    fn default() -> Self {
        Self {
            origin: [-0.5, -0.5, 0.0],
            point1: [0.5, -0.5, 0.0],
            point2: [-0.5, 0.5, 0.0],
            x_resolution: 1,
            y_resolution: 1,
        }
    }
}

/// Generate a rectangular plane as PolyData with normals and texture coordinates.
pub fn plane(params: &PlaneParams) -> PolyData {
    let nx = params.x_resolution.max(1);
    let ny = params.y_resolution.max(1);

    let o = params.origin;
    // Axes
    let ax = [
        params.point1[0] - o[0],
        params.point1[1] - o[1],
        params.point1[2] - o[2],
    ];
    let ay = [
        params.point2[0] - o[0],
        params.point2[1] - o[1],
        params.point2[2] - o[2],
    ];

    // Normal = ax x ay, normalized
    let n = [
        ax[1] * ay[2] - ax[2] * ay[1],
        ax[2] * ay[0] - ax[0] * ay[2],
        ax[0] * ay[1] - ax[1] * ay[0],
    ];
    let nlen = (n[0] * n[0] + n[1] * n[1] + n[2] * n[2]).sqrt();
    let normal = if nlen > 1e-10 {
        [n[0] / nlen, n[1] / nlen, n[2] / nlen]
    } else {
        [0.0, 0.0, 1.0]
    };

    let mut points = Points::new();
    let mut normals = DataArray::<f64>::new("Normals", 3);
    let mut tcoords = DataArray::<f64>::new("TCoords", 2);
    let mut polys = CellArray::new();

    // Generate (nx+1) x (ny+1) points
    for j in 0..=ny {
        let v = j as f64 / ny as f64;
        for i in 0..=nx {
            let u = i as f64 / nx as f64;
            points.push([
                o[0] + u * ax[0] + v * ay[0],
                o[1] + u * ax[1] + v * ay[1],
                o[2] + u * ax[2] + v * ay[2],
            ]);
            normals.push_tuple(&normal);
            tcoords.push_tuple(&[u, v]);
        }
    }

    // Generate quads
    let row_size = nx + 1;
    for j in 0..ny {
        for i in 0..nx {
            let p0 = (j * row_size + i) as i64;
            let p1 = p0 + 1;
            let p2 = p0 + row_size as i64 + 1;
            let p3 = p0 + row_size as i64;
            polys.push_cell(&[p0, p1, p2, p3]);
        }
    }

    let mut pd = PolyData::new();
    pd.points = points;
    pd.polys = polys;
    pd.point_data_mut().add_array(normals.into());
    pd.point_data_mut().set_active_normals("Normals");
    pd.point_data_mut().add_array(tcoords.into());
    pd.point_data_mut().set_active_tcoords("TCoords");
    pd
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn default_plane() {
        let pd = plane(&PlaneParams::default());
        assert_eq!(pd.points.len(), 4); // 2x2 grid
        assert_eq!(pd.polys.num_cells(), 1); // 1 quad
    }

    #[test]
    fn subdivided_plane() {
        let pd = plane(&PlaneParams {
            x_resolution: 3,
            y_resolution: 2,
            ..Default::default()
        });
        assert_eq!(pd.points.len(), 4 * 3); // 4 cols * 3 rows
        assert_eq!(pd.polys.num_cells(), 3 * 2); // 3x2 quads
    }
}
