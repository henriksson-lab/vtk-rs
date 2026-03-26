use vtk_data::{CellArray, ImageData, Points, PolyData};

use crate::marching_cubes::{EDGE_TABLE, TRI_TABLE};

/// Flying Edges 3D — an efficient marching cubes variant.
///
/// Processes the scalar field on an ImageData grid in scanline order,
/// producing an isosurface at the given isovalue. Uses the same lookup
/// tables as marching cubes but processes in a cache-friendly pattern.
pub fn flying_edges_3d(image: &ImageData, scalars: &[f64], isovalue: f64) -> PolyData {
    let dims = image.dimensions();
    if dims[0] < 2 || dims[1] < 2 || dims[2] < 2 {
        return PolyData::new();
    }

    let nx = dims[0];
    let ny = dims[1];
    let nz = dims[2];

    let idx = |i: usize, j: usize, k: usize| -> usize { k * nx * ny + j * nx + i };

    let mut out_points = Points::<f64>::new();
    let mut out_polys = CellArray::new();

    let edge_endpoints: [(usize, usize); 12] = [
        (0, 1), (1, 2), (2, 3), (3, 0),
        (4, 5), (5, 6), (6, 7), (7, 4),
        (0, 4), (1, 5), (2, 6), (3, 7),
    ];

    for k in 0..nz - 1 {
        for j in 0..ny - 1 {
            for i in 0..nx - 1 {
                let v = [
                    scalars[idx(i, j, k)],
                    scalars[idx(i + 1, j, k)],
                    scalars[idx(i + 1, j + 1, k)],
                    scalars[idx(i, j + 1, k)],
                    scalars[idx(i, j, k + 1)],
                    scalars[idx(i + 1, j, k + 1)],
                    scalars[idx(i + 1, j + 1, k + 1)],
                    scalars[idx(i, j + 1, k + 1)],
                ];

                let mut case_idx = 0u8;
                for (bit, val) in v.iter().enumerate() {
                    if *val >= isovalue {
                        case_idx |= 1 << bit;
                    }
                }

                if EDGE_TABLE[case_idx as usize] == 0 {
                    continue;
                }

                let p = [
                    image.point_from_ijk(i, j, k),
                    image.point_from_ijk(i + 1, j, k),
                    image.point_from_ijk(i + 1, j + 1, k),
                    image.point_from_ijk(i, j + 1, k),
                    image.point_from_ijk(i, j, k + 1),
                    image.point_from_ijk(i + 1, j, k + 1),
                    image.point_from_ijk(i + 1, j + 1, k + 1),
                    image.point_from_ijk(i, j + 1, k + 1),
                ];

                let mut edge_verts = [[0.0f64; 3]; 12];
                let edges = EDGE_TABLE[case_idx as usize];
                for e in 0..12 {
                    if edges & (1 << e) != 0 {
                        let (a, b) = edge_endpoints[e];
                        let t = if (v[b] - v[a]).abs() > 1e-20 {
                            (isovalue - v[a]) / (v[b] - v[a])
                        } else {
                            0.5
                        };
                        edge_verts[e] = [
                            p[a][0] + t * (p[b][0] - p[a][0]),
                            p[a][1] + t * (p[b][1] - p[a][1]),
                            p[a][2] + t * (p[b][2] - p[a][2]),
                        ];
                    }
                }

                let tri_row = &TRI_TABLE[case_idx as usize];
                let mut ti = 0;
                while ti < tri_row.len() && tri_row[ti] != -1 {
                    let base = out_points.len() as i64;
                    for off in 0..3 {
                        let edge = tri_row[ti + off] as usize;
                        out_points.push(edge_verts[edge]);
                    }
                    out_polys.push_cell(&[base, base + 1, base + 2]);
                    ti += 3;
                }
            }
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
    use vtk_data::DataSet;

    #[test]
    fn sphere_isosurface() {
        let mut image = ImageData::with_dimensions(10, 10, 10);
        image.set_spacing([0.2, 0.2, 0.2]);
        image.set_origin([-0.9, -0.9, -0.9]);

        let n = image.num_points();
        let scalars: Vec<f64> = (0..n)
            .map(|i| {
                let p = image.point(i);
                p[0] * p[0] + p[1] * p[1] + p[2] * p[2]
            })
            .collect();

        let result = flying_edges_3d(&image, &scalars, 0.5);
        assert!(result.polys.num_cells() > 10, "got {}", result.polys.num_cells());
    }

    #[test]
    fn empty_for_no_crossings() {
        let image = ImageData::with_dimensions(3, 3, 3);
        let scalars = vec![1.0; 27];
        let result = flying_edges_3d(&image, &scalars, 0.5);
        assert_eq!(result.polys.num_cells(), 0);
    }

    #[test]
    fn matches_marching_cubes() {
        let mut image = ImageData::with_dimensions(5, 5, 5);
        image.set_spacing([0.5, 0.5, 0.5]);
        image.set_origin([-1.0, -1.0, -1.0]);

        let n = image.num_points();
        let scalars: Vec<f64> = (0..n)
            .map(|i| {
                let p = image.point(i);
                p[0] * p[0] + p[1] * p[1] + p[2] * p[2]
            })
            .collect();

        let fe = flying_edges_3d(&image, &scalars, 0.5);
        let mc = crate::marching_cubes::marching_cubes(&image, &scalars, 0.5);
        assert_eq!(fe.polys.num_cells(), mc.polys.num_cells());
    }
}
