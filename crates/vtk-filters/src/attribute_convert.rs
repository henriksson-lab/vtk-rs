use vtk_data::{AnyDataArray, DataArray, PolyData};

/// Convert cell data arrays to point data by averaging values from cells
/// that share each point.
///
/// For each point, the output value is the average of values from all
/// cells that contain that point.
pub fn cell_data_to_point_data(input: &PolyData) -> PolyData {
    let n_points = input.points.len();
    let mut pd = input.clone();

    // Clear point data — we'll rebuild from cell data
    // Keep existing point data arrays, add converted ones

    for arr_idx in 0..input.cell_data().num_arrays() {
        let arr = match input.cell_data().get_array_by_index(arr_idx) {
            Some(a) => a,
            None => continue,
        };

        let nc = arr.num_components();
        let mut sums = vec![0.0f64; n_points * nc];
        let mut counts = vec![0usize; n_points];

        let mut cell_idx = 0;
        let process_cells = |cells: &vtk_data::CellArray,
                                  sums: &mut Vec<f64>,
                                  counts: &mut Vec<usize>,
                                  cell_idx: &mut usize,
                                  arr: &AnyDataArray| {
            for cell in cells.iter() {
                if *cell_idx >= arr.num_tuples() {
                    break;
                }
                let mut buf = vec![0.0f64; nc];
                arr.tuple_as_f64(*cell_idx, &mut buf);
                for &pt_id in cell {
                    let pi = pt_id as usize;
                    for c in 0..nc {
                        sums[pi * nc + c] += buf[c];
                    }
                    counts[pi] += 1;
                }
                *cell_idx += 1;
            }
        };

        process_cells(&input.verts, &mut sums, &mut counts, &mut cell_idx, arr);
        process_cells(&input.lines, &mut sums, &mut counts, &mut cell_idx, arr);
        process_cells(&input.polys, &mut sums, &mut counts, &mut cell_idx, arr);
        process_cells(&input.strips, &mut sums, &mut counts, &mut cell_idx, arr);

        // Average
        for i in 0..n_points {
            if counts[i] > 0 {
                for c in 0..nc {
                    sums[i * nc + c] /= counts[i] as f64;
                }
            }
        }

        let out_arr = AnyDataArray::F64(DataArray::from_vec(arr.name(), sums, nc));
        pd.point_data_mut().add_array(out_arr);
    }

    pd
}

/// Convert point data arrays to cell data by averaging point values for
/// each cell.
///
/// For each cell, the output value is the average of the point data values
/// at the cell's vertices.
pub fn point_data_to_cell_data(input: &PolyData) -> PolyData {
    let mut pd = input.clone();

    for arr_idx in 0..input.point_data().num_arrays() {
        let arr = match input.point_data().get_array_by_index(arr_idx) {
            Some(a) => a,
            None => continue,
        };

        let nc = arr.num_components();
        let mut cell_values: Vec<f64> = Vec::new();

        let process_cells =
            |cells: &vtk_data::CellArray, values: &mut Vec<f64>, arr: &AnyDataArray| {
                let mut buf = vec![0.0f64; nc];
                for cell in cells.iter() {
                    let mut avg = vec![0.0f64; nc];
                    for &pt_id in cell {
                        arr.tuple_as_f64(pt_id as usize, &mut buf);
                        for c in 0..nc {
                            avg[c] += buf[c];
                        }
                    }
                    let n = cell.len() as f64;
                    for v in &mut avg {
                        *v /= n;
                    }
                    values.extend_from_slice(&avg);
                }
            };

        process_cells(&input.verts, &mut cell_values, arr);
        process_cells(&input.lines, &mut cell_values, arr);
        process_cells(&input.polys, &mut cell_values, arr);
        process_cells(&input.strips, &mut cell_values, arr);

        let out_arr = AnyDataArray::F64(DataArray::from_vec(arr.name(), cell_values, nc));
        pd.cell_data_mut().add_array(out_arr);
    }

    pd
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn cell_to_point_averaging() {
        let mut pd = PolyData::from_triangles(
            vec![
                [0.0, 0.0, 0.0],
                [1.0, 0.0, 0.0],
                [0.5, 1.0, 0.0],
                [1.5, 1.0, 0.0],
            ],
            vec![[0, 1, 2], [1, 3, 2]],
        );
        // Cell data: cell 0 = 10.0, cell 1 = 20.0
        let cell_scalars = DataArray::from_vec("pressure", vec![10.0, 20.0], 1);
        pd.cell_data_mut().add_array(cell_scalars.into());

        let result = cell_data_to_point_data(&pd);
        let pt_arr = result.point_data().get_array("pressure").unwrap();
        assert_eq!(pt_arr.num_tuples(), 4);

        let mut val = [0.0f64];
        // Point 0: only in cell 0 -> 10.0
        pt_arr.tuple_as_f64(0, &mut val);
        assert!((val[0] - 10.0).abs() < 1e-10);
        // Point 1: in cell 0 (10) and cell 1 (20) -> 15.0
        pt_arr.tuple_as_f64(1, &mut val);
        assert!((val[0] - 15.0).abs() < 1e-10);
        // Point 3: only in cell 1 -> 20.0
        pt_arr.tuple_as_f64(3, &mut val);
        assert!((val[0] - 20.0).abs() < 1e-10);
    }

    #[test]
    fn point_to_cell_averaging() {
        let mut pd = PolyData::from_triangles(
            vec![
                [0.0, 0.0, 0.0],
                [1.0, 0.0, 0.0],
                [0.5, 1.0, 0.0],
            ],
            vec![[0, 1, 2]],
        );
        let pt_scalars = DataArray::from_vec("temp", vec![10.0, 20.0, 30.0], 1);
        pd.point_data_mut().add_array(pt_scalars.into());

        let result = point_data_to_cell_data(&pd);
        let cell_arr = result.cell_data().get_array("temp").unwrap();
        assert_eq!(cell_arr.num_tuples(), 1);

        let mut val = [0.0f64];
        cell_arr.tuple_as_f64(0, &mut val);
        // Average of 10, 20, 30 = 20.0
        assert!((val[0] - 20.0).abs() < 1e-10);
    }

    #[test]
    fn multicomponent_conversion() {
        let mut pd = PolyData::from_triangles(
            vec![
                [0.0, 0.0, 0.0],
                [1.0, 0.0, 0.0],
                [0.5, 1.0, 0.0],
            ],
            vec![[0, 1, 2]],
        );
        let vectors = DataArray::from_vec(
            "velocity",
            vec![1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0],
            3,
        );
        pd.point_data_mut().add_array(vectors.into());

        let result = point_data_to_cell_data(&pd);
        let cell_arr = result.cell_data().get_array("velocity").unwrap();
        assert_eq!(cell_arr.num_components(), 3);

        let mut val = [0.0f64; 3];
        cell_arr.tuple_as_f64(0, &mut val);
        // Average of (1,0,0), (0,1,0), (0,0,1) = (1/3, 1/3, 1/3)
        for v in &val {
            assert!((v - 1.0 / 3.0).abs() < 1e-10);
        }
    }
}
