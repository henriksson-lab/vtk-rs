use vtk_data::{AnyDataArray, DataArray, DataSet, ImageData, PolyData};

/// Resample source data onto the grid points of a target ImageData.
///
/// For each point in the target grid, finds the nearest point in the
/// source PolyData and copies its scalar values. Produces a new
/// ImageData with the resampled scalar arrays.
pub fn resample_with_dataset(source: &PolyData, target: &ImageData) -> ImageData {
    let mut result = target.clone();
    let n_target = target.num_points();
    let n_source = source.points.len();

    if n_source == 0 || n_target == 0 {
        return result;
    }

    for arr_idx in 0..source.point_data().num_arrays() {
        let arr = match source.point_data().get_array_by_index(arr_idx) {
            Some(a) => a,
            None => continue,
        };

        let nc = arr.num_components();
        let mut out_data = vec![0.0f64; n_target * nc];

        for ti in 0..n_target {
            let tp = target.point(ti);

            // Nearest neighbor
            let mut best_dist = f64::MAX;
            let mut best_idx = 0;
            for si in 0..n_source {
                let sp = source.points.get(si);
                let d = (tp[0] - sp[0]) * (tp[0] - sp[0])
                    + (tp[1] - sp[1]) * (tp[1] - sp[1])
                    + (tp[2] - sp[2]) * (tp[2] - sp[2]);
                if d < best_dist {
                    best_dist = d;
                    best_idx = si;
                }
            }

            let mut buf = vec![0.0f64; nc];
            arr.tuple_as_f64(best_idx, &mut buf);
            for c in 0..nc {
                out_data[ti * nc + c] = buf[c];
            }
        }

        let out_arr = AnyDataArray::F64(DataArray::from_vec(arr.name(), out_data, nc));
        let name = out_arr.name().to_string();
        result.point_data_mut().add_array(out_arr);
        if result.point_data().scalars().is_none() {
            result.point_data_mut().set_active_scalars(&name);
        }
    }

    result
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn resample_scalars() {
        let mut source = PolyData::new();
        source.points.push([0.0, 0.0, 0.0]);
        source.points.push([1.0, 0.0, 0.0]);
        let scalars = DataArray::from_vec("temp", vec![10.0, 20.0], 1);
        source.point_data_mut().add_array(scalars.into());

        let target = ImageData::with_dimensions(3, 1, 1);
        // Points at x=0, x=1, x=2
        let result = resample_with_dataset(&source, &target);

        let s = result.point_data().get_array("temp").unwrap();
        assert_eq!(s.num_tuples(), 3);
        let mut val = [0.0f64];
        s.tuple_as_f64(0, &mut val);
        assert!((val[0] - 10.0).abs() < 1e-10); // nearest to (0,0,0)
        s.tuple_as_f64(1, &mut val);
        assert!((val[0] - 20.0).abs() < 1e-10); // nearest to (1,0,0)
    }
}
