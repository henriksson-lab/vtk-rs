use vtk_data::{AnyDataArray, DataArray, PolyData};

/// Interpolate source data at probe point locations.
///
/// For each point in `probe`, finds the nearest point in `source` and
/// copies its scalar data. Uses a brute-force nearest-neighbor search.
pub fn probe(source: &PolyData, probe_points: &PolyData) -> PolyData {
    let mut pd = probe_points.clone();
    let n_source = source.points.len();
    let n_probe = probe_points.points.len();

    if n_source == 0 || n_probe == 0 {
        return pd;
    }

    // For each scalar array in source point data, interpolate at probe points
    for arr_idx in 0..source.point_data().num_arrays() {
        let arr = match source.point_data().get_array_by_index(arr_idx) {
            Some(a) => a,
            None => continue,
        };

        let nc = arr.num_components();
        let mut out_data = vec![0.0f64; n_probe * nc];

        for pi in 0..n_probe {
            let pp = probe_points.points.get(pi);

            // Find nearest source point (brute force)
            let mut best_dist = f64::MAX;
            let mut best_idx = 0;
            for si in 0..n_source {
                let sp = source.points.get(si);
                let d = (pp[0] - sp[0]) * (pp[0] - sp[0])
                    + (pp[1] - sp[1]) * (pp[1] - sp[1])
                    + (pp[2] - sp[2]) * (pp[2] - sp[2]);
                if d < best_dist {
                    best_dist = d;
                    best_idx = si;
                }
            }

            let mut buf = vec![0.0f64; nc];
            arr.tuple_as_f64(best_idx, &mut buf);
            for c in 0..nc {
                out_data[pi * nc + c] = buf[c];
            }
        }

        let out_arr = AnyDataArray::F64(DataArray::from_vec(arr.name(), out_data, nc));
        pd.point_data_mut().add_array(out_arr);
    }

    pd
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn probe_nearest_neighbor() {
        let mut source = PolyData::new();
        source.points.push([0.0, 0.0, 0.0]);
        source.points.push([1.0, 0.0, 0.0]);
        source.points.push([2.0, 0.0, 0.0]);
        let scalars = DataArray::from_vec("temp", vec![10.0, 20.0, 30.0], 1);
        source.point_data_mut().add_array(scalars.into());

        let mut probe_pts = PolyData::new();
        probe_pts.points.push([0.1, 0.0, 0.0]); // nearest to point 0
        probe_pts.points.push([1.6, 0.0, 0.0]); // nearest to point 2

        let result = probe(&source, &probe_pts);
        let arr = result.point_data().get_array("temp").unwrap();

        let mut val = [0.0f64];
        arr.tuple_as_f64(0, &mut val);
        assert!((val[0] - 10.0).abs() < 1e-10);
        arr.tuple_as_f64(1, &mut val);
        assert!((val[0] - 30.0).abs() < 1e-10);
    }

    #[test]
    fn probe_multicomponent() {
        let mut source = PolyData::new();
        source.points.push([0.0, 0.0, 0.0]);
        source.points.push([1.0, 0.0, 0.0]);
        let vecs = DataArray::from_vec("vel", vec![1.0, 0.0, 0.0, 0.0, 1.0, 0.0], 3);
        source.point_data_mut().add_array(vecs.into());

        let mut probe_pts = PolyData::new();
        probe_pts.points.push([0.9, 0.0, 0.0]); // nearest to point 1

        let result = probe(&source, &probe_pts);
        let arr = result.point_data().get_array("vel").unwrap();
        let mut val = [0.0f64; 3];
        arr.tuple_as_f64(0, &mut val);
        assert!((val[1] - 1.0).abs() < 1e-10);
    }
}
