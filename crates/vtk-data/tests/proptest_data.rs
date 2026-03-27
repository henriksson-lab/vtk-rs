use proptest::prelude::*;
use vtk_data::{DataArray, CellArray, Points, PolyData};

proptest! {
    /// DataArray from_vec/num_tuples consistency.
    #[test]
    fn data_array_tuple_count(
        values in prop::collection::vec(any::<f64>(), 0..200),
        nc in 1usize..5,
    ) {
        let trimmed_len = (values.len() / nc) * nc;
        let trimmed: Vec<f64> = values[..trimmed_len].to_vec();
        let arr = DataArray::<f64>::from_vec("test", trimmed.clone(), nc);
        prop_assert_eq!(arr.num_tuples(), trimmed_len / nc);
        prop_assert_eq!(arr.num_components(), nc);
    }

    /// DataArray tuple read-back matches input.
    #[test]
    fn data_array_roundtrip(
        values in prop::collection::vec(-1e6f64..1e6, 3..60),
    ) {
        let nc = 3;
        let trimmed_len = (values.len() / nc) * nc;
        let trimmed: Vec<f64> = values[..trimmed_len].to_vec();
        let arr = DataArray::<f64>::from_vec("test", trimmed.clone(), nc);
        for i in 0..arr.num_tuples() {
            let t = arr.tuple(i);
            for c in 0..nc {
                prop_assert!((t[c] - trimmed[i * nc + c]).abs() < 1e-12);
            }
        }
    }

    /// CellArray push/iter roundtrip.
    #[test]
    fn cell_array_roundtrip(
        cells in prop::collection::vec(
            prop::collection::vec(0i64..100, 1..6),
            1..20
        ),
    ) {
        let mut ca = CellArray::new();
        for cell in &cells {
            ca.push_cell(cell);
        }
        prop_assert_eq!(ca.num_cells(), cells.len());
        for (i, cell) in ca.iter().enumerate() {
            prop_assert_eq!(cell, cells[i].as_slice());
        }
    }

    /// Points push/get roundtrip.
    #[test]
    fn points_roundtrip(
        pts in prop::collection::vec(
            (any::<f64>(), any::<f64>(), any::<f64>()),
            1..50
        ),
    ) {
        let mut points = Points::<f64>::new();
        for &(x, y, z) in &pts {
            if x.is_finite() && y.is_finite() && z.is_finite() {
                points.push([x, y, z]);
            }
        }
        let finite_count = pts.iter().filter(|(x, y, z)| x.is_finite() && y.is_finite() && z.is_finite()).count();
        prop_assert_eq!(points.len(), finite_count);
    }

    /// PolyData from_triangles produces correct cell count.
    #[test]
    fn poly_data_triangle_count(n_tris in 1usize..20) {
        let n_pts = n_tris * 3;
        let mut pts = Vec::new();
        let mut tris = Vec::new();
        for i in 0..n_pts {
            pts.push([i as f64, 0.0, 0.0]);
        }
        for i in 0..n_tris {
            let base = (i * 3) as i64;
            tris.push([base, base + 1, base + 2]);
        }
        let pd = PolyData::from_triangles(pts, tris);
        prop_assert_eq!(pd.polys.num_cells(), n_tris);
        prop_assert_eq!(pd.points.len(), n_pts);
    }

    /// DataArray name is preserved.
    #[test]
    fn data_array_name(name in "[a-zA-Z][a-zA-Z0-9_]{0,20}") {
        let arr = DataArray::<f64>::new(&name, 1);
        prop_assert_eq!(arr.name(), name.as_str());
    }

    /// DataArray statistics are consistent.
    #[test]
    fn data_array_statistics(
        values in prop::collection::vec(-1000.0f64..1000.0, 1..100),
    ) {
        let arr = vtk_data::AnyDataArray::F64(DataArray::from_vec("test", values.clone(), 1));
        let stats = arr.statistics().unwrap();
        prop_assert!(stats.min <= stats.max);
        prop_assert!(stats.mean >= stats.min);
        prop_assert!(stats.mean <= stats.max);
        prop_assert!(stats.variance >= 0.0);
        prop_assert_eq!(stats.count, values.len());
    }

    /// CellArray PartialEq is reflexive.
    #[test]
    fn cell_array_eq(
        cells in prop::collection::vec(
            prop::collection::vec(0i64..50, 1..5),
            1..10
        ),
    ) {
        let mut ca = CellArray::new();
        for cell in &cells {
            ca.push_cell(cell);
        }
        let ca2 = ca.clone();
        prop_assert_eq!(&ca, &ca2);
    }

    /// PolyData append preserves total point count.
    #[test]
    fn poly_data_append_point_count(
        n1 in 1usize..10,
        n2 in 1usize..10,
    ) {
        let pts1: Vec<[f64; 3]> = (0..n1*3).map(|i| [i as f64, 0.0, 0.0]).collect();
        let tris1: Vec<[i64; 3]> = (0..n1).map(|i| [(i*3) as i64, (i*3+1) as i64, (i*3+2) as i64]).collect();
        let pts2: Vec<[f64; 3]> = (0..n2*3).map(|i| [i as f64, 1.0, 0.0]).collect();
        let tris2: Vec<[i64; 3]> = (0..n2).map(|i| [(i*3) as i64, (i*3+1) as i64, (i*3+2) as i64]).collect();
        let pd1 = PolyData::from_triangles(pts1, tris1);
        let pd2 = PolyData::from_triangles(pts2, tris2);
        let mut merged = pd1.clone();
        merged.append(&pd2);
        prop_assert_eq!(merged.points.len(), n1*3 + n2*3);
        prop_assert_eq!(merged.polys.num_cells(), n1 + n2);
    }

    /// DataArray from_fn produces correct length.
    #[test]
    fn data_array_from_fn(count in 0usize..100) {
        let arr = DataArray::<f64>::from_fn("test", count, |i| i as f64);
        prop_assert_eq!(arr.num_tuples(), count);
    }
}
