//! Compile-time assertions that key types are Send + Sync.
//! If any of these types stop being thread-safe, this file will fail to compile.

fn assert_send_sync<T: Send + Sync>() {}

#[test]
fn data_types_are_send_sync() {
    assert_send_sync::<vtk_data::PolyData>();
    assert_send_sync::<vtk_data::ImageData>();
    assert_send_sync::<vtk_data::UnstructuredGrid>();
    assert_send_sync::<vtk_data::RectilinearGrid>();
    assert_send_sync::<vtk_data::StructuredGrid>();
    assert_send_sync::<vtk_data::MultiBlockDataSet>();
    assert_send_sync::<vtk_data::Table>();
    assert_send_sync::<vtk_data::DataArray<f64>>();
    assert_send_sync::<vtk_data::AnyDataArray>();
    assert_send_sync::<vtk_data::CellArray>();
    assert_send_sync::<vtk_data::Points<f64>>();
    assert_send_sync::<vtk_data::FieldData>();
    assert_send_sync::<vtk_data::DataSetAttributes>();
    assert_send_sync::<vtk_data::KdTree>();
    assert_send_sync::<vtk_data::Selection>();
    assert_send_sync::<vtk_data::Graph>();
    assert_send_sync::<vtk_data::Tree>();
    assert_send_sync::<vtk_data::Molecule>();
}

#[test]
fn error_type_is_send_sync() {
    assert_send_sync::<vtk_types::VtkError>();
    assert_send_sync::<vtk_types::BoundingBox>();
    assert_send_sync::<vtk_types::CellType>();
}
