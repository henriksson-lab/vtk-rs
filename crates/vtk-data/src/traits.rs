use vtk_types::BoundingBox;

use crate::{DataSetAttributes, FieldData};

/// Base trait for any data object with field data.
///
/// Analogous to VTK's `vtkDataObject`.
pub trait DataObject {
    fn field_data(&self) -> &FieldData;
    fn field_data_mut(&mut self) -> &mut FieldData;
}

/// Trait for spatial datasets with points and cells.
///
/// Analogous to VTK's `vtkDataSet`.
pub trait DataSet: DataObject {
    fn num_points(&self) -> usize;
    fn num_cells(&self) -> usize;
    fn point(&self, idx: usize) -> [f64; 3];
    fn bounds(&self) -> BoundingBox;
    fn point_data(&self) -> &DataSetAttributes;
    fn point_data_mut(&mut self) -> &mut DataSetAttributes;
    fn cell_data(&self) -> &DataSetAttributes;
    fn cell_data_mut(&mut self) -> &mut DataSetAttributes;
}
