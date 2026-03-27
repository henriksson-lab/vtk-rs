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

    /// Center of the bounding box.
    fn center(&self) -> [f64; 3] {
        self.bounds().center()
    }

    /// Diagonal length of the bounding box.
    fn diagonal(&self) -> f64 {
        self.bounds().diagonal_length()
    }

    /// Whether the dataset has no points.
    fn is_empty(&self) -> bool {
        self.num_points() == 0
    }

    /// Number of point data arrays.
    fn num_point_arrays(&self) -> usize {
        self.point_data().num_arrays()
    }

    /// Number of cell data arrays.
    fn num_cell_arrays(&self) -> usize {
        self.cell_data().num_arrays()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::PolyData;

    #[test]
    fn dataset_default_methods() {
        let pd = PolyData::from_triangles(
            vec![[0.0, 0.0, 0.0], [2.0, 0.0, 0.0], [0.0, 2.0, 0.0]],
            vec![[0, 1, 2]],
        );
        let c = pd.center();
        assert!((c[0] - 1.0).abs() < 1e-10);
        assert!(pd.diagonal() > 0.0);
        assert!(!pd.is_empty());
        assert_eq!(pd.num_point_arrays(), 0);
    }

    #[test]
    fn empty_dataset() {
        let pd = PolyData::new();
        assert!(pd.is_empty());
    }
}
