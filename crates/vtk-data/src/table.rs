use crate::{AnyDataArray, FieldData};
use crate::traits::DataObject;

/// Columnar data table (rows × named columns).
///
/// Analogous to VTK's `vtkTable`. Each column is an `AnyDataArray` where
/// the number of tuples is the number of rows.
#[derive(Debug, Clone, Default)]
pub struct Table {
    columns: Vec<AnyDataArray>,
    field_data: FieldData,
}

impl Table {
    pub fn new() -> Self {
        Self::default()
    }

    /// Add a column. All columns must have the same number of tuples (rows).
    pub fn add_column(&mut self, column: AnyDataArray) {
        if !self.columns.is_empty() {
            assert_eq!(
                column.num_tuples(),
                self.num_rows(),
                "column '{}' has {} rows, expected {}",
                column.name(),
                column.num_tuples(),
                self.num_rows()
            );
        }
        self.columns.push(column);
    }

    /// Number of rows.
    pub fn num_rows(&self) -> usize {
        self.columns.first().map(|c| c.num_tuples()).unwrap_or(0)
    }

    /// Number of columns.
    pub fn num_columns(&self) -> usize {
        self.columns.len()
    }

    /// Get column by index.
    pub fn column(&self, idx: usize) -> Option<&AnyDataArray> {
        self.columns.get(idx)
    }

    /// Get column by name.
    pub fn column_by_name(&self, name: &str) -> Option<&AnyDataArray> {
        self.columns.iter().find(|c| c.name() == name)
    }

    /// Get column names.
    pub fn column_names(&self) -> Vec<&str> {
        self.columns.iter().map(|c| c.name()).collect()
    }

    /// Iterate over columns.
    pub fn columns(&self) -> &[AnyDataArray] {
        &self.columns
    }
}

impl DataObject for Table {
    fn field_data(&self) -> &FieldData {
        &self.field_data
    }

    fn field_data_mut(&mut self) -> &mut FieldData {
        &mut self.field_data
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::DataArray;

    #[test]
    fn basic_table() {
        let mut table = Table::new();
        table.add_column(AnyDataArray::F64(DataArray::from_vec("x", vec![1.0, 2.0, 3.0], 1)));
        table.add_column(AnyDataArray::F64(DataArray::from_vec("y", vec![4.0, 5.0, 6.0], 1)));

        assert_eq!(table.num_rows(), 3);
        assert_eq!(table.num_columns(), 2);
        assert_eq!(table.column_names(), vec!["x", "y"]);
    }

    #[test]
    fn column_lookup() {
        let mut table = Table::new();
        table.add_column(AnyDataArray::I32(DataArray::from_vec("ids", vec![10, 20, 30], 1)));
        table.add_column(AnyDataArray::F64(DataArray::from_vec("values", vec![1.1, 2.2, 3.3], 1)));

        let col = table.column_by_name("values").unwrap();
        assert_eq!(col.num_tuples(), 3);
        let mut buf = [0.0f64];
        col.tuple_as_f64(1, &mut buf);
        assert!((buf[0] - 2.2).abs() < 1e-10);
    }

    #[test]
    fn multicomponent_columns() {
        let mut table = Table::new();
        table.add_column(AnyDataArray::F64(DataArray::from_vec(
            "position",
            vec![0.0, 0.0, 0.0, 1.0, 2.0, 3.0],
            3,
        )));
        assert_eq!(table.num_rows(), 2);
        assert_eq!(table.column(0).unwrap().num_components(), 3);
    }
}
