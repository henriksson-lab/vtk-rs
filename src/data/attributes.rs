use crate::data::{AnyDataArray, FieldData};

/// Field data with "active" attribute designations for scalars, vectors, normals, etc.
///
/// Analogous to VTK's `vtkDataSetAttributes`.
#[derive(Debug, Clone, Default, PartialEq)]
pub struct DataSetAttributes {
    field_data: FieldData,
    active_scalars: Option<usize>,
    active_vectors: Option<usize>,
    active_normals: Option<usize>,
    active_tcoords: Option<usize>,
    active_tensors: Option<usize>,
}

impl DataSetAttributes {
    pub fn new() -> Self {
        Self::default()
    }

    pub fn field_data(&self) -> &FieldData {
        &self.field_data
    }

    pub fn field_data_mut(&mut self) -> &mut FieldData {
        &mut self.field_data
    }

    pub fn add_array(&mut self, array: AnyDataArray) {
        self.field_data.add_array(array);
    }

    pub fn num_arrays(&self) -> usize {
        self.field_data.num_arrays()
    }

    pub fn get_array(&self, name: &str) -> Option<&AnyDataArray> {
        self.field_data.get_array(name)
    }

    pub fn get_array_by_index(&self, idx: usize) -> Option<&AnyDataArray> {
        self.field_data.get_array_by_index(idx)
    }

    // Active attribute setters/getters

    pub fn set_active_scalars(&mut self, name: &str) -> bool {
        if let Some(idx) = self.find_index(name) {
            self.active_scalars = Some(idx);
            true
        } else {
            false
        }
    }

    pub fn set_active_vectors(&mut self, name: &str) -> bool {
        if let Some(idx) = self.find_index(name) {
            self.active_vectors = Some(idx);
            true
        } else {
            false
        }
    }

    pub fn set_active_normals(&mut self, name: &str) -> bool {
        if let Some(idx) = self.find_index(name) {
            self.active_normals = Some(idx);
            true
        } else {
            false
        }
    }

    pub fn set_active_tcoords(&mut self, name: &str) -> bool {
        if let Some(idx) = self.find_index(name) {
            self.active_tcoords = Some(idx);
            true
        } else {
            false
        }
    }

    pub fn set_active_tensors(&mut self, name: &str) -> bool {
        if let Some(idx) = self.find_index(name) {
            self.active_tensors = Some(idx);
            true
        } else {
            false
        }
    }

    pub fn scalars(&self) -> Option<&AnyDataArray> {
        self.active_scalars.and_then(|i| self.field_data.get_array_by_index(i))
    }

    pub fn vectors(&self) -> Option<&AnyDataArray> {
        self.active_vectors.and_then(|i| self.field_data.get_array_by_index(i))
    }

    pub fn normals(&self) -> Option<&AnyDataArray> {
        self.active_normals.and_then(|i| self.field_data.get_array_by_index(i))
    }

    pub fn tcoords(&self) -> Option<&AnyDataArray> {
        self.active_tcoords.and_then(|i| self.field_data.get_array_by_index(i))
    }

    pub fn tensors(&self) -> Option<&AnyDataArray> {
        self.active_tensors.and_then(|i| self.field_data.get_array_by_index(i))
    }

    /// Check if an array with the given name exists.
    pub fn has_array(&self, name: &str) -> bool {
        self.field_data.has_array(name)
    }

    /// Get all array names.
    pub fn array_names(&self) -> Vec<&str> {
        self.field_data.names()
    }

    /// Remove an array by name. Adjusts active attribute indices.
    pub fn remove_array(&mut self, name: &str) -> Option<AnyDataArray> {
        let idx = self.find_index(name);
        let result = self.field_data.remove_array(name);
        if let Some(removed_idx) = idx {
            self.adjust_active_after_remove(removed_idx);
        }
        result
    }

    /// Remove all arrays and clear active attributes.
    pub fn clear(&mut self) {
        self.field_data.clear();
        self.active_scalars = None;
        self.active_vectors = None;
        self.active_normals = None;
        self.active_tcoords = None;
        self.active_tensors = None;
    }

    /// Check if any active attributes are set.
    pub fn has_active_attributes(&self) -> bool {
        self.active_scalars.is_some()
            || self.active_vectors.is_some()
            || self.active_normals.is_some()
            || self.active_tcoords.is_some()
            || self.active_tensors.is_some()
    }

    /// Iterate over all arrays.
    pub fn iter(&self) -> impl Iterator<Item = &AnyDataArray> {
        self.field_data.iter()
    }

    fn find_index(&self, name: &str) -> Option<usize> {
        (0..self.field_data.num_arrays())
            .find(|&i| self.field_data.get_array_by_index(i).map(|a| a.name()) == Some(name))
    }

    fn adjust_active_after_remove(&mut self, removed: usize) {
        fn adjust(active: &mut Option<usize>, removed: usize) {
            *active = match *active {
                Some(i) if i == removed => None,
                Some(i) if i > removed => Some(i - 1),
                other => other,
            };
        }
        adjust(&mut self.active_scalars, removed);
        adjust(&mut self.active_vectors, removed);
        adjust(&mut self.active_normals, removed);
        adjust(&mut self.active_tcoords, removed);
        adjust(&mut self.active_tensors, removed);
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::data::DataArray;

    #[test]
    fn add_and_get() {
        let mut attrs = DataSetAttributes::new();
        attrs.add_array(AnyDataArray::F64(DataArray::from_vec("temp", vec![1.0, 2.0, 3.0], 1)));
        assert_eq!(attrs.num_arrays(), 1);
        assert!(attrs.has_array("temp"));
    }

    #[test]
    fn active_scalars() {
        let mut attrs = DataSetAttributes::new();
        attrs.add_array(AnyDataArray::F64(DataArray::from_vec("temp", vec![1.0], 1)));
        attrs.set_active_scalars("temp");
        assert!(attrs.scalars().is_some());
        assert_eq!(attrs.scalars().unwrap().name(), "temp");
    }

    #[test]
    fn active_normals() {
        let mut attrs = DataSetAttributes::new();
        attrs.add_array(AnyDataArray::F64(DataArray::from_vec("N", vec![0.0, 0.0, 1.0], 3)));
        attrs.set_active_normals("N");
        assert!(attrs.normals().is_some());
    }

    #[test]
    fn remove_adjusts_active() {
        let mut attrs = DataSetAttributes::new();
        attrs.add_array(AnyDataArray::F64(DataArray::from_vec("a", vec![1.0], 1)));
        attrs.add_array(AnyDataArray::F64(DataArray::from_vec("b", vec![2.0], 1)));
        attrs.set_active_scalars("b");
        attrs.remove_array("a");
        assert!(attrs.scalars().is_some());
        assert_eq!(attrs.scalars().unwrap().name(), "b");
    }

    #[test]
    fn remove_active_clears() {
        let mut attrs = DataSetAttributes::new();
        attrs.add_array(AnyDataArray::F64(DataArray::from_vec("x", vec![1.0], 1)));
        attrs.set_active_scalars("x");
        attrs.remove_array("x");
        assert!(attrs.scalars().is_none());
    }

    #[test]
    fn array_names() {
        let mut attrs = DataSetAttributes::new();
        attrs.add_array(AnyDataArray::F64(DataArray::from_vec("a", vec![1.0], 1)));
        attrs.add_array(AnyDataArray::F64(DataArray::from_vec("b", vec![2.0], 1)));
        assert_eq!(attrs.array_names(), vec!["a", "b"]);
    }

    #[test]
    fn has_active() {
        let mut attrs = DataSetAttributes::new();
        assert!(!attrs.has_active_attributes());
        attrs.add_array(AnyDataArray::F64(DataArray::from_vec("s", vec![1.0], 1)));
        attrs.set_active_scalars("s");
        assert!(attrs.has_active_attributes());
    }

    #[test]
    fn iterate() {
        let mut attrs = DataSetAttributes::new();
        attrs.add_array(AnyDataArray::F64(DataArray::from_vec("x", vec![1.0], 1)));
        attrs.add_array(AnyDataArray::F64(DataArray::from_vec("y", vec![2.0], 1)));
        let names: Vec<&str> = attrs.iter().map(|a| a.name()).collect();
        assert_eq!(names, vec!["x", "y"]);
    }

    #[test]
    fn clear() {
        let mut attrs = DataSetAttributes::new();
        attrs.add_array(AnyDataArray::F64(DataArray::from_vec("a", vec![1.0], 1)));
        attrs.set_active_scalars("a");
        assert_eq!(attrs.num_arrays(), 1);
        assert!(attrs.has_active_attributes());

        attrs.clear();
        assert_eq!(attrs.num_arrays(), 0);
        assert!(!attrs.has_active_attributes());
        assert!(attrs.scalars().is_none());
    }
}
