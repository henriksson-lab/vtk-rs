use crate::{AnyDataArray, FieldData};

/// Field data with "active" attribute designations for scalars, vectors, normals, etc.
///
/// Analogous to VTK's `vtkDataSetAttributes`.
#[derive(Debug, Clone, Default)]
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

    fn find_index(&self, name: &str) -> Option<usize> {
        (0..self.field_data.num_arrays())
            .find(|&i| self.field_data.get_array_by_index(i).map(|a| a.name()) == Some(name))
    }
}
