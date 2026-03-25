use crate::AnyDataArray;

/// A named, heterogeneous collection of data arrays.
///
/// Analogous to VTK's `vtkFieldData`.
#[derive(Debug, Clone, Default)]
pub struct FieldData {
    arrays: Vec<AnyDataArray>,
}

impl FieldData {
    pub fn new() -> Self {
        Self { arrays: Vec::new() }
    }

    pub fn add_array(&mut self, array: AnyDataArray) {
        // Replace existing array with the same name.
        if let Some(pos) = self.arrays.iter().position(|a| a.name() == array.name()) {
            self.arrays[pos] = array;
        } else {
            self.arrays.push(array);
        }
    }

    pub fn get_array(&self, name: &str) -> Option<&AnyDataArray> {
        self.arrays.iter().find(|a| a.name() == name)
    }

    pub fn get_array_mut(&mut self, name: &str) -> Option<&mut AnyDataArray> {
        self.arrays.iter_mut().find(|a| a.name() == name)
    }

    pub fn get_array_by_index(&self, idx: usize) -> Option<&AnyDataArray> {
        self.arrays.get(idx)
    }

    pub fn num_arrays(&self) -> usize {
        self.arrays.len()
    }

    pub fn iter(&self) -> impl Iterator<Item = &AnyDataArray> {
        self.arrays.iter()
    }

    pub fn remove_array(&mut self, name: &str) -> Option<AnyDataArray> {
        if let Some(pos) = self.arrays.iter().position(|a| a.name() == name) {
            Some(self.arrays.remove(pos))
        } else {
            None
        }
    }
}
