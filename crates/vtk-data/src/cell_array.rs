/// Storage for cell topology using offsets + connectivity arrays.
///
/// Mirrors VTK's `vtkCellArray` design: an offsets array of length `num_cells + 1`
/// and a connectivity array containing the point indices for all cells concatenated.
///
/// For cell `i`, the point indices are `connectivity[offsets[i]..offsets[i+1]]`.
#[derive(Debug, Clone, Default)]
pub struct CellArray {
    offsets: Vec<i64>,
    connectivity: Vec<i64>,
}

impl CellArray {
    pub fn new() -> Self {
        Self {
            offsets: vec![0],
            connectivity: Vec::new(),
        }
    }

    /// Create a CellArray from raw offsets and connectivity.
    pub fn from_raw(offsets: Vec<i64>, connectivity: Vec<i64>) -> Self {
        assert!(!offsets.is_empty(), "offsets must have at least one element");
        assert_eq!(offsets[0], 0, "first offset must be 0");
        Self {
            offsets,
            connectivity,
        }
    }

    /// Append a cell with the given point indices.
    pub fn push_cell(&mut self, point_ids: &[i64]) {
        self.connectivity.extend_from_slice(point_ids);
        self.offsets.push(self.connectivity.len() as i64);
    }

    pub fn num_cells(&self) -> usize {
        self.offsets.len().saturating_sub(1)
    }

    /// Returns the point indices for cell `idx`.
    pub fn cell(&self, idx: usize) -> &[i64] {
        let start = self.offsets[idx] as usize;
        let end = self.offsets[idx + 1] as usize;
        &self.connectivity[start..end]
    }

    pub fn is_empty(&self) -> bool {
        self.num_cells() == 0
    }

    pub fn offsets(&self) -> &[i64] {
        &self.offsets
    }

    pub fn connectivity(&self) -> &[i64] {
        &self.connectivity
    }

    /// Total number of point index entries across all cells.
    pub fn connectivity_len(&self) -> usize {
        self.connectivity.len()
    }

    pub fn clear(&mut self) {
        self.offsets.clear();
        self.offsets.push(0);
        self.connectivity.clear();
    }

    /// Iterate over cells, yielding a slice of point indices for each.
    pub fn iter(&self) -> CellIter<'_> {
        CellIter {
            cells: self,
            idx: 0,
        }
    }
}

pub struct CellIter<'a> {
    cells: &'a CellArray,
    idx: usize,
}

impl<'a> Iterator for CellIter<'a> {
    type Item = &'a [i64];

    fn next(&mut self) -> Option<Self::Item> {
        if self.idx >= self.cells.num_cells() {
            return None;
        }
        let cell = self.cells.cell(self.idx);
        self.idx += 1;
        Some(cell)
    }

    fn size_hint(&self) -> (usize, Option<usize>) {
        let remaining = self.cells.num_cells() - self.idx;
        (remaining, Some(remaining))
    }
}

impl ExactSizeIterator for CellIter<'_> {}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn push_and_iterate() {
        let mut cells = CellArray::new();
        cells.push_cell(&[0, 1, 2]);
        cells.push_cell(&[2, 3, 4, 5]);
        assert_eq!(cells.num_cells(), 2);
        assert_eq!(cells.cell(0), &[0, 1, 2]);
        assert_eq!(cells.cell(1), &[2, 3, 4, 5]);

        let collected: Vec<&[i64]> = cells.iter().collect();
        assert_eq!(collected.len(), 2);
    }

    #[test]
    fn empty_cell_array() {
        let cells = CellArray::new();
        assert!(cells.is_empty());
        assert_eq!(cells.num_cells(), 0);
        assert_eq!(cells.iter().count(), 0);
    }

    #[test]
    fn from_raw() {
        let cells = CellArray::from_raw(vec![0, 3, 7], vec![0, 1, 2, 3, 4, 5, 6]);
        assert_eq!(cells.num_cells(), 2);
        assert_eq!(cells.cell(0), &[0, 1, 2]);
        assert_eq!(cells.cell(1), &[3, 4, 5, 6]);
    }
}
