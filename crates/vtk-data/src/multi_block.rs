use crate::{FieldData, PolyData, ImageData, UnstructuredGrid, RectilinearGrid, StructuredGrid};
use crate::traits::DataObject;

/// A block in a MultiBlockDataSet, which can hold different dataset types.
#[derive(Debug, Clone)]
pub enum Block {
    PolyData(PolyData),
    ImageData(ImageData),
    UnstructuredGrid(UnstructuredGrid),
    RectilinearGrid(RectilinearGrid),
    StructuredGrid(StructuredGrid),
    MultiBlock(MultiBlockDataSet),
}

/// A composite dataset containing named blocks of heterogeneous data types.
///
/// Analogous to VTK's `vtkMultiBlockDataSet`. Each block can be any dataset
/// type, including nested MultiBlockDataSets.
#[derive(Debug, Clone, Default)]
pub struct MultiBlockDataSet {
    blocks: Vec<(Option<String>, Block)>,
    field_data: FieldData,
}

impl MultiBlockDataSet {
    pub fn new() -> Self {
        Self::default()
    }

    /// Add a named block.
    pub fn add_block(&mut self, name: impl Into<String>, block: Block) {
        self.blocks.push((Some(name.into()), block));
    }

    /// Add an unnamed block.
    pub fn add_unnamed_block(&mut self, block: Block) {
        self.blocks.push((None, block));
    }

    /// Number of blocks.
    pub fn num_blocks(&self) -> usize {
        self.blocks.len()
    }

    /// Get block by index.
    pub fn block(&self, idx: usize) -> Option<&Block> {
        self.blocks.get(idx).map(|(_, b)| b)
    }

    /// Get block name by index.
    pub fn block_name(&self, idx: usize) -> Option<&str> {
        self.blocks.get(idx).and_then(|(n, _)| n.as_deref())
    }

    /// Get block by name (first match).
    pub fn block_by_name(&self, name: &str) -> Option<&Block> {
        self.blocks.iter().find_map(|(n, b)| {
            if n.as_deref() == Some(name) { Some(b) } else { None }
        })
    }

    /// Iterate over all blocks with their optional names.
    pub fn iter(&self) -> impl Iterator<Item = (Option<&str>, &Block)> {
        self.blocks.iter().map(|(n, b)| (n.as_deref(), b))
    }
}

impl DataObject for MultiBlockDataSet {
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

    #[test]
    fn multi_block_basics() {
        let mut mb = MultiBlockDataSet::new();
        let pd = PolyData::from_triangles(
            vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.5, 1.0, 0.0]],
            vec![[0, 1, 2]],
        );
        mb.add_block("mesh", Block::PolyData(pd));
        mb.add_block("image", Block::ImageData(ImageData::with_dimensions(2, 2, 2)));

        assert_eq!(mb.num_blocks(), 2);
        assert_eq!(mb.block_name(0), Some("mesh"));
        assert!(mb.block_by_name("image").is_some());
    }

    #[test]
    fn nested_multi_block() {
        let mut inner = MultiBlockDataSet::new();
        inner.add_block("tri", Block::PolyData(PolyData::new()));

        let mut outer = MultiBlockDataSet::new();
        outer.add_block("sub", Block::MultiBlock(inner));

        assert_eq!(outer.num_blocks(), 1);
        if let Some(Block::MultiBlock(sub)) = outer.block(0) {
            assert_eq!(sub.num_blocks(), 1);
        } else {
            panic!("expected nested MultiBlock");
        }
    }
}
