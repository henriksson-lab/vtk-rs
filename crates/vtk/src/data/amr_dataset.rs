//! Adaptive Mesh Refinement (AMR) dataset.
//!
//! An AMR dataset consists of multiple refinement levels, each containing
//! one or more `ImageData` blocks at that level's resolution.

use crate::data::ImageData;

/// A single refinement level in an AMR hierarchy.
#[derive(Debug, Clone)]
pub struct AMRLevel {
    /// Level index (0 = coarsest).
    pub level: usize,
    /// Blocks at this refinement level.
    pub blocks: Vec<ImageData>,
    /// Grid spacing for this level.
    pub spacing: [f64; 3],
}

/// Adaptive Mesh Refinement dataset with multiple resolution levels.
///
/// Level 0 is the coarsest, and higher levels have finer spacing.
#[derive(Debug, Clone, Default)]
pub struct AMRDataSet {
    /// Refinement levels from coarsest to finest.
    pub levels: Vec<AMRLevel>,
}

impl AMRDataSet {
    /// Create an empty AMR dataset.
    pub fn new() -> Self {
        Self::default()
    }

    /// Add a new refinement level with the given spacing. Returns the level index.
    pub fn add_level(&mut self, spacing: [f64; 3]) -> usize {
        let level = self.levels.len();
        self.levels.push(AMRLevel {
            level,
            blocks: Vec::new(),
            spacing,
        });
        level
    }

    /// Add a block to the given level. Returns the block index within that level.
    pub fn add_block(&mut self, level: usize, block: ImageData) -> usize {
        let idx = self.levels[level].blocks.len();
        self.levels[level].blocks.push(block);
        idx
    }

    /// Number of refinement levels.
    pub fn num_levels(&self) -> usize {
        self.levels.len()
    }

    /// Number of blocks at the given level.
    pub fn num_blocks(&self, level: usize) -> usize {
        self.levels[level].blocks.len()
    }

    /// Get a reference to a specific block.
    pub fn block(&self, level: usize, idx: usize) -> &ImageData {
        &self.levels[level].blocks[idx]
    }

    /// Total number of blocks across all levels.
    pub fn total_blocks(&self) -> usize {
        self.levels.iter().map(|l| l.blocks.len()).sum()
    }

    /// Spacing of the coarsest (level 0) level.
    ///
    /// Returns `[0.0, 0.0, 0.0]` if there are no levels.
    pub fn coarsest_spacing(&self) -> [f64; 3] {
        self.levels.first().map_or([0.0; 3], |l| l.spacing)
    }

    /// Spacing of the finest (highest) level.
    ///
    /// Returns `[0.0, 0.0, 0.0]` if there are no levels.
    pub fn finest_spacing(&self) -> [f64; 3] {
        self.levels.last().map_or([0.0; 3], |l| l.spacing)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn amr_basic() {
        let mut amr = AMRDataSet::new();
        let l0 = amr.add_level([1.0, 1.0, 1.0]);
        let l1 = amr.add_level([0.5, 0.5, 0.5]);

        let mut block0 = ImageData::with_dimensions(10, 10, 1);
        block0.set_spacing([1.0, 1.0, 1.0]);
        amr.add_block(l0, block0);

        let mut block1 = ImageData::with_dimensions(20, 20, 1);
        block1.set_spacing([0.5, 0.5, 0.5]);
        amr.add_block(l1, block1);

        assert_eq!(amr.num_levels(), 2);
        assert_eq!(amr.num_blocks(0), 1);
        assert_eq!(amr.num_blocks(1), 1);
        assert_eq!(amr.total_blocks(), 2);
        assert_eq!(amr.coarsest_spacing(), [1.0, 1.0, 1.0]);
        assert_eq!(amr.finest_spacing(), [0.5, 0.5, 0.5]);
    }

    #[test]
    fn amr_multiple_blocks() {
        let mut amr = AMRDataSet::new();
        let l0 = amr.add_level([2.0, 2.0, 2.0]);
        amr.add_block(l0, ImageData::with_dimensions(5, 5, 5));
        amr.add_block(l0, ImageData::with_dimensions(5, 5, 5));
        amr.add_block(l0, ImageData::with_dimensions(5, 5, 5));

        assert_eq!(amr.num_blocks(0), 3);
        assert_eq!(amr.total_blocks(), 3);
        assert_eq!(amr.block(0, 1).dimensions(), [5, 5, 5]);
    }
}
