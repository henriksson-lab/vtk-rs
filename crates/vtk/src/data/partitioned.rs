//! Partitioned dataset for distributed/decomposed data.
//!
//! Analogous to VTK's `vtkPartitionedDataSet` and `vtkPartitionedDataSetCollection`.

use crate::data::PolyData;

/// A collection of PolyData partitions for distributed or decomposed data.
///
/// Each partition has a name and contains a `PolyData` mesh. Partitions can
/// be merged into a single `PolyData` for rendering or further processing.
#[derive(Debug, Clone, Default)]
pub struct PartitionedDataSet {
    pub partitions: Vec<PolyData>,
    pub partition_names: Vec<String>,
}

impl PartitionedDataSet {
    /// Create an empty partitioned dataset.
    pub fn new() -> Self {
        Self::default()
    }

    /// Add a named partition.
    pub fn add_partition(&mut self, name: &str, data: PolyData) {
        self.partition_names.push(name.to_string());
        self.partitions.push(data);
    }

    /// Get a partition by index.
    pub fn partition(&self, idx: usize) -> Option<&PolyData> {
        self.partitions.get(idx)
    }

    /// Get a partition by name (returns the first match).
    pub fn partition_by_name(&self, name: &str) -> Option<&PolyData> {
        self.partition_names
            .iter()
            .position(|n| n == name)
            .and_then(|idx| self.partitions.get(idx))
    }

    /// Number of partitions.
    pub fn num_partitions(&self) -> usize {
        self.partitions.len()
    }

    /// Merge all partitions into a single PolyData.
    ///
    /// Points and cells from each partition are concatenated, with cell indices
    /// offset appropriately.
    pub fn merge(&self) -> PolyData {
        let mut result = PolyData::new();
        let mut point_offset: i64 = 0;

        for part in &self.partitions {
            // Append points
            for i in 0..part.points.len() {
                result.points.push(part.points.get(i));
            }

            // Append polys with offset
            for ci in 0..part.polys.num_cells() {
                let cell = part.polys.cell(ci);
                let offset_cell: Vec<i64> = cell.iter().map(|&id| id + point_offset).collect();
                result.polys.push_cell(&offset_cell);
            }

            // Append lines with offset
            for ci in 0..part.lines.num_cells() {
                let cell = part.lines.cell(ci);
                let offset_cell: Vec<i64> = cell.iter().map(|&id| id + point_offset).collect();
                result.lines.push_cell(&offset_cell);
            }

            // Append verts with offset
            for ci in 0..part.verts.num_cells() {
                let cell = part.verts.cell(ci);
                let offset_cell: Vec<i64> = cell.iter().map(|&id| id + point_offset).collect();
                result.verts.push_cell(&offset_cell);
            }

            point_offset += part.points.len() as i64;
        }

        result
    }
}

/// A collection of `PartitionedDataSet`s.
///
/// Analogous to VTK's `vtkPartitionedDataSetCollection`.
#[derive(Debug, Clone, Default)]
pub struct PartitionedDataSetCollection {
    pub datasets: Vec<PartitionedDataSet>,
}

impl PartitionedDataSetCollection {
    /// Create an empty collection.
    pub fn new() -> Self {
        Self::default()
    }

    /// Add a partitioned dataset to the collection.
    pub fn add(&mut self, dataset: PartitionedDataSet) {
        self.datasets.push(dataset);
    }

    /// Number of datasets in the collection.
    pub fn num_datasets(&self) -> usize {
        self.datasets.len()
    }

    /// Get a dataset by index.
    pub fn dataset(&self, idx: usize) -> Option<&PartitionedDataSet> {
        self.datasets.get(idx)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn triangle(offset: f64) -> PolyData {
        PolyData::from_triangles(
            vec![
                [offset, 0.0, 0.0],
                [offset + 1.0, 0.0, 0.0],
                [offset + 0.5, 1.0, 0.0],
            ],
            vec![[0, 1, 2]],
        )
    }

    #[test]
    fn add_and_query_partitions() {
        let mut pds = PartitionedDataSet::new();
        pds.add_partition("left", triangle(0.0));
        pds.add_partition("right", triangle(5.0));

        assert_eq!(pds.num_partitions(), 2);
        assert!(pds.partition(0).is_some());
        assert!(pds.partition(2).is_none());
        assert!(pds.partition_by_name("left").is_some());
        assert!(pds.partition_by_name("missing").is_none());
    }

    #[test]
    fn merge_partitions() {
        let mut pds = PartitionedDataSet::new();
        pds.add_partition("a", triangle(0.0));
        pds.add_partition("b", triangle(5.0));

        let merged = pds.merge();
        assert_eq!(merged.points.len(), 6);
        assert_eq!(merged.polys.num_cells(), 2);

        // Second triangle's indices should be offset by 3
        let cell1 = merged.polys.cell(1);
        assert_eq!(cell1, &[3, 4, 5]);
    }

    #[test]
    fn collection() {
        let mut coll = PartitionedDataSetCollection::new();
        let mut pds = PartitionedDataSet::new();
        pds.add_partition("part0", triangle(0.0));
        coll.add(pds);

        assert_eq!(coll.num_datasets(), 1);
        assert_eq!(coll.dataset(0).unwrap().num_partitions(), 1);
    }
}
