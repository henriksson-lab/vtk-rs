use vtk_types::BoundingBox;

use crate::{CellArray, DataSetAttributes, FieldData, Points};
use crate::traits::{DataObject, DataSet};

/// Polygonal mesh consisting of vertices, lines, polygons, and triangle strips.
///
/// Analogous to VTK's `vtkPolyData`. Points are stored explicitly, with four
/// separate `CellArray`s for the four cell categories.
#[derive(Debug, Clone, Default)]
pub struct PolyData {
    pub points: Points<f64>,
    pub verts: CellArray,
    pub lines: CellArray,
    pub polys: CellArray,
    pub strips: CellArray,
    point_data: DataSetAttributes,
    cell_data: DataSetAttributes,
    field_data: FieldData,
}

impl PolyData {
    pub fn new() -> Self {
        Self::default()
    }

    /// Create a PolyData from points and triangle connectivity.
    ///
    /// Each element of `triangles` is `[i0, i1, i2]` indexing into `points`.
    pub fn from_triangles(points: Vec<[f64; 3]>, triangles: Vec<[i64; 3]>) -> Self {
        let pts = Points::from_vec(points);
        let mut polys = CellArray::new();
        for tri in &triangles {
            polys.push_cell(&[tri[0], tri[1], tri[2]]);
        }
        Self {
            points: pts,
            polys,
            ..Default::default()
        }
    }

    /// Total number of cells across all four categories.
    pub fn total_cells(&self) -> usize {
        self.verts.num_cells()
            + self.lines.num_cells()
            + self.polys.num_cells()
            + self.strips.num_cells()
    }

    pub fn point_data(&self) -> &DataSetAttributes {
        &self.point_data
    }

    pub fn point_data_mut(&mut self) -> &mut DataSetAttributes {
        &mut self.point_data
    }

    pub fn cell_data(&self) -> &DataSetAttributes {
        &self.cell_data
    }

    pub fn cell_data_mut(&mut self) -> &mut DataSetAttributes {
        &mut self.cell_data
    }
}

impl DataObject for PolyData {
    fn field_data(&self) -> &FieldData {
        &self.field_data
    }

    fn field_data_mut(&mut self) -> &mut FieldData {
        &mut self.field_data
    }
}

impl DataSet for PolyData {
    fn num_points(&self) -> usize {
        self.points.len()
    }

    fn num_cells(&self) -> usize {
        self.total_cells()
    }

    fn point(&self, idx: usize) -> [f64; 3] {
        self.points.get(idx)
    }

    fn bounds(&self) -> BoundingBox {
        self.points.bounds()
    }

    fn point_data(&self) -> &DataSetAttributes {
        &self.point_data
    }

    fn point_data_mut(&mut self) -> &mut DataSetAttributes {
        &mut self.point_data
    }

    fn cell_data(&self) -> &DataSetAttributes {
        &self.cell_data
    }

    fn cell_data_mut(&mut self) -> &mut DataSetAttributes {
        &mut self.cell_data
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn from_triangles() {
        let pd = PolyData::from_triangles(
            vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.5, 1.0, 0.0]],
            vec![[0, 1, 2]],
        );
        assert_eq!(pd.num_points(), 3);
        assert_eq!(pd.num_cells(), 1);
        assert_eq!(pd.polys.cell(0), &[0, 1, 2]);
    }

    #[test]
    fn bounds() {
        let pd = PolyData::from_triangles(
            vec![[0.0, 0.0, 0.0], [1.0, 2.0, 3.0], [0.5, 1.0, 1.5]],
            vec![[0, 1, 2]],
        );
        let bb = pd.bounds();
        assert_eq!(bb.x_min, 0.0);
        assert_eq!(bb.x_max, 1.0);
        assert_eq!(bb.y_max, 2.0);
        assert_eq!(bb.z_max, 3.0);
    }

    #[test]
    fn empty_poly_data() {
        let pd = PolyData::new();
        assert_eq!(pd.num_points(), 0);
        assert_eq!(pd.num_cells(), 0);
    }
}
