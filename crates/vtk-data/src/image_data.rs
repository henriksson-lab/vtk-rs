use vtk_types::BoundingBox;

use crate::{DataSetAttributes, FieldData};
use crate::traits::{DataObject, DataSet};

/// Regular grid with implicit point coordinates computed from extent, spacing, and origin.
///
/// Analogous to VTK's `vtkImageData`. Points are not stored explicitly — their
/// coordinates are computed as: `point(i,j,k) = origin + spacing * [i, j, k]`.
///
/// The extent is `[x_min, x_max, y_min, y_max, z_min, z_max]` in index space.
///
/// # Examples
///
/// ```
/// use vtk_data::ImageData;
///
/// // Create a 10x10x10 grid
/// let img = ImageData::with_dimensions(10, 10, 10);
/// assert_eq!(img.dimensions(), [10, 10, 10]);
/// ```
#[derive(Debug, Clone)]
pub struct ImageData {
    extent: [i64; 6],
    spacing: [f64; 3],
    origin: [f64; 3],
    point_data: DataSetAttributes,
    cell_data: DataSetAttributes,
    field_data: FieldData,
}

impl Default for ImageData {
    fn default() -> Self {
        Self {
            extent: [0, 0, 0, 0, 0, 0],
            spacing: [1.0, 1.0, 1.0],
            origin: [0.0, 0.0, 0.0],
            point_data: DataSetAttributes::new(),
            cell_data: DataSetAttributes::new(),
            field_data: FieldData::new(),
        }
    }
}

impl ImageData {
    pub fn new() -> Self {
        Self::default()
    }

    /// Create an ImageData with given dimensions (number of points in each direction).
    pub fn with_dimensions(nx: usize, ny: usize, nz: usize) -> Self {
        Self {
            extent: [0, (nx as i64) - 1, 0, (ny as i64) - 1, 0, (nz as i64) - 1],
            ..Default::default()
        }
    }

    pub fn extent(&self) -> [i64; 6] {
        self.extent
    }

    pub fn set_extent(&mut self, extent: [i64; 6]) {
        self.extent = extent;
    }

    pub fn spacing(&self) -> [f64; 3] {
        self.spacing
    }

    pub fn set_spacing(&mut self, spacing: [f64; 3]) {
        self.spacing = spacing;
    }

    pub fn origin(&self) -> [f64; 3] {
        self.origin
    }

    pub fn set_origin(&mut self, origin: [f64; 3]) {
        self.origin = origin;
    }

    /// Number of points in each dimension.
    pub fn dimensions(&self) -> [usize; 3] {
        [
            (self.extent[1] - self.extent[0] + 1).max(0) as usize,
            (self.extent[3] - self.extent[2] + 1).max(0) as usize,
            (self.extent[5] - self.extent[4] + 1).max(0) as usize,
        ]
    }

    /// Compute the world-space position of a point given its (i, j, k) index.
    pub fn point_from_ijk(&self, i: usize, j: usize, k: usize) -> [f64; 3] {
        [
            self.origin[0] + (self.extent[0] as f64 + i as f64) * self.spacing[0],
            self.origin[1] + (self.extent[2] as f64 + j as f64) * self.spacing[1],
            self.origin[2] + (self.extent[4] as f64 + k as f64) * self.spacing[2],
        ]
    }

    /// Convert a flat point index to (i, j, k).
    pub fn ijk_from_index(&self, idx: usize) -> (usize, usize, usize) {
        let dims = self.dimensions();
        let k = idx / (dims[0] * dims[1]);
        let remainder = idx % (dims[0] * dims[1]);
        let j = remainder / dims[0];
        let i = remainder % dims[0];
        (i, j, k)
    }

    /// Get point coordinates by flat index (equivalent to `point_from_ijk(ijk_from_index(idx))`).
    pub fn point_at(&self, idx: usize) -> [f64; 3] {
        let (i, j, k) = self.ijk_from_index(idx);
        self.point_from_ijk(i, j, k)
    }

    /// Convert (i, j, k) to a flat point index.
    pub fn index_from_ijk(&self, i: usize, j: usize, k: usize) -> usize {
        let dims = self.dimensions();
        k * dims[0] * dims[1] + j * dims[0] + i
    }

    /// Get the active scalar value at grid position (i, j, k).
    ///
    /// Returns None if no active scalars are set.
    pub fn scalar_at(&self, i: usize, j: usize, k: usize) -> Option<f64> {
        let scalars = self.point_data.scalars()?;
        let idx = self.index_from_ijk(i, j, k);
        let mut buf = [0.0f64];
        scalars.tuple_as_f64(idx, &mut buf);
        Some(buf[0])
    }

    /// Iterate over all point coordinates in index order.
    pub fn point_positions(&self) -> Vec<[f64; 3]> {
        let n = self.num_points();
        (0..n).map(|idx| self.point_at(idx)).collect()
    }

    /// Number of points.
    pub fn num_points(&self) -> usize {
        let d = self.dimensions();
        d[0] * d[1] * d[2]
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

    /// Builder: set spacing.
    pub fn with_spacing(mut self, spacing: [f64; 3]) -> Self {
        self.spacing = spacing;
        self
    }

    /// Builder: set origin.
    pub fn with_origin(mut self, origin: [f64; 3]) -> Self {
        self.origin = origin;
        self
    }

    /// Create an ImageData with a scalar field generated from a function.
    ///
    /// The function receives `(x, y, z)` world coordinates and returns a scalar value.
    /// The result is stored as point data with the given name.
    ///
    /// # Examples
    ///
    /// ```
    /// use vtk_data::ImageData;
    ///
    /// // Create a sphere distance field
    /// let img = ImageData::from_function(
    ///     [10, 10, 10], [0.1, 0.1, 0.1], [-0.5, -0.5, -0.5],
    ///     "distance",
    ///     |x, y, z| (x*x + y*y + z*z).sqrt(),
    /// );
    /// assert_eq!(img.dimensions(), [10, 10, 10]);
    /// ```
    pub fn from_function(
        dims: [usize; 3],
        spacing: [f64; 3],
        origin: [f64; 3],
        name: &str,
        f: impl Fn(f64, f64, f64) -> f64,
    ) -> Self {
        let mut img = Self::with_dimensions(dims[0], dims[1], dims[2]);
        img.set_spacing(spacing);
        img.set_origin(origin);

        let mut values = Vec::with_capacity(dims[0] * dims[1] * dims[2]);
        for k in 0..dims[2] {
            for j in 0..dims[1] {
                for i in 0..dims[0] {
                    let x = origin[0] + i as f64 * spacing[0];
                    let y = origin[1] + j as f64 * spacing[1];
                    let z = origin[2] + k as f64 * spacing[2];
                    values.push(f(x, y, z));
                }
            }
        }

        let arr = crate::DataArray::from_vec(name, values, 1);
        img.point_data.add_array(crate::AnyDataArray::F64(arr));
        img.point_data.set_active_scalars(name);
        img
    }

    /// Builder: add a point data array.
    pub fn with_point_array(mut self, array: crate::AnyDataArray) -> Self {
        let name = array.name().to_string();
        self.point_data.add_array(array);
        if self.point_data.scalars().is_none() {
            self.point_data.set_active_scalars(&name);
        }
        self
    }
}

impl DataObject for ImageData {
    fn field_data(&self) -> &FieldData {
        &self.field_data
    }

    fn field_data_mut(&mut self) -> &mut FieldData {
        &mut self.field_data
    }
}

impl DataSet for ImageData {
    fn num_points(&self) -> usize {
        let d = self.dimensions();
        d[0] * d[1] * d[2]
    }

    fn num_cells(&self) -> usize {
        let d = self.dimensions();
        let cx = d[0].saturating_sub(1);
        let cy = d[1].saturating_sub(1);
        let cz = d[2].saturating_sub(1);
        // For lower-dimensional data, treat missing dims as 1 cell
        match (cx, cy, cz) {
            (0, _, _) | (_, 0, _) | (_, _, 0) => cx.max(1) * cy.max(1) * cz.max(1),
            _ => cx * cy * cz,
        }
    }

    fn point(&self, idx: usize) -> [f64; 3] {
        let (i, j, k) = self.ijk_from_index(idx);
        self.point_from_ijk(i, j, k)
    }

    fn bounds(&self) -> BoundingBox {
        let d = self.dimensions();
        if d[0] == 0 || d[1] == 0 || d[2] == 0 {
            return BoundingBox::empty();
        }
        let p0 = self.point_from_ijk(0, 0, 0);
        let p1 = self.point_from_ijk(d[0] - 1, d[1] - 1, d[2] - 1);
        let mut bb = BoundingBox::empty();
        bb.expand(p0);
        bb.expand(p1);
        bb
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

impl std::fmt::Display for ImageData {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let d = self.dimensions();
        write!(f, "ImageData: {}x{}x{}, spacing {:?}, origin {:?}, {} point arrays",
            d[0], d[1], d[2], self.spacing(), self.origin(), self.point_data.num_arrays())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn dimensions_and_points() {
        let img = ImageData::with_dimensions(3, 4, 5);
        assert_eq!(img.dimensions(), [3, 4, 5]);
        assert_eq!(img.num_points(), 60);
        assert_eq!(img.point_from_ijk(0, 0, 0), [0.0, 0.0, 0.0]);
        assert_eq!(img.point_from_ijk(2, 3, 4), [2.0, 3.0, 4.0]);
    }

    #[test]
    fn custom_spacing_and_origin() {
        let mut img = ImageData::with_dimensions(2, 2, 2);
        img.set_spacing([0.5, 0.5, 0.5]);
        img.set_origin([1.0, 2.0, 3.0]);
        assert_eq!(img.point_from_ijk(1, 1, 1), [1.5, 2.5, 3.5]);
    }

    #[test]
    fn index_roundtrip() {
        let img = ImageData::with_dimensions(3, 4, 5);
        for idx in 0..60 {
            let (i, j, k) = img.ijk_from_index(idx);
            assert_eq!(img.index_from_ijk(i, j, k), idx);
        }
    }

    #[test]
    fn bounds_computation() {
        let mut img = ImageData::with_dimensions(10, 10, 10);
        img.set_spacing([0.1, 0.1, 0.1]);
        let bb = img.bounds();
        assert!((bb.x_max - 0.9).abs() < 1e-10);
        assert!((bb.y_max - 0.9).abs() < 1e-10);
    }

    #[test]
    fn num_cells() {
        let img = ImageData::with_dimensions(3, 4, 5);
        assert_eq!(img.num_cells(), 2 * 3 * 4); // 24 cells
    }

    #[test]
    fn from_function() {
        let img = ImageData::from_function(
            [5, 5, 5], [0.2, 0.2, 0.2], [0.0, 0.0, 0.0],
            "dist",
            |x, y, z| (x*x + y*y + z*z).sqrt(),
        );
        assert_eq!(img.dimensions(), [5, 5, 5]);
        assert!(img.point_data().scalars().is_some());
        assert_eq!(img.point_data().scalars().unwrap().name(), "dist");
        assert_eq!(img.point_data().scalars().unwrap().num_tuples(), 125);
    }

    #[test]
    fn builder_chain() {
        let img = ImageData::with_dimensions(3, 3, 3)
            .with_spacing([0.5, 0.5, 0.5])
            .with_origin([1.0, 2.0, 3.0]);
        assert_eq!(img.spacing(), [0.5, 0.5, 0.5]);
        assert_eq!(img.origin(), [1.0, 2.0, 3.0]);
    }

    #[test]
    fn scalar_at() {
        let img = ImageData::from_function(
            [3, 3, 3], [1.0, 1.0, 1.0], [0.0, 0.0, 0.0],
            "val", |x, _y, _z| x,
        );
        // At (0,0,0) x=0, at (2,0,0) x=2
        assert!((img.scalar_at(0, 0, 0).unwrap()).abs() < 1e-10);
        assert!((img.scalar_at(2, 0, 0).unwrap() - 2.0).abs() < 1e-10);
    }

    #[test]
    fn point_positions() {
        let img = ImageData::with_dimensions(2, 2, 1);
        let positions = img.point_positions();
        assert_eq!(positions.len(), 4);
        assert_eq!(positions[0], [0.0, 0.0, 0.0]);
        assert_eq!(positions[1], [1.0, 0.0, 0.0]);
    }
}
