use vtk_types::{BoundingBox, Scalar};

use crate::DataArray;

/// A collection of 3D points, backed by a 3-component `DataArray`.
///
/// Default precision is `f64`, matching VTK convention.
#[derive(Debug, Clone)]
pub struct Points<T: Scalar = f64> {
    data: DataArray<T>,
}

impl<T: Scalar> Points<T> {
    pub fn new() -> Self {
        Self {
            data: DataArray::new("Points", 3),
        }
    }

    pub fn with_capacity(_capacity: usize) -> Self {
        // DataArray doesn't expose reserve yet, so just create empty for now.
        Self::new()
    }

    pub fn from_vec(points: Vec<[T; 3]>) -> Self {
        let flat: Vec<T> = points.into_iter().flatten().collect();
        Self {
            data: DataArray::from_vec("Points", flat, 3),
        }
    }

    pub fn push(&mut self, point: [T; 3]) {
        self.data.push_tuple(&point);
    }

    pub fn get(&self, idx: usize) -> [T; 3] {
        let t = self.data.tuple(idx);
        [t[0], t[1], t[2]]
    }

    pub fn set(&mut self, idx: usize, point: [T; 3]) {
        let t = self.data.tuple_mut(idx);
        t[0] = point[0];
        t[1] = point[1];
        t[2] = point[2];
    }

    pub fn len(&self) -> usize {
        self.data.num_tuples()
    }

    pub fn is_empty(&self) -> bool {
        self.data.is_empty()
    }

    pub fn as_data_array(&self) -> &DataArray<T> {
        &self.data
    }

    pub fn into_data_array(self) -> DataArray<T> {
        self.data
    }

    pub fn bounds(&self) -> BoundingBox {
        let mut bb = BoundingBox::empty();
        for i in 0..self.len() {
            let p = self.get(i);
            bb.expand([p[0].to_f64(), p[1].to_f64(), p[2].to_f64()]);
        }
        bb
    }
}

impl<T: Scalar> Default for Points<T> {
    fn default() -> Self {
        Self::new()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn push_and_get() {
        let mut pts = Points::<f64>::new();
        pts.push([1.0, 2.0, 3.0]);
        pts.push([4.0, 5.0, 6.0]);
        assert_eq!(pts.len(), 2);
        assert_eq!(pts.get(0), [1.0, 2.0, 3.0]);
        assert_eq!(pts.get(1), [4.0, 5.0, 6.0]);
    }

    #[test]
    fn from_vec_points() {
        let pts = Points::from_vec(vec![[0.0, 0.0, 0.0], [1.0, 1.0, 1.0]]);
        assert_eq!(pts.len(), 2);
        assert_eq!(pts.get(1), [1.0, 1.0, 1.0]);
    }

    #[test]
    fn bounds_computation() {
        let pts = Points::from_vec(vec![[0.0, -1.0, 2.0], [3.0, 4.0, -5.0]]);
        let bb = pts.bounds();
        assert_eq!(bb.x_min, 0.0);
        assert_eq!(bb.x_max, 3.0);
        assert_eq!(bb.y_min, -1.0);
        assert_eq!(bb.y_max, 4.0);
        assert_eq!(bb.z_min, -5.0);
        assert_eq!(bb.z_max, 2.0);
    }
}
