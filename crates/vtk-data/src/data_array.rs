use vtk_types::{Scalar, ScalarType};

/// A contiguous array of tuples, where each tuple has `num_components` values.
///
/// This is the fundamental data container, analogous to VTK's `vtkDataArray`.
#[derive(Debug, Clone)]
pub struct DataArray<T: Scalar> {
    data: Vec<T>,
    num_components: usize,
    name: String,
}

impl<T: Scalar> DataArray<T> {
    pub fn new(name: impl Into<String>, num_components: usize) -> Self {
        assert!(num_components > 0, "num_components must be > 0");
        Self {
            data: Vec::new(),
            num_components,
            name: name.into(),
        }
    }

    pub fn from_vec(name: impl Into<String>, data: Vec<T>, num_components: usize) -> Self {
        assert!(num_components > 0, "num_components must be > 0");
        assert!(
            data.len().is_multiple_of(num_components),
            "data length {} is not divisible by num_components {}",
            data.len(),
            num_components
        );
        Self {
            data,
            num_components,
            name: name.into(),
        }
    }

    pub fn name(&self) -> &str {
        &self.name
    }

    pub fn set_name(&mut self, name: impl Into<String>) {
        self.name = name.into();
    }

    pub fn scalar_type(&self) -> ScalarType {
        T::TYPE_ID
    }

    pub fn num_components(&self) -> usize {
        self.num_components
    }

    pub fn num_tuples(&self) -> usize {
        self.data.len() / self.num_components
    }

    pub fn tuple(&self, idx: usize) -> &[T] {
        let start = idx * self.num_components;
        &self.data[start..start + self.num_components]
    }

    pub fn tuple_mut(&mut self, idx: usize) -> &mut [T] {
        let start = idx * self.num_components;
        &mut self.data[start..start + self.num_components]
    }

    pub fn push_tuple(&mut self, values: &[T]) {
        assert_eq!(
            values.len(),
            self.num_components,
            "expected {} components, got {}",
            self.num_components,
            values.len()
        );
        self.data.extend_from_slice(values);
    }

    pub fn as_slice(&self) -> &[T] {
        &self.data
    }

    pub fn as_mut_slice(&mut self) -> &mut [T] {
        &mut self.data
    }

    pub fn into_vec(self) -> Vec<T> {
        self.data
    }

    pub fn is_empty(&self) -> bool {
        self.data.is_empty()
    }

    pub fn clear(&mut self) {
        self.data.clear();
    }
}

/// Type-erased data array that can hold any scalar type.
///
/// Uses an enum rather than trait objects for exhaustive matching and zero vtable overhead.
#[derive(Debug, Clone)]
pub enum AnyDataArray {
    F32(DataArray<f32>),
    F64(DataArray<f64>),
    I8(DataArray<i8>),
    I16(DataArray<i16>),
    I32(DataArray<i32>),
    I64(DataArray<i64>),
    U8(DataArray<u8>),
    U16(DataArray<u16>),
    U32(DataArray<u32>),
    U64(DataArray<u64>),
}

impl AnyDataArray {
    pub fn name(&self) -> &str {
        match self {
            AnyDataArray::F32(a) => a.name(),
            AnyDataArray::F64(a) => a.name(),
            AnyDataArray::I8(a) => a.name(),
            AnyDataArray::I16(a) => a.name(),
            AnyDataArray::I32(a) => a.name(),
            AnyDataArray::I64(a) => a.name(),
            AnyDataArray::U8(a) => a.name(),
            AnyDataArray::U16(a) => a.name(),
            AnyDataArray::U32(a) => a.name(),
            AnyDataArray::U64(a) => a.name(),
        }
    }

    pub fn scalar_type(&self) -> ScalarType {
        match self {
            AnyDataArray::F32(_) => ScalarType::F32,
            AnyDataArray::F64(_) => ScalarType::F64,
            AnyDataArray::I8(_) => ScalarType::I8,
            AnyDataArray::I16(_) => ScalarType::I16,
            AnyDataArray::I32(_) => ScalarType::I32,
            AnyDataArray::I64(_) => ScalarType::I64,
            AnyDataArray::U8(_) => ScalarType::U8,
            AnyDataArray::U16(_) => ScalarType::U16,
            AnyDataArray::U32(_) => ScalarType::U32,
            AnyDataArray::U64(_) => ScalarType::U64,
        }
    }

    pub fn num_components(&self) -> usize {
        match self {
            AnyDataArray::F32(a) => a.num_components(),
            AnyDataArray::F64(a) => a.num_components(),
            AnyDataArray::I8(a) => a.num_components(),
            AnyDataArray::I16(a) => a.num_components(),
            AnyDataArray::I32(a) => a.num_components(),
            AnyDataArray::I64(a) => a.num_components(),
            AnyDataArray::U8(a) => a.num_components(),
            AnyDataArray::U16(a) => a.num_components(),
            AnyDataArray::U32(a) => a.num_components(),
            AnyDataArray::U64(a) => a.num_components(),
        }
    }

    pub fn num_tuples(&self) -> usize {
        match self {
            AnyDataArray::F32(a) => a.num_tuples(),
            AnyDataArray::F64(a) => a.num_tuples(),
            AnyDataArray::I8(a) => a.num_tuples(),
            AnyDataArray::I16(a) => a.num_tuples(),
            AnyDataArray::I32(a) => a.num_tuples(),
            AnyDataArray::I64(a) => a.num_tuples(),
            AnyDataArray::U8(a) => a.num_tuples(),
            AnyDataArray::U16(a) => a.num_tuples(),
            AnyDataArray::U32(a) => a.num_tuples(),
            AnyDataArray::U64(a) => a.num_tuples(),
        }
    }

    /// Get tuple values as f64 (type-erased access).
    pub fn tuple_as_f64(&self, idx: usize, out: &mut [f64]) {
        macro_rules! convert {
            ($arr:expr) => {{
                let t = $arr.tuple(idx);
                for (o, v) in out.iter_mut().zip(t.iter()) {
                    *o = v.to_f64();
                }
            }};
        }
        match self {
            AnyDataArray::F32(a) => convert!(a),
            AnyDataArray::F64(a) => convert!(a),
            AnyDataArray::I8(a) => convert!(a),
            AnyDataArray::I16(a) => convert!(a),
            AnyDataArray::I32(a) => convert!(a),
            AnyDataArray::I64(a) => convert!(a),
            AnyDataArray::U8(a) => convert!(a),
            AnyDataArray::U16(a) => convert!(a),
            AnyDataArray::U32(a) => convert!(a),
            AnyDataArray::U64(a) => convert!(a),
        }
    }
}

// Specific From impls for each scalar type.
macro_rules! impl_from_data_array {
    ($ty:ty, $variant:ident) => {
        impl From<DataArray<$ty>> for AnyDataArray {
            fn from(a: DataArray<$ty>) -> Self {
                AnyDataArray::$variant(a)
            }
        }
    };
}

impl_from_data_array!(f32, F32);
impl_from_data_array!(f64, F64);
impl_from_data_array!(i8, I8);
impl_from_data_array!(i16, I16);
impl_from_data_array!(i32, I32);
impl_from_data_array!(i64, I64);
impl_from_data_array!(u8, U8);
impl_from_data_array!(u16, U16);
impl_from_data_array!(u32, U32);
impl_from_data_array!(u64, U64);

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn data_array_basics() {
        let mut arr = DataArray::<f64>::new("test", 3);
        arr.push_tuple(&[1.0, 2.0, 3.0]);
        arr.push_tuple(&[4.0, 5.0, 6.0]);
        assert_eq!(arr.num_tuples(), 2);
        assert_eq!(arr.num_components(), 3);
        assert_eq!(arr.tuple(0), &[1.0, 2.0, 3.0]);
        assert_eq!(arr.tuple(1), &[4.0, 5.0, 6.0]);
    }

    #[test]
    fn data_array_from_vec() {
        let arr = DataArray::from_vec("coords", vec![1.0f32, 2.0, 3.0, 4.0, 5.0, 6.0], 3);
        assert_eq!(arr.num_tuples(), 2);
        assert_eq!(arr.scalar_type(), ScalarType::F32);
    }

    #[test]
    fn any_data_array_type_erased() {
        let arr = DataArray::from_vec("scalars", vec![10i32, 20, 30], 1);
        let any: AnyDataArray = arr.into();
        assert_eq!(any.name(), "scalars");
        assert_eq!(any.scalar_type(), ScalarType::I32);
        assert_eq!(any.num_tuples(), 3);

        let mut val = [0.0f64];
        any.tuple_as_f64(1, &mut val);
        assert_eq!(val[0], 20.0);
    }
}
