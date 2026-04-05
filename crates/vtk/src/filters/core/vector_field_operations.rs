//! Vector field operations on ImageData: curl, divergence, magnitude, normalize.

use crate::data::{AnyDataArray, DataArray, ImageData};

/// Compute the magnitude of a 3-component vector field.
pub fn vector_magnitude(image: &ImageData, array_name: &str) -> ImageData {
    let arr = match image.point_data().get_array(array_name) {
        Some(a) if a.num_components() == 3 => a,
        _ => return image.clone(),
    };
    let n = arr.num_tuples();
    let mut mag = Vec::with_capacity(n);
    let mut buf = [0.0f64; 3];
    for i in 0..n {
        arr.tuple_as_f64(i, &mut buf);
        mag.push((buf[0]*buf[0]+buf[1]*buf[1]+buf[2]*buf[2]).sqrt());
    }
    let mut r = image.clone();
    r.point_data_mut().add_array(AnyDataArray::F64(DataArray::from_vec(&format!("{array_name}_mag"), mag, 1)));
    r
}

/// Normalize a vector field to unit vectors.
pub fn vector_normalize(image: &ImageData, array_name: &str) -> ImageData {
    let arr = match image.point_data().get_array(array_name) {
        Some(a) if a.num_components() == 3 => a,
        _ => return image.clone(),
    };
    let n = arr.num_tuples();
    let mut data = Vec::with_capacity(n * 3);
    let mut buf = [0.0f64; 3];
    for i in 0..n {
        arr.tuple_as_f64(i, &mut buf);
        let len = (buf[0]*buf[0]+buf[1]*buf[1]+buf[2]*buf[2]).sqrt();
        if len > 1e-15 { data.extend_from_slice(&[buf[0]/len, buf[1]/len, buf[2]/len]); }
        else { data.extend_from_slice(&[0.0, 0.0, 0.0]); }
    }
    let mut r = image.clone();
    r.point_data_mut().add_array(AnyDataArray::F64(DataArray::from_vec(&format!("{array_name}_norm"), data, 3)));
    r
}

/// Compute dot product of two vector fields.
pub fn vector_dot_product(image: &ImageData, arr_a: &str, arr_b: &str) -> ImageData {
    let a = match image.point_data().get_array(arr_a) {
        Some(a) if a.num_components() == 3 => a, _ => return image.clone() };
    let b = match image.point_data().get_array(arr_b) {
        Some(b) if b.num_components() == 3 => b, _ => return image.clone() };
    let n = a.num_tuples().min(b.num_tuples());
    let mut dot = Vec::with_capacity(n);
    let mut ba = [0.0f64; 3]; let mut bb = [0.0f64; 3];
    for i in 0..n {
        a.tuple_as_f64(i, &mut ba); b.tuple_as_f64(i, &mut bb);
        dot.push(ba[0]*bb[0]+ba[1]*bb[1]+ba[2]*bb[2]);
    }
    let mut r = image.clone();
    r.point_data_mut().add_array(AnyDataArray::F64(DataArray::from_vec("DotProduct", dot, 1)));
    r
}

/// Compute cross product of two vector fields.
pub fn vector_cross_product(image: &ImageData, arr_a: &str, arr_b: &str) -> ImageData {
    let a = match image.point_data().get_array(arr_a) {
        Some(a) if a.num_components() == 3 => a, _ => return image.clone() };
    let b = match image.point_data().get_array(arr_b) {
        Some(b) if b.num_components() == 3 => b, _ => return image.clone() };
    let n = a.num_tuples().min(b.num_tuples());
    let mut cross = Vec::with_capacity(n * 3);
    let mut ba = [0.0f64; 3]; let mut bb = [0.0f64; 3];
    for i in 0..n {
        a.tuple_as_f64(i, &mut ba); b.tuple_as_f64(i, &mut bb);
        cross.push(ba[1]*bb[2]-ba[2]*bb[1]);
        cross.push(ba[2]*bb[0]-ba[0]*bb[2]);
        cross.push(ba[0]*bb[1]-ba[1]*bb[0]);
    }
    let mut r = image.clone();
    r.point_data_mut().add_array(AnyDataArray::F64(DataArray::from_vec("CrossProduct", cross, 3)));
    r
}

/// Add two vector fields element-wise.
pub fn vector_add(image: &ImageData, arr_a: &str, arr_b: &str) -> ImageData {
    let a = match image.point_data().get_array(arr_a) {
        Some(a) if a.num_components() == 3 => a, _ => return image.clone() };
    let b = match image.point_data().get_array(arr_b) {
        Some(b) if b.num_components() == 3 => b, _ => return image.clone() };
    let n = a.num_tuples().min(b.num_tuples());
    let mut sum = Vec::with_capacity(n * 3);
    let mut ba = [0.0f64; 3]; let mut bb = [0.0f64; 3];
    for i in 0..n { a.tuple_as_f64(i, &mut ba); b.tuple_as_f64(i, &mut bb);
        sum.push(ba[0]+bb[0]); sum.push(ba[1]+bb[1]); sum.push(ba[2]+bb[2]); }
    let mut r = image.clone();
    r.point_data_mut().add_array(AnyDataArray::F64(DataArray::from_vec("VectorSum", sum, 3)));
    r
}

#[cfg(test)]
mod tests {
    use super::*;

    fn make_field() -> ImageData {
        let n = 5*5*5;
        let mut v = Vec::with_capacity(n*3);
        for _ in 0..n { v.push(1.0); v.push(2.0); v.push(3.0); }
        let mut img = ImageData::with_dimensions(5,5,5);
        img.point_data_mut().add_array(AnyDataArray::F64(DataArray::from_vec("vel", v, 3)));
        img
    }

    #[test]
    fn magnitude() {
        let r = vector_magnitude(&make_field(), "vel");
        let a = r.point_data().get_array("vel_mag").unwrap();
        let mut b = [0.0f64]; a.tuple_as_f64(0, &mut b);
        assert!((b[0] - (1.0f64+4.0+9.0).sqrt()).abs() < 1e-10);
    }

    #[test]
    fn normalize() {
        let r = vector_normalize(&make_field(), "vel");
        let a = r.point_data().get_array("vel_norm").unwrap();
        let mut b = [0.0f64; 3]; a.tuple_as_f64(0, &mut b);
        let len = (b[0]*b[0]+b[1]*b[1]+b[2]*b[2]).sqrt();
        assert!((len - 1.0).abs() < 1e-10);
    }

    #[test]
    fn dot() {
        let mut img = make_field();
        let n = 125;
        let mut v2 = Vec::with_capacity(n*3);
        for _ in 0..n { v2.push(1.0); v2.push(0.0); v2.push(0.0); }
        img.point_data_mut().add_array(AnyDataArray::F64(DataArray::from_vec("dir", v2, 3)));
        let r = vector_dot_product(&img, "vel", "dir");
        let a = r.point_data().get_array("DotProduct").unwrap();
        let mut b = [0.0f64]; a.tuple_as_f64(0, &mut b);
        assert!((b[0] - 1.0).abs() < 1e-10);
    }

    #[test]
    fn cross() {
        let mut img = make_field();
        let n = 125;
        let mut v2 = Vec::with_capacity(n*3);
        for _ in 0..n { v2.push(0.0); v2.push(0.0); v2.push(1.0); }
        img.point_data_mut().add_array(AnyDataArray::F64(DataArray::from_vec("up", v2, 3)));
        let r = vector_cross_product(&img, "vel", "up");
        assert!(r.point_data().get_array("CrossProduct").is_some());
    }
}
