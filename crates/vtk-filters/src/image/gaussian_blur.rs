use vtk_data::{AnyDataArray, DataArray, ImageData};

/// Apply a Gaussian blur to a named scalar array on an ImageData.
///
/// Uses a separable 1D Gaussian kernel along each axis. The kernel radius
/// is derived from `sigma` as `ceil(3 * sigma)`. Adds a "Blurred" array
/// to the output point data.
pub fn gaussian_blur(input: &ImageData, scalars: &str, sigma: f64) -> ImageData {
    let arr = match input.point_data().get_array(scalars) {
        Some(a) => a,
        None => return input.clone(),
    };

    let dims = input.dimensions();
    let nx: usize = dims[0];
    let ny: usize = dims[1];
    let nz: usize = dims[2];
    let n: usize = nx * ny * nz;

    // Read source values
    let mut values: Vec<f64> = vec![0.0; n];
    let mut buf = [0.0f64];
    for i in 0..n {
        arr.tuple_as_f64(i, &mut buf);
        values[i] = buf[0];
    }

    // Build 1D Gaussian kernel
    let sigma: f64 = sigma.max(0.1);
    let r: usize = (3.0 * sigma).ceil() as usize;
    let r: usize = r.max(1);
    let kernel_size: usize = 2 * r + 1;
    let mut kernel: Vec<f64> = vec![0.0; kernel_size];
    let mut ksum: f64 = 0.0;
    for i in 0..kernel_size {
        let x: f64 = i as f64 - r as f64;
        let w: f64 = (-x * x / (2.0 * sigma * sigma)).exp();
        kernel[i] = w;
        ksum += w;
    }
    for w in &mut kernel {
        *w /= ksum;
    }

    // Separable: X pass
    let mut tmp: Vec<f64> = vec![0.0; n];
    for k in 0..nz {
        for j in 0..ny {
            for i in 0..nx {
                let mut acc: f64 = 0.0;
                for ki in 0..kernel_size {
                    let ii: usize =
                        (i as i64 + ki as i64 - r as i64).clamp(0, nx as i64 - 1) as usize;
                    acc += values[k * ny * nx + j * nx + ii] * kernel[ki];
                }
                tmp[k * ny * nx + j * nx + i] = acc;
            }
        }
    }

    // Y pass
    let mut tmp2: Vec<f64> = vec![0.0; n];
    for k in 0..nz {
        for j in 0..ny {
            for i in 0..nx {
                let mut acc: f64 = 0.0;
                for ki in 0..kernel_size {
                    let jj: usize =
                        (j as i64 + ki as i64 - r as i64).clamp(0, ny as i64 - 1) as usize;
                    acc += tmp[k * ny * nx + jj * nx + i] * kernel[ki];
                }
                tmp2[k * ny * nx + j * nx + i] = acc;
            }
        }
    }

    // Z pass
    let mut result: Vec<f64> = vec![0.0; n];
    for k in 0..nz {
        for j in 0..ny {
            for i in 0..nx {
                let mut acc: f64 = 0.0;
                for ki in 0..kernel_size {
                    let kk: usize =
                        (k as i64 + ki as i64 - r as i64).clamp(0, nz as i64 - 1) as usize;
                    acc += tmp2[kk * ny * nx + j * nx + i] * kernel[ki];
                }
                result[k * ny * nx + j * nx + i] = acc;
            }
        }
    }

    let mut output = input.clone();
    output
        .point_data_mut()
        .add_array(AnyDataArray::F64(DataArray::from_vec("Blurred", result, 1)));
    output
}

#[cfg(test)]
mod tests {
    use super::*;
    use vtk_data::DataSet;

    #[test]
    fn blur_preserves_dimensions() {
        let mut img = ImageData::with_dimensions(5, 5, 5);
        let n: usize = img.num_points();
        let scalars: Vec<f64> = vec![1.0; n];
        img.point_data_mut()
            .add_array(AnyDataArray::F64(DataArray::from_vec("S", scalars, 1)));

        let out = gaussian_blur(&img, "S", 1.0);
        assert_eq!(out.dimensions(), [5, 5, 5]);
        let blurred = out.point_data().get_array("Blurred").unwrap();
        assert_eq!(blurred.num_tuples(), n);
    }

    #[test]
    fn uniform_field_stays_uniform() {
        let mut img = ImageData::with_dimensions(7, 7, 7);
        let n: usize = img.num_points();
        let val: f64 = 3.5;
        let scalars: Vec<f64> = vec![val; n];
        img.point_data_mut()
            .add_array(AnyDataArray::F64(DataArray::from_vec("S", scalars, 1)));

        let out = gaussian_blur(&img, "S", 1.5);
        let blurred = out.point_data().get_array("Blurred").unwrap();
        let mut buf = [0.0f64];
        // Check center voxel — should still be ~3.5
        let center: usize = n / 2;
        blurred.tuple_as_f64(center, &mut buf);
        assert!((buf[0] - val).abs() < 1e-10);
    }

    #[test]
    fn blur_reduces_peak() {
        let mut img = ImageData::with_dimensions(11, 11, 11);
        let n: usize = img.num_points();
        let mut scalars: Vec<f64> = vec![0.0; n];
        // Single hot voxel in the center
        let cx: usize = 5;
        let cy: usize = 5;
        let cz: usize = 5;
        let center_idx: usize = cz * 11 * 11 + cy * 11 + cx;
        scalars[center_idx] = 100.0;
        img.point_data_mut()
            .add_array(AnyDataArray::F64(DataArray::from_vec("S", scalars, 1)));

        let out = gaussian_blur(&img, "S", 2.0);
        let blurred = out.point_data().get_array("Blurred").unwrap();
        let mut buf = [0.0f64];
        blurred.tuple_as_f64(center_idx, &mut buf);
        // Peak should be reduced significantly
        assert!(buf[0] < 50.0);
        assert!(buf[0] > 0.0);
    }
}
