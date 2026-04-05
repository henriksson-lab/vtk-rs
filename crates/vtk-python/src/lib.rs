use pyo3::prelude::*;
use pyo3::exceptions::PyValueError;

use vtk::data::{DataArray, PolyData};

/// A Python-accessible PolyData mesh.
#[pyclass(name = "PolyData")]
#[derive(Clone)]
struct PyPolyData {
    inner: PolyData,
}

#[pymethods]
impl PyPolyData {
    #[new]
    fn new() -> Self {
        Self { inner: PolyData::new() }
    }

    /// Create from flat point and triangle index arrays.
    #[staticmethod]
    fn from_triangles(points: Vec<f64>, triangles: Vec<i64>) -> PyResult<Self> {
        if points.len() % 3 != 0 {
            return Err(PyValueError::new_err("points length must be divisible by 3"));
        }
        if triangles.len() % 3 != 0 {
            return Err(PyValueError::new_err("triangles length must be divisible by 3"));
        }
        let pts: Vec<[f64; 3]> = points.chunks_exact(3)
            .map(|c| [c[0], c[1], c[2]])
            .collect();
        let tris: Vec<[i64; 3]> = triangles.chunks_exact(3)
            .map(|c| [c[0], c[1], c[2]])
            .collect();
        Ok(Self { inner: PolyData::from_triangles(pts, tris) })
    }

    #[getter]
    fn num_points(&self) -> usize {
        self.inner.points.len()
    }

    #[getter]
    fn num_cells(&self) -> usize {
        self.inner.polys.num_cells()
    }

    fn get_points(&self) -> Vec<f64> {
        let mut result = Vec::with_capacity(self.inner.points.len() * 3);
        for i in 0..self.inner.points.len() {
            let p = self.inner.points.get(i);
            result.extend_from_slice(&p);
        }
        result
    }

    fn add_point_scalars(&mut self, name: &str, values: Vec<f64>) -> PyResult<()> {
        if values.len() != self.inner.points.len() {
            return Err(PyValueError::new_err(format!(
                "expected {} values, got {}", self.inner.points.len(), values.len()
            )));
        }
        let arr = DataArray::<f64>::from_vec(name, values, 1);
        self.inner.point_data_mut().add_array(arr.into());
        self.inner.point_data_mut().set_active_scalars(name);
        Ok(())
    }

    fn get_point_scalars(&self, name: &str) -> PyResult<Vec<f64>> {
        let arr = self.inner.point_data().get_array(name)
            .ok_or_else(|| PyValueError::new_err(format!("array '{name}' not found")))?;
        let mut result = Vec::with_capacity(arr.num_tuples());
        let mut buf = [0.0f64];
        for i in 0..arr.num_tuples() {
            arr.tuple_as_f64(i, &mut buf);
            result.push(buf[0]);
        }
        Ok(result)
    }

    fn compute_normals(&self) -> Self {
        Self { inner: vtk::filters::core::normals::compute_normals(&self.inner) }
    }

    fn triangulate(&self) -> Self {
        Self { inner: vtk::filters::core::triangulate::triangulate(&self.inner) }
    }

    fn elevation_z(&self) -> Self {
        Self { inner: vtk::filters::core::elevation::elevation_z(&self.inner) }
    }

    fn decimate(&self, target_reduction: f64) -> Self {
        Self { inner: vtk::filters::core::decimate::decimate(&self.inner, target_reduction) }
    }

    fn write_vtk(&self, path: &str) -> PyResult<()> {
        let p = std::path::Path::new(path);
        vtk::io::legacy::LegacyWriter::ascii()
            .write_poly_data(p, &self.inner)
            .map_err(|e| PyValueError::new_err(format!("{e}")))
    }

    #[staticmethod]
    fn read_vtk(path: &str) -> PyResult<Self> {
        let p = std::path::Path::new(path);
        let pd = vtk::io::legacy::LegacyReader::read_poly_data(p)
            .map_err(|e| PyValueError::new_err(format!("{e}")))?;
        Ok(Self { inner: pd })
    }

    fn write_stl(&self, path: &str) -> PyResult<()> {
        let p = std::path::Path::new(path);
        vtk::io::stl::StlWriter::binary()
            .write(p, &self.inner)
            .map_err(|e| PyValueError::new_err(format!("{e}")))
    }

    #[staticmethod]
    fn read_stl(path: &str) -> PyResult<Self> {
        let p = std::path::Path::new(path);
        let pd = vtk::io::stl::StlReader::read(p)
            .map_err(|e| PyValueError::new_err(format!("{e}")))?;
        Ok(Self { inner: pd })
    }

    fn __repr__(&self) -> String {
        format!("PolyData(points={}, cells={})", self.num_points(), self.num_cells())
    }
}

#[pyfunction]
fn sphere(theta_resolution: usize, phi_resolution: usize) -> PyPolyData {
    use vtk::filters::core::sources::sphere::{sphere as make_sphere, SphereParams};
    PyPolyData {
        inner: make_sphere(&SphereParams {
            theta_resolution,
            phi_resolution,
            ..SphereParams::default()
        }),
    }
}

#[pyfunction]
fn cube() -> PyPolyData {
    use vtk::filters::core::sources::cube::{cube as make_cube, CubeParams};
    PyPolyData { inner: make_cube(&CubeParams::default()) }
}

#[pyfunction]
fn cone() -> PyPolyData {
    use vtk::filters::core::sources::cone::{cone as make_cone, ConeParams};
    PyPolyData { inner: make_cone(&ConeParams::default()) }
}

#[pymodule]
fn vtk_python(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<PyPolyData>()?;
    m.add_function(wrap_pyfunction!(sphere, m)?)?;
    m.add_function(wrap_pyfunction!(cube, m)?)?;
    m.add_function(wrap_pyfunction!(cone, m)?)?;
    Ok(())
}
