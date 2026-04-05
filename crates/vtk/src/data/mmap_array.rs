//! File-backed data array loading and binary I/O utilities.
//!
//! Provides `MmapDataArray` for lazily loading binary f64 data from files,
//! `MmapPointCloud` for loading XYZ point clouds, and raw binary
//! read/write helpers for `DataArray<f64>`.

use std::fs::File;
use std::io::{Read, Seek, SeekFrom, Write};
use std::path::{Path, PathBuf};

use crate::data::{CellArray, DataArray, Points, PolyData};

/// Descriptor for a binary f64 array stored in a file.
///
/// The file is not opened until [`load`](MmapDataArray::load) is called,
/// enabling lazy / deferred loading of large datasets.
#[derive(Debug, Clone)]
pub struct MmapDataArray {
    /// Path to the binary file.
    pub path: PathBuf,
    /// Byte offset where the array data begins.
    pub offset: u64,
    /// Number of tuples in the array.
    pub num_tuples: usize,
    /// Number of f64 components per tuple.
    pub num_components: usize,
}

impl MmapDataArray {
    /// Create a new descriptor.
    pub fn new(path: impl Into<PathBuf>, offset: u64, num_tuples: usize, num_components: usize) -> Self {
        Self { path: path.into(), offset, num_tuples, num_components }
    }

    /// Read the data from the file and return a `DataArray<f64>`.
    pub fn load(&self) -> std::io::Result<DataArray<f64>> {
        let total = self.num_tuples * self.num_components;
        let byte_len = total * 8;
        let mut file = File::open(&self.path)?;
        file.seek(SeekFrom::Start(self.offset))?;
        let mut buf = vec![0u8; byte_len];
        file.read_exact(&mut buf)?;
        let data: Vec<f64> = buf
            .chunks_exact(8)
            .map(|c| f64::from_le_bytes(c.try_into().unwrap()))
            .collect();
        Ok(DataArray::from_vec("mmap_array", data, self.num_components))
    }
}

/// Descriptor for a binary XYZ point cloud stored as packed f64 triples.
#[derive(Debug, Clone)]
pub struct MmapPointCloud {
    /// Path to the binary file.
    pub path: PathBuf,
    /// Byte offset where the point data begins.
    pub offset: u64,
    /// Number of points (each point is 3 consecutive f64 values).
    pub num_points: usize,
}

impl MmapPointCloud {
    pub fn new(path: impl Into<PathBuf>, offset: u64, num_points: usize) -> Self {
        Self { path: path.into(), offset, num_points }
    }

    /// Load the points from the file and return a `PolyData` with vertex cells.
    pub fn load(&self) -> std::io::Result<PolyData> {
        let byte_len = self.num_points * 3 * 8;
        let mut file = File::open(&self.path)?;
        file.seek(SeekFrom::Start(self.offset))?;
        let mut buf = vec![0u8; byte_len];
        file.read_exact(&mut buf)?;

        let mut points = Points::<f64>::new();
        let mut verts = CellArray::new();
        for i in 0..self.num_points {
            let base = i * 3 * 8;
            let x = f64::from_le_bytes(buf[base..base + 8].try_into().unwrap());
            let y = f64::from_le_bytes(buf[base + 8..base + 16].try_into().unwrap());
            let z = f64::from_le_bytes(buf[base + 16..base + 24].try_into().unwrap());
            points.push([x, y, z]);
            verts.push_cell(&[i as i64]);
        }

        let mut pd = PolyData::new();
        pd.points = points;
        pd.verts = verts;
        Ok(pd)
    }
}

/// Write a `DataArray<f64>` as raw little-endian binary.
pub fn write_binary_array(path: &Path, array: &DataArray<f64>) -> std::io::Result<()> {
    let mut file = File::create(path)?;
    for &v in array.as_slice() {
        file.write_all(&v.to_le_bytes())?;
    }
    file.flush()?;
    Ok(())
}

/// Read a raw little-endian f64 binary file into a `DataArray<f64>`.
pub fn read_binary_array(path: &Path, name: &str, num_components: usize) -> std::io::Result<DataArray<f64>> {
    let mut file = File::open(path)?;
    let mut buf = Vec::new();
    file.read_to_end(&mut buf)?;
    let data: Vec<f64> = buf
        .chunks_exact(8)
        .map(|c| f64::from_le_bytes(c.try_into().unwrap()))
        .collect();
    Ok(DataArray::from_vec(name, data, num_components))
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write as _;

    #[test]
    fn write_read_roundtrip() {
        let dir = std::env::temp_dir().join("vtk_mmap_test_roundtrip");
        let _ = std::fs::create_dir_all(&dir);
        let path = dir.join("array.bin");

        let orig = DataArray::<f64>::from_vec("temp", vec![1.0, 2.0, 3.0, 4.0, 5.0, 6.0], 3);
        write_binary_array(&path, &orig).unwrap();
        let loaded = read_binary_array(&path, "temp", 3).unwrap();
        assert_eq!(loaded.num_tuples(), 2);
        assert_eq!(loaded.num_components(), 3);
        assert_eq!(loaded.as_slice(), orig.as_slice());

        let _ = std::fs::remove_dir_all(&dir);
    }

    #[test]
    fn mmap_data_array_load() {
        let dir = std::env::temp_dir().join("vtk_mmap_test_load");
        let _ = std::fs::create_dir_all(&dir);
        let path = dir.join("data.bin");

        // Write header bytes then array data
        let mut f = File::create(&path).unwrap();
        f.write_all(&[0u8; 16]).unwrap(); // 16-byte header
        for v in &[10.0f64, 20.0, 30.0] {
            f.write_all(&v.to_le_bytes()).unwrap();
        }
        f.flush().unwrap();

        let desc = MmapDataArray::new(&path, 16, 3, 1);
        let arr = desc.load().unwrap();
        assert_eq!(arr.num_tuples(), 3);
        assert_eq!(arr.as_slice(), &[10.0, 20.0, 30.0]);

        let _ = std::fs::remove_dir_all(&dir);
    }

    #[test]
    fn mmap_point_cloud_load() {
        let dir = std::env::temp_dir().join("vtk_mmap_test_pc");
        let _ = std::fs::create_dir_all(&dir);
        let path = dir.join("points.bin");

        let mut f = File::create(&path).unwrap();
        for v in &[1.0f64, 2.0, 3.0, 4.0, 5.0, 6.0] {
            f.write_all(&v.to_le_bytes()).unwrap();
        }
        f.flush().unwrap();

        let pc = MmapPointCloud::new(&path, 0, 2);
        let pd = pc.load().unwrap();
        assert_eq!(pd.points.len(), 2);
        assert_eq!(pd.verts.num_cells(), 2);
        let p0 = pd.points.get(0);
        assert!((p0[0] - 1.0).abs() < 1e-12);
        assert!((p0[1] - 2.0).abs() < 1e-12);
        assert!((p0[2] - 3.0).abs() < 1e-12);

        let _ = std::fs::remove_dir_all(&dir);
    }
}
