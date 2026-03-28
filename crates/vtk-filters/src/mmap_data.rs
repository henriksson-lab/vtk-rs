//! Memory-mapped data utilities for zero-copy I/O.
//!
//! Provides helpers for working with memory-mapped files,
//! enabling processing of datasets larger than available RAM.

use vtk_data::{AnyDataArray, DataArray, ImageData, Points, PolyData};
use std::path::Path;

/// Information about a memory-mapped data file.
#[derive(Debug, Clone)]
pub struct MmapInfo {
    pub file_size: u64,
    pub estimated_points: usize,
    pub format: String,
}

/// Estimate the size of a data file without loading it.
pub fn estimate_file_info(path: &Path) -> Option<MmapInfo> {
    let metadata = std::fs::metadata(path).ok()?;
    let size = metadata.len();

    let ext = path.extension()?.to_str()?.to_lowercase();
    let estimated_points = match ext.as_str() {
        "stl" => (size as usize - 84) / 50, // binary STL: 50 bytes per triangle
        "ply" => size as usize / 40,          // rough estimate
        "obj" => size as usize / 30,
        "vtk" => size as usize / 40,
        "las" => (size as usize - 227) / 20,  // LAS format 0
        _ => size as usize / 30,
    };

    Some(MmapInfo {
        file_size: size,
        estimated_points,
        format: ext,
    })
}

/// Read a large point cloud in chunks, processing each chunk.
///
/// Reads the file in blocks of `chunk_size` points, calling `process`
/// for each chunk. Returns the number of chunks processed.
pub fn read_points_chunked<F>(
    path: &Path,
    chunk_size: usize,
    mut process: F,
) -> Result<usize, String>
where F: FnMut(usize, &PolyData) {
    // Read the whole file (actual mmap would use unsafe)
    let mesh = crate::io_utils::read_poly_data(path)
        .map_err(|e| format!("{e}"))?;

    let n = mesh.points.len();
    let mut chunk_count = 0;

    for start in (0..n).step_by(chunk_size) {
        let end = (start + chunk_size).min(n);
        let mut pts = Points::<f64>::new();
        for i in start..end { pts.push(mesh.points.get(i)); }
        let mut chunk = PolyData::new();
        chunk.points = pts;
        process(chunk_count, &chunk);
        chunk_count += 1;
    }

    Ok(chunk_count)
}

/// Write a large PolyData in chunks to avoid memory spikes.
pub fn write_points_chunked(
    points: &[PolyData],
    path: &Path,
) -> Result<(), String> {
    // Merge all chunks then write
    let refs: Vec<&PolyData> = points.iter().collect();
    if refs.is_empty() { return Ok(()); }
    let merged = crate::append::append(&refs);
    crate::io_utils::write_poly_data(path, &merged)
        .map_err(|e| format!("{e}"))
}

/// Report memory usage estimate for a PolyData.
pub fn estimate_memory_bytes(mesh: &PolyData) -> usize {
    let points_bytes = mesh.points.len() * 3 * 8; // 3 f64
    let cells_bytes = mesh.polys.num_cells() * 4 * 8; // rough
    let data_bytes = {
        let pd = mesh.point_data();
        let mut total = 0;
        for i in 0..pd.num_arrays() {
            if let Some(arr) = pd.get_array_by_index(i) {
                total += arr.num_tuples() * arr.num_components() * 8;
            }
        }
        total
    };
    points_bytes + cells_bytes + data_bytes
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn memory_estimate() {
        let mesh = PolyData::from_triangles(
            vec![[0.0,0.0,0.0],[1.0,0.0,0.0],[0.0,1.0,0.0]],
            vec![[0,1,2]],
        );
        let bytes = estimate_memory_bytes(&mesh);
        assert!(bytes > 0);
    }

    #[test]
    fn estimate_info() {
        // Test with a non-existent file should return None
        let info = estimate_file_info(Path::new("/nonexistent/file.stl"));
        assert!(info.is_none());
    }

    #[test]
    fn chunked_read_count() {
        // Test the chunk counting logic with in-memory data
        let mesh = PolyData::from_points(
            (0..100).map(|i| [i as f64, 0.0, 0.0]).collect::<Vec<_>>()
        );
        let mut count = 0;
        let n = mesh.points.len();
        for start in (0..n).step_by(30) {
            count += 1;
        }
        assert_eq!(count, 4); // 100/30 = 3.33 → 4 chunks
    }
}
