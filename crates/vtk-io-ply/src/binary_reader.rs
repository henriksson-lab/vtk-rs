use std::io::BufRead;
use std::path::Path;

use vtk_data::{CellArray, Points, PolyData};
use vtk_types::VtkError;

/// Reader for Stanford PLY format (binary little-endian).
pub struct PlyBinaryReader;

impl PlyBinaryReader {
    pub fn read(path: &Path) -> Result<PolyData, VtkError> {
        let file = std::fs::File::open(path)?;
        let reader = std::io::BufReader::new(file);
        Self::read_from(reader)
    }

    pub fn read_from<R: BufRead>(mut reader: R) -> Result<PolyData, VtkError> {
        let mut n_vertices = 0usize;
        let mut n_faces = 0usize;
        let mut is_binary_le = false;

        // Parse ASCII header
        loop {
            let mut line = String::new();
            reader.read_line(&mut line).map_err(VtkError::Io)?;
            let trimmed = line.trim();

            if trimmed == "end_header" {
                break;
            }

            let parts: Vec<&str> = trimmed.split_whitespace().collect();
            if parts.is_empty() {
                continue;
            }

            match parts[0] {
                "format" => {
                    if parts.get(1) == Some(&"binary_little_endian") {
                        is_binary_le = true;
                    }
                }
                "element" => {
                    if parts.len() >= 3 {
                        let count: usize = parts[2].parse().unwrap_or(0);
                        if parts[1] == "vertex" {
                            n_vertices = count;
                        } else if parts[1] == "face" {
                            n_faces = count;
                        }
                    }
                }
                _ => {}
            }
        }

        if !is_binary_le {
            return Err(VtkError::Unsupported(
                "only binary_little_endian PLY supported by PlyBinaryReader".into(),
            ));
        }

        // Read binary vertex data (3 × f32 per vertex)
        let mut points = Points::<f64>::new();
        let mut vert_buf = vec![0u8; n_vertices * 12];
        reader.read_exact(&mut vert_buf).map_err(VtkError::Io)?;
        for i in 0..n_vertices {
            let base = i * 12;
            let x = f32::from_le_bytes(vert_buf[base..base + 4].try_into().unwrap()) as f64;
            let y = f32::from_le_bytes(vert_buf[base + 4..base + 8].try_into().unwrap()) as f64;
            let z = f32::from_le_bytes(vert_buf[base + 8..base + 12].try_into().unwrap()) as f64;
            points.push([x, y, z]);
        }

        // Read binary face data (u8 count + count × i32)
        let mut polys = CellArray::new();
        for _ in 0..n_faces {
            let mut count_buf = [0u8; 1];
            reader.read_exact(&mut count_buf).map_err(VtkError::Io)?;
            let n = count_buf[0] as usize;
            let mut face_buf = vec![0u8; n * 4];
            reader.read_exact(&mut face_buf).map_err(VtkError::Io)?;
            let ids: Vec<i64> = face_buf
                .chunks_exact(4)
                .map(|c| i32::from_le_bytes(c.try_into().unwrap()) as i64)
                .collect();
            polys.push_cell(&ids);
        }

        let mut pd = PolyData::new();
        pd.points = points;
        pd.polys = polys;
        Ok(pd)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::PlyBinaryWriter;

    #[test]
    fn roundtrip_binary_ply() {
        let pd = PolyData::from_triangles(
            vec![[1.0, 2.0, 3.0], [4.0, 5.0, 6.0], [7.0, 8.0, 9.0]],
            vec![[0, 1, 2]],
        );

        let mut buf = Vec::new();
        PlyBinaryWriter::write_to(&mut buf, &pd).unwrap();

        let reader = std::io::BufReader::new(&buf[..]);
        let result = PlyBinaryReader::read_from(reader).unwrap();

        assert_eq!(result.points.len(), 3);
        assert_eq!(result.polys.num_cells(), 1);
        assert_eq!(result.polys.cell(0), &[0, 1, 2]);

        let p0 = result.points.get(0);
        assert!((p0[0] - 1.0).abs() < 1e-4); // f32 precision
        assert!((p0[1] - 2.0).abs() < 1e-4);
    }

    #[test]
    fn roundtrip_binary_quad() {
        let mut pd = PolyData::new();
        pd.points.push([0.0, 0.0, 0.0]);
        pd.points.push([1.0, 0.0, 0.0]);
        pd.points.push([1.0, 1.0, 0.0]);
        pd.points.push([0.0, 1.0, 0.0]);
        pd.polys.push_cell(&[0, 1, 2, 3]);

        let mut buf = Vec::new();
        PlyBinaryWriter::write_to(&mut buf, &pd).unwrap();

        let reader = std::io::BufReader::new(&buf[..]);
        let result = PlyBinaryReader::read_from(reader).unwrap();

        assert_eq!(result.polys.num_cells(), 1);
        assert_eq!(result.polys.cell(0), &[0, 1, 2, 3]);
    }
}
