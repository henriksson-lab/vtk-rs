use std::io::BufRead;
use std::path::Path;

use crate::data::{CellArray, DataArray, Points, PolyData};
use crate::types::VtkError;

/// Reader for Wavefront OBJ format.
pub struct ObjReader;

impl ObjReader {
    pub fn read(path: &Path) -> Result<PolyData, VtkError> {
        let file = std::fs::File::open(path)?;
        let reader = std::io::BufReader::new(file);
        Self::read_from(reader)
    }

    pub fn read_from<R: BufRead>(reader: R) -> Result<PolyData, VtkError> {
        let mut points = Points::<f64>::new();
        let mut normals = DataArray::<f64>::new("Normals", 3);
        let mut tcoords = DataArray::<f64>::new("TCoords", 2);
        let mut polys = CellArray::new();

        let mut has_normals = false;
        let mut has_tcoords = false;

        for line in reader.lines() {
            let line = line.map_err(VtkError::Io)?;
            let trimmed = line.trim();

            if trimmed.is_empty() || trimmed.starts_with('#') {
                continue;
            }

            let mut parts = trimmed.split_whitespace();
            let Some(keyword) = parts.next() else {
                continue;
            };

            match keyword {
                "v" => {
                    let coords: Vec<f64> = parts.filter_map(|s| s.parse().ok()).collect();
                    if coords.len() >= 3 {
                        points.push([coords[0], coords[1], coords[2]]);
                    }
                }
                "vn" => {
                    let coords: Vec<f64> = parts.filter_map(|s| s.parse().ok()).collect();
                    if coords.len() >= 3 {
                        normals.push_tuple(&[coords[0], coords[1], coords[2]]);
                        has_normals = true;
                    }
                }
                "vt" => {
                    let coords: Vec<f64> = parts.filter_map(|s| s.parse().ok()).collect();
                    if coords.len() >= 2 {
                        tcoords.push_tuple(&[coords[0], coords[1]]);
                    } else if coords.len() == 1 {
                        tcoords.push_tuple(&[coords[0], 0.0]);
                    }
                    has_tcoords = true;
                }
                "f" => {
                    let indices: Vec<i64> = parts
                        .filter_map(|s| {
                            // Parse "v", "v/vt", "v/vt/vn", or "v//vn"
                            let idx_str = s.split('/').next()?;
                            let idx: i64 = idx_str.parse().ok()?;
                            Some(idx - 1) // Convert 1-based to 0-based
                        })
                        .collect();
                    if indices.len() >= 3 {
                        polys.push_cell(&indices);
                    }
                }
                _ => {
                    // Skip: mtllib, usemtl, g, s, o, etc.
                }
            }
        }

        let mut pd = PolyData::new();
        pd.points = points;
        pd.polys = polys;

        if has_normals && normals.num_tuples() > 0 {
            pd.point_data_mut().add_array(normals.into());
            pd.point_data_mut().set_active_normals("Normals");
        }
        if has_tcoords && tcoords.num_tuples() > 0 {
            pd.point_data_mut().add_array(tcoords.into());
            pd.point_data_mut().set_active_tcoords("TCoords");
        }

        Ok(pd)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::io::obj::ObjWriter;

    #[test]
    fn roundtrip() {
        let pd = PolyData::from_triangles(
            vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]],
            vec![[0, 1, 2]],
        );

        let mut buf = Vec::new();
        ObjWriter::write_to(&mut buf, &pd).unwrap();

        let reader = std::io::BufReader::new(&buf[..]);
        let result = ObjReader::read_from(reader).unwrap();

        assert_eq!(result.points.len(), 3);
        assert_eq!(result.polys.num_cells(), 1);
        assert_eq!(result.polys.cell(0), &[0, 1, 2]);
    }

    #[test]
    fn parse_complex_faces() {
        let obj = b"# test\nv 0 0 0\nv 1 0 0\nv 1 1 0\nv 0 1 0\nf 1/1/1 2/2/2 3/3/3 4/4/4\n";
        let reader = std::io::BufReader::new(&obj[..]);
        let result = ObjReader::read_from(reader).unwrap();

        assert_eq!(result.points.len(), 4);
        assert_eq!(result.polys.num_cells(), 1);
        assert_eq!(result.polys.cell(0), &[0, 1, 2, 3]);
    }
}
