use std::io::Write;
use std::path::Path;

use vtk_data::PolyData;
use vtk_types::VtkError;

/// Writer for Wavefront OBJ format.
pub struct ObjWriter;

impl ObjWriter {
    pub fn write(path: &Path, data: &PolyData) -> Result<(), VtkError> {
        let file = std::fs::File::create(path)?;
        let mut w = std::io::BufWriter::new(file);
        Self::write_to(&mut w, data)
    }

    pub fn write_to<W: Write>(w: &mut W, data: &PolyData) -> Result<(), VtkError> {
        writeln!(w, "# vtk-rs OBJ export")?;

        // Vertices
        for i in 0..data.points.len() {
            let p = data.points.get(i);
            writeln!(w, "v {} {} {}", p[0], p[1], p[2])?;
        }

        // Normals (if present)
        let has_normals = data.point_data().normals().is_some();
        if let Some(normals) = data.point_data().normals() {
            let mut buf = [0.0f64; 3];
            for i in 0..normals.num_tuples() {
                normals.tuple_as_f64(i, &mut buf);
                writeln!(w, "vn {} {} {}", buf[0], buf[1], buf[2])?;
            }
        }

        // Texture coordinates (if present)
        let has_tcoords = data.point_data().tcoords().is_some();
        if let Some(tcoords) = data.point_data().tcoords() {
            let nc = tcoords.num_components();
            let mut buf = vec![0.0f64; nc];
            for i in 0..tcoords.num_tuples() {
                tcoords.tuple_as_f64(i, &mut buf);
                if nc >= 2 {
                    writeln!(w, "vt {} {}", buf[0], buf[1])?;
                } else {
                    writeln!(w, "vt {}", buf[0])?;
                }
            }
        }

        // Faces (OBJ uses 1-based indices)
        for cell in data.polys.iter() {
            write!(w, "f")?;
            for &id in cell {
                let idx = id + 1; // 1-based
                if has_normals && has_tcoords {
                    write!(w, " {}/{}/{}", idx, idx, idx)?;
                } else if has_normals {
                    write!(w, " {}//{}", idx, idx)?;
                } else if has_tcoords {
                    write!(w, " {}/{}", idx, idx)?;
                } else {
                    write!(w, " {}", idx)?;
                }
            }
            writeln!(w)?;
        }

        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn write_triangle() {
        let pd = PolyData::from_triangles(
            vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]],
            vec![[0, 1, 2]],
        );
        let mut buf = Vec::new();
        ObjWriter::write_to(&mut buf, &pd).unwrap();
        let output = String::from_utf8(buf).unwrap();

        assert!(output.contains("v 0 0 0"));
        assert!(output.contains("v 1 0 0"));
        assert!(output.contains("f 1 2 3"));
    }
}
