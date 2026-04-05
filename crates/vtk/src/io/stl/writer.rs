use std::io::Write;
use std::path::Path;

use crate::data::PolyData;
use crate::types::VtkError;

/// File type for STL format.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum StlFormat {
    Ascii,
    Binary,
}

/// Writer for STL (stereolithography) format.
///
/// STL only supports triangle meshes. Non-triangle polygons are fan-triangulated
/// on export. Only polygon cells are exported (vertices, lines, strips are ignored).
pub struct StlWriter {
    pub format: StlFormat,
    pub solid_name: String,
}

impl Default for StlWriter {
    fn default() -> Self {
        Self {
            format: StlFormat::Binary,
            solid_name: "vtk-rs".to_string(),
        }
    }
}

impl StlWriter {
    pub fn ascii() -> Self {
        Self {
            format: StlFormat::Ascii,
            ..Default::default()
        }
    }

    pub fn binary() -> Self {
        Self::default()
    }

    pub fn write(&self, path: &Path, data: &PolyData) -> Result<(), VtkError> {
        let file = std::fs::File::create(path)?;
        let mut w = std::io::BufWriter::new(file);
        self.write_to(&mut w, data)
    }

    pub fn write_to<W: Write>(&self, w: &mut W, data: &PolyData) -> Result<(), VtkError> {
        match self.format {
            StlFormat::Ascii => self.write_ascii(w, data),
            StlFormat::Binary => self.write_binary(w, data),
        }
    }

    fn write_ascii<W: Write>(&self, w: &mut W, data: &PolyData) -> Result<(), VtkError> {
        writeln!(w, "solid {}", self.solid_name)?;

        for cell in data.polys.iter() {
            if cell.len() < 3 {
                continue;
            }
            let p0 = data.points.get(cell[0] as usize);
            for i in 1..cell.len() - 1 {
                let p1 = data.points.get(cell[i] as usize);
                let p2 = data.points.get(cell[i + 1] as usize);
                let n = triangle_normal(p0, p1, p2);
                writeln!(w, "  facet normal {} {} {}", n[0], n[1], n[2])?;
                writeln!(w, "    outer loop")?;
                writeln!(w, "      vertex {} {} {}", p0[0], p0[1], p0[2])?;
                writeln!(w, "      vertex {} {} {}", p1[0], p1[1], p1[2])?;
                writeln!(w, "      vertex {} {} {}", p2[0], p2[1], p2[2])?;
                writeln!(w, "    endloop")?;
                writeln!(w, "  endfacet")?;
            }
        }

        writeln!(w, "endsolid {}", self.solid_name)?;
        Ok(())
    }

    fn write_binary<W: Write>(&self, w: &mut W, data: &PolyData) -> Result<(), VtkError> {
        // 80-byte header
        let mut header = [0u8; 80];
        let name_bytes = self.solid_name.as_bytes();
        let len = name_bytes.len().min(80);
        header[..len].copy_from_slice(&name_bytes[..len]);
        w.write_all(&header)?;

        // Count total triangles (fan-triangulating larger polygons)
        let num_triangles: u32 = data
            .polys
            .iter()
            .map(|cell| {
                if cell.len() >= 3 {
                    (cell.len() - 2) as u32
                } else {
                    0
                }
            })
            .sum();
        w.write_all(&num_triangles.to_le_bytes())?;

        // Write triangles
        for cell in data.polys.iter() {
            if cell.len() < 3 {
                continue;
            }
            let p0 = data.points.get(cell[0] as usize);
            for i in 1..cell.len() - 1 {
                let p1 = data.points.get(cell[i] as usize);
                let p2 = data.points.get(cell[i + 1] as usize);
                let n = triangle_normal(p0, p1, p2);

                // Normal (3 x f32 LE)
                for &v in &n {
                    w.write_all(&(v as f32).to_le_bytes())?;
                }
                // Vertices (3 x 3 x f32 LE)
                for p in &[p0, p1, p2] {
                    for &v in p {
                        w.write_all(&(v as f32).to_le_bytes())?;
                    }
                }
                // Attribute byte count
                w.write_all(&0u16.to_le_bytes())?;
            }
        }

        Ok(())
    }
}

fn triangle_normal(p0: [f64; 3], p1: [f64; 3], p2: [f64; 3]) -> [f64; 3] {
    let e1 = [p1[0] - p0[0], p1[1] - p0[1], p1[2] - p0[2]];
    let e2 = [p2[0] - p0[0], p2[1] - p0[1], p2[2] - p0[2]];
    let n = [
        e1[1] * e2[2] - e1[2] * e2[1],
        e1[2] * e2[0] - e1[0] * e2[2],
        e1[0] * e2[1] - e1[1] * e2[0],
    ];
    let len = (n[0] * n[0] + n[1] * n[1] + n[2] * n[2]).sqrt();
    if len > 1e-10 {
        [n[0] / len, n[1] / len, n[2] / len]
    } else {
        [0.0, 0.0, 0.0]
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn write_ascii_triangle() {
        let pd = PolyData::from_triangles(
            vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]],
            vec![[0, 1, 2]],
        );
        let writer = StlWriter::ascii();
        let mut buf = Vec::new();
        writer.write_to(&mut buf, &pd).unwrap();
        let output = String::from_utf8(buf).unwrap();
        assert!(output.starts_with("solid vtk-rs"));
        assert!(output.contains("facet normal"));
        assert!(output.contains("vertex"));
        assert!(output.contains("endsolid"));
    }

    #[test]
    fn write_binary_triangle() {
        let pd = PolyData::from_triangles(
            vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]],
            vec![[0, 1, 2]],
        );
        let writer = StlWriter::binary();
        let mut buf = Vec::new();
        writer.write_to(&mut buf, &pd).unwrap();
        // 80 header + 4 count + 50 per triangle = 134
        assert_eq!(buf.len(), 134);
    }
}
