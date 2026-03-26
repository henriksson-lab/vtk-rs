use std::io::BufRead;
use std::path::Path;

use vtk_data::{AnyDataArray, DataArray, Points, StructuredGrid};
use vtk_types::VtkError;

/// Reader for VTK XML StructuredGrid format (.vts).
pub struct VtsReader;

impl VtsReader {
    pub fn read(path: &Path) -> Result<StructuredGrid, VtkError> {
        let file = std::fs::File::open(path)?;
        let reader = std::io::BufReader::new(file);
        Self::read_from(reader)
    }

    pub fn read_from<R: BufRead>(reader: R) -> Result<StructuredGrid, VtkError> {
        let content: String = reader
            .lines()
            .collect::<Result<Vec<_>, _>>()
            .map_err(VtkError::Io)?
            .join("\n");

        // Parse extent from StructuredGrid tag
        let mut dims = [1usize; 3];
        if let Some(sg_tag) = find_tag(&content, "StructuredGrid") {
            if let Some(ext_str) = extract_attr(&sg_tag, "WholeExtent") {
                let vals: Vec<i64> = ext_str.split_whitespace().filter_map(|s| s.parse().ok()).collect();
                if vals.len() == 6 {
                    dims[0] = (vals[1] - vals[0] + 1).max(1) as usize;
                    dims[1] = (vals[3] - vals[2] + 1).max(1) as usize;
                    dims[2] = (vals[5] - vals[4] + 1).max(1) as usize;
                }
            }
        }

        // Parse Points
        let mut points = Points::new();
        if let Some(pts_data) = extract_data_array_content(&content, "Points") {
            let values: Vec<f64> = pts_data.split_whitespace().filter_map(|s| s.parse().ok()).collect();
            for chunk in values.chunks(3) {
                if chunk.len() == 3 {
                    points.push([chunk[0], chunk[1], chunk[2]]);
                }
            }
        }

        let expected = dims[0] * dims[1] * dims[2];
        if points.len() != expected {
            return Err(VtkError::Parse(format!(
                "expected {} points for dims {:?}, got {}",
                expected, dims, points.len()
            )));
        }

        let mut grid = StructuredGrid::from_dimensions_and_points(dims, points);

        // Parse PointData
        if let Some(pd_section) = extract_section(&content, "PointData") {
            parse_attribute_arrays(&pd_section, grid.point_data_mut())?;
        }

        // Parse CellData
        if let Some(cd_section) = extract_section(&content, "CellData") {
            parse_attribute_arrays(&cd_section, grid.cell_data_mut())?;
        }

        Ok(grid)
    }
}

fn find_tag(content: &str, tag: &str) -> Option<String> {
    let pat = format!("<{}", tag);
    let start = content.find(&pat)?;
    let end = content[start..].find('>')?;
    Some(content[start..start + end + 1].to_string())
}

fn extract_section(content: &str, tag: &str) -> Option<String> {
    let open = format!("<{}", tag);
    let close = format!("</{}>", tag);
    let start = content.find(&open)?;
    let te = content[start..].find('>')?;
    let cs = start + te + 1;
    let end = content[cs..].find(&close)?;
    Some(content[cs..cs + end].to_string())
}

fn extract_data_array_content(content: &str, parent: &str) -> Option<String> {
    let section = extract_section(content, parent)?;
    let start = section.find("<DataArray")?;
    let te = section[start..].find('>')?;
    let cs = start + te + 1;
    let end = section[cs..].find("</DataArray>")?;
    Some(section[cs..cs + end].trim().to_string())
}

fn extract_attr(tag: &str, name: &str) -> Option<String> {
    let pat = format!("{}=\"", name);
    let start = tag.find(&pat)?;
    let vs = start + pat.len();
    let end = tag[vs..].find('"')?;
    Some(tag[vs..vs + end].to_string())
}

fn parse_attribute_arrays(section: &str, attrs: &mut vtk_data::DataSetAttributes) -> Result<(), VtkError> {
    let mut sp = 0;
    while let Some(ds) = section[sp..].find("<DataArray") {
        let abs = sp + ds;
        let te = section[abs..].find('>').ok_or_else(|| VtkError::Parse("unclosed tag".into()))?;
        let tag = &section[abs..abs + te + 1];
        let cs = abs + te + 1;
        let ce = section[cs..].find("</DataArray>").ok_or_else(|| VtkError::Parse("missing close".into()))?;
        let content = section[cs..cs + ce].trim();

        let name = extract_attr(tag, "Name").unwrap_or_else(|| "data".to_string());
        let nc: usize = extract_attr(tag, "NumberOfComponents").and_then(|s| s.parse().ok()).unwrap_or(1);
        let values: Vec<f64> = content.split_whitespace().filter_map(|s| s.parse().ok()).collect();

        let arr = AnyDataArray::F64(DataArray::from_vec(&name, values, nc));
        let arr_name = arr.name().to_string();
        attrs.add_array(arr);
        if attrs.scalars().is_none() { attrs.set_active_scalars(&arr_name); }

        sp = cs + ce + "</DataArray>".len();
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::VtsWriter;
    use vtk_data::{DataArray as DA, DataSet};

    #[test]
    fn roundtrip_vts() {
        let mut pts = Points::new();
        for j in 0..2 {
            for i in 0..3 {
                pts.push([i as f64, j as f64, 0.0]);
            }
        }
        let mut grid = StructuredGrid::from_dimensions_and_points([3, 2, 1], pts);
        let scalars: Vec<f64> = (0..6).map(|i| i as f64).collect();
        grid.point_data_mut().add_array(DA::from_vec("idx", scalars, 1).into());
        grid.point_data_mut().set_active_scalars("idx");

        let mut buf = Vec::new();
        VtsWriter::write_to(&mut buf, &grid).unwrap();

        let reader = std::io::BufReader::new(&buf[..]);
        let result = VtsReader::read_from(reader).unwrap();

        assert_eq!(result.dimensions(), [3, 2, 1]);
        assert_eq!(result.num_points(), 6);
        let s = result.point_data().scalars().unwrap();
        assert_eq!(s.num_tuples(), 6);
    }
}
