use std::io::BufRead;
use std::path::Path;

use vtk_data::{AnyDataArray, CellArray, DataArray, PolyData, Points};
use vtk_types::VtkError;

/// Reader for VTK XML PolyData format (.vtp).
///
/// Reads ASCII VTK XML files. Binary/appended data is not yet supported.
pub struct VtpReader;

impl VtpReader {
    pub fn read(path: &Path) -> Result<PolyData, VtkError> {
        let file = std::fs::File::open(path)?;
        let reader = std::io::BufReader::new(file);
        Self::read_from(reader)
    }

    pub fn read_from<R: BufRead>(reader: R) -> Result<PolyData, VtkError> {
        let content: String = reader
            .lines()
            .collect::<Result<Vec<_>, _>>()
            .map_err(VtkError::Io)?
            .join("\n");

        // Simple XML parsing — we handle the subset VTK XML uses
        let mut pd = PolyData::new();

        // Extract Points
        if let Some(points_data) = extract_data_array_content(&content, "Points") {
            pd.points = parse_points(&points_data)?;
        }

        // Extract cell sections
        if let Some(polys_section) = extract_section(&content, "Polys") {
            pd.polys = parse_cell_section(&polys_section)?;
        }
        if let Some(lines_section) = extract_section(&content, "Lines") {
            pd.lines = parse_cell_section(&lines_section)?;
        }
        if let Some(verts_section) = extract_section(&content, "Verts") {
            pd.verts = parse_cell_section(&verts_section)?;
        }
        if let Some(strips_section) = extract_section(&content, "Strips") {
            pd.strips = parse_cell_section(&strips_section)?;
        }

        // Extract PointData
        if let Some(pd_section) = extract_section(&content, "PointData") {
            parse_attribute_arrays(&pd_section, pd.point_data_mut())?;
        }

        // Extract CellData
        if let Some(cd_section) = extract_section(&content, "CellData") {
            parse_attribute_arrays(&cd_section, pd.cell_data_mut())?;
        }

        Ok(pd)
    }
}

/// Extract content between <tag> and </tag>.
fn extract_section(content: &str, tag: &str) -> Option<String> {
    // Handle both <Tag> and <Tag attrs...>
    let open_pattern = format!("<{}", tag);
    let close_pattern = format!("</{}>", tag);

    let start = content.find(&open_pattern)?;
    let after_open = &content[start..];
    let tag_end = after_open.find('>')?;
    let content_start = start + tag_end + 1;

    let end = content[content_start..].find(&close_pattern)?;
    Some(content[content_start..content_start + end].to_string())
}

/// Extract the first DataArray content inside a section (for Points).
fn extract_data_array_content(content: &str, parent_tag: &str) -> Option<String> {
    let section = extract_section(content, parent_tag)?;
    extract_inner_data_array_values(&section)
}

fn extract_inner_data_array_values(section: &str) -> Option<String> {
    let start = section.find("<DataArray")?;
    let after = &section[start..];
    let tag_end = after.find('>')?;
    let content_start = start + tag_end + 1;
    let end = section[content_start..].find("</DataArray>")?;
    Some(section[content_start..content_start + end].trim().to_string())
}

/// Parse attribute for a DataArray tag.
fn extract_attr(tag: &str, attr_name: &str) -> Option<String> {
    let pattern = format!("{}=\"", attr_name);
    let start = tag.find(&pattern)?;
    let value_start = start + pattern.len();
    let end = tag[value_start..].find('"')?;
    Some(tag[value_start..value_start + end].to_string())
}

fn parse_points(data: &str) -> Result<Points<f64>, VtkError> {
    let values: Vec<f64> = data
        .split_whitespace()
        .filter_map(|s| s.parse().ok())
        .collect();
    let mut pts = Points::new();
    for chunk in values.chunks(3) {
        if chunk.len() == 3 {
            pts.push([chunk[0], chunk[1], chunk[2]]);
        }
    }
    Ok(pts)
}

fn parse_cell_section(section: &str) -> Result<CellArray, VtkError> {
    // Find connectivity and offsets DataArrays
    let mut connectivity_values: Vec<i64> = Vec::new();
    let mut offsets_values: Vec<i64> = Vec::new();

    // Find all DataArray tags in this section
    let mut search_pos = 0;
    while let Some(da_start) = section[search_pos..].find("<DataArray") {
        let abs_start = search_pos + da_start;
        let tag_end = section[abs_start..]
            .find('>')
            .ok_or_else(|| VtkError::Parse("unclosed DataArray tag".into()))?;
        let tag = &section[abs_start..abs_start + tag_end + 1];

        let content_start = abs_start + tag_end + 1;
        let content_end = section[content_start..]
            .find("</DataArray>")
            .ok_or_else(|| VtkError::Parse("missing </DataArray>".into()))?;
        let content = section[content_start..content_start + content_end].trim();

        let name = extract_attr(tag, "Name").unwrap_or_default();
        let values: Vec<i64> = content
            .split_whitespace()
            .filter_map(|s| s.parse().ok())
            .collect();

        match name.as_str() {
            "connectivity" => connectivity_values = values,
            "offsets" => offsets_values = values,
            _ => {}
        }

        search_pos = content_start + content_end + "</DataArray>".len();
    }

    // Build CellArray from connectivity + offsets
    let mut cells = CellArray::new();
    let mut prev_offset = 0;
    for &offset in &offsets_values {
        let start = prev_offset as usize;
        let end = offset as usize;
        if end <= connectivity_values.len() && start < end {
            cells.push_cell(&connectivity_values[start..end]);
        }
        prev_offset = offset;
    }

    Ok(cells)
}

fn parse_attribute_arrays(
    section: &str,
    attrs: &mut vtk_data::DataSetAttributes,
) -> Result<(), VtkError> {
    let mut search_pos = 0;
    while let Some(da_start) = section[search_pos..].find("<DataArray") {
        let abs_start = search_pos + da_start;
        let tag_end = section[abs_start..]
            .find('>')
            .ok_or_else(|| VtkError::Parse("unclosed DataArray tag".into()))?;
        let tag = &section[abs_start..abs_start + tag_end + 1];

        let content_start = abs_start + tag_end + 1;
        let content_end = section[content_start..]
            .find("</DataArray>")
            .ok_or_else(|| VtkError::Parse("missing </DataArray>".into()))?;
        let content = section[content_start..content_start + content_end].trim();

        let name = extract_attr(tag, "Name").unwrap_or_else(|| "data".to_string());
        let type_str = extract_attr(tag, "type").unwrap_or_else(|| "Float64".to_string());
        let nc: usize = extract_attr(tag, "NumberOfComponents")
            .and_then(|s| s.parse().ok())
            .unwrap_or(1);

        let values: Vec<f64> = content
            .split_whitespace()
            .filter_map(|s| s.parse().ok())
            .collect();

        let arr = match type_str.as_str() {
            "Float32" => {
                let data: Vec<f32> = values.iter().map(|&v| v as f32).collect();
                AnyDataArray::F32(DataArray::from_vec(&name, data, nc))
            }
            "Int32" => {
                let data: Vec<i32> = values.iter().map(|&v| v as i32).collect();
                AnyDataArray::I32(DataArray::from_vec(&name, data, nc))
            }
            "Int64" => {
                let data: Vec<i64> = values.iter().map(|&v| v as i64).collect();
                AnyDataArray::I64(DataArray::from_vec(&name, data, nc))
            }
            "UInt8" => {
                let data: Vec<u8> = values.iter().map(|&v| v as u8).collect();
                AnyDataArray::U8(DataArray::from_vec(&name, data, nc))
            }
            _ => AnyDataArray::F64(DataArray::from_vec(&name, values, nc)),
        };
        let arr_name = arr.name().to_string();
        attrs.add_array(arr);
        if attrs.scalars().is_none() {
            attrs.set_active_scalars(&arr_name);
        }

        search_pos = content_start + content_end + "</DataArray>".len();
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::VtpWriter;
    use vtk_data::DataArray as DA;

    #[test]
    fn roundtrip_vtp_triangle() {
        let pd = PolyData::from_triangles(
            vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]],
            vec![[0, 1, 2]],
        );

        let mut buf = Vec::new();
        VtpWriter::write_to(&mut buf, &pd).unwrap();

        let reader = std::io::BufReader::new(&buf[..]);
        let result = VtpReader::read_from(reader).unwrap();

        assert_eq!(result.points.len(), 3);
        assert_eq!(result.polys.num_cells(), 1);
        assert_eq!(result.polys.cell(0), &[0, 1, 2]);
    }

    #[test]
    fn roundtrip_vtp_with_scalars() {
        let mut pd = PolyData::from_triangles(
            vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]],
            vec![[0, 1, 2]],
        );
        let scalars = DA::from_vec("temperature", vec![10.0f64, 20.0, 30.0], 1);
        pd.point_data_mut().add_array(scalars.into());
        pd.point_data_mut().set_active_scalars("temperature");

        let mut buf = Vec::new();
        VtpWriter::write_to(&mut buf, &pd).unwrap();

        let reader = std::io::BufReader::new(&buf[..]);
        let result = VtpReader::read_from(reader).unwrap();

        let s = result.point_data().scalars().unwrap();
        assert_eq!(s.name(), "temperature");
        assert_eq!(s.num_tuples(), 3);
    }

    #[test]
    fn roundtrip_vtp_quad() {
        let pd = PolyData::from_triangles(
            vec![
                [0.0, 0.0, 0.0],
                [1.0, 0.0, 0.0],
                [1.0, 1.0, 0.0],
                [0.0, 1.0, 0.0],
                [2.0, 0.0, 0.0],
                [2.0, 1.0, 0.0],
            ],
            vec![[0, 1, 2], [0, 2, 3], [1, 4, 5], [1, 5, 2]],
        );

        let mut buf = Vec::new();
        VtpWriter::write_to(&mut buf, &pd).unwrap();

        let reader = std::io::BufReader::new(&buf[..]);
        let result = VtpReader::read_from(reader).unwrap();

        assert_eq!(result.points.len(), 6);
        assert_eq!(result.polys.num_cells(), 4);
    }
}
