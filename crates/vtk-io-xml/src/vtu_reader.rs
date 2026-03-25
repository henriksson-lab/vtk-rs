use std::io::BufRead;
use std::path::Path;

use vtk_data::{AnyDataArray, DataArray, UnstructuredGrid};
use vtk_types::{CellType, VtkError};

/// Reader for VTK XML UnstructuredGrid format (.vtu).
///
/// Reads ASCII VTK XML files. Binary/appended data is not yet supported.
pub struct VtuReader;

impl VtuReader {
    pub fn read(path: &Path) -> Result<UnstructuredGrid, VtkError> {
        let file = std::fs::File::open(path)?;
        let reader = std::io::BufReader::new(file);
        Self::read_from(reader)
    }

    pub fn read_from<R: BufRead>(reader: R) -> Result<UnstructuredGrid, VtkError> {
        let content: String = reader
            .lines()
            .collect::<Result<Vec<_>, _>>()
            .map_err(VtkError::Io)?
            .join("\n");

        let mut grid = UnstructuredGrid::new();

        // Extract Points
        if let Some(points_data) = extract_data_array_content(&content, "Points") {
            let values: Vec<f64> = points_data
                .split_whitespace()
                .filter_map(|s| s.parse().ok())
                .collect();
            for chunk in values.chunks(3) {
                if chunk.len() == 3 {
                    grid.points.push([chunk[0], chunk[1], chunk[2]]);
                }
            }
        }

        // Extract Cells
        if let Some(cells_section) = extract_section(&content, "Cells") {
            let mut connectivity: Vec<i64> = Vec::new();
            let mut offsets: Vec<i64> = Vec::new();
            let mut types: Vec<u8> = Vec::new();

            let mut search_pos = 0;
            while let Some(da_start) = cells_section[search_pos..].find("<DataArray") {
                let abs_start = search_pos + da_start;
                let tag_end = cells_section[abs_start..]
                    .find('>')
                    .ok_or_else(|| VtkError::Parse("unclosed DataArray tag".into()))?;
                let tag = &cells_section[abs_start..abs_start + tag_end + 1];

                let content_start = abs_start + tag_end + 1;
                let content_end = cells_section[content_start..]
                    .find("</DataArray>")
                    .ok_or_else(|| VtkError::Parse("missing </DataArray>".into()))?;
                let da_content = cells_section[content_start..content_start + content_end].trim();

                let name = extract_attr(tag, "Name").unwrap_or_default();
                match name.as_str() {
                    "connectivity" => {
                        connectivity = da_content
                            .split_whitespace()
                            .filter_map(|s| s.parse().ok())
                            .collect();
                    }
                    "offsets" => {
                        offsets = da_content
                            .split_whitespace()
                            .filter_map(|s| s.parse().ok())
                            .collect();
                    }
                    "types" => {
                        types = da_content
                            .split_whitespace()
                            .filter_map(|s| s.parse().ok())
                            .collect();
                    }
                    _ => {}
                }

                search_pos = content_start + content_end + "</DataArray>".len();
            }

            // Build cells from connectivity + offsets + types
            let mut prev_offset: usize = 0;
            for (i, &offset) in offsets.iter().enumerate() {
                let end = offset as usize;
                if end <= connectivity.len() && prev_offset < end {
                    let cell_pts = &connectivity[prev_offset..end];
                    let ct = types
                        .get(i)
                        .and_then(|&v| CellType::from_u8(v))
                        .unwrap_or(CellType::Triangle);
                    grid.push_cell(ct, cell_pts);
                }
                prev_offset = end;
            }
        }

        // Extract PointData
        if let Some(pd_section) = extract_section(&content, "PointData") {
            parse_attribute_arrays(&pd_section, grid.point_data_mut())?;
        }

        // Extract CellData
        if let Some(cd_section) = extract_section(&content, "CellData") {
            parse_attribute_arrays(&cd_section, grid.cell_data_mut())?;
        }

        Ok(grid)
    }
}

fn extract_section(content: &str, tag: &str) -> Option<String> {
    let open_pattern = format!("<{}", tag);
    let close_pattern = format!("</{}>", tag);

    let start = content.find(&open_pattern)?;
    let after_open = &content[start..];
    let tag_end = after_open.find('>')?;
    let content_start = start + tag_end + 1;

    let end = content[content_start..].find(&close_pattern)?;
    Some(content[content_start..content_start + end].to_string())
}

fn extract_data_array_content(content: &str, parent_tag: &str) -> Option<String> {
    let section = extract_section(content, parent_tag)?;
    let start = section.find("<DataArray")?;
    let after = &section[start..];
    let tag_end = after.find('>')?;
    let content_start = start + tag_end + 1;
    let end = section[content_start..].find("</DataArray>")?;
    Some(section[content_start..content_start + end].trim().to_string())
}

fn extract_attr(tag: &str, attr_name: &str) -> Option<String> {
    let pattern = format!("{}=\"", attr_name);
    let start = tag.find(&pattern)?;
    let value_start = start + pattern.len();
    let end = tag[value_start..].find('"')?;
    Some(tag[value_start..value_start + end].to_string())
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
    use crate::VtuWriter;
    use vtk_data::{DataArray as DA, DataSet};
    use vtk_types::CellType;

    #[test]
    fn roundtrip_vtu_tetra() {
        let mut grid = UnstructuredGrid::new();
        grid.points.push([0.0, 0.0, 0.0]);
        grid.points.push([1.0, 0.0, 0.0]);
        grid.points.push([0.5, 1.0, 0.0]);
        grid.points.push([0.5, 0.5, 1.0]);
        grid.push_cell(CellType::Tetra, &[0, 1, 2, 3]);

        let mut buf = Vec::new();
        VtuWriter::write_to(&mut buf, &grid).unwrap();

        let reader = std::io::BufReader::new(&buf[..]);
        let result = VtuReader::read_from(reader).unwrap();

        assert_eq!(result.num_points(), 4);
        assert_eq!(result.num_cells(), 1);
        assert_eq!(result.cell_type(0), CellType::Tetra);
        assert_eq!(result.cell_points(0), &[0, 1, 2, 3]);
    }

    #[test]
    fn roundtrip_vtu_mixed() {
        let mut grid = UnstructuredGrid::new();
        grid.points.push([0.0, 0.0, 0.0]);
        grid.points.push([1.0, 0.0, 0.0]);
        grid.points.push([0.5, 1.0, 0.0]);
        grid.points.push([0.5, 0.5, 1.0]);
        grid.points.push([2.0, 0.0, 0.0]);

        grid.push_cell(CellType::Tetra, &[0, 1, 2, 3]);
        grid.push_cell(CellType::Triangle, &[1, 4, 2]);

        let mut buf = Vec::new();
        VtuWriter::write_to(&mut buf, &grid).unwrap();

        let reader = std::io::BufReader::new(&buf[..]);
        let result = VtuReader::read_from(reader).unwrap();

        assert_eq!(result.num_cells(), 2);
        assert_eq!(result.cell_type(0), CellType::Tetra);
        assert_eq!(result.cell_type(1), CellType::Triangle);
    }

    #[test]
    fn roundtrip_vtu_with_scalars() {
        let mut grid = UnstructuredGrid::new();
        grid.points.push([0.0, 0.0, 0.0]);
        grid.points.push([1.0, 0.0, 0.0]);
        grid.points.push([0.5, 1.0, 0.0]);
        grid.push_cell(CellType::Triangle, &[0, 1, 2]);

        let scalars = DA::from_vec("temp", vec![10.0, 20.0, 30.0], 1);
        grid.point_data_mut().add_array(scalars.into());
        grid.point_data_mut().set_active_scalars("temp");

        let mut buf = Vec::new();
        VtuWriter::write_to(&mut buf, &grid).unwrap();

        let reader = std::io::BufReader::new(&buf[..]);
        let result = VtuReader::read_from(reader).unwrap();

        let s = result.point_data().scalars().unwrap();
        assert_eq!(s.name(), "temp");
        assert_eq!(s.num_tuples(), 3);
    }
}
