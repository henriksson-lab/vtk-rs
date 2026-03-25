use std::io::BufRead;
use std::path::Path;

use vtk_data::{AnyDataArray, DataArray, ImageData};
use vtk_types::VtkError;

/// Reader for VTK XML ImageData format (.vti).
pub struct VtiReader;

impl VtiReader {
    pub fn read(path: &Path) -> Result<ImageData, VtkError> {
        let file = std::fs::File::open(path)?;
        let reader = std::io::BufReader::new(file);
        Self::read_from(reader)
    }

    pub fn read_from<R: BufRead>(reader: R) -> Result<ImageData, VtkError> {
        let content: String = reader
            .lines()
            .collect::<Result<Vec<_>, _>>()
            .map_err(VtkError::Io)?
            .join("\n");

        let mut image = ImageData::new();

        // Parse ImageData attributes from the tag
        if let Some(id_tag) = find_tag(&content, "ImageData") {
            if let Some(extent_str) = extract_attr(&id_tag, "WholeExtent") {
                let vals: Vec<i64> = extent_str
                    .split_whitespace()
                    .filter_map(|s| s.parse().ok())
                    .collect();
                if vals.len() == 6 {
                    image.set_extent([vals[0], vals[1], vals[2], vals[3], vals[4], vals[5]]);
                }
            }
            if let Some(origin_str) = extract_attr(&id_tag, "Origin") {
                let vals: Vec<f64> = origin_str
                    .split_whitespace()
                    .filter_map(|s| s.parse().ok())
                    .collect();
                if vals.len() == 3 {
                    image.set_origin([vals[0], vals[1], vals[2]]);
                }
            }
            if let Some(spacing_str) = extract_attr(&id_tag, "Spacing") {
                let vals: Vec<f64> = spacing_str
                    .split_whitespace()
                    .filter_map(|s| s.parse().ok())
                    .collect();
                if vals.len() == 3 {
                    image.set_spacing([vals[0], vals[1], vals[2]]);
                }
            }
        }

        // Parse PointData
        if let Some(pd_section) = extract_section(&content, "PointData") {
            parse_attribute_arrays(&pd_section, image.point_data_mut())?;
        }

        // Parse CellData
        if let Some(cd_section) = extract_section(&content, "CellData") {
            parse_attribute_arrays(&cd_section, image.cell_data_mut())?;
        }

        Ok(image)
    }
}

fn find_tag(content: &str, tag: &str) -> Option<String> {
    let pattern = format!("<{}", tag);
    let start = content.find(&pattern)?;
    let end = content[start..].find('>')?;
    Some(content[start..start + end + 1].to_string())
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
    use crate::VtiWriter;
    use vtk_data::{DataArray as DA, DataSet};

    #[test]
    fn roundtrip_vti() {
        let mut img = ImageData::with_dimensions(3, 4, 5);
        img.set_spacing([0.5, 0.5, 0.5]);
        img.set_origin([1.0, 2.0, 3.0]);

        let n = img.num_points();
        let scalars: Vec<f64> = (0..n).map(|i| i as f64 * 0.1).collect();
        let arr = DA::from_vec("density", scalars, 1);
        img.point_data_mut().add_array(arr.into());
        img.point_data_mut().set_active_scalars("density");

        let mut buf = Vec::new();
        VtiWriter::write_to(&mut buf, &img).unwrap();

        let reader = std::io::BufReader::new(&buf[..]);
        let result = VtiReader::read_from(reader).unwrap();

        assert_eq!(result.dimensions(), [3, 4, 5]);
        assert_eq!(result.spacing(), [0.5, 0.5, 0.5]);
        assert_eq!(result.origin(), [1.0, 2.0, 3.0]);

        let s = result.point_data().scalars().unwrap();
        assert_eq!(s.num_tuples(), 60);
        let mut val = [0.0f64];
        s.tuple_as_f64(10, &mut val);
        assert!((val[0] - 1.0).abs() < 1e-6);
    }

    #[test]
    fn roundtrip_vti_no_data() {
        let mut img = ImageData::with_dimensions(2, 2, 2);
        img.set_spacing([1.0, 2.0, 3.0]);

        let mut buf = Vec::new();
        VtiWriter::write_to(&mut buf, &img).unwrap();

        let reader = std::io::BufReader::new(&buf[..]);
        let result = VtiReader::read_from(reader).unwrap();

        assert_eq!(result.dimensions(), [2, 2, 2]);
        assert_eq!(result.spacing(), [1.0, 2.0, 3.0]);
    }
}
