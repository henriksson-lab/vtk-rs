use std::io::BufRead;
use std::path::Path;

use vtk_data::{AnyDataArray, DataArray, RectilinearGrid};
use vtk_types::VtkError;

/// Reader for VTK XML RectilinearGrid format (.vtr).
pub struct VtrReader;

impl VtrReader {
    pub fn read(path: &Path) -> Result<RectilinearGrid, VtkError> {
        let file = std::fs::File::open(path)?;
        let reader = std::io::BufReader::new(file);
        Self::read_from(reader)
    }

    pub fn read_from<R: BufRead>(reader: R) -> Result<RectilinearGrid, VtkError> {
        let content: String = reader
            .lines()
            .collect::<Result<Vec<_>, _>>()
            .map_err(VtkError::Io)?
            .join("\n");

        let mut x_coords = vec![0.0];
        let mut y_coords = vec![0.0];
        let mut z_coords = vec![0.0];

        // Parse Coordinates section
        if let Some(coords_section) = extract_section(&content, "Coordinates") {
            let mut search_pos = 0;
            while let Some(da_start) = coords_section[search_pos..].find("<DataArray") {
                let abs_start = search_pos + da_start;
                let tag_end = coords_section[abs_start..].find('>').ok_or_else(|| VtkError::Parse("unclosed tag".into()))?;
                let tag = &coords_section[abs_start..abs_start + tag_end + 1];
                let content_start = abs_start + tag_end + 1;
                let content_end = coords_section[content_start..].find("</DataArray>").ok_or_else(|| VtkError::Parse("missing close".into()))?;
                let da_content = coords_section[content_start..content_start + content_end].trim();

                let name = extract_attr(tag, "Name").unwrap_or_default();
                let values: Vec<f64> = da_content.split_whitespace().filter_map(|s| s.parse().ok()).collect();

                match name.as_str() {
                    "x" => x_coords = values,
                    "y" => y_coords = values,
                    "z" => z_coords = values,
                    _ => {}
                }

                search_pos = content_start + content_end + "</DataArray>".len();
            }
        }

        if x_coords.is_empty() { x_coords = vec![0.0]; }
        if y_coords.is_empty() { y_coords = vec![0.0]; }
        if z_coords.is_empty() { z_coords = vec![0.0]; }

        let mut grid = RectilinearGrid::from_coords(x_coords, y_coords, z_coords);

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

fn extract_section(content: &str, tag: &str) -> Option<String> {
    let open = format!("<{}", tag);
    let close = format!("</{}>", tag);
    let start = content.find(&open)?;
    let tag_end = content[start..].find('>')?;
    let cs = start + tag_end + 1;
    let end = content[cs..].find(&close)?;
    Some(content[cs..cs + end].to_string())
}

fn extract_attr(tag: &str, name: &str) -> Option<String> {
    let pat = format!("{}=\"", name);
    let start = tag.find(&pat)?;
    let vs = start + pat.len();
    let end = tag[vs..].find('"')?;
    Some(tag[vs..vs + end].to_string())
}

fn parse_attribute_arrays(section: &str, attrs: &mut vtk_data::DataSetAttributes) -> Result<(), VtkError> {
    let mut search_pos = 0;
    while let Some(da_start) = section[search_pos..].find("<DataArray") {
        let abs_start = search_pos + da_start;
        let tag_end = section[abs_start..].find('>').ok_or_else(|| VtkError::Parse("unclosed tag".into()))?;
        let tag = &section[abs_start..abs_start + tag_end + 1];
        let cs = abs_start + tag_end + 1;
        let ce = section[cs..].find("</DataArray>").ok_or_else(|| VtkError::Parse("missing close".into()))?;
        let content = section[cs..cs + ce].trim();

        let name = extract_attr(tag, "Name").unwrap_or_else(|| "data".to_string());
        let nc: usize = extract_attr(tag, "NumberOfComponents").and_then(|s| s.parse().ok()).unwrap_or(1);
        let values: Vec<f64> = content.split_whitespace().filter_map(|s| s.parse().ok()).collect();

        let arr = AnyDataArray::F64(DataArray::from_vec(&name, values, nc));
        let arr_name = arr.name().to_string();
        attrs.add_array(arr);
        if attrs.scalars().is_none() {
            attrs.set_active_scalars(&arr_name);
        }

        search_pos = cs + ce + "</DataArray>".len();
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::VtrWriter;
    use vtk_data::{DataArray as DA, DataSet};

    #[test]
    fn roundtrip_vtr() {
        let mut grid = RectilinearGrid::from_coords(
            vec![0.0, 1.0, 3.0],
            vec![0.0, 2.0],
            vec![0.0, 5.0],
        );
        let n = grid.num_points();
        let scalars: Vec<f64> = (0..n).map(|i| i as f64).collect();
        let arr = DA::from_vec("idx", scalars, 1);
        grid.point_data_mut().add_array(arr.into());
        grid.point_data_mut().set_active_scalars("idx");

        let mut buf = Vec::new();
        VtrWriter::write_to(&mut buf, &grid).unwrap();

        let reader = std::io::BufReader::new(&buf[..]);
        let result = VtrReader::read_from(reader).unwrap();

        assert_eq!(result.dimensions(), [3, 2, 2]);
        assert_eq!(result.x_coords(), &[0.0, 1.0, 3.0]);
        assert_eq!(result.y_coords(), &[0.0, 2.0]);
        assert_eq!(result.z_coords(), &[0.0, 5.0]);

        let s = result.point_data().scalars().unwrap();
        assert_eq!(s.num_tuples(), 12);
    }
}
