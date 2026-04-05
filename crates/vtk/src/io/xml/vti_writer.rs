use std::io::Write;
use std::path::Path;

use crate::data::{AnyDataArray, DataSetAttributes, ImageData};
use crate::types::VtkError;

/// Writer for VTK XML ImageData format (.vti).
pub struct VtiWriter;

impl VtiWriter {
    pub fn write(path: &Path, data: &ImageData) -> Result<(), VtkError> {
        let file = std::fs::File::create(path)?;
        let mut w = std::io::BufWriter::new(file);
        Self::write_to(&mut w, data)
    }

    pub fn write_to<W: Write>(w: &mut W, data: &ImageData) -> Result<(), VtkError> {
        let ext = data.extent();
        let spacing = data.spacing();
        let origin = data.origin();

        writeln!(w, "<?xml version=\"1.0\"?>")?;
        writeln!(
            w,
            "<VTKFile type=\"ImageData\" version=\"1.0\" byte_order=\"LittleEndian\">"
        )?;
        writeln!(
            w,
            "  <ImageData WholeExtent=\"{} {} {} {} {} {}\" Origin=\"{} {} {}\" Spacing=\"{} {} {}\">",
            ext[0], ext[1], ext[2], ext[3], ext[4], ext[5],
            origin[0], origin[1], origin[2],
            spacing[0], spacing[1], spacing[2],
        )?;
        writeln!(
            w,
            "    <Piece Extent=\"{} {} {} {} {} {}\">",
            ext[0], ext[1], ext[2], ext[3], ext[4], ext[5],
        )?;

        // Point data
        if data.point_data().num_arrays() > 0 {
            write_data_section(w, "PointData", data.point_data())?;
        }

        // Cell data
        if data.cell_data().num_arrays() > 0 {
            write_data_section(w, "CellData", data.cell_data())?;
        }

        writeln!(w, "    </Piece>")?;
        writeln!(w, "  </ImageData>")?;
        writeln!(w, "</VTKFile>")?;

        Ok(())
    }
}

fn write_data_section<W: Write>(
    w: &mut W,
    section: &str,
    attrs: &DataSetAttributes,
) -> Result<(), VtkError> {
    let scalars_name = attrs.scalars().map(|a| a.name().to_string());
    let mut attrs_str = String::new();
    if let Some(ref name) = scalars_name {
        attrs_str.push_str(&format!(" Scalars=\"{}\"", name));
    }

    writeln!(w, "      <{}{}>", section, attrs_str)?;
    for i in 0..attrs.num_arrays() {
        if let Some(arr) = attrs.get_array_by_index(i) {
            write_any_data_array(w, arr)?;
        }
    }
    writeln!(w, "      </{}>", section)?;
    Ok(())
}

fn write_any_data_array<W: Write>(w: &mut W, arr: &AnyDataArray) -> Result<(), VtkError> {
    let type_name = match arr.scalar_type() {
        crate::types::ScalarType::F32 => "Float32",
        crate::types::ScalarType::F64 => "Float64",
        crate::types::ScalarType::I8 => "Int8",
        crate::types::ScalarType::I16 => "Int16",
        crate::types::ScalarType::I32 => "Int32",
        crate::types::ScalarType::I64 => "Int64",
        crate::types::ScalarType::U8 => "UInt8",
        crate::types::ScalarType::U16 => "UInt16",
        crate::types::ScalarType::U32 => "UInt32",
        crate::types::ScalarType::U64 => "UInt64",
    };

    writeln!(
        w,
        "        <DataArray type=\"{}\" Name=\"{}\" NumberOfComponents=\"{}\" format=\"ascii\">",
        type_name,
        arr.name(),
        arr.num_components()
    )?;
    write!(w, "          ")?;
    let nt = arr.num_tuples();
    let nc = arr.num_components();
    let mut buf = vec![0.0f64; nc];
    for i in 0..nt {
        arr.tuple_as_f64(i, &mut buf);
        for v in &buf {
            write!(w, "{} ", v)?;
        }
    }
    writeln!(w)?;
    writeln!(w, "        </DataArray>")?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::data::{DataArray, DataSet, ImageData};

    #[test]
    fn write_simple_vti() {
        let mut img = ImageData::with_dimensions(3, 4, 5);
        img.set_spacing([0.5, 0.5, 0.5]);
        img.set_origin([1.0, 2.0, 3.0]);

        let n = img.num_points();
        let scalars: Vec<f64> = (0..n).map(|i| i as f64 * 0.1).collect();
        let arr = DataArray::from_vec("density", scalars, 1);
        img.point_data_mut().add_array(arr.into());
        img.point_data_mut().set_active_scalars("density");

        let mut buf = Vec::new();
        VtiWriter::write_to(&mut buf, &img).unwrap();
        let output = String::from_utf8(buf).unwrap();

        assert!(output.contains("<VTKFile type=\"ImageData\""));
        assert!(output.contains("WholeExtent=\"0 2 0 3 0 4\""));
        assert!(output.contains("Origin=\"1 2 3\""));
        assert!(output.contains("Spacing=\"0.5 0.5 0.5\""));
        assert!(output.contains("Name=\"density\""));
    }
}
