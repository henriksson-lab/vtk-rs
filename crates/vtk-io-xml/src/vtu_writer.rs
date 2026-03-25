use std::io::Write;
use std::path::Path;

use vtk_data::{AnyDataArray, DataSetAttributes, UnstructuredGrid};
use vtk_types::VtkError;

/// Writer for VTK XML UnstructuredGrid format (.vtu).
///
/// Produces ASCII XML files compatible with ParaView and other VTK-based tools.
pub struct VtuWriter;

impl VtuWriter {
    pub fn write(path: &Path, grid: &UnstructuredGrid) -> Result<(), VtkError> {
        let file = std::fs::File::create(path)?;
        let mut w = std::io::BufWriter::new(file);
        Self::write_to(&mut w, grid)
    }

    pub fn write_to<W: Write>(w: &mut W, grid: &UnstructuredGrid) -> Result<(), VtkError> {
        writeln!(w, "<?xml version=\"1.0\"?>")?;
        writeln!(
            w,
            "<VTKFile type=\"UnstructuredGrid\" version=\"1.0\" byte_order=\"LittleEndian\">"
        )?;
        writeln!(w, "  <UnstructuredGrid>")?;

        let n_points = grid.points.len();
        let n_cells = grid.cells().num_cells();

        writeln!(
            w,
            "    <Piece NumberOfPoints=\"{}\" NumberOfCells=\"{}\">",
            n_points, n_cells
        )?;

        // Points
        writeln!(w, "      <Points>")?;
        writeln!(
            w,
            "        <DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">"
        )?;
        write!(w, "          ")?;
        for i in 0..n_points {
            let p = grid.points.get(i);
            write!(w, "{} {} {} ", p[0], p[1], p[2])?;
        }
        writeln!(w)?;
        writeln!(w, "        </DataArray>")?;
        writeln!(w, "      </Points>")?;

        // Cells
        writeln!(w, "      <Cells>")?;

        // Connectivity
        writeln!(
            w,
            "        <DataArray type=\"Int64\" Name=\"connectivity\" format=\"ascii\">"
        )?;
        write!(w, "          ")?;
        for i in 0..n_cells {
            for &id in grid.cell_points(i) {
                write!(w, "{} ", id)?;
            }
        }
        writeln!(w)?;
        writeln!(w, "        </DataArray>")?;

        // Offsets
        writeln!(
            w,
            "        <DataArray type=\"Int64\" Name=\"offsets\" format=\"ascii\">"
        )?;
        write!(w, "          ")?;
        let mut offset: i64 = 0;
        for i in 0..n_cells {
            offset += grid.cell_points(i).len() as i64;
            write!(w, "{} ", offset)?;
        }
        writeln!(w)?;
        writeln!(w, "        </DataArray>")?;

        // Types
        writeln!(
            w,
            "        <DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">"
        )?;
        write!(w, "          ")?;
        for i in 0..n_cells {
            write!(w, "{} ", grid.cell_type(i) as u8)?;
        }
        writeln!(w)?;
        writeln!(w, "        </DataArray>")?;

        writeln!(w, "      </Cells>")?;

        // Point data
        if grid.point_data().num_arrays() > 0 {
            write_data_section(w, "PointData", grid.point_data())?;
        }

        // Cell data
        if grid.cell_data().num_arrays() > 0 {
            write_data_section(w, "CellData", grid.cell_data())?;
        }

        writeln!(w, "    </Piece>")?;
        writeln!(w, "  </UnstructuredGrid>")?;
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
        vtk_types::ScalarType::F32 => "Float32",
        vtk_types::ScalarType::F64 => "Float64",
        vtk_types::ScalarType::I8 => "Int8",
        vtk_types::ScalarType::I16 => "Int16",
        vtk_types::ScalarType::I32 => "Int32",
        vtk_types::ScalarType::I64 => "Int64",
        vtk_types::ScalarType::U8 => "UInt8",
        vtk_types::ScalarType::U16 => "UInt16",
        vtk_types::ScalarType::U32 => "UInt32",
        vtk_types::ScalarType::U64 => "UInt64",
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
    use vtk_data::{DataArray, DataSet, UnstructuredGrid};
    use vtk_types::CellType;

    #[test]
    fn write_single_tetra_vtu() {
        let mut grid = UnstructuredGrid::new();
        grid.points.push([0.0, 0.0, 0.0]);
        grid.points.push([1.0, 0.0, 0.0]);
        grid.points.push([0.5, 1.0, 0.0]);
        grid.points.push([0.5, 0.5, 1.0]);
        grid.push_cell(CellType::Tetra, &[0, 1, 2, 3]);

        let mut buf = Vec::new();
        VtuWriter::write_to(&mut buf, &grid).unwrap();
        let output = String::from_utf8(buf).unwrap();

        assert!(output.contains("<VTKFile type=\"UnstructuredGrid\""));
        assert!(output.contains("NumberOfPoints=\"4\""));
        assert!(output.contains("NumberOfCells=\"1\""));
        assert!(output.contains("Name=\"connectivity\""));
        assert!(output.contains("Name=\"offsets\""));
        assert!(output.contains("Name=\"types\""));
    }

    #[test]
    fn write_with_scalars_vtu() {
        let mut grid = UnstructuredGrid::new();
        grid.points.push([0.0, 0.0, 0.0]);
        grid.points.push([1.0, 0.0, 0.0]);
        grid.points.push([0.5, 1.0, 0.0]);
        grid.points.push([0.5, 0.5, 1.0]);
        grid.push_cell(CellType::Tetra, &[0, 1, 2, 3]);

        let scalars = DataArray::from_vec("temp", vec![10.0, 20.0, 30.0, 40.0], 1);
        grid.point_data_mut().add_array(scalars.into());
        grid.point_data_mut().set_active_scalars("temp");

        let mut buf = Vec::new();
        VtuWriter::write_to(&mut buf, &grid).unwrap();
        let output = String::from_utf8(buf).unwrap();

        assert!(output.contains("<PointData Scalars=\"temp\">"));
        assert!(output.contains("Name=\"temp\""));
    }
}
