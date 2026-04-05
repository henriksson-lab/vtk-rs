use std::path::Path;

use vtk_data::PolyData;
use vtk_types::VtkError;

/// Read a PolyData from a file, auto-detecting the format by extension.
///
/// Supported extensions: `.vtk`, `.vtp`, `.stl`, `.obj`, `.ply`, `.glb`, `.off`
pub fn read_poly_data(path: &Path) -> Result<PolyData, VtkError> {
    let ext = path.extension()
        .and_then(|e| e.to_str())
        .map(|e| e.to_lowercase())
        .unwrap_or_default();

    match ext.as_str() {
        "vtk" => vtk_io_legacy::LegacyReader::read_poly_data(path),
        "vtp" => vtk_io_xml::VtpReader::read(path),
        "stl" => vtk_io_stl::StlReader::read(path),
        "obj" => vtk_io_obj::ObjReader::read(path),
        "ply" => vtk_io_ply::PlyReader::read(path),
        "glb" => vtk_io_gltf::GlbReader::read(path),
        "off" => vtk_io_off::read_off_file(path).map_err(|e| VtkError::Parse(e)),
        _ => Err(VtkError::Unsupported(format!("unknown file extension: .{ext}"))),
    }
}

/// Write a PolyData to a file, auto-detecting the format by extension.
///
/// Supported extensions: `.vtk`, `.vtp`, `.stl`, `.obj`, `.ply`, `.glb`, `.off`
pub fn write_poly_data(path: &Path, data: &PolyData) -> Result<(), VtkError> {
    let ext = path.extension()
        .and_then(|e| e.to_str())
        .map(|e| e.to_lowercase())
        .unwrap_or_default();

    match ext.as_str() {
        "vtk" => vtk_io_legacy::LegacyWriter::ascii().write_poly_data(path, data),
        "vtp" => vtk_io_xml::VtpWriter::write(path, data),
        "stl" => vtk_io_stl::StlWriter::binary().write(path, data),
        "obj" => vtk_io_obj::ObjWriter::write(path, data),
        "ply" => vtk_io_ply::PlyWriter::write(path, data),
        "glb" => vtk_io_gltf::GlbWriter::write(path, data),
        "off" => vtk_io_off::write_off_file(data, path).map_err(|e| VtkError::Io(std::io::Error::new(std::io::ErrorKind::Other, e))),
        _ => Err(VtkError::Unsupported(format!("unknown file extension: .{ext}"))),
    }
}

/// List of supported file extensions for reading.
pub fn supported_read_extensions() -> &'static [&'static str] {
    &["vtk", "vtp", "stl", "obj", "ply", "glb", "off"]
}

/// List of supported file extensions for writing.
pub fn supported_write_extensions() -> &'static [&'static str] {
    &["vtk", "vtp", "stl", "obj", "ply", "glb", "off"]
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn roundtrip_auto_vtk() {
        let pd = PolyData::from_triangles(
            vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]],
            vec![[0, 1, 2]],
        );
        let dir = std::env::temp_dir().join("vtk_auto_io_test2");
        let _ = std::fs::remove_dir_all(&dir);
        std::fs::create_dir_all(&dir).unwrap();

        let path = dir.join("test.vtk");
        write_poly_data(&path, &pd).unwrap();
        let result = read_poly_data(&path).unwrap();
        assert_eq!(result.points.len(), 3);

        let _ = std::fs::remove_dir_all(&dir);
    }

    #[test]
    fn unsupported_extension() {
        let result = read_poly_data(Path::new("test.xyz"));
        assert!(result.is_err());
    }
}
