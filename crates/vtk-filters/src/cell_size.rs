use vtk_data::{AnyDataArray, DataArray, PolyData};

/// Compute the area of each polygon cell and add it as cell data.
///
/// Adds a "CellSize" scalar array to the cell data. For triangles,
/// this is the exact area; for general polygons, it uses fan triangulation.
pub fn cell_size(input: &PolyData) -> PolyData {
    let mut areas = Vec::new();

    for cell in input.polys.iter() {
        if cell.len() < 3 {
            areas.push(0.0);
            continue;
        }

        let p0 = input.points.get(cell[0] as usize);
        let mut total_area = 0.0;

        for i in 1..cell.len() - 1 {
            let p1 = input.points.get(cell[i] as usize);
            let p2 = input.points.get(cell[i + 1] as usize);
            let e1 = [p1[0] - p0[0], p1[1] - p0[1], p1[2] - p0[2]];
            let e2 = [p2[0] - p0[0], p2[1] - p0[1], p2[2] - p0[2]];
            let cross = [
                e1[1] * e2[2] - e1[2] * e2[1],
                e1[2] * e2[0] - e1[0] * e2[2],
                e1[0] * e2[1] - e1[1] * e2[0],
            ];
            total_area += 0.5 * (cross[0] * cross[0] + cross[1] * cross[1] + cross[2] * cross[2]).sqrt();
        }

        areas.push(total_area);
    }

    let mut pd = input.clone();
    pd.cell_data_mut().add_array(AnyDataArray::F64(
        DataArray::from_vec("CellSize", areas, 1),
    ));
    pd
}

/// Compute the edge length of each line cell.
pub fn cell_size_lines(input: &PolyData) -> PolyData {
    let mut lengths = Vec::new();

    for cell in input.lines.iter() {
        let mut total_len = 0.0;
        for i in 0..cell.len().saturating_sub(1) {
            let p0 = input.points.get(cell[i] as usize);
            let p1 = input.points.get(cell[i + 1] as usize);
            let d = [p1[0] - p0[0], p1[1] - p0[1], p1[2] - p0[2]];
            total_len += (d[0] * d[0] + d[1] * d[1] + d[2] * d[2]).sqrt();
        }
        lengths.push(total_len);
    }

    let mut pd = input.clone();
    pd.cell_data_mut().add_array(AnyDataArray::F64(
        DataArray::from_vec("CellSize", lengths, 1),
    ));
    pd
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn triangle_area() {
        let pd = PolyData::from_triangles(
            vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]],
            vec![[0, 1, 2]],
        );
        let result = cell_size(&pd);
        let arr = result.cell_data().get_array("CellSize").unwrap();
        let mut val = [0.0f64];
        arr.tuple_as_f64(0, &mut val);
        assert!((val[0] - 0.5).abs() < 1e-10);
    }

    #[test]
    fn quad_area() {
        let mut pd = PolyData::new();
        pd.points.push([0.0, 0.0, 0.0]);
        pd.points.push([2.0, 0.0, 0.0]);
        pd.points.push([2.0, 3.0, 0.0]);
        pd.points.push([0.0, 3.0, 0.0]);
        pd.polys.push_cell(&[0, 1, 2, 3]);
        let result = cell_size(&pd);
        let arr = result.cell_data().get_array("CellSize").unwrap();
        let mut val = [0.0f64];
        arr.tuple_as_f64(0, &mut val);
        assert!((val[0] - 6.0).abs() < 1e-10);
    }

    #[test]
    fn line_length() {
        let mut pd = PolyData::new();
        pd.points.push([0.0, 0.0, 0.0]);
        pd.points.push([3.0, 4.0, 0.0]);
        pd.lines.push_cell(&[0, 1]);
        let result = cell_size_lines(&pd);
        let arr = result.cell_data().get_array("CellSize").unwrap();
        let mut val = [0.0f64];
        arr.tuple_as_f64(0, &mut val);
        assert!((val[0] - 5.0).abs() < 1e-10);
    }
}
