use vtk_data::{CellArray, Points, PolyData};

/// Separate all cells so they don't share any vertices.
///
/// Each cell gets its own copy of its vertices. Useful for flat shading
/// or for exploded views. Opposite of `clean` which merges shared vertices.
pub fn separate_cells(input: &PolyData) -> PolyData {
    let mut out_points = Points::<f64>::new();
    let mut out_polys = CellArray::new();
    let mut out_lines = CellArray::new();

    for cell in input.polys.iter() {
        let base = out_points.len() as i64;
        let mut ids = Vec::with_capacity(cell.len());
        for &pid in cell.iter() {
            out_points.push(input.points.get(pid as usize));
            ids.push(base + ids.len() as i64);
        }
        out_polys.push_cell(&ids);
    }

    for cell in input.lines.iter() {
        let base = out_points.len() as i64;
        let mut ids = Vec::with_capacity(cell.len());
        for &pid in cell.iter() {
            out_points.push(input.points.get(pid as usize));
            ids.push(base + ids.len() as i64);
        }
        out_lines.push_cell(&ids);
    }

    let mut pd = PolyData::new();
    pd.points = out_points;
    pd.polys = out_polys;
    pd.lines = out_lines;
    pd
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn shared_vertices_duplicated() {
        let mut pd = PolyData::new();
        pd.points.push([0.0, 0.0, 0.0]); // shared by both tris
        pd.points.push([1.0, 0.0, 0.0]); // shared
        pd.points.push([0.5, 1.0, 0.0]);
        pd.points.push([1.5, 1.0, 0.0]);
        pd.polys.push_cell(&[0, 1, 2]);
        pd.polys.push_cell(&[0, 1, 3]);

        let result = separate_cells(&pd);
        assert_eq!(result.points.len(), 6); // 3 + 3, no sharing
        assert_eq!(result.polys.num_cells(), 2);
    }

    #[test]
    fn single_cell_unchanged_count() {
        let mut pd = PolyData::new();
        pd.points.push([0.0, 0.0, 0.0]);
        pd.points.push([1.0, 0.0, 0.0]);
        pd.points.push([0.5, 1.0, 0.0]);
        pd.polys.push_cell(&[0, 1, 2]);

        let result = separate_cells(&pd);
        assert_eq!(result.points.len(), 3);
        assert_eq!(result.polys.num_cells(), 1);
    }

    #[test]
    fn empty_input() {
        let pd = PolyData::new();
        let result = separate_cells(&pd);
        assert_eq!(result.polys.num_cells(), 0);
    }
}
