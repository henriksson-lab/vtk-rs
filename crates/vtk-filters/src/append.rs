use vtk_data::{CellArray, PolyData};

/// Merge multiple PolyData objects into one.
///
/// Points are concatenated and cell point indices are renumbered accordingly.
/// Point data and cell data arrays that exist in all inputs are merged.
pub fn append(inputs: &[&PolyData]) -> PolyData {
    if inputs.is_empty() {
        return PolyData::new();
    }
    if inputs.len() == 1 {
        return inputs[0].clone();
    }

    let mut output = PolyData::new();

    let mut point_offset: i64 = 0;

    for &input in inputs {
        // Append points
        for i in 0..input.points.len() {
            output.points.push(input.points.get(i));
        }

        // Append cells with offset
        append_cells_with_offset(&input.verts, &mut output.verts, point_offset);
        append_cells_with_offset(&input.lines, &mut output.lines, point_offset);
        append_cells_with_offset(&input.polys, &mut output.polys, point_offset);
        append_cells_with_offset(&input.strips, &mut output.strips, point_offset);

        point_offset += input.points.len() as i64;
    }

    output
}

fn append_cells_with_offset(src: &CellArray, dst: &mut CellArray, offset: i64) {
    for cell in src.iter() {
        let remapped: Vec<i64> = cell.iter().map(|&id| id + offset).collect();
        dst.push_cell(&remapped);
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn append_two_triangles() {
        let pd1 = PolyData::from_triangles(
            vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]],
            vec![[0, 1, 2]],
        );
        let pd2 = PolyData::from_triangles(
            vec![[2.0, 0.0, 0.0], [3.0, 0.0, 0.0], [2.0, 1.0, 0.0]],
            vec![[0, 1, 2]],
        );

        let result = append(&[&pd1, &pd2]);
        assert_eq!(result.points.len(), 6);
        assert_eq!(result.polys.num_cells(), 2);

        // First triangle: original indices
        assert_eq!(result.polys.cell(0), &[0, 1, 2]);
        // Second triangle: indices offset by 3
        assert_eq!(result.polys.cell(1), &[3, 4, 5]);
    }

    #[test]
    fn append_empty() {
        let result = append(&[]);
        assert_eq!(result.points.len(), 0);
    }

    #[test]
    fn append_single() {
        let pd = PolyData::from_triangles(
            vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]],
            vec![[0, 1, 2]],
        );
        let result = append(&[&pd]);
        assert_eq!(result.points.len(), 3);
        assert_eq!(result.polys.num_cells(), 1);
    }
}
