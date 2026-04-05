use crate::data::{CellArray, PolyData};

/// Convert all polygons to triangles using fan triangulation.
///
/// Triangles pass through unchanged. Quads and larger polygons are decomposed
/// into triangle fans from the first vertex. Triangle strips are also decomposed.
pub fn triangulate(input: &PolyData) -> PolyData {
    // Fast path: if already all triangles and no strips, just clone
    let has_strips = !input.strips.is_empty();
    let all_tri = input.polys.is_empty() || input.polys.iter().all(|c| c.len() == 3);
    if all_tri && !has_strips {
        return input.clone();
    }

    let mut output = input.clone();

    // Triangulate polygons
    if !input.polys.is_empty() {
        output.polys = triangulate_cells(&input.polys);
    }

    // Decompose triangle strips
    if has_strips {
        let tri_from_strips = decompose_strips(&input.strips);
        // Append strip-derived triangles to polys
        for cell in tri_from_strips.iter() {
            output.polys.push_cell(cell);
        }
        output.strips = CellArray::new();
    }

    output
}

fn triangulate_cells(polys: &CellArray) -> CellArray {
    // Fast path: if all cells are already triangles, return a clone directly.
    // This matches VTK's TriangleFilter behavior which is a no-op on triangle meshes.
    let all_triangles = polys.iter().all(|cell| cell.len() == 3);
    if all_triangles {
        return polys.clone();
    }

    let mut out = CellArray::new();

    for cell in polys.iter() {
        if cell.len() < 3 {
            // Degenerate, skip
            continue;
        }
        if cell.len() == 3 {
            // Already a triangle
            out.push_cell(cell);
        } else {
            // Fan triangulation from vertex 0
            for i in 1..cell.len() - 1 {
                out.push_cell(&[cell[0], cell[i], cell[i + 1]]);
            }
        }
    }

    out
}

fn decompose_strips(strips: &CellArray) -> CellArray {
    let mut out = CellArray::new();

    for strip in strips.iter() {
        if strip.len() < 3 {
            continue;
        }
        for i in 0..strip.len() - 2 {
            if i % 2 == 0 {
                out.push_cell(&[strip[i], strip[i + 1], strip[i + 2]]);
            } else {
                // Flip winding for odd triangles to maintain consistent orientation
                out.push_cell(&[strip[i + 1], strip[i], strip[i + 2]]);
            }
        }
    }

    out
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn triangle_passthrough() {
        let pd = PolyData::from_triangles(
            vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]],
            vec![[0, 1, 2]],
        );
        let result = triangulate(&pd);
        assert_eq!(result.polys.num_cells(), 1);
        assert_eq!(result.polys.cell(0), &[0, 1, 2]);
    }

    #[test]
    fn quad_to_triangles() {
        let mut pd = PolyData::new();
        pd.points = crate::data::Points::from_vec(vec![
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [1.0, 1.0, 0.0],
            [0.0, 1.0, 0.0],
        ]);
        pd.polys.push_cell(&[0, 1, 2, 3]);

        let result = triangulate(&pd);
        assert_eq!(result.polys.num_cells(), 2);
        assert_eq!(result.polys.cell(0), &[0, 1, 2]);
        assert_eq!(result.polys.cell(1), &[0, 2, 3]);
    }

    #[test]
    fn strip_decomposition() {
        let mut pd = PolyData::new();
        pd.points = crate::data::Points::from_vec(vec![
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.5, 1.0, 0.0],
            [1.5, 1.0, 0.0],
        ]);
        pd.strips.push_cell(&[0, 1, 2, 3]);

        let result = triangulate(&pd);
        assert!(result.strips.is_empty());
        assert_eq!(result.polys.num_cells(), 2);
    }
}
