use crate::data::{CellArray, Points, PolyData};

/// Extract connected components from a PolyData mesh.
///
/// Returns a Vec of PolyData, one per connected component, sorted by size (largest first).
pub fn extract_components(input: &PolyData) -> Vec<PolyData> {
    let n = input.points.len();
    if n == 0 {
        return Vec::new();
    }

    // Union-find
    let mut parent: Vec<usize> = (0..n).collect();
    let mut rank = vec![0usize; n];

    // Union all points in each cell
    for cell in input.polys.iter() {
        if cell.len() < 2 {
            continue;
        }
        let first = cell[0] as usize;
        for &id in &cell[1..] {
            union(&mut parent, &mut rank, first, id as usize);
        }
    }
    // Also for other cell types
    for cells in [&input.verts, &input.lines, &input.strips] {
        for cell in cells.iter() {
            if cell.len() < 2 {
                continue;
            }
            let first = cell[0] as usize;
            for &id in &cell[1..] {
                union(&mut parent, &mut rank, first, id as usize);
            }
        }
    }

    // Group cells by root
    let mut component_cells: std::collections::HashMap<usize, Vec<Vec<i64>>> =
        std::collections::HashMap::new();

    for cell in input.polys.iter() {
        if cell.is_empty() {
            continue;
        }
        let root = find(&mut parent, cell[0] as usize);
        component_cells
            .entry(root)
            .or_default()
            .push(cell.to_vec());
    }

    // Sort by size (largest first)
    let mut components: Vec<(usize, Vec<Vec<i64>>)> = component_cells.into_iter().collect();
    components.sort_by(|a, b| b.1.len().cmp(&a.1.len()));

    // Build output PolyData for each component
    components
        .into_iter()
        .map(|(_, cells)| {
            build_component(input, &cells)
        })
        .collect()
}

/// Extract only the largest connected component.
pub fn extract_largest_component(input: &PolyData) -> PolyData {
    let mut components = extract_components(input);
    if components.is_empty() {
        PolyData::new()
    } else {
        components.remove(0)
    }
}

fn build_component(input: &PolyData, cells: &[Vec<i64>]) -> PolyData {
    let mut point_map = std::collections::HashMap::new();
    let mut new_points = Points::<f64>::new();
    let mut polys = CellArray::new();

    for cell in cells {
        let remapped: Vec<i64> = cell
            .iter()
            .map(|&id| {
                let old = id as usize;
                *point_map.entry(old).or_insert_with(|| {
                    let idx = new_points.len();
                    new_points.push(input.points.get(old));
                    idx
                }) as i64
            })
            .collect();
        polys.push_cell(&remapped);
    }

    let mut output = PolyData::new();
    output.points = new_points;
    output.polys = polys;
    output
}

fn find(parent: &mut [usize], x: usize) -> usize {
    if parent[x] != x {
        parent[x] = find(parent, parent[x]);
    }
    parent[x]
}

fn union(parent: &mut [usize], rank: &mut [usize], a: usize, b: usize) {
    let ra = find(parent, a);
    let rb = find(parent, b);
    if ra == rb {
        return;
    }
    if rank[ra] < rank[rb] {
        parent[ra] = rb;
    } else if rank[ra] > rank[rb] {
        parent[rb] = ra;
    } else {
        parent[rb] = ra;
        rank[ra] += 1;
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn two_disconnected_triangles() {
        let pd = PolyData::from_triangles(
            vec![
                [0.0, 0.0, 0.0],
                [1.0, 0.0, 0.0],
                [0.0, 1.0, 0.0],
                // Separate triangle
                [10.0, 0.0, 0.0],
                [11.0, 0.0, 0.0],
                [10.0, 1.0, 0.0],
            ],
            vec![[0, 1, 2], [3, 4, 5]],
        );

        let components = extract_components(&pd);
        assert_eq!(components.len(), 2);
        assert_eq!(components[0].polys.num_cells(), 1);
        assert_eq!(components[1].polys.num_cells(), 1);
    }

    #[test]
    fn single_connected_mesh() {
        let pd = PolyData::from_triangles(
            vec![
                [0.0, 0.0, 0.0],
                [1.0, 0.0, 0.0],
                [0.5, 1.0, 0.0],
                [1.5, 1.0, 0.0],
            ],
            vec![[0, 1, 2], [1, 3, 2]],
        );

        let components = extract_components(&pd);
        assert_eq!(components.len(), 1);
        assert_eq!(components[0].polys.num_cells(), 2);
    }

    #[test]
    fn largest_component() {
        let pd = PolyData::from_triangles(
            vec![
                [0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.5, 1.0, 0.0], [1.5, 1.0, 0.0],
                // Separate single triangle
                [10.0, 0.0, 0.0], [11.0, 0.0, 0.0], [10.0, 1.0, 0.0],
            ],
            vec![[0, 1, 2], [1, 3, 2], [4, 5, 6]],
        );

        let largest = extract_largest_component(&pd);
        assert_eq!(largest.polys.num_cells(), 2);
    }
}
