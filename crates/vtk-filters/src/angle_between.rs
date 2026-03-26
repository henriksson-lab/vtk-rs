use vtk_data::{AnyDataArray, DataArray, PolyData};

/// Compute dihedral angles between adjacent triangles.
///
/// For each edge shared by two triangles, computes the angle between
/// their face normals. Adds a "DihedralAngle" cell data array with
/// the minimum dihedral angle per cell (in degrees).
pub fn dihedral_angles(input: &PolyData) -> PolyData {
    use std::collections::HashMap;

    let cells: Vec<Vec<i64>> = input.polys.iter().map(|c| c.to_vec()).collect();
    let n_cells = cells.len();

    // Compute face normals
    let normals: Vec<[f64; 3]> = cells.iter().map(|cell| {
        if cell.len() < 3 { return [0.0, 0.0, 0.0]; }
        let v0 = input.points.get(cell[0] as usize);
        let v1 = input.points.get(cell[1] as usize);
        let v2 = input.points.get(cell[2] as usize);
        let e1 = [v1[0]-v0[0], v1[1]-v0[1], v1[2]-v0[2]];
        let e2 = [v2[0]-v0[0], v2[1]-v0[1], v2[2]-v0[2]];
        let n = [e1[1]*e2[2]-e1[2]*e2[1], e1[2]*e2[0]-e1[0]*e2[2], e1[0]*e2[1]-e1[1]*e2[0]];
        let len = (n[0]*n[0]+n[1]*n[1]+n[2]*n[2]).sqrt();
        if len > 1e-15 { [n[0]/len, n[1]/len, n[2]/len] } else { [0.0; 3] }
    }).collect();

    // Build edge adjacency
    let mut edge_faces: HashMap<(i64,i64), Vec<usize>> = HashMap::new();
    for (fi, cell) in cells.iter().enumerate() {
        for i in 0..cell.len() {
            let a = cell[i]; let b = cell[(i+1)%cell.len()];
            let key = if a < b { (a,b) } else { (b,a) };
            edge_faces.entry(key).or_default().push(fi);
        }
    }

    // For each cell, find minimum dihedral angle with neighbors
    let mut min_angles = vec![180.0f64; n_cells];
    for faces in edge_faces.values() {
        if faces.len() == 2 {
            let na = normals[faces[0]];
            let nb = normals[faces[1]];
            let dot = (na[0]*nb[0] + na[1]*nb[1] + na[2]*nb[2]).clamp(-1.0, 1.0);
            let angle = dot.acos().to_degrees();
            min_angles[faces[0]] = min_angles[faces[0]].min(angle);
            min_angles[faces[1]] = min_angles[faces[1]].min(angle);
        }
    }

    let mut pd = input.clone();
    pd.cell_data_mut().add_array(AnyDataArray::F64(
        DataArray::from_vec("DihedralAngle", min_angles, 1),
    ));
    pd
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn flat_surface_zero_angle() {
        let mut pd = PolyData::new();
        pd.points.push([0.0, 0.0, 0.0]);
        pd.points.push([1.0, 0.0, 0.0]);
        pd.points.push([1.0, 1.0, 0.0]);
        pd.points.push([0.0, 1.0, 0.0]);
        pd.polys.push_cell(&[0, 1, 2]);
        pd.polys.push_cell(&[0, 2, 3]);

        let result = dihedral_angles(&pd);
        let arr = result.cell_data().get_array("DihedralAngle").unwrap();
        let mut buf = [0.0f64];
        arr.tuple_as_f64(0, &mut buf);
        assert!(buf[0] < 1.0); // nearly 0 degrees for coplanar
    }

    #[test]
    fn right_angle() {
        let mut pd = PolyData::new();
        pd.points.push([0.0, 0.0, 0.0]);
        pd.points.push([1.0, 0.0, 0.0]);
        pd.points.push([0.5, 1.0, 0.0]); // XY plane
        pd.points.push([0.5, 0.0, 1.0]); // XZ plane
        pd.polys.push_cell(&[0, 1, 2]);
        pd.polys.push_cell(&[0, 1, 3]);

        let result = dihedral_angles(&pd);
        let arr = result.cell_data().get_array("DihedralAngle").unwrap();
        let mut buf = [0.0f64];
        arr.tuple_as_f64(0, &mut buf);
        assert!((buf[0] - 90.0).abs() < 5.0);
    }

    #[test]
    fn single_triangle() {
        let mut pd = PolyData::new();
        pd.points.push([0.0, 0.0, 0.0]);
        pd.points.push([1.0, 0.0, 0.0]);
        pd.points.push([0.0, 1.0, 0.0]);
        pd.polys.push_cell(&[0, 1, 2]);

        let result = dihedral_angles(&pd);
        let arr = result.cell_data().get_array("DihedralAngle").unwrap();
        let mut buf = [0.0f64];
        arr.tuple_as_f64(0, &mut buf);
        assert_eq!(buf[0], 180.0); // no neighbor -> default max
    }
}
