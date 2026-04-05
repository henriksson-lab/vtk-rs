use crate::data::{CellArray, Points, PolyData};

/// Generate points at the centroid of each cell. Pre-sized buffers.
pub fn cell_centers(input: &PolyData) -> PolyData {
    let nc = input.polys.num_cells() + input.lines.num_cells();
    let mut pts = Vec::with_capacity(nc * 3);
    let mut vert_off = Vec::with_capacity(nc + 1);
    let mut vert_conn = Vec::with_capacity(nc);
    vert_off.push(0i64);

    let mut idx = 0i64;
    for cells in [&input.polys, &input.lines] {
        for ci in 0..cells.num_cells() {
            let cell = cells.cell(ci);
            if cell.is_empty() { continue; }
            let n = cell.len() as f64;
            let (mut cx, mut cy, mut cz) = (0.0, 0.0, 0.0);
            for &id in cell {
                let p = input.points.get(id as usize);
                cx += p[0]; cy += p[1]; cz += p[2];
            }
            pts.push(cx/n); pts.push(cy/n); pts.push(cz/n);
            vert_conn.push(idx);
            vert_off.push(idx + 1);
            idx += 1;
        }
    }

    let mut pd = PolyData::new();
    pd.points = Points::from_flat_vec(pts);
    pd.verts = CellArray::from_raw(vert_off, vert_conn);
    pd
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn centers_of_triangles() {
        let pd = PolyData::from_triangles(
            vec![[0.0,0.0,0.0],[3.0,0.0,0.0],[0.0,3.0,0.0],[3.0,3.0,0.0]],
            vec![[0,1,2],[1,3,2]]);
        let r = cell_centers(&pd);
        assert_eq!(r.points.len(), 2);
        let c = r.points.get(0);
        assert!((c[0] - 1.0).abs() < 1e-10);
    }
}
