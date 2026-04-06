use rayon::prelude::*;
use crate::data::{CellArray, DataArray, ImageData, Points, PolyData};
use crate::filters::core::marching_cubes::{EDGE_TABLE, TRI_TABLE};

/// Flying Edges 3D — scanline-optimized isosurface with row trimming.
///
/// Key optimizations: x-edge classification with trim ranges, edge caching,
/// edge-case-to-MC-case conversion via inline bit mapping.
pub fn flying_edges_3d(image: &ImageData, scalars: &[f64], isovalue: f64) -> PolyData {
    let dims = image.dimensions();
    if dims[0] < 2 || dims[1] < 2 || dims[2] < 2 { return PolyData::new(); }

    let (nx, ny, nz) = (dims[0], dims[1], dims[2]);
    let sp = image.spacing();
    let org = image.origin();
    let nxy = nx * ny;
    let nxm1 = nx - 1;

    // Pass 1: classify x-edges + trim ranges
    let n_rows = ny * nz;
    let mut x_cases: Vec<u8> = vec![0; nxm1 * n_rows];
    let mut trim_min: Vec<usize> = vec![nxm1; n_rows];
    let mut trim_max: Vec<usize> = vec![0; n_rows];

    for k in 0..nz {
        for j in 0..ny {
            let row = k * ny + j;
            let rs = k * nxy + j * nx;
            let cs = row * nxm1;
            let mut s1 = scalars[rs];
            for i in 0..nxm1 {
                let s0 = s1;
                s1 = unsafe { *scalars.get_unchecked(rs + i + 1) };
                let mut ec = 0u8;
                if s0 >= isovalue { ec |= 1; }
                if s1 >= isovalue { ec |= 2; }
                unsafe { *x_cases.get_unchecked_mut(cs + i) = ec; }
                if ec == 1 || ec == 2 {
                    if i < trim_min[row] { trim_min[row] = i; }
                    trim_max[row] = i + 1;
                }
            }
        }
    }

    // Corner offsets (MC ordering)
    const C: [[usize;3];8] = [[0,0,0],[1,0,0],[1,1,0],[0,1,0],[0,0,1],[1,0,1],[1,1,1],[0,1,1]];
    const ELUT: [(usize,usize,u8,usize,usize); 12] = [
        (0,1,0,0,0),(1,2,2,1,0),(2,3,0,0,1),(3,0,2,0,0),
        (4,5,1,0,0),(5,6,3,1,0),(6,7,1,0,1),(7,4,3,0,0),
        (0,4,4,0,0),(1,5,4,1,0),(2,6,4,1,1),(3,7,4,0,1),
    ];

    // Pass 2: count points and triangles per cell for pre-allocation
    let mut est_pts: usize = 0;
    let mut est_tris: usize = 0;
    for k in 0..nz-1 {
        for j in 0..ny-1 {
            let rows = [k*ny+j, k*ny+j+1, (k+1)*ny+j, (k+1)*ny+j+1];
            let mut xl = nxm1; let mut xr = 0;
            for &r in &rows { xl = xl.min(trim_min[r]); xr = xr.max(trim_max[r]); }
            if xl >= xr { continue; }
            let xc = [rows[0]*nxm1, rows[1]*nxm1, rows[2]*nxm1, rows[3]*nxm1];
            for i in xl..xr {
                let e00 = unsafe { *x_cases.get_unchecked(xc[0]+i) };
                let e10 = unsafe { *x_cases.get_unchecked(xc[1]+i) };
                let e01 = unsafe { *x_cases.get_unchecked(xc[2]+i) };
                let e11 = unsafe { *x_cases.get_unchecked(xc[3]+i) };
                let ci = (e00 & 1) | (e00 & 2)
                    | ((e10 & 2) << 1) | ((e10 & 1) << 3)
                    | ((e01 & 1) << 4) | ((e01 & 2) << 4)
                    | ((e11 & 2) << 5) | ((e11 & 1) << 7);
                let edges = EDGE_TABLE[ci as usize];
                if edges == 0 { continue; }
                est_pts += (edges as u32).count_ones() as usize; // upper bound
                let tr = &TRI_TABLE[ci as usize];
                let mut t = 0;
                while t < 15 && tr[t] != -1 { t += 3; }
                est_tris += t / 3;
            }
        }
    }

    let mut pts_flat: Vec<f64> = Vec::with_capacity(est_pts * 3);
    let mut conn: Vec<i64> = Vec::with_capacity(est_tris * 3);
    let mut n_pts: i64 = 0;

    let ss = nx * ny;
    let mut xe0: Vec<i64> = vec![-1; ss];
    let mut xe1: Vec<i64> = vec![-1; ss];
    let mut ye0: Vec<i64> = vec![-1; ss];
    let mut ye1: Vec<i64> = vec![-1; ss];
    let mut ze:  Vec<i64> = vec![-1; ss];

    for k in 0..nz-1 {
        for j in 0..ny-1 {
            let rows = [k*ny+j, k*ny+j+1, (k+1)*ny+j, (k+1)*ny+j+1];
            let mut xl = nxm1; let mut xr = 0;
            for &r in &rows { xl = xl.min(trim_min[r]); xr = xr.max(trim_max[r]); }
            if xl >= xr { continue; }

            let xc = [rows[0]*nxm1, rows[1]*nxm1, rows[2]*nxm1, rows[3]*nxm1];

            for i in xl..xr {
                // Build MC case index from 4 x-edge cases
                let e00 = unsafe { *x_cases.get_unchecked(xc[0]+i) };
                let e10 = unsafe { *x_cases.get_unchecked(xc[1]+i) };
                let e01 = unsafe { *x_cases.get_unchecked(xc[2]+i) };
                let e11 = unsafe { *x_cases.get_unchecked(xc[3]+i) };

                let ci = (e00 & 1) | (e00 & 2)
                    | ((e10 & 2) << 1) | ((e10 & 1) << 3)
                    | ((e01 & 1) << 4) | ((e01 & 2) << 4)
                    | ((e11 & 2) << 5) | ((e11 & 1) << 7);

                let edges = EDGE_TABLE[ci as usize];
                if edges == 0 { continue; }

                let base = k * nxy + j * nx + i;
                let v = unsafe {[
                    *scalars.get_unchecked(base),
                    *scalars.get_unchecked(base + 1),
                    *scalars.get_unchecked(base + nx + 1),
                    *scalars.get_unchecked(base + nx),
                    *scalars.get_unchecked(base + nxy),
                    *scalars.get_unchecked(base + nxy + 1),
                    *scalars.get_unchecked(base + nxy + nx + 1),
                    *scalars.get_unchecked(base + nxy + nx),
                ]};

                let mut ep = [0i64; 12];
                for e in 0..12u32 {
                    if edges & (1 << e) == 0 { continue; }
                    let (c0, c1, cid, di, dj) = ELUT[e as usize];
                    let sf = (j+dj)*nx + (i+di);
                    let cv = unsafe { match cid {
                        0 => xe0.get_unchecked_mut(sf),
                        1 => xe1.get_unchecked_mut(sf),
                        2 => ye0.get_unchecked_mut(sf),
                        3 => ye1.get_unchecked_mut(sf),
                        _ => ze.get_unchecked_mut(sf),
                    }};
                    if *cv >= 0 {
                        ep[e as usize] = *cv;
                    } else {
                        let d = v[c1] - v[c0];
                        let t = if d.abs() > 1e-30 { (isovalue-v[c0])/d } else { 0.5 };
                        let (g0, g1) = (C[c0], C[c1]);
                        let idx = n_pts;
                        n_pts += 1;
                        pts_flat.push(org[0]+((i+g0[0]) as f64+t*(g1[0] as f64-g0[0] as f64))*sp[0]);
                        pts_flat.push(org[1]+((j+g0[1]) as f64+t*(g1[1] as f64-g0[1] as f64))*sp[1]);
                        pts_flat.push(org[2]+((k+g0[2]) as f64+t*(g1[2] as f64-g0[2] as f64))*sp[2]);
                        *cv = idx;
                        ep[e as usize] = idx;
                    }
                }

                let tr = &TRI_TABLE[ci as usize];
                let mut t = 0;
                while t < 15 && tr[t] != -1 {
                    conn.push(ep[tr[t] as usize]);
                    conn.push(ep[tr[t+1] as usize]);
                    conn.push(ep[tr[t+2] as usize]);
                    t += 3;
                }
            }
        }
        std::mem::swap(&mut xe0, &mut xe1);
        std::mem::swap(&mut ye0, &mut ye1);
        for v in xe1.iter_mut() { *v = -1; }
        for v in ye1.iter_mut() { *v = -1; }
        for v in ze.iter_mut()  { *v = -1; }
    }

    if conn.is_empty() { return PolyData::new(); }
    let points = Points::from_flat_vec(pts_flat);
    let nt = conn.len() / 3;
    let offsets: Vec<i64> = (0..=nt).map(|i| (i*3) as i64).collect();
    let polys = CellArray::from_raw(offsets, conn);
    let mut pd = PolyData::new(); pd.points = points; pd.polys = polys; pd
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn sphere_isosurface() {
        let img = ImageData::with_dimensions(32, 32, 32);
        let mut v = Vec::with_capacity(32*32*32);
        for k in 0..32 { for j in 0..32 { for i in 0..32 {
            let (x,y,z) = ((i as f64-16.0),(j as f64-16.0),(k as f64-16.0));
            v.push(x*x+y*y+z*z);
        }}}
        let r = flying_edges_3d(&img, &v, 100.0);
        assert!(r.polys.num_cells() > 100);
    }

    #[test]
    fn empty_field() {
        let img = ImageData::with_dimensions(4,4,4);
        let r = flying_edges_3d(&img, &vec![0.0;64], 1.0);
        assert_eq!(r.polys.num_cells(), 0);
    }

    #[test]
    fn deduplicates() {
        let img = ImageData::with_dimensions(8,8,8);
        let mut v = Vec::with_capacity(512);
        for k in 0..8{for j in 0..8{for i in 0..8{
            v.push(((i as f64-4.0).powi(2)+(j as f64-4.0).powi(2)+(k as f64-4.0).powi(2)));
        }}}
        let r = flying_edges_3d(&img, &v, 5.0);
        assert!(r.points.len() < r.polys.num_cells()*3);
    }
}
