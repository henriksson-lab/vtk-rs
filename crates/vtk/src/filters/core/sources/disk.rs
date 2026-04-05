use std::f64::consts::PI;

use crate::data::{CellArray, DataArray, Points, PolyData};

/// Parameters for generating a disk (annulus).
pub struct DiskParams {
    /// Inner radius. 0.0 for a filled disk. Default: 0.0
    pub inner_radius: f64,
    /// Outer radius. Default: 0.5
    pub outer_radius: f64,
    /// Number of sides around the circumference. Default: 32
    pub circumferential_resolution: usize,
    /// Number of rings between inner and outer radius. Default: 1
    pub radial_resolution: usize,
    /// Center of the disk. Default: [0, 0, 0]
    pub center: [f64; 3],
}

impl Default for DiskParams {
    fn default() -> Self {
        Self {
            inner_radius: 0.0,
            outer_radius: 0.5,
            circumferential_resolution: 32,
            radial_resolution: 1,
            center: [0.0, 0.0, 0.0],
        }
    }
}

/// Generate a flat disk (or annulus) in the XY plane as PolyData.
pub fn disk(params: &DiskParams) -> PolyData {
    let n_sides = params.circumferential_resolution.max(3);
    let n_rings = params.radial_resolution.max(1);
    let r_inner = params.inner_radius.max(0.0);
    let r_outer = params.outer_radius.max(r_inner + 1e-10);
    let cx = params.center[0];
    let cy = params.center[1];
    let cz = params.center[2];

    let mut points = Points::new();
    let mut normals = DataArray::<f64>::new("Normals", 3);
    let mut polys = CellArray::new();

    let normal = [0.0, 0.0, 1.0];

    if r_inner < 1e-10 {
        // Filled disk: center point + rings
        points.push([cx, cy, cz]);
        normals.push_tuple(&normal);

        for ring in 1..=n_rings {
            let r = r_outer * ring as f64 / n_rings as f64;
            for s in 0..n_sides {
                let angle = 2.0 * PI * s as f64 / n_sides as f64;
                points.push([cx + r * angle.cos(), cy + r * angle.sin(), cz]);
                normals.push_tuple(&normal);
            }
        }

        // Inner fan triangles (center to first ring)
        for s in 0..n_sides {
            let sn = (s + 1) % n_sides;
            polys.push_cell(&[0, (1 + s) as i64, (1 + sn) as i64]);
        }

        // Quads between rings
        for ring in 0..n_rings - 1 {
            let base_inner = 1 + ring * n_sides;
            let base_outer = 1 + (ring + 1) * n_sides;
            for s in 0..n_sides {
                let sn = (s + 1) % n_sides;
                polys.push_cell(&[
                    (base_inner + s) as i64,
                    (base_outer + s) as i64,
                    (base_outer + sn) as i64,
                    (base_inner + sn) as i64,
                ]);
            }
        }
    } else {
        // Annulus: concentric rings, no center point
        for ring in 0..=n_rings {
            let t = ring as f64 / n_rings as f64;
            let r = r_inner + t * (r_outer - r_inner);
            for s in 0..n_sides {
                let angle = 2.0 * PI * s as f64 / n_sides as f64;
                points.push([cx + r * angle.cos(), cy + r * angle.sin(), cz]);
                normals.push_tuple(&normal);
            }
        }

        // Quads between rings
        for ring in 0..n_rings {
            let base_inner = ring * n_sides;
            let base_outer = (ring + 1) * n_sides;
            for s in 0..n_sides {
                let sn = (s + 1) % n_sides;
                polys.push_cell(&[
                    (base_inner + s) as i64,
                    (base_outer + s) as i64,
                    (base_outer + sn) as i64,
                    (base_inner + sn) as i64,
                ]);
            }
        }
    }

    let mut pd = PolyData::new();
    pd.points = points;
    pd.polys = polys;
    pd.point_data_mut().add_array(normals.into());
    pd.point_data_mut().set_active_normals("Normals");
    pd
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn default_disk() {
        let pd = disk(&DiskParams::default());
        // Center + 32 rim points
        assert_eq!(pd.points.len(), 33);
        // 32 fan triangles
        assert_eq!(pd.polys.num_cells(), 32);
    }

    #[test]
    fn annulus() {
        let pd = disk(&DiskParams {
            inner_radius: 0.25,
            outer_radius: 0.5,
            circumferential_resolution: 8,
            radial_resolution: 2,
            ..Default::default()
        });
        // 3 rings * 8 points = 24
        assert_eq!(pd.points.len(), 24);
        // 2 ring bands * 8 quads = 16
        assert_eq!(pd.polys.num_cells(), 16);
    }

    #[test]
    fn multi_ring_disk() {
        let pd = disk(&DiskParams {
            inner_radius: 0.0,
            outer_radius: 1.0,
            circumferential_resolution: 6,
            radial_resolution: 3,
            ..Default::default()
        });
        // Center + 3 rings * 6 = 19 points
        assert_eq!(pd.points.len(), 19);
        // 6 inner triangles + 2 ring bands * 6 quads = 18
        assert_eq!(pd.polys.num_cells(), 18);
    }
}
