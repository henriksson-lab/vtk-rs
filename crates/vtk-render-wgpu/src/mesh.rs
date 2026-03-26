use bytemuck::{Pod, Zeroable};
use vtk_data::PolyData;
use vtk_render::Coloring;

/// GPU-ready vertex with position, normal, and color.
#[repr(C)]
#[derive(Debug, Clone, Copy, Pod, Zeroable)]
pub struct Vertex {
    pub position: [f32; 3],
    pub normal: [f32; 3],
    pub color: [f32; 3],
}

impl Vertex {
    pub fn layout() -> wgpu::VertexBufferLayout<'static> {
        wgpu::VertexBufferLayout {
            array_stride: std::mem::size_of::<Vertex>() as wgpu::BufferAddress,
            step_mode: wgpu::VertexStepMode::Vertex,
            attributes: &[
                // position
                wgpu::VertexAttribute {
                    offset: 0,
                    shader_location: 0,
                    format: wgpu::VertexFormat::Float32x3,
                },
                // normal
                wgpu::VertexAttribute {
                    offset: 12,
                    shader_location: 1,
                    format: wgpu::VertexFormat::Float32x3,
                },
                // color
                wgpu::VertexAttribute {
                    offset: 24,
                    shader_location: 2,
                    format: wgpu::VertexFormat::Float32x3,
                },
            ],
        }
    }
}

/// Convert PolyData triangles to GPU vertex and index buffers.
///
/// Uses smooth normals from point data if available, otherwise computes flat normals.
/// Supports solid color and scalar color mapping via the `Coloring` enum.
pub fn poly_data_to_mesh(poly_data: &PolyData, coloring: &Coloring) -> (Vec<Vertex>, Vec<u32>) {
    let mut vertices = Vec::new();
    let mut indices = Vec::new();

    let point_colors = resolve_colors(poly_data, coloring);
    let point_normals = extract_point_normals(poly_data);

    for cell in poly_data.polys.iter() {
        if cell.len() < 3 {
            continue;
        }

        let p0 = get_point(poly_data, cell[0] as usize);
        for i in 1..cell.len() - 1 {
            let p1 = get_point(poly_data, cell[i] as usize);
            let p2 = get_point(poly_data, cell[i + 1] as usize);

            // Use smooth normals if available, otherwise flat
            let (n0, n1, n2) = if let Some(ref pn) = point_normals {
                (pn[cell[0] as usize], pn[cell[i] as usize], pn[cell[i + 1] as usize])
            } else {
                let e1 = sub(p1, p0);
                let e2 = sub(p2, p0);
                let flat = normalize(cross(e1, e2));
                (flat, flat, flat)
            };

            let base = vertices.len() as u32;

            vertices.push(Vertex { position: p0, normal: n0, color: point_colors[cell[0] as usize] });
            vertices.push(Vertex { position: p1, normal: n1, color: point_colors[cell[i] as usize] });
            vertices.push(Vertex { position: p2, normal: n2, color: point_colors[cell[i + 1] as usize] });

            indices.push(base);
            indices.push(base + 1);
            indices.push(base + 2);
        }
    }

    (vertices, indices)
}

fn extract_point_normals(poly_data: &PolyData) -> Option<Vec<[f32; 3]>> {
    let normals = poly_data.point_data().normals()?;
    if normals.num_components() != 3 {
        return None;
    }
    let nt = normals.num_tuples();
    let mut result = Vec::with_capacity(nt);
    let mut buf = [0.0f64; 3];
    for i in 0..nt {
        normals.tuple_as_f64(i, &mut buf);
        result.push([buf[0] as f32, buf[1] as f32, buf[2] as f32]);
    }
    Some(result)
}

/// Resolve per-vertex colors for a PolyData based on the Coloring mode.
pub fn resolve_colors_pub(poly_data: &PolyData, coloring: &Coloring) -> Vec<[f32; 3]> {
    resolve_colors(poly_data, coloring)
}

fn resolve_colors(poly_data: &PolyData, coloring: &Coloring) -> Vec<[f32; 3]> {
    let n = poly_data.points.len();

    match coloring {
        Coloring::Solid(color) => vec![*color; n],
        Coloring::ScalarMap { color_map, range } => {
            // Try to get active scalars from point data
            let scalars = poly_data.point_data().scalars();
            let Some(scalars) = scalars else {
                return vec![[1.0, 1.0, 1.0]; n];
            };

            let nt = scalars.num_tuples().min(n);

            // Compute range if not provided
            let (smin, smax) = if let Some([lo, hi]) = range {
                (*lo, *hi)
            } else {
                let mut lo = f64::INFINITY;
                let mut hi = f64::NEG_INFINITY;
                let mut buf = [0.0f64];
                for i in 0..nt {
                    scalars.tuple_as_f64(i, &mut buf);
                    lo = lo.min(buf[0]);
                    hi = hi.max(buf[0]);
                }
                (lo, hi)
            };

            let mut colors = vec![[1.0f32, 1.0, 1.0]; n];
            let mut buf = [0.0f64];
            for (i, color) in colors.iter_mut().enumerate().take(nt) {
                scalars.tuple_as_f64(i, &mut buf);
                *color = color_map.map_value(buf[0], smin, smax);
            }
            colors
        }
    }
}

fn get_point(pd: &PolyData, idx: usize) -> [f32; 3] {
    let p = pd.points.get(idx);
    [p[0] as f32, p[1] as f32, p[2] as f32]
}

fn sub(a: [f32; 3], b: [f32; 3]) -> [f32; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

fn cross(a: [f32; 3], b: [f32; 3]) -> [f32; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

fn normalize(v: [f32; 3]) -> [f32; 3] {
    let len = (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt();
    if len < 1e-10 {
        [0.0, 0.0, 1.0]
    } else {
        [v[0] / len, v[1] / len, v[2] / len]
    }
}
