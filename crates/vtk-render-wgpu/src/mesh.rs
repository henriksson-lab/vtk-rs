use bytemuck::{Pod, Zeroable};
use vtk_data::PolyData;

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
/// Currently supports only triangle polygons. Each triangle gets a flat normal.
pub fn poly_data_to_mesh(poly_data: &PolyData, default_color: [f32; 3]) -> (Vec<Vertex>, Vec<u32>) {
    let mut vertices = Vec::new();
    let mut indices = Vec::new();

    // Extract per-vertex colors from point data if available
    let point_colors = extract_point_colors(poly_data);

    for cell in poly_data.polys.iter() {
        if cell.len() < 3 {
            continue;
        }

        // Fan triangulation for polygons with > 3 vertices
        let p0 = get_point(poly_data, cell[0] as usize);
        for i in 1..cell.len() - 1 {
            let p1 = get_point(poly_data, cell[i] as usize);
            let p2 = get_point(poly_data, cell[i + 1] as usize);

            // Flat normal
            let e1 = sub(p1, p0);
            let e2 = sub(p2, p0);
            let normal = normalize(cross(e1, e2));

            let base = vertices.len() as u32;

            let c0 = point_colors.as_ref().map(|c| c[cell[0] as usize]).unwrap_or(default_color);
            let c1 = point_colors.as_ref().map(|c| c[cell[i] as usize]).unwrap_or(default_color);
            let c2 = point_colors.as_ref().map(|c| c[cell[i + 1] as usize]).unwrap_or(default_color);

            vertices.push(Vertex { position: p0, normal, color: c0 });
            vertices.push(Vertex { position: p1, normal, color: c1 });
            vertices.push(Vertex { position: p2, normal, color: c2 });

            indices.push(base);
            indices.push(base + 1);
            indices.push(base + 2);
        }
    }

    (vertices, indices)
}

fn extract_point_colors(poly_data: &PolyData) -> Option<Vec<[f32; 3]>> {
    let scalars = poly_data.point_data().scalars()?;
    let nc = scalars.num_components();
    let nt = scalars.num_tuples();

    if nc < 3 {
        return None;
    }

    let mut colors = Vec::with_capacity(nt);
    let mut buf = vec![0.0f64; nc];
    for i in 0..nt {
        scalars.tuple_as_f64(i, &mut buf);
        colors.push([buf[0] as f32, buf[1] as f32, buf[2] as f32]);
    }
    Some(colors)
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
