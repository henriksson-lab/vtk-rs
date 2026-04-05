// Compute per-face normals from triangle vertex positions.
// Input: positions array (3 floats per vertex), indices array (3 per triangle)
// Output: normals array (3 floats per triangle)

@group(0) @binding(0) var<storage, read> positions: array<f32>;
@group(0) @binding(1) var<storage, read> indices: array<u32>;
@group(0) @binding(2) var<storage, read_write> normals: array<f32>;
@group(0) @binding(3) var<uniform> params: vec4<f32>; // x=num_triangles

@compute @workgroup_size(256)
fn main(@builtin(global_invocation_id) gid: vec3<u32>) {
    let tri_idx = gid.x;
    let num_tris = u32(params.x);
    if tri_idx >= num_tris { return; }

    let i0 = indices[tri_idx * 3u];
    let i1 = indices[tri_idx * 3u + 1u];
    let i2 = indices[tri_idx * 3u + 2u];

    let p0 = vec3<f32>(positions[i0 * 3u], positions[i0 * 3u + 1u], positions[i0 * 3u + 2u]);
    let p1 = vec3<f32>(positions[i1 * 3u], positions[i1 * 3u + 1u], positions[i1 * 3u + 2u]);
    let p2 = vec3<f32>(positions[i2 * 3u], positions[i2 * 3u + 1u], positions[i2 * 3u + 2u]);

    let e1 = p1 - p0;
    let e2 = p2 - p0;
    let n = cross(e1, e2);
    let len = length(n);
    var normal: vec3<f32>;
    if len > 1e-20 {
        normal = n / len;
    } else {
        normal = vec3<f32>(0.0, 0.0, 1.0);
    }

    normals[tri_idx * 3u] = normal.x;
    normals[tri_idx * 3u + 1u] = normal.y;
    normals[tri_idx * 3u + 2u] = normal.z;
}
