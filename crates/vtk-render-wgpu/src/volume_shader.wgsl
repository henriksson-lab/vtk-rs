// Volume ray marching shader.
// Renders a 3D texture as a volume using front-to-back compositing.

struct VolumeUniforms {
    mvp: mat4x4<f32>,
    camera_pos: vec3<f32>,
    num_steps: f32,
    volume_min: vec3<f32>,
    opacity_scale: f32,
    volume_max: vec3<f32>,
    _pad: f32,
};

@group(0) @binding(0) var<uniform> params: VolumeUniforms;
@group(0) @binding(1) var volume_tex: texture_3d<f32>;
@group(0) @binding(2) var volume_sampler: sampler;
@group(0) @binding(3) var lut_tex: texture_1d<f32>;
@group(0) @binding(4) var lut_sampler: sampler;

struct VertexOutput {
    @builtin(position) position: vec4<f32>,
    @location(0) world_pos: vec3<f32>,
};

// Proxy cube vertices (8 corners, 36 indices for 12 triangles)
@vertex
fn vs_volume(@builtin(vertex_index) idx: u32) -> VertexOutput {
    // Generate cube vertices procedurally
    let corners = array<vec3<f32>, 8>(
        vec3<f32>(0.0, 0.0, 0.0), vec3<f32>(1.0, 0.0, 0.0),
        vec3<f32>(1.0, 1.0, 0.0), vec3<f32>(0.0, 1.0, 0.0),
        vec3<f32>(0.0, 0.0, 1.0), vec3<f32>(1.0, 0.0, 1.0),
        vec3<f32>(1.0, 1.0, 1.0), vec3<f32>(0.0, 1.0, 1.0),
    );
    let indices = array<u32, 36>(
        0u,1u,2u, 0u,2u,3u, // -Z
        4u,6u,5u, 4u,7u,6u, // +Z
        0u,4u,5u, 0u,5u,1u, // -Y
        3u,2u,6u, 3u,6u,7u, // +Y
        0u,3u,7u, 0u,7u,4u, // -X
        1u,5u,6u, 1u,6u,2u, // +X
    );

    let vi = indices[idx];
    let unit_pos = corners[vi];
    let world_pos = params.volume_min + unit_pos * (params.volume_max - params.volume_min);

    var out: VertexOutput;
    out.position = params.mvp * vec4<f32>(world_pos, 1.0);
    out.world_pos = world_pos;
    return out;
}

@fragment
fn fs_volume(in: VertexOutput) -> @location(0) vec4<f32> {
    let ray_origin = params.camera_pos;
    let ray_dir = normalize(in.world_pos - params.camera_pos);

    // Ray-AABB intersection
    let inv_dir = 1.0 / ray_dir;
    let t1 = (params.volume_min - ray_origin) * inv_dir;
    let t2 = (params.volume_max - ray_origin) * inv_dir;
    let tmin_v = min(t1, t2);
    let tmax_v = max(t1, t2);
    let t_enter = max(max(tmin_v.x, tmin_v.y), tmin_v.z);
    let t_exit = min(min(tmax_v.x, tmax_v.y), tmax_v.z);

    if t_enter > t_exit || t_exit < 0.0 {
        discard;
    }

    let t_start = max(t_enter, 0.0);
    let steps = u32(params.num_steps);
    let step_size = (t_exit - t_start) / f32(steps);
    let volume_size = params.volume_max - params.volume_min;

    var color = vec3<f32>(0.0);
    var alpha = 0.0;

    for (var i = 0u; i < steps; i = i + 1u) {
        if alpha >= 0.99 { break; }

        let t = t_start + (f32(i) + 0.5) * step_size;
        let pos = ray_origin + t * ray_dir;

        // Normalize to [0,1] texture coordinates
        let uvw = (pos - params.volume_min) / volume_size;

        if any(uvw < vec3<f32>(0.0)) || any(uvw > vec3<f32>(1.0)) {
            continue;
        }

        // Sample 3D volume
        let scalar = textureSample(volume_tex, volume_sampler, uvw).r;

        // Look up transfer function
        let rgba = textureSample(lut_tex, lut_sampler, scalar);

        let sample_alpha = rgba.a * params.opacity_scale * step_size * 100.0;
        let sa = clamp(sample_alpha, 0.0, 1.0);

        // Front-to-back compositing
        color += (1.0 - alpha) * sa * rgba.rgb;
        alpha += (1.0 - alpha) * sa;
    }

    if alpha < 0.001 {
        discard;
    }

    return vec4<f32>(color, alpha);
}
