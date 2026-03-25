struct Uniforms {
    mvp: mat4x4<f32>,
    model: mat4x4<f32>,
    light_dir: vec3<f32>,
    _pad: f32,
    ambient: f32,
    _pad2: vec3<f32>,
};

@group(0) @binding(0)
var<uniform> uniforms: Uniforms;

struct VertexInput {
    @location(0) position: vec3<f32>,
    @location(1) normal: vec3<f32>,
    @location(2) color: vec3<f32>,
};

struct VertexOutput {
    @builtin(position) clip_position: vec4<f32>,
    @location(0) world_normal: vec3<f32>,
    @location(1) color: vec3<f32>,
};

@vertex
fn vs_main(in: VertexInput) -> VertexOutput {
    var out: VertexOutput;
    out.clip_position = uniforms.mvp * vec4<f32>(in.position, 1.0);
    // Transform normal by model matrix (assuming no non-uniform scale)
    out.world_normal = (uniforms.model * vec4<f32>(in.normal, 0.0)).xyz;
    out.color = in.color;
    return out;
}

@fragment
fn fs_main(in: VertexOutput) -> @location(0) vec4<f32> {
    let n = normalize(in.world_normal);
    let light = normalize(uniforms.light_dir);
    let diffuse = max(dot(n, light), 0.0);
    let intensity = uniforms.ambient + (1.0 - uniforms.ambient) * diffuse;
    return vec4<f32>(in.color * intensity, 1.0);
}
