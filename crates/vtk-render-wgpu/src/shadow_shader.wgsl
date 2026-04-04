struct ShadowUniforms {
    light_vp: mat4x4<f32>,
    model: mat4x4<f32>,
};

@group(0) @binding(0)
var<uniform> uniforms: ShadowUniforms;

struct VertexInput {
    @location(0) position: vec3<f32>,
};

@vertex
fn vs_shadow(in: VertexInput) -> @builtin(position) vec4<f32> {
    return uniforms.light_vp * vec4<f32>(in.position, 1.0);
}
