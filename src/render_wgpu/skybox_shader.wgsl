// Skybox gradient shader: renders a vertical color gradient as background.

struct SkyboxUniforms {
    bottom: vec3<f32>,
    _pad0: f32,
    horizon: vec3<f32>,
    _pad1: f32,
    top: vec3<f32>,
    mode: f32, // 0=solid, 1=gradient, 2=three_stop
};

@group(0) @binding(0)
var<uniform> sky: SkyboxUniforms;

struct SkyOutput {
    @builtin(position) position: vec4<f32>,
    @location(0) uv: vec2<f32>,
};

@vertex
fn vs_skybox(@builtin(vertex_index) idx: u32) -> SkyOutput {
    var out: SkyOutput;
    let x = f32(i32(idx) / 2) * 4.0 - 1.0;
    let y = f32(i32(idx) % 2) * 4.0 - 1.0;
    out.position = vec4<f32>(x, y, 0.9999, 1.0);
    out.uv = vec2<f32>((x + 1.0) * 0.5, (y + 1.0) * 0.5);
    return out;
}

@fragment
fn fs_skybox(in: SkyOutput) -> @location(0) vec4<f32> {
    let t = in.uv.y; // 0=bottom, 1=top

    var color: vec3<f32>;

    if sky.mode < 0.5 {
        // Solid
        color = sky.bottom;
    } else if sky.mode < 1.5 {
        // Gradient
        color = mix(sky.bottom, sky.top, t);
    } else {
        // Three-stop
        if t < 0.5 {
            color = mix(sky.bottom, sky.horizon, t * 2.0);
        } else {
            color = mix(sky.horizon, sky.top, (t - 0.5) * 2.0);
        }
    }

    return vec4<f32>(color, 1.0);
}
