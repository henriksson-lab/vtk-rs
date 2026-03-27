// Silhouette post-process: depth-based edge detection.
// Samples the depth buffer at neighboring texels and draws edges
// where depth discontinuities are found.

struct SilhouetteUniforms {
    outline_color: vec3<f32>,
    depth_threshold: f32,
    texel_size: vec2<f32>,
    _pad: vec2<f32>,
};

@group(0) @binding(0)
var depth_tex: texture_2d<f32>;
@group(0) @binding(1)
var depth_sampler: sampler;
@group(0) @binding(2)
var<uniform> params: SilhouetteUniforms;

struct FullscreenOutput {
    @builtin(position) position: vec4<f32>,
    @location(0) uv: vec2<f32>,
};

@vertex
fn vs_fullscreen(@builtin(vertex_index) idx: u32) -> FullscreenOutput {
    // Full-screen triangle
    var out: FullscreenOutput;
    let x = f32(i32(idx) / 2) * 4.0 - 1.0;
    let y = f32(i32(idx) % 2) * 4.0 - 1.0;
    out.position = vec4<f32>(x, y, 0.0, 1.0);
    out.uv = vec2<f32>((x + 1.0) * 0.5, (1.0 - y) * 0.5);
    return out;
}

@fragment
fn fs_silhouette(in: FullscreenOutput) -> @location(0) vec4<f32> {
    let center = textureSample(depth_tex, depth_sampler, in.uv).r;

    // Skip background (depth == 1.0)
    if center >= 0.9999 {
        discard;
    }

    let tx = params.texel_size.x;
    let ty = params.texel_size.y;

    let d_left = textureSample(depth_tex, depth_sampler, in.uv + vec2<f32>(-tx, 0.0)).r;
    let d_right = textureSample(depth_tex, depth_sampler, in.uv + vec2<f32>(tx, 0.0)).r;
    let d_up = textureSample(depth_tex, depth_sampler, in.uv + vec2<f32>(0.0, -ty)).r;
    let d_down = textureSample(depth_tex, depth_sampler, in.uv + vec2<f32>(0.0, ty)).r;

    let dx = abs(d_right - d_left);
    let dy = abs(d_down - d_up);
    let edge = max(dx, dy);

    if edge > params.depth_threshold {
        return vec4<f32>(params.outline_color, 1.0);
    }

    discard;
}
