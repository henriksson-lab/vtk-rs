// Depth-of-field post-processing shader.
//
// Three passes:
// 1. CoC computation: sample depth buffer, compute circle of confusion
// 2. Blur: variable-radius blur based on CoC
// 3. Composite: blend blurred result with original

struct DofUniforms {
    focal_distance: f32,
    aperture: f32,
    max_blur: f32,
    texel_size_x: f32,
    texel_size_y: f32,
    near: f32,
    far: f32,
    horizontal: f32,
};

@group(0) @binding(0) var<uniform> u: DofUniforms;
@group(0) @binding(1) var color_tex: texture_2d<f32>;
@group(0) @binding(2) var depth_tex: texture_2d<f32>;
@group(0) @binding(3) var tex_sampler: sampler;

struct VertexOutput {
    @builtin(position) position: vec4<f32>,
    @location(0) uv: vec2<f32>,
};

// Full-screen triangle
@vertex
fn vs_main(@builtin(vertex_index) idx: u32) -> VertexOutput {
    var out: VertexOutput;
    let x = f32(i32(idx & 1u)) * 4.0 - 1.0;
    let y = f32(i32(idx >> 1u)) * 4.0 - 1.0;
    out.position = vec4<f32>(x, y, 0.0, 1.0);
    out.uv = vec2<f32>((x + 1.0) * 0.5, (1.0 - y) * 0.5);
    return out;
}

fn linearize_depth(d: f32) -> f32 {
    let z_ndc = d * 2.0 - 1.0;
    return (2.0 * u.near * u.far) / (u.far + u.near - z_ndc * (u.far - u.near));
}

// Pass 1: compute CoC and store in alpha channel
@fragment
fn fs_coc(in: VertexOutput) -> @location(0) vec4<f32> {
    let color = textureSample(color_tex, tex_sampler, in.uv);
    let depth_raw = textureSample(depth_tex, tex_sampler, in.uv).r;
    let depth = linearize_depth(depth_raw);
    let coc = clamp(u.aperture * abs(depth - u.focal_distance) / max(depth, 0.001), 0.0, u.max_blur);
    return vec4<f32>(color.rgb, coc);
}

// Pass 2: separable Gaussian blur weighted by CoC
@fragment
fn fs_blur(in: VertexOutput) -> @location(0) vec4<f32> {
    let center = textureSample(color_tex, tex_sampler, in.uv);
    let coc = center.a;
    let radius = coc;

    if (radius < 0.5) {
        return center;
    }

    var dir: vec2<f32>;
    if (u.horizontal > 0.5) {
        dir = vec2<f32>(u.texel_size_x, 0.0);
    } else {
        dir = vec2<f32>(0.0, u.texel_size_y);
    }

    var result = center.rgb;
    var weight_sum = 1.0;
    let steps = i32(min(radius, 10.0));

    for (var i = 1; i <= steps; i = i + 1) {
        let offset = dir * f32(i);
        let s1 = textureSample(color_tex, tex_sampler, in.uv + offset);
        let s2 = textureSample(color_tex, tex_sampler, in.uv - offset);
        let w = 1.0 - f32(i) / (f32(steps) + 1.0);
        result = result + s1.rgb * w + s2.rgb * w;
        weight_sum = weight_sum + 2.0 * w;
    }

    return vec4<f32>(result / weight_sum, coc);
}

// Pass 3: composite (just output the blurred result)
@fragment
fn fs_composite(in: VertexOutput) -> @location(0) vec4<f32> {
    let color = textureSample(color_tex, tex_sampler, in.uv);
    return vec4<f32>(color.rgb, 1.0);
}
