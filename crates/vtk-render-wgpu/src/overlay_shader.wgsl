// 2D overlay shader for scalar bar and UI elements.
// Positions are in NDC [0,1] mapped to clip space [-1,1].

struct OverlayVertex {
    @location(0) position: vec2<f32>,
    @location(1) color: vec4<f32>,
};

struct OverlayOutput {
    @builtin(position) clip_position: vec4<f32>,
    @location(0) color: vec4<f32>,
};

@vertex
fn vs_overlay(in: OverlayVertex) -> OverlayOutput {
    var out: OverlayOutput;
    // Convert from [0,1] NDC to clip space [-1,1], with Y flipped for screen coords
    out.clip_position = vec4<f32>(
        in.position.x * 2.0 - 1.0,
        in.position.y * 2.0 - 1.0,
        0.0,
        1.0
    );
    out.color = in.color;
    return out;
}

@fragment
fn fs_overlay(in: OverlayOutput) -> @location(0) vec4<f32> {
    return in.color;
}
