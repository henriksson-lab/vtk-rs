// GPU picking shader: renders actor and cell IDs as colors.
// Each pixel encodes [actor_id, cell_id_high, cell_id_low, 255].

struct PickUniforms {
    mvp: mat4x4<f32>,
    actor_id: f32,
    _pad: vec3<f32>,
};

@group(0) @binding(0)
var<uniform> pick: PickUniforms;

struct PickVertexInput {
    @location(0) position: vec3<f32>,
    @location(1) normal: vec3<f32>,
    @location(2) color: vec3<f32>,
};

struct PickVertexOutput {
    @builtin(position) position: vec4<f32>,
    @location(0) @interpolate(flat) cell_id: u32,
};

@vertex
fn vs_pick(in: PickVertexInput, @builtin(vertex_index) vid: u32) -> PickVertexOutput {
    var out: PickVertexOutput;
    out.position = pick.mvp * vec4<f32>(in.position, 1.0);
    out.cell_id = vid / 3u; // triangle index (assuming triangulated mesh)
    return out;
}

@fragment
fn fs_pick(in: PickVertexOutput) -> @location(0) vec4<f32> {
    let actor = u32(pick.actor_id);
    let cell = in.cell_id;
    // Encode IDs as normalized color: R = actor, G = cell_high, B = cell_low
    return vec4<f32>(
        f32(actor) / 255.0,
        f32((cell >> 8u) & 0xFFu) / 255.0,
        f32(cell & 0xFFu) / 255.0,
        1.0
    );
}
