use std::sync::Arc;

use bytemuck::{Pod, Zeroable};
use glam::Mat4;
use wgpu::util::DeviceExt;
use winit::window::Window;

use vtk_render::{Light, LightType, Material, Renderer, Representation, Scene};
use vtk_types::VtkError;

use crate::mesh::{self, Vertex};
use crate::overlay::OverlayPipeline;
use crate::skybox_pass::SkyboxPass;
use crate::volume_pass::VolumePass;
use crate::wireframe;

const MAX_LIGHTS: usize = 8;
const MSAA_SAMPLE_COUNT: u32 = 4;

#[repr(C)]
#[derive(Debug, Clone, Copy, Pod, Zeroable)]
struct GpuLight {
    light_type: f32,
    intensity: f32,
    cone_angle: f32,
    exponent: f32,
    position: [f32; 3],
    _pad0: f32,
    direction: [f32; 3],
    _pad1: f32,
    color: [f32; 3],
    _pad2: f32,
}

#[repr(C)]
#[derive(Debug, Clone, Copy, Pod, Zeroable)]
struct Uniforms {
    mvp: [[f32; 4]; 4],
    model: [[f32; 4]; 4],
    camera_pos: [f32; 3],
    opacity: f32,
    mat_ambient: f32,
    mat_diffuse: f32,
    mat_specular: f32,
    mat_specular_power: f32,
    specular_color: [f32; 3],
    use_lighting: f32,
    num_lights: f32,
    metallic: f32,
    roughness: f32,
    use_pbr: f32,
    flat_shading: f32,
    num_clip_planes: f32,
    fog_enabled: f32,
    fog_mode: f32,
    clip_planes: [[f32; 4]; 6],
    fog_color: [f32; 3],
    fog_near: f32,
    fog_far: f32,
    fog_density: f32,
    _fog_pad: [f32; 2],
    lights: [GpuLight; MAX_LIGHTS],
}

struct GpuMesh {
    vertex_buffer: wgpu::Buffer,
    index_buffer: wgpu::Buffer,
    num_indices: u32,
}

pub struct WgpuRenderer {
    device: wgpu::Device,
    queue: wgpu::Queue,
    surface: wgpu::Surface<'static>,
    surface_config: wgpu::SurfaceConfiguration,
    surface_format: wgpu::TextureFormat,
    pipeline_triangles: wgpu::RenderPipeline,
    pipeline_lines: wgpu::RenderPipeline,
    pipeline_points: wgpu::RenderPipeline,
    pipeline_triangles_blend: wgpu::RenderPipeline,
    pipeline_lines_blend: wgpu::RenderPipeline,
    pipeline_points_blend: wgpu::RenderPipeline,
    depth_texture: wgpu::TextureView,
    msaa_texture: wgpu::TextureView,
    uniform_buffer: wgpu::Buffer,
    bind_group: wgpu::BindGroup,
    overlay_pipeline: OverlayPipeline,
    pipeline_triangles_cull: wgpu::RenderPipeline,
    pipeline_triangles_blend_cull: wgpu::RenderPipeline,
    volume_pass: VolumePass,
    skybox_pass: SkyboxPass,
    pipeline_lines_no_msaa: wgpu::RenderPipeline,
    width: u32,
    height: u32,
}

impl WgpuRenderer {
    pub async fn new(window: Arc<Window>) -> Result<Self, VtkError> {
        let size = window.inner_size();
        let width = size.width.max(1);
        let height = size.height.max(1);

        let instance = wgpu::Instance::default();
        let surface = instance
            .create_surface(window)
            .map_err(|e| VtkError::InvalidData(format!("surface creation failed: {e}")))?;

        let adapter = instance
            .request_adapter(&wgpu::RequestAdapterOptions {
                power_preference: wgpu::PowerPreference::default(),
                compatible_surface: Some(&surface),
                force_fallback_adapter: false,
            })
            .await
            .ok_or_else(|| VtkError::InvalidData("no suitable GPU adapter found".into()))?;

        let (device, queue) = adapter
            .request_device(&wgpu::DeviceDescriptor::default(), None)
            .await
            .map_err(|e| VtkError::InvalidData(format!("device request failed: {e}")))?;

        let surface_caps = surface.get_capabilities(&adapter);
        let surface_format = surface_caps
            .formats
            .iter()
            .find(|f| f.is_srgb())
            .copied()
            .unwrap_or(surface_caps.formats[0]);

        let surface_config = wgpu::SurfaceConfiguration {
            usage: wgpu::TextureUsages::RENDER_ATTACHMENT,
            format: surface_format,
            width,
            height,
            present_mode: wgpu::PresentMode::Fifo,
            alpha_mode: surface_caps.alpha_modes[0],
            view_formats: vec![],
            desired_maximum_frame_latency: 2,
        };
        surface.configure(&device, &surface_config);

        let shader = device.create_shader_module(wgpu::ShaderModuleDescriptor {
            label: Some("vtk shader"),
            source: wgpu::ShaderSource::Wgsl(include_str!("shader.wgsl").into()),
        });

        let uniform_buffer = device.create_buffer(&wgpu::BufferDescriptor {
            label: Some("uniforms"),
            size: std::mem::size_of::<Uniforms>() as u64,
            usage: wgpu::BufferUsages::UNIFORM | wgpu::BufferUsages::COPY_DST,
            mapped_at_creation: false,
        });

        let bind_group_layout = device.create_bind_group_layout(&wgpu::BindGroupLayoutDescriptor {
            label: Some("bind group layout"),
            entries: &[wgpu::BindGroupLayoutEntry {
                binding: 0,
                visibility: wgpu::ShaderStages::VERTEX | wgpu::ShaderStages::FRAGMENT,
                ty: wgpu::BindingType::Buffer {
                    ty: wgpu::BufferBindingType::Uniform,
                    has_dynamic_offset: false,
                    min_binding_size: None,
                },
                count: None,
            }],
        });

        let bind_group = device.create_bind_group(&wgpu::BindGroupDescriptor {
            label: Some("bind group"),
            layout: &bind_group_layout,
            entries: &[wgpu::BindGroupEntry {
                binding: 0,
                resource: uniform_buffer.as_entire_binding(),
            }],
        });

        let pipeline_layout = device.create_pipeline_layout(&wgpu::PipelineLayoutDescriptor {
            label: Some("pipeline layout"),
            bind_group_layouts: &[&bind_group_layout],
            push_constant_ranges: &[],
        });

        let pipeline_triangles =
            create_pipeline_with_cull(&device, &pipeline_layout, &shader, surface_format, wgpu::PrimitiveTopology::TriangleList, false, MSAA_SAMPLE_COUNT, None);
        let pipeline_lines =
            create_pipeline_with_cull(&device, &pipeline_layout, &shader, surface_format, wgpu::PrimitiveTopology::LineList, false, MSAA_SAMPLE_COUNT, None);
        let pipeline_points =
            create_pipeline_with_cull(&device, &pipeline_layout, &shader, surface_format, wgpu::PrimitiveTopology::PointList, false, MSAA_SAMPLE_COUNT, None);
        let pipeline_triangles_blend =
            create_pipeline_with_cull(&device, &pipeline_layout, &shader, surface_format, wgpu::PrimitiveTopology::TriangleList, true, MSAA_SAMPLE_COUNT, None);
        let pipeline_lines_blend =
            create_pipeline_with_cull(&device, &pipeline_layout, &shader, surface_format, wgpu::PrimitiveTopology::LineList, true, MSAA_SAMPLE_COUNT, None);
        let pipeline_points_blend =
            create_pipeline_with_cull(&device, &pipeline_layout, &shader, surface_format, wgpu::PrimitiveTopology::PointList, true, MSAA_SAMPLE_COUNT, None);
        let pipeline_triangles_cull =
            create_pipeline_with_cull(&device, &pipeline_layout, &shader, surface_format, wgpu::PrimitiveTopology::TriangleList, false, MSAA_SAMPLE_COUNT, Some(wgpu::Face::Back));
        let pipeline_triangles_blend_cull =
            create_pipeline_with_cull(&device, &pipeline_layout, &shader, surface_format, wgpu::PrimitiveTopology::TriangleList, true, MSAA_SAMPLE_COUNT, Some(wgpu::Face::Back));

        let depth_texture = create_depth_texture(&device, width, height, MSAA_SAMPLE_COUNT);
        let msaa_texture = create_msaa_texture(&device, surface_format, width, height);
        let overlay_pipeline = OverlayPipeline::new(&device, surface_format);
        let volume_pass = VolumePass::new(&device, surface_format);
        let skybox_pass = SkyboxPass::new(&device, surface_format);
        let pipeline_lines_no_msaa =
            create_pipeline_with_cull(&device, &pipeline_layout, &shader, surface_format, wgpu::PrimitiveTopology::LineList, false, 1, None);

        Ok(Self {
            device,
            queue,
            surface,
            surface_config,
            surface_format,
            pipeline_triangles,
            pipeline_lines,
            pipeline_points,
            pipeline_triangles_blend,
            pipeline_lines_blend,
            pipeline_points_blend,
            depth_texture,
            msaa_texture,
            uniform_buffer,
            bind_group,
            overlay_pipeline,
            volume_pass,
            skybox_pass,
            pipeline_triangles_cull,
            pipeline_triangles_blend_cull,
            pipeline_lines_no_msaa,
            width,
            height,
        })
    }

    fn prepare_actor_mesh(&self, actor: &vtk_render::Actor, camera_distance: f64) -> Option<GpuMesh> {
        // Select LOD level if available
        let data = if let Some(ref lod) = actor.lod {
            lod.select(camera_distance)
        } else {
            &actor.data
        };
        let (vertices, indices) = match actor.representation {
            Representation::Surface => mesh::poly_data_to_mesh(data, &actor.coloring),
            Representation::Wireframe => wireframe::poly_data_to_wireframe(data, &actor.coloring),
            Representation::Points => poly_data_to_points(data, &actor.coloring),
        };
        if indices.is_empty() {
            return None;
        }
        Some(upload_mesh(&self.device, &vertices, &indices))
    }

    fn prepare_edge_mesh(&self, actor: &vtk_render::Actor) -> Option<GpuMesh> {
        let coloring = vtk_render::Coloring::Solid(actor.material.edge_color);
        let (vertices, indices) = wireframe::poly_data_to_wireframe(&actor.data, &coloring);
        if indices.is_empty() {
            return None;
        }
        Some(upload_mesh(&self.device, &vertices, &indices))
    }

    fn select_pipeline(&self, repr: Representation, translucent: bool, backface_cull: bool) -> &wgpu::RenderPipeline {
        match (repr, translucent, backface_cull) {
            (Representation::Surface, false, true) => &self.pipeline_triangles_cull,
            (Representation::Surface, true, true) => &self.pipeline_triangles_blend_cull,
            (Representation::Surface, false, _) => &self.pipeline_triangles,
            (Representation::Wireframe, false, _) => &self.pipeline_lines,
            (Representation::Points, false, _) => &self.pipeline_points,
            (Representation::Surface, true, _) => &self.pipeline_triangles_blend,
            (Representation::Wireframe, true, _) => &self.pipeline_lines_blend,
            (Representation::Points, true, _) => &self.pipeline_points_blend,
        }
    }

    fn write_uniforms(
        &self,
        scene: &Scene,
        material: &Material,
        opacity: f32,
        use_lighting: bool,
        position: [f64; 3],
        scale: f64,
    ) {
        let aspect = self.width as f64 / self.height as f64;
        let view_mat = scene.camera.view_matrix();
        let proj_mat = scene.camera.projection_matrix(aspect);
        let vp = proj_mat * view_mat;

        // Compute model matrix from actor position and scale
        let model = Mat4::from_scale_rotation_translation(
            glam::Vec3::splat(scale as f32),
            glam::Quat::IDENTITY,
            glam::Vec3::new(position[0] as f32, position[1] as f32, position[2] as f32),
        );
        let mvp_f32 = Mat4::from_cols_array(&vp.to_cols_array().map(|v| v as f32)) * model;
        let cam_pos = scene.camera.position;

        let mut gpu_lights = [GpuLight::zeroed(); MAX_LIGHTS];
        let num_lights = scene.lights.len().min(MAX_LIGHTS);
        for (i, light) in scene.lights.iter().take(MAX_LIGHTS).enumerate() {
            if !light.enabled {
                continue;
            }
            gpu_lights[i] = light_to_gpu(light);
        }

        let uniforms = Uniforms {
            mvp: mvp_f32.to_cols_array_2d(),
            model: model.to_cols_array_2d(),
            camera_pos: [cam_pos.x as f32, cam_pos.y as f32, cam_pos.z as f32],
            opacity,
            mat_ambient: material.ambient as f32,
            mat_diffuse: material.diffuse as f32,
            mat_specular: material.specular as f32,
            mat_specular_power: material.specular_power as f32,
            specular_color: material.specular_color,
            use_lighting: if use_lighting { 1.0 } else { 0.0 },
            num_lights: num_lights as f32,
            metallic: material.metallic as f32,
            roughness: material.roughness as f32,
            use_pbr: if material.pbr { 1.0 } else { 0.0 },
            flat_shading: if material.flat_shading { 1.0 } else { 0.0 },
            num_clip_planes: scene.clip_planes.iter().filter(|c| c.enabled).count() as f32,
            fog_enabled: if scene.fog.enabled { 1.0 } else { 0.0 },
            fog_mode: match scene.fog.mode {
                vtk_render::FogMode::Linear => 0.0,
                vtk_render::FogMode::Exponential => 1.0,
                vtk_render::FogMode::ExponentialSquared => 2.0,
            },
            clip_planes: {
                let mut cp = [[0.0f32; 4]; 6];
                for (i, plane) in scene.clip_planes.iter().filter(|c| c.enabled).take(6).enumerate() {
                    cp[i] = [
                        plane.normal[0] as f32,
                        plane.normal[1] as f32,
                        plane.normal[2] as f32,
                        plane.distance as f32,
                    ];
                }
                cp
            },
            fog_color: scene.fog.color,
            fog_near: scene.fog.near as f32,
            fog_far: scene.fog.far as f32,
            fog_density: scene.fog.density as f32,
            _fog_pad: [0.0; 2],
            lights: gpu_lights,
        };

        self.queue
            .write_buffer(&self.uniform_buffer, 0, bytemuck::bytes_of(&uniforms));
    }

    /// Render the 3D scene with MSAA, resolving to `resolve_target`.
    fn render_3d_msaa(
        &self,
        scene: &Scene,
        msaa_view: &wgpu::TextureView,
        resolve_target: &wgpu::TextureView,
        depth_view: &wgpu::TextureView,
        encoder: &mut wgpu::CommandEncoder,
    ) {
        struct ActorDraw {
            mesh: GpuMesh,
            edge_mesh: Option<GpuMesh>,
            repr: Representation,
            opacity: f32,
            material: Material,
            position: [f64; 3],
            scale: f64,
        }

        let mut opaque_draws = Vec::new();
        let mut translucent_draws = Vec::new();

        let cam_dist = (scene.camera.position - scene.camera.focal_point).length();

        for actor in &scene.actors {
            if !actor.visible {
                continue;
            }
            let Some(gpu_mesh) = self.prepare_actor_mesh(actor, cam_dist) else {
                continue;
            };
            let edge_mesh = if actor.material.edge_visibility
                && actor.representation == Representation::Surface
            {
                self.prepare_edge_mesh(actor)
            } else {
                None
            };
            let draw = ActorDraw {
                mesh: gpu_mesh,
                edge_mesh,
                repr: actor.representation,
                opacity: actor.opacity,
                material: actor.material.clone(),
                position: actor.position,
                scale: actor.scale,
            };
            if actor.opacity < 1.0 {
                translucent_draws.push(draw);
            } else {
                opaque_draws.push(draw);
            }
        }

        {
            let mut pass = encoder.begin_render_pass(&wgpu::RenderPassDescriptor {
                label: Some("render pass"),
                color_attachments: &[Some(wgpu::RenderPassColorAttachment {
                    view: msaa_view,
                    resolve_target: Some(resolve_target),
                    ops: wgpu::Operations {
                        load: wgpu::LoadOp::Clear(wgpu::Color {
                            r: scene.background[0] as f64,
                            g: scene.background[1] as f64,
                            b: scene.background[2] as f64,
                            a: scene.background[3] as f64,
                        }),
                        store: wgpu::StoreOp::Store,
                    },
                })],
                depth_stencil_attachment: Some(wgpu::RenderPassDepthStencilAttachment {
                    view: depth_view,
                    depth_ops: Some(wgpu::Operations {
                        load: wgpu::LoadOp::Clear(1.0),
                        store: wgpu::StoreOp::Store,
                    }),
                    stencil_ops: None,
                }),
                ..Default::default()
            });

            // Opaque pass
            for draw in &opaque_draws {
                let use_lighting = draw.repr == Representation::Surface;
                self.write_uniforms(scene, &draw.material, draw.opacity, use_lighting, draw.position, draw.scale);
                let cull = draw.material.backface_culling;
                pass.set_pipeline(self.select_pipeline(draw.repr, false, cull));
                pass.set_bind_group(0, &self.bind_group, &[]);
                pass.set_vertex_buffer(0, draw.mesh.vertex_buffer.slice(..));
                pass.set_index_buffer(draw.mesh.index_buffer.slice(..), wgpu::IndexFormat::Uint32);
                pass.draw_indexed(0..draw.mesh.num_indices, 0, 0..1);

                if let Some(ref edge) = draw.edge_mesh {
                    self.write_uniforms(scene, &draw.material, 1.0, false, draw.position, draw.scale);
                    pass.set_pipeline(self.select_pipeline(Representation::Wireframe, false, false));
                    pass.set_bind_group(0, &self.bind_group, &[]);
                    pass.set_vertex_buffer(0, edge.vertex_buffer.slice(..));
                    pass.set_index_buffer(edge.index_buffer.slice(..), wgpu::IndexFormat::Uint32);
                    pass.draw_indexed(0..edge.num_indices, 0, 0..1);
                }
            }

            // Translucent pass
            for draw in &translucent_draws {
                let use_lighting = draw.repr == Representation::Surface;
                self.write_uniforms(scene, &draw.material, draw.opacity, use_lighting, draw.position, draw.scale);
                let cull = draw.material.backface_culling;
                pass.set_pipeline(self.select_pipeline(draw.repr, true, cull));
                pass.set_bind_group(0, &self.bind_group, &[]);
                pass.set_vertex_buffer(0, draw.mesh.vertex_buffer.slice(..));
                pass.set_index_buffer(draw.mesh.index_buffer.slice(..), wgpu::IndexFormat::Uint32);
                pass.draw_indexed(0..draw.mesh.num_indices, 0, 0..1);

                if let Some(ref edge) = draw.edge_mesh {
                    self.write_uniforms(scene, &draw.material, draw.opacity, false, draw.position, draw.scale);
                    pass.set_pipeline(self.select_pipeline(Representation::Wireframe, true, false));
                    pass.set_bind_group(0, &self.bind_group, &[]);
                    pass.set_vertex_buffer(0, edge.vertex_buffer.slice(..));
                    pass.set_index_buffer(edge.index_buffer.slice(..), wgpu::IndexFormat::Uint32);
                    pass.draw_indexed(0..edge.num_indices, 0, 0..1);
                }
            }
        }
    }

    fn render_to_view(
        &self,
        scene: &Scene,
        msaa_view: &wgpu::TextureView,
        resolve_target: &wgpu::TextureView,
        depth_view: &wgpu::TextureView,
    ) -> Result<(), VtkError> {
        let mut encoder = self.device.create_command_encoder(&wgpu::CommandEncoderDescriptor {
            label: Some("render encoder"),
        });

        // Skybox (render gradient background before 3D scene)
        if !matches!(scene.skybox, vtk_render::Skybox::Solid(_)) {
            self.skybox_pass.render(&self.queue, &mut encoder, resolve_target, &scene.skybox);
        }

        // 3D scene with MSAA
        self.render_3d_msaa(scene, msaa_view, resolve_target, depth_view, &mut encoder);

        // Silhouette edges (rendered as overlay lines after MSAA resolve)
        if scene.silhouette.enabled {
            let cam_pos = scene.camera.position;
            let cam = [cam_pos.x, cam_pos.y, cam_pos.z];
            let mut all_sil_verts = Vec::new();
            let mut all_sil_idxs: Vec<u32> = Vec::new();
            for actor in &scene.actors {
                if actor.representation != Representation::Surface {
                    continue;
                }
                let (verts, idxs) = crate::silhouette_pass::extract_silhouette_edges(
                    &actor.data,
                    cam,
                    &scene.silhouette,
                );
                let base = all_sil_verts.len() as u32;
                all_sil_verts.extend_from_slice(&verts);
                all_sil_idxs.extend(idxs.iter().map(|&i| i + base));
            }
            if !all_sil_idxs.is_empty() {
                let sil_mesh = upload_mesh(&self.device, &all_sil_verts, &all_sil_idxs);
                let sil_material = Material::default();
                self.write_uniforms(scene, &sil_material, 1.0, false, [0.0, 0.0, 0.0], 1.0);
                let mut pass = encoder.begin_render_pass(&wgpu::RenderPassDescriptor {
                    label: Some("silhouette pass"),
                    color_attachments: &[Some(wgpu::RenderPassColorAttachment {
                        view: resolve_target,
                        resolve_target: None,
                        ops: wgpu::Operations {
                            load: wgpu::LoadOp::Load,
                            store: wgpu::StoreOp::Store,
                        },
                    })],
                    depth_stencil_attachment: None,
                    ..Default::default()
                });
                pass.set_pipeline(&self.pipeline_lines_no_msaa);
                pass.set_bind_group(0, &self.bind_group, &[]);
                pass.set_vertex_buffer(0, sil_mesh.vertex_buffer.slice(..));
                pass.set_index_buffer(sil_mesh.index_buffer.slice(..), wgpu::IndexFormat::Uint32);
                pass.draw_indexed(0..sil_mesh.num_indices, 0, 0..1);
            }
        }

        // Volume rendering pass (after surface, before overlays)
        if !scene.volumes.is_empty() {
            let aspect = self.width as f64 / self.height as f64;
            let view_mat = scene.camera.view_matrix();
            let proj_mat = scene.camera.projection_matrix(aspect);
            let vp = proj_mat * view_mat;
            let mvp = Mat4::from_cols_array(&vp.to_cols_array().map(|v| v as f32));
            let cam = scene.camera.position;
            let cam_pos = [cam.x as f32, cam.y as f32, cam.z as f32];

            for volume in &scene.volumes {
                self.volume_pass.render(
                    &self.device, &self.queue, &mut encoder,
                    resolve_target, volume, mvp, cam_pos,
                );
            }
        }

        // 2D overlay (scalar bars) — rendered after resolve, directly on resolve target
        self.overlay_pipeline.render_scalar_bars(
            &self.device,
            &mut encoder,
            resolve_target,
            &scene.scalar_bars,
        );

        // Axes widget overlay
        if let Some(ref axes) = scene.axes_widget {
            let view_mat = scene.camera.view_matrix();
            let cols = view_mat.to_cols_array();
            let view_arr = [
                [cols[0], cols[1], cols[2], cols[3]],
                [cols[4], cols[5], cols[6], cols[7]],
                [cols[8], cols[9], cols[10], cols[11]],
                [cols[12], cols[13], cols[14], cols[15]],
            ];
            self.overlay_pipeline.render_axes_widget(
                &self.device,
                &mut encoder,
                resolve_target,
                axes,
                &view_arr,
            );
        }

        // Annotations overlay (labels, rulers, protractors)
        if !scene.annotations.is_empty() {
            let aspect = self.width as f64 / self.height as f64;
            let view_mat = scene.camera.view_matrix();
            let proj_mat = scene.camera.projection_matrix(aspect);
            let vp = proj_mat * view_mat;

            let mut vertices = Vec::new();
            let mut indices: Vec<u32> = Vec::new();

            // Project 3D positions to NDC and render as text
            for label in &scene.annotations.labels {
                let ndc = project_to_ndc(&vp, label.position);
                if let Some((nx, ny)) = ndc {
                    let color = [label.color[0], label.color[1], label.color[2], 1.0];
                    crate::bitmap_font::render_text(
                        &mut vertices, &mut indices,
                        &label.text, nx, ny, label.scale, color,
                    );
                }
            }

            for ruler in &scene.annotations.rulers {
                let start_ndc = project_to_ndc(&vp, ruler.start);
                let end_ndc = project_to_ndc(&vp, ruler.end);
                if let (Some((sx, sy)), Some((ex, ey))) = (start_ndc, end_ndc) {
                    let color = [ruler.color[0], ruler.color[1], ruler.color[2], 1.0];
                    // Draw line as thin quad
                    let dx = ex - sx;
                    let dy = ey - sy;
                    let len = (dx * dx + dy * dy).sqrt();
                    if len > 0.001 {
                        let nx = -dy / len * 0.002;
                        let ny = dx / len * 0.002;
                        let base = vertices.len() as u32;
                        vertices.push(crate::overlay::OverlayVertex { position: [sx - nx, sy - ny], color });
                        vertices.push(crate::overlay::OverlayVertex { position: [sx + nx, sy + ny], color });
                        vertices.push(crate::overlay::OverlayVertex { position: [ex + nx, ey + ny], color });
                        vertices.push(crate::overlay::OverlayVertex { position: [ex - nx, ey - ny], color });
                        indices.extend_from_slice(&[base, base+1, base+2, base, base+2, base+3]);
                    }

                    if ruler.show_label {
                        let mid_ndc = project_to_ndc(&vp, ruler.midpoint());
                        if let Some((mx, my)) = mid_ndc {
                            let text = format!("{:.2}", ruler.distance());
                            crate::bitmap_font::render_text(
                                &mut vertices, &mut indices,
                                &text, mx + 0.01, my + 0.01, 0.012, color,
                            );
                        }
                    }
                }
            }

            if !vertices.is_empty() {
                use wgpu::util::DeviceExt;
                let vb = self.device.create_buffer_init(&wgpu::util::BufferInitDescriptor {
                    label: Some("annotation vb"),
                    contents: bytemuck::cast_slice(&vertices),
                    usage: wgpu::BufferUsages::VERTEX,
                });
                let ib = self.device.create_buffer_init(&wgpu::util::BufferInitDescriptor {
                    label: Some("annotation ib"),
                    contents: bytemuck::cast_slice(&indices),
                    usage: wgpu::BufferUsages::INDEX,
                });
                let mut pass = encoder.begin_render_pass(&wgpu::RenderPassDescriptor {
                    label: Some("annotation pass"),
                    color_attachments: &[Some(wgpu::RenderPassColorAttachment {
                        view: resolve_target,
                        resolve_target: None,
                        ops: wgpu::Operations {
                            load: wgpu::LoadOp::Load,
                            store: wgpu::StoreOp::Store,
                        },
                    })],
                    depth_stencil_attachment: None,
                    ..Default::default()
                });
                pass.set_pipeline(&self.overlay_pipeline.pipeline());
                pass.set_vertex_buffer(0, vb.slice(..));
                pass.set_index_buffer(ib.slice(..), wgpu::IndexFormat::Uint32);
                pass.draw_indexed(0..indices.len() as u32, 0, 0..1);
            }
        }

        self.queue.submit(std::iter::once(encoder.finish()));
        Ok(())
    }
}

impl Renderer for WgpuRenderer {
    fn render(&mut self, scene: &Scene) -> Result<(), VtkError> {
        let output = self
            .surface
            .get_current_texture()
            .map_err(|e| VtkError::InvalidData(format!("surface texture error: {e}")))?;
        let resolve_view = output.texture.create_view(&Default::default());
        self.render_to_view(scene, &self.msaa_texture, &resolve_view, &self.depth_texture)?;
        output.present();
        Ok(())
    }

    fn resize(&mut self, width: u32, height: u32) {
        if width > 0 && height > 0 {
            self.width = width;
            self.height = height;
            self.surface_config.width = width;
            self.surface_config.height = height;
            self.surface.configure(&self.device, &self.surface_config);
            self.depth_texture = create_depth_texture(&self.device, width, height, MSAA_SAMPLE_COUNT);
            self.msaa_texture = create_msaa_texture(&self.device, self.surface_format, width, height);
        }
    }

    fn render_to_image(
        &mut self,
        scene: &Scene,
        width: u32,
        height: u32,
    ) -> Result<Vec<u8>, VtkError> {
        let saved_w = self.width;
        let saved_h = self.height;
        self.width = width;
        self.height = height;

        let color_texture = self.device.create_texture(&wgpu::TextureDescriptor {
            label: Some("offscreen color"),
            size: wgpu::Extent3d { width, height, depth_or_array_layers: 1 },
            mip_level_count: 1,
            sample_count: 1,
            dimension: wgpu::TextureDimension::D2,
            format: self.surface_format,
            usage: wgpu::TextureUsages::RENDER_ATTACHMENT | wgpu::TextureUsages::COPY_SRC,
            view_formats: &[],
        });
        let color_view = color_texture.create_view(&Default::default());
        let depth_view = create_depth_texture(&self.device, width, height, MSAA_SAMPLE_COUNT);
        let msaa_view = create_msaa_texture(&self.device, self.surface_format, width, height);

        self.render_to_view(scene, &msaa_view, &color_view, &depth_view)?;

        let bytes_per_pixel = 4u32;
        let unpadded_bytes_per_row = width * bytes_per_pixel;
        let padded_bytes_per_row = (unpadded_bytes_per_row + 255) & !255;

        let staging_buffer = self.device.create_buffer(&wgpu::BufferDescriptor {
            label: Some("offscreen staging"),
            size: (padded_bytes_per_row * height) as u64,
            usage: wgpu::BufferUsages::COPY_DST | wgpu::BufferUsages::MAP_READ,
            mapped_at_creation: false,
        });

        let mut encoder = self.device.create_command_encoder(&wgpu::CommandEncoderDescriptor {
            label: Some("copy encoder"),
        });
        encoder.copy_texture_to_buffer(
            wgpu::TexelCopyTextureInfo {
                texture: &color_texture,
                mip_level: 0,
                origin: wgpu::Origin3d::ZERO,
                aspect: wgpu::TextureAspect::All,
            },
            wgpu::TexelCopyBufferInfo {
                buffer: &staging_buffer,
                layout: wgpu::TexelCopyBufferLayout {
                    offset: 0,
                    bytes_per_row: Some(padded_bytes_per_row),
                    rows_per_image: Some(height),
                },
            },
            wgpu::Extent3d { width, height, depth_or_array_layers: 1 },
        );
        self.queue.submit(std::iter::once(encoder.finish()));

        let buffer_slice = staging_buffer.slice(..);
        let (sender, receiver) = std::sync::mpsc::channel();
        buffer_slice.map_async(wgpu::MapMode::Read, move |result| {
            let _ = sender.send(result);
        });
        self.device.poll(wgpu::Maintain::Wait);

        receiver
            .recv()
            .map_err(|_| VtkError::InvalidData("buffer map recv failed".into()))?
            .map_err(|e| VtkError::InvalidData(format!("buffer map failed: {e}")))?;

        let data = buffer_slice.get_mapped_range();
        let mut result = Vec::with_capacity((width * height * bytes_per_pixel) as usize);
        for row in 0..height {
            let start = (row * padded_bytes_per_row) as usize;
            let end = start + unpadded_bytes_per_row as usize;
            result.extend_from_slice(&data[start..end]);
        }
        drop(data);
        staging_buffer.unmap();

        self.width = saved_w;
        self.height = saved_h;
        Ok(result)
    }
}

fn light_to_gpu(light: &Light) -> GpuLight {
    let (lt, cone_angle, exponent) = match light.light_type {
        LightType::Directional => (0.0, 0.0, 0.0),
        LightType::Point => (1.0, 0.0, 0.0),
        LightType::Spot { cone_angle, exponent } => (2.0, cone_angle as f32, exponent as f32),
        LightType::Ambient => (3.0, 0.0, 0.0),
    };
    GpuLight {
        light_type: lt,
        intensity: light.intensity as f32,
        cone_angle,
        exponent,
        position: [light.position[0] as f32, light.position[1] as f32, light.position[2] as f32],
        _pad0: 0.0,
        direction: [light.direction[0] as f32, light.direction[1] as f32, light.direction[2] as f32],
        _pad1: 0.0,
        color: light.color,
        _pad2: 0.0,
    }
}

fn upload_mesh(device: &wgpu::Device, vertices: &[Vertex], indices: &[u32]) -> GpuMesh {
    let vertex_buffer = device.create_buffer_init(&wgpu::util::BufferInitDescriptor {
        label: Some("vertex buffer"),
        contents: bytemuck::cast_slice(vertices),
        usage: wgpu::BufferUsages::VERTEX,
    });
    let index_buffer = device.create_buffer_init(&wgpu::util::BufferInitDescriptor {
        label: Some("index buffer"),
        contents: bytemuck::cast_slice(indices),
        usage: wgpu::BufferUsages::INDEX,
    });
    GpuMesh {
        vertex_buffer,
        index_buffer,
        num_indices: indices.len() as u32,
    }
}

fn create_pipeline_with_cull(
    device: &wgpu::Device,
    layout: &wgpu::PipelineLayout,
    shader: &wgpu::ShaderModule,
    format: wgpu::TextureFormat,
    topology: wgpu::PrimitiveTopology,
    alpha_blend: bool,
    sample_count: u32,
    cull_mode: Option<wgpu::Face>,
) -> wgpu::RenderPipeline {
    let blend = if alpha_blend {
        Some(wgpu::BlendState::ALPHA_BLENDING)
    } else {
        Some(wgpu::BlendState::REPLACE)
    };
    device.create_render_pipeline(&wgpu::RenderPipelineDescriptor {
        label: Some("render pipeline"),
        layout: Some(layout),
        vertex: wgpu::VertexState {
            module: shader,
            entry_point: Some("vs_main"),
            buffers: &[Vertex::layout()],
            compilation_options: Default::default(),
        },
        fragment: Some(wgpu::FragmentState {
            module: shader,
            entry_point: Some("fs_main"),
            targets: &[Some(wgpu::ColorTargetState { format, blend, write_mask: wgpu::ColorWrites::ALL })],
            compilation_options: Default::default(),
        }),
        primitive: wgpu::PrimitiveState {
            topology,
            front_face: wgpu::FrontFace::Ccw,
            cull_mode,
            ..Default::default()
        },
        depth_stencil: Some(wgpu::DepthStencilState {
            format: wgpu::TextureFormat::Depth32Float,
            depth_write_enabled: !alpha_blend,
            depth_compare: wgpu::CompareFunction::Less,
            stencil: Default::default(),
            bias: Default::default(),
        }),
        multisample: wgpu::MultisampleState {
            count: sample_count,
            mask: !0,
            alpha_to_coverage_enabled: false,
        },
        multiview: None,
        cache: None,
    })
}

fn create_depth_texture(device: &wgpu::Device, width: u32, height: u32, sample_count: u32) -> wgpu::TextureView {
    let texture = device.create_texture(&wgpu::TextureDescriptor {
        label: Some("depth texture"),
        size: wgpu::Extent3d { width, height, depth_or_array_layers: 1 },
        mip_level_count: 1,
        sample_count,
        dimension: wgpu::TextureDimension::D2,
        format: wgpu::TextureFormat::Depth32Float,
        usage: wgpu::TextureUsages::RENDER_ATTACHMENT,
        view_formats: &[],
    });
    texture.create_view(&Default::default())
}

fn create_msaa_texture(device: &wgpu::Device, format: wgpu::TextureFormat, width: u32, height: u32) -> wgpu::TextureView {
    let texture = device.create_texture(&wgpu::TextureDescriptor {
        label: Some("msaa texture"),
        size: wgpu::Extent3d { width, height, depth_or_array_layers: 1 },
        mip_level_count: 1,
        sample_count: MSAA_SAMPLE_COUNT,
        dimension: wgpu::TextureDimension::D2,
        format,
        usage: wgpu::TextureUsages::RENDER_ATTACHMENT,
        view_formats: &[],
    });
    texture.create_view(&Default::default())
}

fn poly_data_to_points(
    poly_data: &vtk_data::PolyData,
    coloring: &vtk_render::Coloring,
) -> (Vec<Vertex>, Vec<u32>) {
    let point_colors = mesh::resolve_colors_pub(poly_data, coloring);
    let n = poly_data.points.len();
    let mut vertices = Vec::with_capacity(n);
    let mut indices = Vec::with_capacity(n);
    for i in 0..n {
        let p = poly_data.points.get(i);
        vertices.push(Vertex {
            position: [p[0] as f32, p[1] as f32, p[2] as f32],
            normal: [0.0, 0.0, 1.0],
            color: point_colors[i],
        });
        indices.push(i as u32);
    }
    (vertices, indices)
}

/// Project a 3D world position to NDC [0,1] screen coordinates.
/// Returns None if behind camera.
fn project_to_ndc(vp: &glam::DMat4, pos: [f64; 3]) -> Option<(f32, f32)> {
    let p = *vp * glam::DVec4::new(pos[0], pos[1], pos[2], 1.0);
    if p.w.abs() < 1e-10 || p.z / p.w < -1.0 { return None; }
    let nx = ((p.x / p.w + 1.0) * 0.5) as f32;
    let ny = ((p.y / p.w + 1.0) * 0.5) as f32;
    Some((nx, ny))
}
