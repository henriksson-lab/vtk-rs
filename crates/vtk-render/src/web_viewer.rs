//! Web viewer scaffolding for embedding scene data in a standalone HTML page.
//!
//! Generates an HTML page with a `<canvas>` element and a JavaScript stub
//! that could be wired up to a WASM + WebGPU renderer.  The scene
//! configuration is embedded as JSON inside the page.

use crate::Scene;

/// Configuration for the web viewer page.
#[derive(Debug, Clone)]
pub struct WebViewerConfig {
    /// TCP port (informational — no server is started).
    pub port: u16,
    /// Page title.
    pub title: String,
    /// Canvas width in CSS pixels.
    pub width: u32,
    /// Canvas height in CSS pixels.
    pub height: u32,
}

impl Default for WebViewerConfig {
    fn default() -> Self {
        Self {
            port: 8080,
            title: "vtk-rs Web Viewer".to_string(),
            width: 800,
            height: 600,
        }
    }
}

/// Generate a standalone HTML page with embedded scene JSON and a canvas.
///
/// The `scene_json` string is placed verbatim into a `<script>` block as a
/// JavaScript variable.
pub fn generate_html(scene_json: &str) -> String {
    let config = WebViewerConfig::default();
    generate_html_with_config(scene_json, &config)
}

/// Generate HTML using a custom [`WebViewerConfig`].
pub fn generate_html_with_config(scene_json: &str, config: &WebViewerConfig) -> String {
    format!(
        r#"<!DOCTYPE html>
<html>
<head>
<meta charset="utf-8">
<title>{title}</title>
<style>
body {{ margin: 0; background: #222; display: flex; justify-content: center; align-items: center; height: 100vh; }}
canvas {{ border: 1px solid #555; }}
</style>
</head>
<body>
<canvas id="vtk-canvas" width="{w}" height="{h}"></canvas>
<script>
const SCENE_DATA = {json};

// WebGPU initialization stub
async function initWebGPU() {{
    if (!navigator.gpu) {{
        console.error("WebGPU not supported");
        return;
    }}
    const adapter = await navigator.gpu.requestAdapter();
    const device = await adapter.requestDevice();
    const canvas = document.getElementById("vtk-canvas");
    const context = canvas.getContext("webgpu");
    // TODO: configure swap chain and render pipeline using SCENE_DATA
    console.log("WebGPU device acquired, scene has", SCENE_DATA.num_actors, "actors");
}}
initWebGPU();
</script>
</body>
</html>"#,
        title = config.title,
        w = config.width,
        h = config.height,
        json = scene_json,
    )
}

/// Generate an HTML page with the scene embedded as JSON.
///
/// Uses [`scene_to_json_string`](crate::scene_json::scene_to_json_string) to
/// serialise the scene configuration.
pub fn scene_to_embedded_html(scene: &Scene) -> String {
    let json = crate::scene_json::scene_to_json_string(scene);
    generate_html(&json)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::Scene;

    #[test]
    fn html_contains_canvas() {
        let html = generate_html(r#"{"num_actors": 0}"#);
        assert!(html.contains("<canvas"));
        assert!(html.contains("vtk-canvas"));
        assert!(html.contains("WebGPU"));
    }

    #[test]
    fn scene_json_embedded() {
        let scene = Scene::new();
        let html = scene_to_embedded_html(&scene);
        assert!(html.contains("<canvas"));
        // The JSON should contain camera info
        assert!(html.contains("camera"));
        assert!(html.contains("SCENE_DATA"));
    }
}
