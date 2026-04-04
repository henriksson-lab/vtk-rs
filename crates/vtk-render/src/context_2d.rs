//! 2D rendering context for overlay drawing.
//!
//! Accumulates 2D drawing commands (lines, rectangles, circles, text) as
//! vertex/index buffers suitable for the overlay pipeline.

/// A vertex in the 2D context: position (x, y) and RGBA color.
pub type Vertex2D = (f32, f32, [f32; 4]);

/// 2D drawing context that generates triangles for overlay rendering.
///
/// Accumulates geometry from drawing commands and exposes raw vertex/index
/// buffers for GPU submission.
#[derive(Debug, Clone)]
pub struct Context2D {
    verts: Vec<Vertex2D>,
    idxs: Vec<u32>,
}

impl Default for Context2D {
    fn default() -> Self {
        Self::new()
    }
}

impl Context2D {
    /// Create an empty 2D context.
    pub fn new() -> Self {
        Self {
            verts: Vec::new(),
            idxs: Vec::new(),
        }
    }

    /// Draw a line from (x0, y0) to (x1, y1) with given color and width.
    ///
    /// The line is rendered as a thin quad (two triangles).
    pub fn line(&mut self, x0: f32, y0: f32, x1: f32, y1: f32, color: [f32; 4], width: f32) {
        let dx = x1 - x0;
        let dy = y1 - y0;
        let len = (dx * dx + dy * dy).sqrt();
        if len < 1e-8 {
            return;
        }
        let half = width * 0.5;
        // Perpendicular direction
        let nx = -dy / len * half;
        let ny = dx / len * half;

        let base = self.verts.len() as u32;
        self.verts.push((x0 + nx, y0 + ny, color));
        self.verts.push((x0 - nx, y0 - ny, color));
        self.verts.push((x1 - nx, y1 - ny, color));
        self.verts.push((x1 + nx, y1 + ny, color));

        self.idxs.extend_from_slice(&[base, base + 1, base + 2]);
        self.idxs.extend_from_slice(&[base, base + 2, base + 3]);
    }

    /// Draw a rectangle outline.
    pub fn rect(&mut self, x: f32, y: f32, w: f32, h: f32, color: [f32; 4]) {
        let lw = 1.0;
        self.line(x, y, x + w, y, color, lw);
        self.line(x + w, y, x + w, y + h, color, lw);
        self.line(x + w, y + h, x, y + h, color, lw);
        self.line(x, y + h, x, y, color, lw);
    }

    /// Draw a filled rectangle.
    pub fn filled_rect(&mut self, x: f32, y: f32, w: f32, h: f32, color: [f32; 4]) {
        let base = self.verts.len() as u32;
        self.verts.push((x, y, color));
        self.verts.push((x + w, y, color));
        self.verts.push((x + w, y + h, color));
        self.verts.push((x, y + h, color));

        self.idxs.extend_from_slice(&[base, base + 1, base + 2]);
        self.idxs.extend_from_slice(&[base, base + 2, base + 3]);
    }

    /// Draw a circle outline with the given number of segments.
    pub fn circle(&mut self, cx: f32, cy: f32, r: f32, color: [f32; 4], segments: u32) {
        let segs = segments.max(3);
        for i in 0..segs {
            let a0 = 2.0 * std::f32::consts::PI * i as f32 / segs as f32;
            let a1 = 2.0 * std::f32::consts::PI * ((i + 1) % segs) as f32 / segs as f32;
            let x0 = cx + r * a0.cos();
            let y0 = cy + r * a0.sin();
            let x1 = cx + r * a1.cos();
            let y1 = cy + r * a1.sin();
            self.line(x0, y0, x1, y1, color, 1.0);
        }
    }

    /// Draw a filled circle with the given number of segments.
    pub fn filled_circle(&mut self, cx: f32, cy: f32, r: f32, color: [f32; 4], segments: u32) {
        let segs = segments.max(3);
        let center_idx = self.verts.len() as u32;
        self.verts.push((cx, cy, color));

        for i in 0..segs {
            let a = 2.0 * std::f32::consts::PI * i as f32 / segs as f32;
            self.verts.push((cx + r * a.cos(), cy + r * a.sin(), color));
        }

        for i in 0..segs {
            let i0 = center_idx + 1 + i;
            let i1 = center_idx + 1 + (i + 1) % segs;
            self.idxs.extend_from_slice(&[center_idx, i0, i1]);
        }
    }

    /// Draw text at (x, y) with given scale and color.
    ///
    /// Uses a simple bitmap font approach — each character is a small filled rect.
    /// This is a placeholder that generates quads at character positions.
    pub fn text(&mut self, x: f32, y: f32, text: &str, scale: f32, color: [f32; 4]) {
        let char_w = 8.0 * scale;
        let char_h = 12.0 * scale;
        for (i, _ch) in text.chars().enumerate() {
            let cx = x + i as f32 * char_w;
            self.filled_rect(cx, y, char_w * 0.8, char_h, color);
        }
    }

    /// Get the accumulated vertices.
    pub fn vertices(&self) -> &[Vertex2D] {
        &self.verts
    }

    /// Get the accumulated indices.
    pub fn indices(&self) -> &[u32] {
        &self.idxs
    }

    /// Clear all accumulated geometry.
    pub fn clear(&mut self) {
        self.verts.clear();
        self.idxs.clear();
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_filled_rect() {
        let mut ctx = Context2D::new();
        ctx.filled_rect(0.0, 0.0, 10.0, 5.0, [1.0, 0.0, 0.0, 1.0]);
        assert_eq!(ctx.vertices().len(), 4);
        assert_eq!(ctx.indices().len(), 6);
    }

    #[test]
    fn test_line() {
        let mut ctx = Context2D::new();
        ctx.line(0.0, 0.0, 10.0, 0.0, [1.0, 1.0, 1.0, 1.0], 2.0);
        assert_eq!(ctx.vertices().len(), 4);
        assert_eq!(ctx.indices().len(), 6);
    }

    #[test]
    fn test_clear() {
        let mut ctx = Context2D::new();
        ctx.filled_rect(0.0, 0.0, 10.0, 5.0, [1.0, 0.0, 0.0, 1.0]);
        ctx.filled_circle(5.0, 5.0, 3.0, [0.0, 1.0, 0.0, 1.0], 16);
        assert!(ctx.vertices().len() > 0);
        ctx.clear();
        assert_eq!(ctx.vertices().len(), 0);
        assert_eq!(ctx.indices().len(), 0);
    }
}
