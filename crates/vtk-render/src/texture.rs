/// A 2D texture for mapping onto surfaces.
///
/// Stores RGBA pixel data. Used with `Coloring::Texture` to apply
/// image textures to actors that have UV (texture coordinate) data.
#[derive(Debug, Clone)]
pub struct Texture {
    /// RGBA pixel data (4 bytes per pixel, row-major, top-to-bottom).
    pub data: Vec<u8>,
    /// Width in pixels.
    pub width: u32,
    /// Height in pixels.
    pub height: u32,
}

impl Texture {
    /// Create a texture from RGBA data.
    pub fn from_rgba(data: Vec<u8>, width: u32, height: u32) -> Self {
        assert_eq!(data.len(), (width * height * 4) as usize);
        Self { data, width, height }
    }

    /// Create a 1x1 solid color texture.
    pub fn solid(r: u8, g: u8, b: u8) -> Self {
        Self {
            data: vec![r, g, b, 255],
            width: 1,
            height: 1,
        }
    }

    /// Create a checkerboard texture for testing.
    pub fn checkerboard(size: u32, tile_size: u32) -> Self {
        let mut data = Vec::with_capacity((size * size * 4) as usize);
        for y in 0..size {
            for x in 0..size {
                let is_white = ((x / tile_size) + (y / tile_size)) % 2 == 0;
                let v = if is_white { 220u8 } else { 40u8 };
                data.extend_from_slice(&[v, v, v, 255]);
            }
        }
        Self { data, width: size, height: size }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn solid_texture() {
        let t = Texture::solid(255, 0, 0);
        assert_eq!(t.width, 1);
        assert_eq!(t.height, 1);
        assert_eq!(t.data, vec![255, 0, 0, 255]);
    }

    #[test]
    fn checkerboard_texture() {
        let t = Texture::checkerboard(8, 4);
        assert_eq!(t.width, 8);
        assert_eq!(t.height, 8);
        assert_eq!(t.data.len(), 8 * 8 * 4);
    }

    #[test]
    fn from_rgba() {
        let data = vec![0u8; 4 * 2 * 2];
        let t = Texture::from_rgba(data, 2, 2);
        assert_eq!(t.width, 2);
    }
}
