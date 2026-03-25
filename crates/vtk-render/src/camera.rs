use glam::{DMat4, DVec3};

/// Virtual camera for 3D viewing.
#[derive(Debug, Clone)]
pub struct Camera {
    pub position: DVec3,
    pub focal_point: DVec3,
    pub view_up: DVec3,
    pub fov: f64,
    pub near_clip: f64,
    pub far_clip: f64,
}

impl Default for Camera {
    fn default() -> Self {
        Self {
            position: DVec3::new(0.0, 0.0, 1.0),
            focal_point: DVec3::ZERO,
            view_up: DVec3::Y,
            fov: 30.0,
            near_clip: 0.01,
            far_clip: 1000.0,
        }
    }
}

impl Camera {
    pub fn new() -> Self {
        Self::default()
    }

    /// View direction (normalized).
    pub fn direction(&self) -> DVec3 {
        (self.focal_point - self.position).normalize()
    }

    /// View matrix (world → camera).
    pub fn view_matrix(&self) -> DMat4 {
        DMat4::look_at_rh(self.position, self.focal_point, self.view_up)
    }

    /// Perspective projection matrix.
    pub fn projection_matrix(&self, aspect_ratio: f64) -> DMat4 {
        let fov_rad = self.fov.to_radians();
        DMat4::perspective_rh(fov_rad, aspect_ratio, self.near_clip, self.far_clip)
    }

    /// Position camera to view the given bounding box.
    pub fn reset_to_bounds(&mut self, center: [f64; 3], diagonal: f64) {
        let c = DVec3::from_array(center);
        let distance = diagonal / (2.0 * (self.fov.to_radians() / 2.0).tan());
        self.focal_point = c;
        self.position = c + DVec3::new(0.0, 0.0, distance);
        self.near_clip = distance * 0.01;
        self.far_clip = distance * 100.0;
    }

    /// Orbit around the focal point by azimuth (degrees) and elevation (degrees).
    pub fn orbit(&mut self, azimuth_deg: f64, elevation_deg: f64) {
        let offset = self.position - self.focal_point;
        let distance = offset.length();

        // Current spherical angles
        let theta = offset.z.atan2(offset.x) + azimuth_deg.to_radians();
        let phi = (offset.y / distance).asin() + elevation_deg.to_radians();
        let phi = phi.clamp(-std::f64::consts::FRAC_PI_2 + 0.01, std::f64::consts::FRAC_PI_2 - 0.01);

        self.position = self.focal_point
            + DVec3::new(
                distance * phi.cos() * theta.cos(),
                distance * phi.sin(),
                distance * phi.cos() * theta.sin(),
            );
    }

    /// Dolly (move along view direction).
    pub fn dolly(&mut self, factor: f64) {
        let dir = self.direction();
        let distance = (self.focal_point - self.position).length();
        let delta = dir * distance * (1.0 - 1.0 / factor);
        self.position += delta;
    }
}
