//! GPU compute filters for vtk-rs.
//!
//! Uses wgpu compute shaders for data-parallel operations on large datasets.
//! Provides a `GpuContext` for device management and compute filters that
//! operate on `DataArray` via GPU storage buffers.
//!
//! # Example
//!
//! ```no_run
//! use vtk_filters_gpu::{GpuContext, gpu_scale_array};
//! use vtk_data::DataArray;
//!
//! let ctx = GpuContext::new().unwrap();
//! let input = DataArray::from_vec("data", vec![1.0f32, 2.0, 3.0, 4.0], 1);
//! let output = gpu_scale_array(&ctx, &input, 2.0);
//! ```

mod context;
mod compute_scale;
mod compute_map;
mod compute_reduce;
mod compute_normals;

pub use context::GpuContext;
pub use compute_scale::gpu_scale_array;
pub use compute_map::gpu_map_array;
pub use compute_reduce::gpu_reduce_sum;
pub use compute_normals::gpu_compute_normals;
