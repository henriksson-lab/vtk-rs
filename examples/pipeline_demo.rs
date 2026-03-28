//! Pipeline demo: demonstrates the filter pipeline, I/O, and offline processing.
//!
//! This example:
//! 1. Generates a sphere source
//! 2. Builds a processing pipeline (smooth → normals → elevation)
//! 3. Writes the result to multiple formats
//! 4. Prints mesh statistics

use vtk_data::{AnyDataArray, DataArray, DataSet, PolyData};
use vtk_filters::pipeline::Pipeline;
use vtk_filters::sources::sphere::{sphere, SphereParams};
use vtk_filters::topology::analyze_topology;

fn main() {
    println!("vtk-rs Pipeline Demo");
    println!("====================\n");

    // 1. Generate a sphere
    let src = sphere(&SphereParams {
        theta_resolution: 32,
        phi_resolution: 32,
        ..Default::default()
    });
    println!("Source: {src}");

    // 2. Build a processing pipeline
    let mut pipe = Pipeline::new(src)
        .with_normals()
        .with_elevation_z()
        .with_decimate(0.5);

    println!("Pipeline stages: {:?}", pipe.stage_names());

    // 3. Get the output (lazy evaluation)
    let result = pipe.output();
    println!("Output: {result}");

    // 4. Topology analysis
    let topo = analyze_topology(result);
    println!("\nTopology:");
    println!("  Points:        {}", topo.num_points);
    println!("  Edges:         {}", topo.num_edges);
    println!("  Faces:         {}", topo.num_faces);
    println!("  Boundary edges:{}", topo.num_boundary_edges);
    println!("  Euler:         {}", topo.euler_characteristic);
    println!("  Components:    {}", topo.num_components);
    println!("  Manifold:      {}", topo.is_manifold);
    println!("  Triangle mesh: {}", topo.is_triangle_mesh);
    if let Some(g) = topo.genus {
        println!("  Genus:         {g}");
    }

    // 5. Data array statistics
    if let Some(scalars) = result.point_data().scalars() {
        if let Some(stats) = scalars.statistics() {
            println!("\nScalar statistics ({}):", scalars.name());
            println!("  Range: [{:.3}, {:.3}]", stats.min, stats.max);
            println!("  Mean:  {:.3}", stats.mean);
            println!("  Std:   {:.3}", stats.std_dev());
        }
    }

    // 6. Write to multiple formats
    let dir = std::env::temp_dir().join("vtk_pipeline_demo");
    let _ = std::fs::create_dir_all(&dir);

    let formats = ["vtk", "vtp", "stl", "obj", "ply", "glb"];
    println!("\nWriting to:");
    for ext in &formats {
        let path = dir.join(format!("sphere.{ext}"));
        match vtk_filters::io_utils::write_poly_data(&path, result) {
            Ok(()) => println!("  {} ✓", path.display()),
            Err(e) => println!("  {} ✗ {e}", path.display()),
        }
    }

    // 7. Verify roundtrip
    let vtk_path = dir.join("sphere.vtk");
    match vtk_filters::io_utils::read_poly_data(&vtk_path) {
        Ok(loaded) => {
            println!("\nRoundtrip verification:");
            println!("  Original: {} points", result.points.len());
            println!("  Loaded:   {} points", loaded.points.len());
            println!("  Match:    {}", result.approx_eq(&loaded, 1e-6));
        }
        Err(e) => println!("\nRoundtrip failed: {e}"),
    }

    let _ = std::fs::remove_dir_all(&dir);
    println!("\nDone.");
}
