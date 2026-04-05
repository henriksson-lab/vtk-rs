//! Mesh info tool: read any supported mesh format and print comprehensive analysis.
//!
//! Usage: cargo run --example mesh_info -- path/to/mesh.stl

use vtk_data::DataSet;

fn main() {
    let args: Vec<String> = std::env::args().collect();

    if args.len() < 2 {
        // Demo mode: generate a sphere and analyze it
        println!("vtk-rs mesh info tool");
        println!("Usage: cargo run --example mesh_info -- <file.vtk|.vtp|.stl|.obj|.ply|.glb>");
        println!();
        println!("Demo mode: analyzing generated sphere...");
        println!();

        let pd = vtk_filters::sources::sphere::sphere(
            &vtk_filters::sources::sphere::SphereParams::default(),
        );
        let pd = vtk_filters::normals::compute_normals(&pd);
        let pd = vtk_filters::elevation::elevation_z(&pd);
        analyze(&pd);
        return;
    }

    let path = std::path::Path::new(&args[1]);
    match vtk::io::read_poly_data(path) {
        Ok(pd) => {
            println!("File: {}", path.display());
            println!();
            analyze(&pd);
        }
        Err(e) => {
            eprintln!("Error reading {}: {e}", path.display());
            std::process::exit(1);
        }
    }
}

fn analyze(pd: &vtk_data::PolyData) {
    println!("=== Geometry ===");
    println!("  {pd}");
    println!("  Triangles: {} / {} polys", pd.num_triangles(), pd.num_polys());
    println!("  All triangles: {}", pd.is_all_triangles());
    println!("  Edges: {}", pd.num_edges());
    println!();

    println!("=== Bounds ===");
    let bb = pd.bounds();
    println!("  X: [{:.4}, {:.4}] (size: {:.4})", bb.x_min, bb.x_max, bb.size()[0]);
    println!("  Y: [{:.4}, {:.4}] (size: {:.4})", bb.y_min, bb.y_max, bb.size()[1]);
    println!("  Z: [{:.4}, {:.4}] (size: {:.4})", bb.z_min, bb.z_max, bb.size()[2]);
    println!("  Center: {:?}", bb.center());
    println!("  Diagonal: {:.4}", bb.diagonal_length());
    println!();

    println!("=== Topology ===");
    let topo = vtk_filters::topology::analyze_topology(pd);
    println!("  Vertices: {}", topo.num_points);
    println!("  Edges: {}", topo.num_edges);
    println!("  Faces: {}", topo.num_faces);
    println!("  Boundary edges: {}", topo.num_boundary_edges);
    println!("  Non-manifold edges: {}", topo.num_non_manifold_edges);
    println!("  Euler characteristic: {}", topo.euler_characteristic);
    println!("  Components: {}", topo.num_components);
    println!("  Manifold: {}", topo.is_manifold);
    println!("  Closed: {}", topo.is_closed);
    if let Some(g) = topo.genus {
        println!("  Genus: {g}");
    }
    println!();

    println!("=== Measurements ===");
    let m = vtk_render::measurement::measure(pd);
    println!("  Surface area: {:.6}", m.surface_area);
    println!("  Edge lengths: min={:.4}, max={:.4}, avg={:.4}", m.min_edge_length, m.max_edge_length, m.avg_edge_length);
    println!("  Total edge length: {:.4}", m.total_edge_length);
    println!();

    println!("=== Point Data ===");
    if pd.point_data().num_arrays() == 0 {
        println!("  (none)");
    }
    for i in 0..pd.point_data().num_arrays() {
        if let Some(arr) = pd.point_data().get_array_by_index(i) {
            print!("  {arr}");
            if let Some(stats) = arr.statistics() {
                print!(" range=[{:.4}, {:.4}] mean={:.4} std={:.4}", stats.min, stats.max, stats.mean, stats.std_dev());
            }
            println!();
        }
    }
    println!();

    println!("=== Cell Data ===");
    if pd.cell_data().num_arrays() == 0 {
        println!("  (none)");
    }
    for i in 0..pd.cell_data().num_arrays() {
        if let Some(arr) = pd.cell_data().get_array_by_index(i) {
            println!("  {arr}");
        }
    }
}
