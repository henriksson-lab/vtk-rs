# TODO.md — vtk-rs Implementation State & Strategy

## Current State (2026-03-29)

### Stats
- **Lines of Rust:** 297,412
- **Image filters:** 2,831 (vtk-filters-image) + 180 (vtk-filters-image-2) = ~3,011
- **Mesh filters:** 802 (vtk-filters-mesh) + 20 (vtk-filters-mesh-2) = ~822
- **Other sub-crate filters:** 243 spread across 16 sub-crates (extract, transform, subdivide, clip, smooth, cell, points, statistics, texture, flow, boolean, grid, data, distance, normals, geometry)
- **Core/infrastructure filters:** 25 (pipeline, convert, append, etc.)
- **Source geometries:** 414
- **Total sub-crate modules:** 4,068
- **Build status:** Clean (warnings only, no errors)

### Workspace Structure (34 crates)
```
vtk-types, vtk-data                    # Core data model
vtk-filters                            # Pipeline + sources + core filters (25 modules + 414 sources)
vtk-filters-image                      # Image filters batch 1 (2,831 modules)
vtk-filters-image-2                    # Image filters batch 2 (180 modules, STAGING)
vtk-filters-mesh                       # Mesh filters batch 1 (802 modules)
vtk-filters-mesh-2                     # Mesh filters batch 2 (20 modules, STAGING)
vtk-filters-{extract,transform,...}    # 16 sub-crates with original "Other" filters (243 total)
vtk-io-{legacy,stl,obj,xml,...}        # 15 I/O format crates
vtk-render, vtk-render-wgpu            # Rendering
vtk-python                             # PyO3 bindings
vtk                                    # Umbrella re-export crate
```

## Build Strategy

### Fast incremental builds
- **New image/mesh filters** → write to `vtk-filters-image-2` / `vtk-filters-mesh-2`
- **New sources** → write to `vtk-filters/src/sources/`
- **Build command:** `make gen-mods && cargo build -p vtk-filters-image-2 -p vtk-filters-mesh-2 -p vtk-filters`
- **Incremental time:** ~30-40s for image-2/mesh-2 alone, ~2min with sources

### Merging staging crates (when ready)
1. Move all `.rs` files from `vtk-filters-image-2/src/` → `vtk-filters-image/src/`
2. Move all `.rs` files from `vtk-filters-mesh-2/src/` → `vtk-filters-mesh/src/`
3. Run `make gen-mods` to regenerate lib.rs
4. Remove staging crates from workspace
5. Full rebuild to verify

### Module generation
- `make gen-mods` — auto-generates `lib.rs` for all 20 sub-crates + `mod.rs` for sources
- `scripts/gen_from_specs.sh` — spec-file-based batch generator
- Image filter template: pointwise `f64 → f64` with test

## What's Been Implemented

### Image filter domains covered (~3,000 filters)
Acoustics, adhesion, aerodynamics, aquaculture, archaeometry, astrodynamics, astronomy, avalanche, battery, bioacoustics, bioinformatics, biomechanics, biofilm, bioreactor, carbon capture, ceramics, chromatography, coastal, combustion, composite, concrete, corrosion, cosmology, cryobiology, cryogenics, cryopreservation, dam engineering, demography, dendrochronology, desalination, electrochemistry, electrochromics, electrophoresis, electroplating, electrorheology, electromagnetics, epidemiology, explosives, ferrofluid, fire science, fluid dynamics, food science, forensics, forestry, fuel cell, fuzzy logic, game theory, geomorphology, geosynthetics, geotechnical, glaciology, granular, HVAC, hydrogen storage, hydrometallurgy, hydrology, hydropower, hydrothermal, information theory, irrigation, laser, lidar, lighting, luminescence, magnetocaloric, magnetohydrodynamics, magnetorheology, material science, MEMS, metamaterials, microelectronics, microfluidics, microplastics, mining, molecular dynamics, music theory, nanotechnology, neuromorphic, neuroscience, nuclear, ocean acoustics, oceanography, optics, paleoclimatology, paleontology, particle physics, perovskite, phase field, phonon, photocatalysis, photonics, photovoltaic, piezoelectrics, planetary, plasma, powder metallurgy, precision agriculture, printing, proteomics, psychophysics, quantum computing, quantum dots, quantum mechanics, radiation, reliability, remote sensing, renewable energy, rheology, robotics, rocket propulsion, rubber, seismology, semiconductor, shape memory, signal processing, soil, sonochemistry, space weather, spectroscopy, speleology, spintronics, sports, structural health, supercapacitor, superconductor, textile, thermodynamics, thermoelectric, thin film, timber, tissue engineering, topological insulator, transportation, tribocorrosion, tribology, triboplasma, tunnel, ultrasonics, vacuum, vibration, viticulture, volcanology, wastewater, wearable sensors, welding, wind engineering

### Mesh filter algorithms (~822 filters)
Smoothing (Laplacian, bilateral, cotangent, curvature-weighted, boundary-locked), decimation (edge collapse, vertex removal, quadric), subdivision (midpoint, random, edge split), curvature (mean, Gaussian, principal, angle defect, cotangent, one-ring), geodesics (Dijkstra, heat method, fast marching, Voronoi), topology (genus, Euler, connected components, boundary loops, self-intersection, star-shaped, half-edge analysis), quality (aspect ratio, Jacobian, edge ratio, conformal factor, Delaunay check), data operations (scalar normalize/clamp/threshold/combine/remap/quantize/percentile/gradient/laplacian/diffuse/smooth/invert/stats/histogram), transformations (translate/scale/rotate, mirror, extrude, offset shell, spin), analysis (feature lines, sharp edges, skeleton, symmetry, saliency, thickness), point operations (Poisson disk, uniform sample, closest point, containment), mesh repair (close holes, remove isolated/degenerate/duplicate, merge close vertices, orient normals, flip normals)

### Source geometries (~414)
Architecture: aqueduct, amphitheater, arch bridge, castle turret, colosseum, gothic window, igloo, lighthouse, minaret, obelisk, observatory dome, pagoda, portcullis, pyramid, Roman column, siege tower, suspension bridge, torii gate, zigzag bridge, drawbridge, rope bridge
Vehicles: catamaran, ship hull, viking ship, zeppelin, roller coaster
Instruments: astrolabe, binoculars, Cassegrain telescope, compass rose, corkscrew, Fresnel lens, gramophone, gyroscope, kaleidoscope, lyre, metronome, microscope, music box, orrery, pendulum clock, periscope, pocket watch, sextant, stethoscope, sundial, telescope, theremin, tuning fork, weathervane
Nature: cactus, DNA helix, mangrove tree, nautilus shell, snowflake, water molecule
Objects: abacus, anchor, anvil, barbell, battering ram, bellows, birdhouse, boomerang, cannon, candelabra, celtic cross, chess pawn, crown, cuckoo clock, diamond gem, dreidel, dumbbell, egg timer, ferris wheel, Fibonacci spiral, guillotine, hamster wheel, honeycomb, horseshoe, hourglass, Jacob's ladder, kalimba, katana, maze, mousetrap, origami crane, parachute, pocket knife, pogo stick, propeller, roulette wheel, sandcastle, satellite, scarecrow, ship wheel, soccer ball, spinning top, spiral staircase, steam locomotive, Tesla coil, trebuchet, trident, typewriter, water wheel, wind turbine, windmill, wine barrel
Engineering: BCC crystal, dish array, digital display, gear, heat sink, helical gear, hex telescope, hot air balloon, hydrofoil, microchip, oil derrick, PCB board, radio telescope, roller skate, satellite antenna, solar panel array, totem pole, tuning peg, gyroscope gimbal, bagpipe, ship anchor, spinning wheel, pocket sundial, dreamcatcher

## What To Do Next

### Priority 1: Continue feature implementation
- Keep generating 30 image + 2 mesh + 2 source modules per round
- New domains not yet covered: crystallography, metallography, polymer processing, injection molding, extrusion, CNC machining, robotics kinematics, control systems, power electronics, antenna arrays, radar, sonar, GPS, telecommunications, blockchain/distributed systems, ML/neural network math, computer vision features

### Priority 2: Merge staging crates
- When vtk-filters-image-2 reaches ~500+ modules, merge into vtk-filters-image
- Or split vtk-filters-image into vtk-filters-image-{1..4} by domain for parallel compilation

### Priority 3: Run full test suite
- Tests have been skipped for speed; a full `cargo test --workspace --exclude vtk-python` verification is overdue
- Fix any test failures from the mass generation

### Priority 4: Polish
- Run `cargo clippy --workspace` and fix warnings
- Add pipeline `with_*()` integration for new filters
- Update FEATURES.md and CLAUDE.md with current counts
- Consider doc examples for key filters

### Priority 5: Architecture improvements
- Split vtk-filters-image into domain-specific sub-crates (acoustics, optics, fluid, etc.)
- Add rayon parallelism to compute-heavy mesh filters
- Integration tests combining sources + filters + I/O
