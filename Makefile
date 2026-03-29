all:
	cargo run

gitaddall:
	git add src crates examples tests

loc:
	find src tests crates examples -name '*.rs' | xargs wc -l

loc_orig:
	find VTK -name '*.cpp' | xargs wc -l
	find VTK -name '*.h' | xargs wc -l
	find VTK -name '*.c' | xargs wc -l


# Auto-generate lib.rs/mod.rs for all sub-crates and sources
SUBCRATES = vtk-filters-image vtk-filters-mesh vtk-filters-extract vtk-filters-transform \
	vtk-filters-subdivide vtk-filters-clip vtk-filters-smooth vtk-filters-cell \
	vtk-filters-points vtk-filters-statistics vtk-filters-texture vtk-filters-flow \
	vtk-filters-boolean vtk-filters-grid vtk-filters-data vtk-filters-distance \
	vtk-filters-normals vtk-filters-geometry vtk-filters-image-2 vtk-filters-mesh-2

gen-mods:
	@for crate in $(SUBCRATES); do \
		dir="crates/$$crate/src"; \
		desc=$$(echo "$$crate" | sed 's/vtk-filters-//' | tr '-' ' '); \
		echo "//! $$desc filters." > "$$dir/lib.rs"; \
		echo '' >> "$$dir/lib.rs"; \
		for f in "$$dir"/*.rs; do \
			b=$$(basename "$$f" .rs); \
			if [ "$$b" != "lib" ]; then \
				echo "pub mod $$b;" >> "$$dir/lib.rs"; \
			fi; \
		done; \
	done
	@echo '//! Geometry sources.' > crates/vtk-filters/src/sources/mod.rs
	@echo '' >> crates/vtk-filters/src/sources/mod.rs
	@for f in crates/vtk-filters/src/sources/*.rs; do \
		b=$$(basename "$$f" .rs); \
		if [ "$$b" != "mod" ]; then \
			echo "pub mod $$b; pub use $$b::*;" >> crates/vtk-filters/src/sources/mod.rs; \
		fi; \
	done
	@echo "Generated lib.rs for $(words $(SUBCRATES)) sub-crates + mod.rs for sources/"

# Print project stats
stats:
	@echo "=== vtk-rs stats ==="
	@printf "Lines of Rust: "; find crates examples -name '*.rs' | xargs cat | wc -l | tr -d ' '
	@printf "Image filters: "; ls crates/vtk-filters-image/src/*.rs 2>/dev/null | grep -v lib.rs | wc -l | tr -d ' '
	@printf "Mesh filters:  "; ls crates/vtk-filters-mesh/src/*.rs 2>/dev/null | grep -v lib.rs | wc -l | tr -d ' '
	@total=0; for crate in $(SUBCRATES); do \
		n=$$(ls "crates/$$crate/src/"*.rs 2>/dev/null | grep -v lib.rs | wc -l); \
		total=$$((total + n)); \
	done; printf "Sub-crate filters: $$total\n"
	@printf "Core filters:  "; grep '^pub mod' crates/vtk-filters/src/lib.rs | wc -l | tr -d ' '
	@printf "Sources:       "; ls crates/vtk-filters/src/sources/*.rs 2>/dev/null | grep -v mod.rs | wc -l | tr -d ' '
