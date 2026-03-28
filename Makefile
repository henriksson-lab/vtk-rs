all:
	cargo run

gitaddall:
	git add src crates examples tests

loc:
	find src tests crates examples -name '*.rs' | xargs wc -l

# Auto-generate mod.rs for image/ and mesh/ subdirectories
gen-mods:
	@echo '//! Image processing filters.' > crates/vtk-filters/src/image/mod.rs
	@echo '' >> crates/vtk-filters/src/image/mod.rs
	@for f in crates/vtk-filters/src/image/*.rs; do \
		b=$$(basename "$$f" .rs); \
		if [ "$$b" != "mod" ]; then \
			echo "pub mod $$b;" >> crates/vtk-filters/src/image/mod.rs; \
		fi; \
	done
	@echo '//! Mesh processing filters.' > crates/vtk-filters/src/mesh/mod.rs
	@echo '' >> crates/vtk-filters/src/mesh/mod.rs
	@for f in crates/vtk-filters/src/mesh/*.rs; do \
		b=$$(basename "$$f" .rs); \
		if [ "$$b" != "mod" ]; then \
			echo "pub mod $$b;" >> crates/vtk-filters/src/mesh/mod.rs; \
		fi; \
	done
	@echo '//! Geometry sources.' > crates/vtk-filters/src/sources/mod.rs
	@echo '' >> crates/vtk-filters/src/sources/mod.rs
	@for f in crates/vtk-filters/src/sources/*.rs; do \
		b=$$(basename "$$f" .rs); \
		if [ "$$b" != "mod" ]; then \
			echo "pub mod $$b;" >> crates/vtk-filters/src/sources/mod.rs; \
		fi; \
	done
	@echo "Generated mod.rs for image/ mesh/ sources/"

# Print project stats
stats:
	@echo "=== vtk-rs stats ==="
	@printf "Lines of Rust: "; find crates examples -name '*.rs' | xargs cat | wc -l | tr -d ' '
	@printf "Image filters: "; ls crates/vtk-filters/src/image/*.rs 2>/dev/null | grep -v mod.rs | wc -l | tr -d ' '
	@printf "Mesh filters:  "; ls crates/vtk-filters/src/mesh/*.rs 2>/dev/null | grep -v mod.rs | wc -l | tr -d ' '
	@printf "Other filters: "; grep '^pub mod ' crates/vtk-filters/src/lib.rs | grep -v -E '^pub mod (sources|image|mesh);' | wc -l | tr -d ' '
	@printf "Sources:       "; ls crates/vtk-filters/src/sources/*.rs 2>/dev/null | grep -v mod.rs | wc -l | tr -d ' '
	@printf "Tests:         "; cargo test --workspace --exclude vtk-python -- --list 2>/dev/null | grep ': test$$' | wc -l | tr -d ' '
