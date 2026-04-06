# NOTES.md

Developer notes on known performance gaps and architectural decisions.

## append: 9x slower than VTK C++

**Our time:** 0.53ms for 3 spheres (32×32). **VTK C++:** 0.06ms.

**Root cause:** Rust ownership semantics require deep-copying all point and cell data when merging meshes. Our `append()` does:
1. Iterate all input points → push each to output Vec (memcpy per point)
2. Iterate all input cell connectivity → copy with offset adjustment
3. Allocate new offsets/connectivity Vecs

VTK C++ uses `vtkDataArray::InsertTuples()` which does bulk `memcpy` of entire typed arrays (one system call per array). VTK also pre-allocates exact output sizes via `SetNumberOfTuples()` then writes via raw pointers — no bounds checks, no Vec growth.

**What we already do:**
- Pre-sized `Vec::with_capacity` for points and connectivity
- Direct `CellArray::from_raw` construction (no per-cell push)

**What we could do:**
- Use `unsafe { ptr::copy_nonoverlapping }` to bulk-copy point data instead of per-element push. This requires the internal `DataArray` storage to expose a raw slice.
- Add `Points::extend_from_points(&other)` that does a single `extend_from_slice` on the flat backing buffer.
- Add `CellArray::append_with_offset(&other, offset)` that bulk-copies connectivity with SIMD-friendly offset addition.

**Future: generic storage backend (zero-copy append)**

The fundamental fix is making `DataArray<T>` generic over its storage, so it can hold either owned data or shared references:

```rust
trait ArrayStorage<T>: AsRef<[T]> {
    fn to_owned(&self) -> Vec<T>;
    fn len(&self) -> usize;
}

impl<T: Clone> ArrayStorage<T> for Vec<T> { ... }       // owned, mutable
impl<T> ArrayStorage<T> for Arc<Vec<T>> { ... }          // shared, zero-copy clone

struct DataArray<T, S: ArrayStorage<T> = Vec<T>> {
    name: String,
    data: S,
    num_components: usize,
}
```

With this, `append()` on `Arc`-backed arrays just clones the `Arc` (one atomic increment, zero data copy). Mutation triggers copy-on-write via `make_owned()`.

A simpler variant that avoids the generic parameter on every type:

```rust
enum PointsStorage {
    Owned(Vec<f64>),
    Shared(Arc<Vec<f64>>, usize, usize),  // data, byte_start, byte_len
}

impl Points {
    fn as_slice(&self) -> &[f64] { ... }  // works for both variants
    fn make_mut(&mut self) -> &mut [f64] {
        // Cow: clone shared data on first write
    }
}
```

This approach:
- **No API change** for read-only consumers (they use `as_slice()` / `get()`)
- **Zero-copy append** when inputs aren't modified afterward
- **Copy-on-write** when a filter needs to mutate the data
- **Thread-safe** — `Arc` is `Send + Sync`
- Matches VTK's `ShallowCopy` semantics without sacrificing Rust's safety guarantees

**Effort:** Medium refactor. `DataArray` internals change, but the public API (`get()`, `tuple()`, `as_slice()`) stays the same. Filters that mutate would call `make_mut()`. The `Points`, `CellArray`, and `DataSetAttributes` wrappers need updating.

**Quick win (no API change):**
Add `Points::extend_from_slice` and `CellArray::append_with_offset` that bulk-copy via `extend_from_slice` on the raw backing `Vec`. Gets within 2-3x of VTK without architectural changes.

**Impact:** Only affects batch mesh merging and operations that concatenate data (append, merge, glyph). Individual filter operations are unaffected. For large meshes the gap narrows since bulk copy dominates over per-element overhead.

## triangulate: no-op in VTK C++ on triangle meshes

**Our time:** ~0.2ms for 128×128 sphere. **VTK C++:** 0.06ms.

**Root cause:** VTK's `vtkTriangleFilter` on an already-triangulated mesh performs a shallow copy (ref-count increment on internal data arrays). Our `triangulate()` detects the all-triangles case and early-exits, but must still `input.clone()` which deep-copies all point and cell data.

This is the same fundamental limitation as `append`: Rust ownership semantics require deep copying where VTK C++ uses `vtkDataObject::ShallowCopy()` (increments `vtkDataArray` reference counts, zero data copy).

**Affected operations (general pattern):** Any filter that returns an unmodified or minimally-modified copy of its input will be slower than VTK C++ by the cost of deep-cloning. This includes:
- `triangulate` on already-triangulated input
- `append` (bulk data merging)
- Any pass-through / identity filter
- Filters that modify a small subset of data but must clone the rest

**Mitigation:** The future generic storage backend (`Arc`-backed `DataArray`, see append section above) would make these shallow copies near-free.

## transform_large: 3.4x slower than VTK C++

**Our time:** ~1ms for 128×128 sphere. **VTK C++:** 0.3ms.

**Root cause:** Same deep-clone pattern. `transform()` clones the entire PolyData, then mutates points in-place. VTK's `vtkTransformPolyDataFilter` shallow-copies the input and only allocates a new point array. The transform itself (matrix multiply per point) is O(n) and fast — the clone is the bottleneck.

**Also applies to:** `warp_by_scalar`, `elevation`, and any filter that clones input then modifies a subset of data.

## Faster than VTK C++: windowed_sinc_smooth (0.34x)

**Our time:** 0.46ms for 20 iterations on 32×32 sphere. **VTK C++:** 1.38ms.

**Technique:** Replaced `HashMap<usize, HashSet<usize>>` adjacency with CSR-like flat arrays (`Vec<u32>` with offset/data layout). Replaced `HashSet<usize>` boundary detection with sorted adjacency + edge counting. Split x/y/z coordinates into separate arrays (SoA layout) for cache-line-friendly iteration during the smoothing loop. Used `unsafe get_unchecked` in the inner neighbor accumulation loop.

**Why it's faster:** VTK's `vtkWindowedSincPolyDataFilter` uses `vtkPolyData::GetPointCells` + `vtkCell::GetPointIds` per vertex per iteration, which involves virtual dispatch and linked-list traversal. Our CSR adjacency is a flat contiguous array — the inner loop reads sequential memory with no indirection.

## Faster than VTK C++: topology_analysis (0.45x)

**Our time:** 3.78ms for 128×128 sphere. **VTK C++:** 8.45ms.

**Technique:** Used packed `u64` edge keys (`(a << 32) | b`) instead of `(usize, usize)` tuple keys — halves hash map key size from 16 bytes to 8 bytes. Union-find with rank compression using `Vec<u32>` (compact) + `Vec<bool>` used-set instead of `HashSet<usize>`. Raw connectivity access via `offsets()`/`connectivity()` slices instead of cell iterator.

## Faster than VTK C++: poly_data_distance (0.04x)

**Our time:** 1.37ms. **VTK C++:** 30.46ms. Already faster — our spatial bucketing is more efficient for this workload.

## Faster than VTK C++: signed_distance_32 (0.06x)

**Our time:** 11.14ms. **VTK C++:** 188.81ms. Our Curless-Levoy with uniform grid bucketing is significantly faster than VTK's approach for this grid size.

## Faster than VTK C++: reflect (0.07x)

**Our time:** 0.30ms. **VTK C++:** 4.13ms. Our flat-slice point transformation avoids VTK's pipeline overhead.

## Faster than VTK C++: gradient (0.41x)

**Our time:** 0.64ms for 32×32 sphere. **VTK C++:** 1.55ms.

**Technique:** Replaced `HashMap<usize, HashSet<usize>>` adjacency with CSR (Compressed Sparse Row) flat arrays. Built via counting pass → offset prefix sum → fill pass → per-vertex sort+dedup. Inner loop uses flat `pts[]` slice with `get_unchecked` for neighbor coordinates. Symmetric ATA matrix only computes upper triangle, then mirrors. All per-point allocations eliminated.

## Faster than VTK C++: validate (0.51x)

**Our time:** 0.40ms for 32×32 sphere. **VTK C++:** 0.79ms.

**Technique:** Replaced spatial hash (HashMap with 3x3x3 neighborhood scan — 27 lookups per point) with sorting-based duplicate detection. Sort point indices by quantized (x,y,z) coordinates at 1e-10 resolution, then scan adjacent pairs. O(n log n) vs O(n) theoretical, but much faster in practice due to cache-friendly sequential access and zero hash overhead.

## Faster than VTK C++: orient_normals (0.54x)

**Our time:** 0.67ms for 32×32 sphere. **VTK C++:** 1.23ms.

**Technique:** Replaced `cells.iter().map(|c| c.to_vec()).collect()` (copies all cell data) with raw offsets/connectivity access. Edge adjacency uses packed u64 keys with (cell_idx, forward) packed into single u64 values. BFS state uses `Vec<u8>` instead of `Vec<Option<bool>>`. Output built with `extend_from_slice` for kept cells, in-place reverse for flipped cells.

## Faster than VTK C++: quadric_clustering (0.25x / 0.61x)

**Our time:** 0.23ms (small) / 3.95ms (large). **VTK C++:** 0.95ms / 6.45ms.

**Technique:** Replaced HashMap-based bin accumulation with flat arrays sized to `divisions^3`. Bin assignment uses precomputed `inv_dx/dy/dz` multipliers. Point output built via flat `out_flat` Vec instead of `Points::push()`. Cell remapping uses raw offsets/connectivity with reusable `mapped` buffer — no per-cell Vec allocation. Degenerate check for triangles is an inline 3-way equality test instead of sort+dedup.

## Faster than VTK C++: cell_quality (0.69x)

**Our time:** 0.30ms. **VTK C++:** 0.44ms.

**Technique:** Replaced cell iterator + per-cell `Vec<[f64;3]>` allocation with raw offsets/connectivity access and a reusable point buffer. Flat `pts[]` slice eliminates per-vertex `points.get()` overhead.

## Faster than VTK C++: connectivity_large (0.28x)

**Our time:** 4.38ms for 128×128 sphere. **VTK C++:** 15.86ms.

**Technique:** Replaced `cell.to_vec()` per-cell allocation with cell index tracking. `build_component_fast()` uses flat `Vec<i64>` point remapping (vec[-1] sentinel) instead of HashMap, and `extend_from_slice` for point data. Union-find uses raw offsets/connectivity instead of cell iterator. Component grouping stores cell indices, not cell data.

## Faster than VTK C++: curvatures_large (0.47x)

**Our time:** 5.83ms for 128×128 sphere. **VTK C++:** 12.45ms.

**Technique:** Switched from cell iterator + `points.get()` to raw offsets/connectivity + flat `pts` slice access. Reading `pts[b0], pts[b0+1], pts[b0+2]` directly avoids the intermediate `[f64;3]` copy and bounds checks per vertex. With 3 vertices per triangle and ~32K triangles, this eliminates ~96K function calls.

## Faster than VTK C++: probe_100 (0.98x)

**Our time:** 0.53ms. **VTK C++:** 0.54ms (effectively matched).

**Technique:** Pre-compute nearest source index using flat slice access for both source and probe points. Separate the nearest-neighbor search from the data copy, avoiding redundant searches when multiple arrays exist. Inner loop uses direct `src_pts[sb]` access instead of `source.points.get(si)`.

## signed_distance_32: 3.6x slower than VTK C++

**Our time:** ~680ms for 32³ grid on 32×32 sphere. **VTK C++:** 189ms.

**Root cause:** Our nearest-triangle search is O(grid_points × num_triangles) brute-force. VTK uses `vtkSignedDistance` with a spatial locator (octree/kd-tree) for O(log n) per query. We added AABB culling on the ray-cast which improved from 6x to 3.6x, but the distance query is still brute-force.

**What would help:** A BVH (bounding volume hierarchy) over triangles for O(log n) nearest-triangle queries. This is a significant implementation effort but would bring us to parity.
