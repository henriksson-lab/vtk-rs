#!/bin/bash
# Universal module generator from spec files.
#
# SPEC FILE FORMAT:
# Lines starting with # are comments
# Empty lines are ignored
#
# For pointwise image filters (use :: as delimiter to avoid | conflicts):
#   img <name> :: <description> :: <rust_expression_using_buf[0]>
#
# For mesh filter files (inline Rust):
#   mesh <name>
#   <full rust source code>
#   ---
#
# For source files (inline Rust):
#   src <name>
#   <full rust source code>
#   ---
#
# Usage:
#   bash scripts/gen_from_specs.sh specs/round42.txt
#   # Then: make gen-mods && cargo test ...

set -e
SPECFILE="${1:?Usage: $0 <specfile>}"
IMG="crates/vtk-filters-image/src"
MESH="crates/vtk-filters-mesh/src"
SRC="crates/vtk-filters/src/sources"
COUNT=0
SKIP=0

gen_img() {
    local name="$1" desc="$2" op="$3"
    [ -f "$IMG/${name}.rs" ] && { SKIP=$((SKIP+1)); return; }
    cat > "$IMG/${name}.rs" << IEOF
//! ${desc}
use vtk_data::{AnyDataArray, DataArray, ImageData};
pub fn ${name}(input: &ImageData, scalars: &str) -> ImageData {
    let arr = match input.point_data().get_array(scalars) { Some(a) if a.num_components()==1=>a, _=>return input.clone() };
    let n = arr.num_tuples(); let mut buf = [0.0f64];
    let data: Vec<f64> = (0..n).map(|i|{arr.tuple_as_f64(i,&mut buf);${op}}).collect();
    let dims = input.dimensions();
    ImageData::with_dimensions(dims[0],dims[1],dims[2]).with_spacing(input.spacing()).with_origin(input.origin())
        .with_point_array(AnyDataArray::F64(DataArray::from_vec(scalars,data,1)))
}
#[cfg(test)] mod tests { use super::*;
    #[test] fn test() { let img=ImageData::from_function([5,5,1],[1.0,1.0,1.0],[0.0,0.0,0.0],"v",|x,_,_|x+1.0);
        let r=${name}(&img,"v"); assert_eq!(r.dimensions(),[5,5,1]); } }
IEOF
    COUNT=$((COUNT+1))
}

gen_block() {
    local dir="$1" name="$2" content="$3"
    [ -f "${dir}/${name}.rs" ] && { SKIP=$((SKIP+1)); return; }
    echo "$content" > "${dir}/${name}.rs"
    COUNT=$((COUNT+1))
}

MODE=""
BLOCK_NAME=""
BLOCK_DIR=""
BLOCK_CONTENT=""

while IFS= read -r line || [ -n "$line" ]; do
    # Skip comments and empty lines (unless inside a block)
    if [ -z "$MODE" ]; then
        [[ "$line" =~ ^#.*$ ]] && continue
        [[ -z "$line" ]] && continue
    fi

    if [ "$MODE" = "block" ]; then
        if [ "$line" = "---" ]; then
            gen_block "$BLOCK_DIR" "$BLOCK_NAME" "$BLOCK_CONTENT"
            MODE=""
            BLOCK_CONTENT=""
        else
            BLOCK_CONTENT="${BLOCK_CONTENT}${line}
"
        fi
        continue
    fi

    # Parse line type
    if [[ "$line" =~ ^img\ (.+)::(.+)::(.+)$ ]]; then
        local_name=$(echo "${BASH_REMATCH[1]}" | xargs)
        local_desc=$(echo "${BASH_REMATCH[2]}" | xargs)
        local_op=$(echo "${BASH_REMATCH[3]}" | xargs)
        gen_img "$local_name" "$local_desc" "$local_op"
    elif [[ "$line" =~ ^mesh\ (.+)$ ]]; then
        BLOCK_NAME=$(echo "${BASH_REMATCH[1]}" | xargs)
        BLOCK_DIR="$MESH"
        MODE="block"
        BLOCK_CONTENT=""
    elif [[ "$line" =~ ^src\ (.+)$ ]]; then
        BLOCK_NAME=$(echo "${BASH_REMATCH[1]}" | xargs)
        BLOCK_DIR="$SRC"
        MODE="block"
        BLOCK_CONTENT=""
    fi
done < "$SPECFILE"

echo "Generated $COUNT new modules ($SKIP skipped)"
echo "Run: make gen-mods && cargo test -p vtk-filters-image --lib && cargo test -p vtk-filters-mesh --lib && cargo test -p vtk-filters --lib"
