#!/bin/bash
# Fast batch module generator. Reads specs from stdin or file.
# Format per line:
#   img NAME DESCRIPTION OPERATION
#   mesh NAME  (followed by heredoc on stdin)
#   src NAME   (followed by heredoc on stdin)
#   gi NAME DESCRIPTION OPERATION  (shorthand for pointwise image filter)
#
# Usage: bash scripts/gen_round.sh < specs.txt
# Or:    echo "gi image_foo Desc op" | bash scripts/gen_round.sh

IMG="crates/vtk-filters-image/src"
MESH="crates/vtk-filters-mesh/src"
SRC="crates/vtk-filters/src/sources"
COUNT=0

gi() {
    local n=$1 d=$2 o=$3
    [ -f "$IMG/${n}.rs" ] && return
    cat > "$IMG/${n}.rs" << IEOF
//! ${d}
use vtk_data::{AnyDataArray, DataArray, ImageData};
pub fn ${n}(input: &ImageData, scalars: &str) -> ImageData {
    let arr = match input.point_data().get_array(scalars) { Some(a) if a.num_components()==1=>a, _=>return input.clone() };
    let n = arr.num_tuples(); let mut buf = [0.0f64];
    let data: Vec<f64> = (0..n).map(|i|{arr.tuple_as_f64(i,&mut buf);${o}}).collect();
    let dims = input.dimensions();
    ImageData::with_dimensions(dims[0],dims[1],dims[2]).with_spacing(input.spacing()).with_origin(input.origin())
        .with_point_array(AnyDataArray::F64(DataArray::from_vec(scalars,data,1)))
}
#[cfg(test)] mod tests { use super::*;
    #[test] fn test() { let img=ImageData::from_function([5,5,1],[1.0,1.0,1.0],[0.0,0.0,0.0],"v",|x,_,_|x+1.0);
        let r=${n}(&img,"v"); assert_eq!(r.dimensions(),[5,5,1]); } }
IEOF
    COUNT=$((COUNT+1))
}

while IFS= read -r line; do
    case "$line" in
        gi\ *) eval "$line" ;;
    esac
done

echo "Generated $COUNT modules"
