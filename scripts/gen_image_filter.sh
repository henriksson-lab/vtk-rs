#!/bin/bash
# Generate image filter modules from a simple spec.
# Usage: ./scripts/gen_image_filter.sh name "description" "operation"
# Example: ./scripts/gen_image_filter.sh sqrt "Square root of pixel values" "buf[0].sqrt()"

NAME=$1
DESC=$2
OP=$3

if [ -z "$NAME" ] || [ -z "$DESC" ] || [ -z "$OP" ]; then
    echo "Usage: $0 <name> <description> <operation_on_buf[0]>"
    exit 1
fi

cat > "crates/vtk-filters-image/src/${NAME}.rs" << ENDOFFILE
//! ${DESC}

use vtk_data::{AnyDataArray, DataArray, ImageData};

/// ${DESC}
pub fn ${NAME}(input: &ImageData, scalars: &str) -> ImageData {
    let arr = match input.point_data().get_array(scalars) {
        Some(a) if a.num_components() == 1 => a,
        _ => return input.clone(),
    };
    let n = arr.num_tuples();
    let mut buf = [0.0f64];
    let data: Vec<f64> = (0..n).map(|i| {
        arr.tuple_as_f64(i, &mut buf);
        ${OP}
    }).collect();
    let dims = input.dimensions();
    ImageData::with_dimensions(dims[0], dims[1], dims[2])
        .with_spacing(input.spacing()).with_origin(input.origin())
        .with_point_array(AnyDataArray::F64(DataArray::from_vec(scalars, data, 1)))
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn test_${NAME}() {
        let img = ImageData::from_function([5,5,1],[1.0,1.0,1.0],[0.0,0.0,0.0],"v",|x,_,_|x+1.0);
        let r = ${NAME}(&img, "v");
        assert_eq!(r.dimensions(), [5, 5, 1]);
    }
}
ENDOFFILE

echo "Created crates/vtk-filters-image/src/${NAME}.rs"
