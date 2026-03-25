/// VTK cell type identifiers.
///
/// Matches the constants from VTK's `vtkCellType.h`.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
#[repr(u8)]
pub enum CellType {
    // Linear cells
    Empty = 0,
    Vertex = 1,
    PolyVertex = 2,
    Line = 3,
    PolyLine = 4,
    Triangle = 5,
    TriangleStrip = 6,
    Polygon = 7,
    Pixel = 8,
    Quad = 9,
    Tetra = 10,
    Voxel = 11,
    Hexahedron = 12,
    Wedge = 13,
    Pyramid = 14,
    PentagonalPrism = 15,
    HexagonalPrism = 16,

    // Quadratic cells
    QuadraticEdge = 21,
    QuadraticTriangle = 22,
    QuadraticQuad = 23,
    QuadraticTetra = 24,
    QuadraticHexahedron = 25,
    QuadraticWedge = 26,
    QuadraticPyramid = 27,

    // Special
    ConvexPointSet = 41,
    Polyhedron = 42,
}

impl CellType {
    /// Number of points for fixed-size linear cell types.
    /// Returns `None` for variable-size types (PolyVertex, PolyLine, Polygon, etc).
    pub fn num_points(&self) -> Option<usize> {
        match self {
            CellType::Vertex => Some(1),
            CellType::Line => Some(2),
            CellType::Triangle => Some(3),
            CellType::Pixel | CellType::Quad | CellType::Tetra => Some(4),
            CellType::Hexahedron | CellType::Voxel => Some(8),
            CellType::Wedge => Some(6),
            CellType::Pyramid => Some(5),
            CellType::PentagonalPrism => Some(10),
            CellType::HexagonalPrism => Some(12),
            CellType::QuadraticEdge => Some(3),
            CellType::QuadraticTriangle => Some(6),
            CellType::QuadraticQuad => Some(8),
            CellType::QuadraticTetra => Some(10),
            CellType::QuadraticHexahedron => Some(20),
            CellType::QuadraticWedge => Some(15),
            CellType::QuadraticPyramid => Some(13),
            _ => None,
        }
    }

    /// Topological dimension of this cell type.
    pub fn dimension(&self) -> usize {
        match self {
            CellType::Empty => 0,
            CellType::Vertex | CellType::PolyVertex => 0,
            CellType::Line | CellType::PolyLine | CellType::QuadraticEdge => 1,
            CellType::Triangle
            | CellType::TriangleStrip
            | CellType::Polygon
            | CellType::Pixel
            | CellType::Quad
            | CellType::QuadraticTriangle
            | CellType::QuadraticQuad => 2,
            _ => 3,
        }
    }

    /// Create from raw u8 value.
    pub fn from_u8(v: u8) -> Option<CellType> {
        match v {
            0 => Some(CellType::Empty),
            1 => Some(CellType::Vertex),
            2 => Some(CellType::PolyVertex),
            3 => Some(CellType::Line),
            4 => Some(CellType::PolyLine),
            5 => Some(CellType::Triangle),
            6 => Some(CellType::TriangleStrip),
            7 => Some(CellType::Polygon),
            8 => Some(CellType::Pixel),
            9 => Some(CellType::Quad),
            10 => Some(CellType::Tetra),
            11 => Some(CellType::Voxel),
            12 => Some(CellType::Hexahedron),
            13 => Some(CellType::Wedge),
            14 => Some(CellType::Pyramid),
            15 => Some(CellType::PentagonalPrism),
            16 => Some(CellType::HexagonalPrism),
            21 => Some(CellType::QuadraticEdge),
            22 => Some(CellType::QuadraticTriangle),
            23 => Some(CellType::QuadraticQuad),
            24 => Some(CellType::QuadraticTetra),
            25 => Some(CellType::QuadraticHexahedron),
            26 => Some(CellType::QuadraticWedge),
            27 => Some(CellType::QuadraticPyramid),
            41 => Some(CellType::ConvexPointSet),
            42 => Some(CellType::Polyhedron),
            _ => None,
        }
    }
}
