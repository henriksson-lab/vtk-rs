use crate::data::{DataSetAttributes, FieldData, Points};

/// An atom in a molecule.
#[derive(Debug, Clone)]
pub struct Atom {
    /// Atomic number (1=H, 6=C, 7=N, 8=O, etc.).
    pub atomic_number: u16,
    /// 3D position.
    pub position: [f64; 3],
}

/// A bond between two atoms.
#[derive(Debug, Clone, Copy)]
pub struct Bond {
    /// Index of the first atom.
    pub atom1: usize,
    /// Index of the second atom.
    pub atom2: usize,
    /// Bond order (1=single, 2=double, 3=triple).
    pub order: u16,
}

/// A molecular data structure with atoms and bonds.
///
/// Analogous to VTK's `vtkMolecule`. Stores atomic positions, atomic numbers,
/// and bond connectivity. Supports per-atom and per-bond data arrays.
#[derive(Debug, Clone)]
pub struct Molecule {
    atoms: Vec<Atom>,
    bonds: Vec<Bond>,
    /// Per-atom data attributes.
    atom_data: DataSetAttributes,
    /// Per-bond data attributes.
    bond_data: DataSetAttributes,
    /// Field data (metadata).
    field_data: FieldData,
}

impl Molecule {
    /// Create an empty molecule.
    pub fn new() -> Self {
        Self {
            atoms: Vec::new(),
            bonds: Vec::new(),
            atom_data: DataSetAttributes::new(),
            bond_data: DataSetAttributes::new(),
            field_data: FieldData::new(),
        }
    }

    /// Add an atom and return its index.
    pub fn add_atom(&mut self, atomic_number: u16, position: [f64; 3]) -> usize {
        let idx = self.atoms.len();
        self.atoms.push(Atom { atomic_number, position });
        idx
    }

    /// Add a bond between two atoms. Returns the bond index.
    pub fn add_bond(&mut self, atom1: usize, atom2: usize, order: u16) -> usize {
        let idx = self.bonds.len();
        self.bonds.push(Bond { atom1, atom2, order });
        idx
    }

    /// Number of atoms.
    pub fn num_atoms(&self) -> usize {
        self.atoms.len()
    }

    /// Number of bonds.
    pub fn num_bonds(&self) -> usize {
        self.bonds.len()
    }

    /// Get atom by index.
    pub fn atom(&self, idx: usize) -> &Atom {
        &self.atoms[idx]
    }

    /// Get bond by index.
    pub fn bond(&self, idx: usize) -> &Bond {
        &self.bonds[idx]
    }

    /// Iterate over atoms.
    pub fn atoms(&self) -> &[Atom] {
        &self.atoms
    }

    /// Iterate over bonds.
    pub fn bonds(&self) -> &[Bond] {
        &self.bonds
    }

    /// Get atom positions as Points.
    pub fn positions(&self) -> Points<f64> {
        let mut pts = Points::new();
        for atom in &self.atoms {
            pts.push(atom.position);
        }
        pts
    }

    /// Element symbol for an atomic number.
    pub fn element_symbol(atomic_number: u16) -> &'static str {
        match atomic_number {
            1 => "H", 2 => "He", 3 => "Li", 4 => "Be", 5 => "B",
            6 => "C", 7 => "N", 8 => "O", 9 => "F", 10 => "Ne",
            11 => "Na", 12 => "Mg", 13 => "Al", 14 => "Si", 15 => "P",
            16 => "S", 17 => "Cl", 18 => "Ar", 19 => "K", 20 => "Ca",
            26 => "Fe", 29 => "Cu", 30 => "Zn", 35 => "Br", 53 => "I",
            _ => "?",
        }
    }

    /// Approximate covalent radius in Angstroms for common elements.
    pub fn covalent_radius(atomic_number: u16) -> f64 {
        match atomic_number {
            1 => 0.31, 6 => 0.76, 7 => 0.71, 8 => 0.66, 9 => 0.57,
            15 => 1.07, 16 => 1.05, 17 => 1.02, 35 => 1.20, 53 => 1.39,
            _ => 0.8,
        }
    }

    /// CPK color for common elements (RGB 0-1).
    pub fn cpk_color(atomic_number: u16) -> [f32; 3] {
        match atomic_number {
            1 => [1.0, 1.0, 1.0],       // H: white
            6 => [0.2, 0.2, 0.2],       // C: dark gray
            7 => [0.0, 0.0, 1.0],       // N: blue
            8 => [1.0, 0.0, 0.0],       // O: red
            9 => [0.0, 1.0, 0.0],       // F: green
            15 => [1.0, 0.5, 0.0],      // P: orange
            16 => [1.0, 1.0, 0.0],      // S: yellow
            17 => [0.0, 1.0, 0.0],      // Cl: green
            26 => [0.5, 0.3, 0.0],      // Fe: brown
            _ => [0.8, 0.5, 1.0],       // default: pink
        }
    }

    pub fn atom_data(&self) -> &DataSetAttributes {
        &self.atom_data
    }

    pub fn atom_data_mut(&mut self) -> &mut DataSetAttributes {
        &mut self.atom_data
    }

    pub fn bond_data(&self) -> &DataSetAttributes {
        &self.bond_data
    }

    pub fn bond_data_mut(&mut self) -> &mut DataSetAttributes {
        &mut self.bond_data
    }

    pub fn field_data(&self) -> &FieldData {
        &self.field_data
    }

    pub fn field_data_mut(&mut self) -> &mut FieldData {
        &mut self.field_data
    }
}

impl Default for Molecule {
    fn default() -> Self {
        Self::new()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn water_molecule() {
        let mut mol = Molecule::new();
        let o = mol.add_atom(8, [0.0, 0.0, 0.0]);
        let h1 = mol.add_atom(1, [0.757, 0.586, 0.0]);
        let h2 = mol.add_atom(1, [-0.757, 0.586, 0.0]);
        mol.add_bond(o, h1, 1);
        mol.add_bond(o, h2, 1);

        assert_eq!(mol.num_atoms(), 3);
        assert_eq!(mol.num_bonds(), 2);
        assert_eq!(mol.atom(0).atomic_number, 8);
        assert_eq!(Molecule::element_symbol(8), "O");
    }

    #[test]
    fn methane_molecule() {
        let mut mol = Molecule::new();
        let c = mol.add_atom(6, [0.0, 0.0, 0.0]);
        for pos in &[
            [0.63, 0.63, 0.63],
            [-0.63, -0.63, 0.63],
            [-0.63, 0.63, -0.63],
            [0.63, -0.63, -0.63],
        ] {
            let h = mol.add_atom(1, *pos);
            mol.add_bond(c, h, 1);
        }

        assert_eq!(mol.num_atoms(), 5);
        assert_eq!(mol.num_bonds(), 4);
    }

    #[test]
    fn positions() {
        let mut mol = Molecule::new();
        mol.add_atom(6, [1.0, 2.0, 3.0]);
        mol.add_atom(8, [4.0, 5.0, 6.0]);

        let pts = mol.positions();
        assert_eq!(pts.len(), 2);
        assert_eq!(pts.get(0), [1.0, 2.0, 3.0]);
    }

    #[test]
    fn cpk_colors() {
        let c = Molecule::cpk_color(6);
        assert_eq!(c, [0.2, 0.2, 0.2]);
        let o = Molecule::cpk_color(8);
        assert_eq!(o, [1.0, 0.0, 0.0]);
    }

    #[test]
    fn covalent_radii() {
        assert!(Molecule::covalent_radius(1) < Molecule::covalent_radius(6));
    }
}
