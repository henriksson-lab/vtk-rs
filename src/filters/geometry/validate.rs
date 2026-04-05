use crate::data::PolyData;

/// Validation issues found in a PolyData mesh.
#[derive(Debug, Clone, Default)]
pub struct ValidationReport {
    pub degenerate_cells: usize,
    pub out_of_range_indices: usize,
    pub duplicate_points: usize,
    pub empty_cells: usize,
    pub warnings: Vec<String>,
}

impl ValidationReport {
    pub fn is_valid(&self) -> bool {
        self.degenerate_cells == 0
            && self.out_of_range_indices == 0
            && self.empty_cells == 0
    }
}

impl std::fmt::Display for ValidationReport {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        if self.is_valid() {
            write!(f, "Valid")
        } else {
            write!(f, "Invalid: {} degenerate, {} OOB, {} empty",
                self.degenerate_cells, self.out_of_range_indices, self.empty_cells)
        }
    }
}

/// Validate a PolyData mesh for common issues.
pub fn validate(pd: &PolyData) -> ValidationReport {
    let mut report = ValidationReport::default();
    let n_pts = pd.points.len();

    for cell in pd.polys.iter() {
        if cell.is_empty() {
            report.empty_cells += 1;
            continue;
        }
        if cell.len() < 3 {
            report.degenerate_cells += 1;
        }
        for &pid in cell {
            if pid < 0 || (pid as usize) >= n_pts {
                report.out_of_range_indices += 1;
            }
        }
    }

    for cell in pd.lines.iter() {
        if cell.is_empty() {
            report.empty_cells += 1;
        }
        if cell.len() < 2 {
            report.degenerate_cells += 1;
        }
        for &pid in cell {
            if pid < 0 || (pid as usize) >= n_pts {
                report.out_of_range_indices += 1;
            }
        }
    }

    // Check for duplicate points (within tolerance)
    if n_pts <= 10000 {
        let mut dupes = 0;
        for i in 0..n_pts {
            let pi = pd.points.get(i);
            for j in (i + 1)..n_pts {
                let pj = pd.points.get(j);
                let d2 = (pi[0]-pj[0]).powi(2) + (pi[1]-pj[1]).powi(2) + (pi[2]-pj[2]).powi(2);
                if d2 < 1e-20 {
                    dupes += 1;
                    break;
                }
            }
        }
        report.duplicate_points = dupes;
    }

    if report.duplicate_points > 0 {
        report.warnings.push(format!("{} duplicate points detected", report.duplicate_points));
    }

    report
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn valid_mesh() {
        let pd = PolyData::from_triangles(
            vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]],
            vec![[0, 1, 2]],
        );
        let report = validate(&pd);
        assert!(report.is_valid());
    }

    #[test]
    fn duplicate_points() {
        let pd = PolyData::from_triangles(
            vec![[0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [1.0, 0.0, 0.0]],
            vec![[0, 1, 2]],
        );
        let report = validate(&pd);
        assert!(report.duplicate_points > 0);
    }

    #[test]
    fn display() {
        let pd = PolyData::from_triangles(
            vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]],
            vec![[0, 1, 2]],
        );
        let report = validate(&pd);
        assert_eq!(format!("{report}"), "Valid");
    }
}
