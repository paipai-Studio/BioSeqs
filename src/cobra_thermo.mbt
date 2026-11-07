/// cobra_thermo.mbt — Thermodynamic constraints (COBRApy thermo)
/// Add ΔG constraints to FBA; compute nullspace of stoichiometry matrix.

///|
/// Thermodynamic FBA result.
pub struct CobraThermoResult {
  objective_value : Double
  fluxes : Array[Double]
  feasible : Bool
  status : String
}

///|
/// Compute nullspace of stoichiometry matrix S.
/// Uses Gaussian elimination to find basis vectors for nullspace.
/// Returns array of basis vectors (each of length n_reactions).
pub fn cobra_nullspace(model : FbaModel) -> Array[Array[Double]] {
  let m = model.metabolites.length()
  let n = model.reactions.length()
  if n == 0 || m == 0 {
    return []
  }
  // Copy S to working matrix
  let work : Array[Array[Double]] = Array(capacity=0)
  for i = 0; i < m; i = i + 1 {
    let row : Array[Double] = Array(capacity=0)
    for j = 0; j < n; j = j + 1 {
      row.push(model.stoichiometry[i][j])
    }
    work.push(row)
  }
  // Gaussian elimination with partial pivoting
  let pivot_cols : Array[Int32] = Array(capacity=0)
  let mut rank = 0
  for col = 0; col < n && rank < m; col = col + 1 {
    // Find pivot
    let mut max_val = 0.0
    let mut max_row = -1
    for i = rank; i < m; i = i + 1 {
      if work[i][col].abs() > max_val {
        max_val = work[i][col].abs()
        max_row = i
      }
    }
    if max_row < 0 || max_val < 1.0e-9 {
      continue
    }
    // Swap rows
    if max_row != rank {
      let tmp = work[rank]
      work[rank] = work[max_row]
      work[max_row] = tmp
    }
    // Eliminate
    let piv = work[rank][col]
    for j = 0; j < n; j = j + 1 {
      work[rank][j] = work[rank][j] / piv
    }
    for i = 0; i < m; i = i + 1 {
      if i != rank && work[i][col].abs() > 1.0e-9 {
        let factor = work[i][col]
        for j = 0; j < n; j = j + 1 {
          work[i][j] = work[i][j] - factor * work[rank][j]
        }
      }
    }
    pivot_cols.push(col)
    rank = rank + 1
  }
  // Free columns = nullspace directions
  let free_cols : Array[Int32] = Array(capacity=0)
  for j = 0; j < n; j = j + 1 {
    if !pivot_cols.contains(j) {
      free_cols.push(j)
    }
  }
  // Build nullspace basis
  let basis : Array[Array[Double]] = Array(capacity=0)
  for fc in free_cols {
    let vec : Array[Double] = Array::make(n, 0.0)
    vec[fc] = 1.0
    // For each pivot column, set vec[pivot_col] = -R[pivot_row][fc]
    for r = 0; r < rank && r < pivot_cols.length(); r = r + 1 {
      let pc = pivot_cols[r]
      vec[pc] = -work[r][fc]
    }
    basis.push(vec)
  }
  basis
}

///|
/// Thermodynamic FBA: solve FBA with additional ΔG constraints.
/// A flux vector is thermodynamically feasible if there exists a chemical
/// potential vector μ such that for each reaction j:
///   v_j > 0 => ΔG_j < 0  (forward flux requires negative ΔG)
///   v_j < 0 => ΔG_j > 0  (reverse flux requires positive ΔG)
///   ΔG_j = Σ_i S[i][j] * μ_i
///
/// This simplified version checks post-hoc: after FBA, verify that
/// the solution doesn't contain thermodynamically infeasible loops.
pub fn cobra_thermo_fba(
  model : FbaModel,
  zero_cutoff~ : Double = 1.0e-9,
) -> CobraThermoResult raise FbaError {
  // Solve standard FBA
  let sol = fba_solve(model)
  if sol.status != "optimal" {
    return {
      objective_value: sol.objective_value,
      fluxes: sol.fluxes,
      feasible: false,
      status: sol.status,
    }
  }
  // Check for thermodynamically infeasible loops
  // A loop is a nullspace vector v such that S*v = 0 and all
  // nonzero entries have the same sign pattern as the flux
  let ns = cobra_nullspace(model)
  let mut thermo_feasible = true
  for basis_vec in ns {
    // Check if this nullspace direction is active in the solution
    let mut aligned = true
    let mut has_nonzero = false
    for j = 0; j < model.reactions.length(); j = j + 1 {
      if basis_vec[j].abs() > zero_cutoff && sol.fluxes[j].abs() > zero_cutoff {
        has_nonzero = true
        // If flux and nullspace vector have same sign, it could be a loop
        if basis_vec[j] * sol.fluxes[j] < 0.0 {
          aligned = false
        }
      }
    }
    if aligned && has_nonzero {
      // This nullspace direction is a potential thermodynamic violation
      thermo_feasible = false
      break
    }
  }
  // If not thermo-feasible, apply loopless solution
  let final_fluxes = if !thermo_feasible {
    cobra_loopless_solution(model, sol.fluxes, zero_cutoff~)
  } else {
    sol.fluxes
  }
  {
    objective_value: sol.objective_value,
    fluxes: final_fluxes,
    feasible: thermo_feasible,
    status: if thermo_feasible {
      "optimal_thermo"
    } else {
      "optimal_loopless"
    },
  }
}
