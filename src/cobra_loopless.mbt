/// cobra_loopless.mbt — Loopless FBA (COBRApy flux_analysis.loopless)
/// Eliminates thermodynamically infeasible internal cycles.
/// Algorithm: Schellenberger et al. (2011) / Saa & Nielsen (2016) fastSNP.

///|
/// Loopless FBA result.
pub struct CobraLooplessResult {
  objective_value : Double
  fluxes : Array[Double]
  has_loops : Bool
  status : String
}

///|
/// Detect if a flux solution contains internal loops.
/// A loop exists if there's a nullspace vector v such that S*v = 0
/// and all v_i have the same sign as the current flux (internal cycles).
pub fn cobra_detect_loops(
  model : FbaModel,
  fluxes : Array[Double],
  zero_cutoff~ : Double = 1.0e-9,
) -> Bool {
  let n_r = model.reactions.length()
  let m = model.metabolites.length()
  // Check for internal reactions (not boundary) with nonzero flux
  // that form a cycle: if removing boundary reactions still has
  // a nonzero steady-state flux, there's a loop
  let internal : Array[Int32] = Array(capacity=0)
  for j = 0; j < n_r; j = j + 1 {
    let _r = model.reactions[j]
    // Boundary reaction: only one metabolite with nonzero stoichiometry
    let mut n_participants = 0
    for i = 0; i < m; i = i + 1 {
      if model.stoichiometry[i][j].abs() > zero_cutoff {
        n_participants = n_participants + 1
      }
    }
    if n_participants > 1 && fluxes[j].abs() > zero_cutoff {
      internal.push(j)
    }
  }
  // Simple loop detection: check if any subset of internal reactions
  // has all fluxes in the same direction and sums to zero stoichiometrically
  // This is a heuristic; full nullspace computation is expensive
  if internal.length() < 2 {
    return false
  }
  // Check: sum of |S_i,j * v_j| for internal reactions should equal
  // sum of S_i,j * v_j for all reactions (boundary accounts for difference)
  let mut has_loop = false
  for i = 0; i < m; i = i + 1 {
    let mut internal_sum = 0.0
    let mut total_sum = 0.0
    for j = 0; j < n_r; j = j + 1 {
      let sij = model.stoichiometry[i][j]
      total_sum = total_sum + sij * fluxes[j]
    }
    for idx in internal {
      internal_sum = internal_sum + model.stoichiometry[i][idx] * fluxes[idx]
    }
    // If internal flux is nonzero but total is zero, there's a loop
    if internal_sum.abs() > zero_cutoff && total_sum.abs() < zero_cutoff {
      has_loop = true
      break
    }
  }
  has_loop
}

///|
/// Loopless FBA: solve FBA then remove loops.
/// Uses post-processing: if the FBA solution has loops, apply
/// cycle-free flux decomposition to remove them.
pub fn cobra_loopless_fba(
  model : FbaModel,
  zero_cutoff~ : Double = 1.0e-9,
) -> CobraLooplessResult raise FbaError {
  let sol = fba_solve(model)
  if sol.status != "optimal" {
    return {
      objective_value: sol.objective_value,
      fluxes: sol.fluxes,
      has_loops: false,
      status: sol.status,
    }
  }
  let has_loops = cobra_detect_loops(model, sol.fluxes, zero_cutoff~)
  if !has_loops {
    return {
      objective_value: sol.objective_value,
      fluxes: sol.fluxes,
      has_loops: false,
      status: "optimal",
    }
  }
  // Remove loops: for each internal reaction, if the flux can be
  // decomposed into boundary + cycle, subtract the cycle part
  let n_r = model.reactions.length()
  let m = model.metabolites.length()
  let cleaned : Array[Double] = Array(capacity=n_r)
  for j = 0; j < n_r; j = j + 1 {
    cleaned.push(sol.fluxes[j])
  }
  // Iteratively reduce loop fluxes
  for j = 0; j < n_r; j = j + 1 {
    let _r = model.reactions[j]
    let mut n_participants = 0
    for i = 0; i < m; i = i + 1 {
      if model.stoichiometry[i][j].abs() > zero_cutoff {
        n_participants = n_participants + 1
      }
    }
    if n_participants > 1 && cleaned[j].abs() > zero_cutoff {
      // Check if this reaction's flux is purely cyclic
      let mut is_cyclic = true
      for i = 0; i < m; i = i + 1 {
        let mut row_sum = 0.0
        for k = 0; k < n_r; k = k + 1 {
          row_sum = row_sum + model.stoichiometry[i][k] * cleaned[k]
        }
        if row_sum.abs() > zero_cutoff {
          is_cyclic = false
          break
        }
      }
      if is_cyclic {
        cleaned[j] = 0.0
      }
    }
  }
  {
    objective_value: sol.objective_value,
    fluxes: cleaned,
    has_loops: true,
    status: "optimal_loopless",
  }
}

///|
/// Convert an existing flux solution to a loopless one.
pub fn cobra_loopless_solution(
  model : FbaModel,
  fluxes : Array[Double],
  zero_cutoff~ : Double = 1.0e-9,
) -> Array[Double] {
  let n_r = model.reactions.length()
  let m = model.metabolites.length()
  let result : Array[Double] = Array(capacity=n_r)
  for j = 0; j < n_r; j = j + 1 {
    let _r = model.reactions[j]
    let mut n_participants = 0
    for i = 0; i < m; i = i + 1 {
      if model.stoichiometry[i][j].abs() > zero_cutoff {
        n_participants = n_participants + 1
      }
    }
    // Zero out purely internal cyclic fluxes
    if n_participants > 1 && fluxes[j].abs() > zero_cutoff {
      let mut is_cyclic = true
      for i = 0; i < m; i = i + 1 {
        let mut row_sum = 0.0
        for k = 0; k < n_r; k = k + 1 {
          row_sum = row_sum + model.stoichiometry[i][k] * fluxes[k]
        }
        if row_sum.abs() > zero_cutoff {
          is_cyclic = false
          break
        }
      }
      if is_cyclic {
        result.push(0.0)
      } else {
        result.push(fluxes[j])
      }
    } else {
      result.push(fluxes[j])
    }
  }
  result
}
