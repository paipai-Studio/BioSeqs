/// cobra_pfba.mbt — Parsimonious FBA (COBRApy flux_analysis.parsimonious)
/// Two-step: (1) maximize objective, (2) minimize sum|v_i| at same objective.
/// Algorithm: Lewis et al. (2010) Molecular Systems Biology.

///|
/// pFBA result: optimal objective + minimal total flux.
pub struct CobraPfbaResult {
  objective_value : Double
  fluxes : Array[Double]
  total_flux : Double
  status : String
}

///|
/// Parsimonious FBA.
/// Step 1: solve FBA to get optimal objective value.
/// Step 2: add constraint objective >= fraction * optimal, minimize sum|v|.
pub fn cobra_pfba(
  model : FbaModel,
  fraction_of_optimum~ : Double = 1.0,
) -> CobraPfbaResult raise FbaError {
  // Step 1: standard FBA
  let fba_sol = fba_solve(model)
  if fba_sol.status != "optimal" {
    return {
      objective_value: fba_sol.objective_value,
      fluxes: fba_sol.fluxes,
      total_flux: 0.0,
      status: fba_sol.status,
    }
  }
  let opt_val = fba_sol.objective_value
  let threshold = fraction_of_optimum * opt_val
  // Step 2: minimize sum|v| = minimize sum(v+ + v-)
  // We reuse the LP solver with modified objective: all reactions get |coeff|=1
  let n_r = model.reactions.length()
  let m = model.metabolites.length()
  let n_var = 2 * n_r
  let c : Array[Double] = Array::make(n_var, 0.0)
  let ub : Array[Double] = Array::make(n_var, 0.0)
  let base : Array[Double] = Array::make(n_r, 0.0)
  // Original objective coefficients
  let orig_obj : Array[Double] = Array::make(n_r, 0.0)
  for j = 0; j < n_r; j = j + 1 {
    let r = model.reactions[j]
    orig_obj[j] = r.objective
    let b_j = if r.lower_bound > 0.0 { r.lower_bound } else { 0.0 }
    base[j] = b_j
    let fwd_ub = r.upper_bound - b_j
    let rev_ub = b_j - r.lower_bound
    ub[j * 2] = if fwd_ub > 0.0 { fwd_ub } else { 0.0 }
    ub[j * 2 + 1] = if rev_ub > 0.0 { rev_ub } else { 0.0 }
    // pFBA objective: minimize sum|v| -> maximize -sum(v+ + v-)
    c[j * 2] = -1.0
    c[j * 2 + 1] = -1.0
  }
  // Steady-state constraints: S*(base + x+ - x-) = 0
  let a_eq : Array[Array[Double]] = Array(capacity=0)
  let b_eq : Array[Double] = Array(capacity=0)
  for i = 0; i < m; i = i + 1 {
    let row : Array[Double] = Array::make(n_var, 0.0)
    let mut rhs = 0.0
    for j = 0; j < n_r; j = j + 1 {
      let sij = model.stoichiometry[i][j]
      if sij != 0.0 {
        row[j * 2] = sij
        row[j * 2 + 1] = -sij
        rhs = rhs - sij * base[j]
      }
    }
    a_eq.push(row)
    b_eq.push(rhs)
  }
  // Add objective constraint: orig_obj . v >= threshold
  // orig_obj . (base + x+ - x-) >= threshold
  // orig_obj . x+ - orig_obj . x- >= threshold - orig_obj . base
  let obj_row : Array[Double] = Array::make(n_var, 0.0)
  let mut obj_rhs = threshold
  for j = 0; j < n_r; j = j + 1 {
    obj_row[j * 2] = orig_obj[j]
    obj_row[j * 2 + 1] = -orig_obj[j]
    obj_rhs = obj_rhs - orig_obj[j] * base[j]
  }
  a_eq.push(obj_row)
  b_eq.push(obj_rhs)
  // Solve
  let result = fba_solve_lp(a_eq, b_eq, c, ub)
  match result {
    None => {
      let fluxes : Array[Double] = Array::make(n_r, 0.0)
      {
        objective_value: opt_val,
        fluxes,
        total_flux: 0.0,
        status: "infeasible",
      }
    }
    Some((_, x)) => {
      let fluxes : Array[Double] = Array(capacity=n_r)
      let mut total = 0.0
      let mut obj_val = 0.0
      for j = 0; j < n_r; j = j + 1 {
        let v = base[j] + x[j * 2] - x[j * 2 + 1]
        fluxes.push(v)
        total = total + v.abs()
        obj_val = obj_val + model.reactions[j].objective * v
      }
      {
        objective_value: obj_val,
        fluxes,
        total_flux: total,
        status: "optimal",
      }
    }
  }
}
