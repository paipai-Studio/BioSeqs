/// cobra_room.mbt — Regulatory On/Off Minimization (COBRApy flux_analysis.room)
/// Minimizes number of reaction state changes from reference solution.
/// Linear relaxation: minimize sum of on/off deviation weights.
/// Algorithm: Shlomi et al. (2005) PNAS.

///|
/// ROOM result.
pub struct CobraRoomResult {
  objective_value : Double
  fluxes : Array[Double]
  num_changes : Int32
  status : String
}

///|
/// ROOM — Regulatory On/Off Minimization.
/// Uses linear relaxation: for each reaction, if flux changes from reference,
/// add a penalty. delta and epsilon control the tolerance.
pub fn cobra_room(
  model : FbaModel,
  reference : Array[Double],
  delta~ : Double = 0.03,
  epsilon~ : Double = 0.001,
) -> CobraRoomResult {
  let n_r = model.reactions.length()
  let m = model.metabolites.length()
  // Linear ROOM: minimize sum(y_j) where y_j = 1 if |v_j - v_wt_j| > delta * |v_wt_j| + epsilon
  // Approximation: minimize sum(w_j * |v_j - v_wt_j|) with adaptive weights
  let n_var = 2 * n_r + n_r // x+, x-, y (on/off indicator)
  let c : Array[Double] = Array::make(n_var, 0.0)
  let ub : Array[Double] = Array::make(n_var, 0.0)
  let base : Array[Double] = Array::make(n_r, 0.0)
  for j = 0; j < n_r; j = j + 1 {
    let r = model.reactions[j]
    let b_j = if r.lower_bound > 0.0 { r.lower_bound } else { 0.0 }
    base[j] = b_j
    let fwd_ub = r.upper_bound - b_j
    let rev_ub = b_j - r.lower_bound
    ub[j * 2] = if fwd_ub > 0.0 { fwd_ub } else { 0.0 }
    ub[j * 2 + 1] = if rev_ub > 0.0 { rev_ub } else { 0.0 }
    ub[2 * n_r + j] = 1.0 // y_j in [0, 1]
    // Minimize sum(y_j) -> maximize -sum(y_j)
    c[2 * n_r + j] = -1.0
  }
  // Steady-state constraints
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
  // Deviation constraints: |v_j - v_wt_j| <= delta * |v_wt_j| + epsilon + M * y_j
  // where M is a big-M constant
  let big_m = 1000.0
  for j = 0; j < n_r; j = j + 1 {
    let tol = delta * reference[j].abs() + epsilon
    let diff = reference[j] - base[j]
    // v_j - v_wt_j <= tol + M * y_j
    // (base + x+ - x-) - v_wt <= tol + M * y
    // x+ - x- - M*y <= tol + v_wt - base
    // x+ - x- - M*y = tol + diff (with slack)
    let row1 : Array[Double] = Array::make(n_var, 0.0)
    row1[j * 2] = 1.0
    row1[j * 2 + 1] = -1.0
    row1[2 * n_r + j] = -big_m
    a_eq.push(row1)
    b_eq.push(tol + diff)
  }
  // Solve
  let result = fba_solve_lp(a_eq, b_eq, c, ub)
  match result {
    None => {
      let fluxes : Array[Double] = Array::make(n_r, 0.0)
      { objective_value: 0.0, fluxes, num_changes: 0, status: "infeasible", }
    }
    Some((_, x)) => {
      let fluxes : Array[Double] = Array(capacity=n_r)
      let mut changes = 0
      let mut obj_val = 0.0
      for j = 0; j < n_r; j = j + 1 {
        let v = base[j] + x[j * 2] - x[j * 2 + 1]
        fluxes.push(v)
        let tol = delta * reference[j].abs() + epsilon
        if (v - reference[j]).abs() > tol {
          changes = changes + 1
        }
        obj_val = obj_val + model.reactions[j].objective * v
      }
      {
        objective_value: obj_val,
        fluxes,
        num_changes: changes,
        status: "optimal",
      }
    }
  }
}
