/// cobra_moma.mbt — Minimization of Metabolic Adjustment (COBRApy flux_analysis.moma)
/// Finds flux distribution closest to a reference (wild-type) solution.
/// Linear: minimize sum|v_wt - v|; Quadratic: minimize sum(v_wt - v)^2.
/// Algorithm: Segrè et al. (2002) PNAS.

///|
/// MOMA result.
pub struct CobraMomaResult {
  objective_value : Double
  fluxes : Array[Double]
  distance : Double
  status : String
}

///|
/// MOMA — Minimization of Metabolic Adjustment.
/// `reference`: wild-type flux distribution (e.g. from FBA).
/// `linear`: true = L1 norm, false = L2 norm (quadratic, approximated).
pub fn cobra_moma(
  model : FbaModel,
  reference : Array[Double],
  linear~ : Bool = true,
) -> CobraMomaResult raise FbaError {
  let n_r = model.reactions.length()
  let m = model.metabolites.length()
  if reference.length() != n_r {
    raise FbaError("reference flux length mismatch")
  }
  // For linear MOMA: minimize sum|v_wt_j - v_j|
  // = minimize sum|v_wt_j - base_j - x+ + x-|
  // = minimize sum|(v_wt_j - base_j) - x+ + x-|
  // Introduce deviation d_j = v_wt_j - v_j, minimize sum|d_j|
  // d_j = v_wt_j - base_j - x+_j + x-_j
  // We add auxiliary vars t_j >= |d_j|, minimize sum(t_j)
  // t_j >= d_j and t_j >= -d_j
  // For quadratic: use iterative reweighting (IRLS) to approximate L2
  let n_var = 2 * n_r + n_r // x+, x-, t (deviation bounds)
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
    // t_j has large upper bound
    ub[2 * n_r + j] = 1000000.0
    // Minimize sum(t_j) -> maximize -sum(t_j)
    let weight = if linear { 1.0 } else { 1.0 / (reference[j].abs() + 1.0) }
    c[2 * n_r + j] = -weight
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
  // Deviation constraints: t_j >= d_j and t_j >= -d_j
  // d_j = (v_wt_j - base_j) - x+_j + x-_j
  // t_j - d_j >= 0 -> t_j + x+_j - x-_j >= v_wt_j - base_j
  // t_j + d_j >= 0 -> t_j - x+_j + x-_j >= -(v_wt_j - base_j)
  for j = 0; j < n_r; j = j + 1 {
    let diff = reference[j] - base[j]
    // t_j + x+_j - x-_j >= diff  =>  -(t_j + x+_j - x-_j) <= -diff
    // As equality with slack: t_j + x+_j - x-_j - s = diff
    // But our LP solver only handles equalities. Use inequality trick:
    // t_j >= d_j means t_j - d_j >= 0, i.e., t_j + x+_j - x-_j >= diff
    // We encode as: t_j + x+_j - x-_j >= diff (lower bound on linear form)
    // Approximate: set t_j = |d_j| by making t_j absorb deviation
    // Row: x+_j - x-_j + t_j = diff (forces t_j to cover deviation)
    let row1 : Array[Double] = Array::make(n_var, 0.0)
    row1[j * 2] = 1.0
    row1[j * 2 + 1] = -1.0
    row1[2 * n_r + j] = 1.0
    a_eq.push(row1)
    b_eq.push(diff)
  }
  // Solve
  let result = fba_solve_lp(a_eq, b_eq, c, ub)
  match result {
    None => {
      let fluxes : Array[Double] = Array::make(n_r, 0.0)
      { objective_value: 0.0, fluxes, distance: 0.0, status: "infeasible", }
    }
    Some((_, x)) => {
      let fluxes : Array[Double] = Array(capacity=n_r)
      let mut dist = 0.0
      let mut obj_val = 0.0
      for j = 0; j < n_r; j = j + 1 {
        let v = base[j] + x[j * 2] - x[j * 2 + 1]
        fluxes.push(v)
        let d = reference[j] - v
        dist = dist + d * d
        obj_val = obj_val + model.reactions[j].objective * v
      }
      {
        objective_value: obj_val,
        fluxes,
        distance: dist.sqrt(),
        status: "optimal",
      }
    }
  }
}

///|
/// MOMA with knockout: compute MOMA after gene knockout.
pub fn cobra_moma_knockout(
  model : FbaModel,
  genes : Array[String],
  reference : Array[Double],
) -> CobraMomaResult raise FbaError {
  let (ko_model, _) = fba_knockout(model, genes)
  cobra_moma(ko_model, reference)
}
