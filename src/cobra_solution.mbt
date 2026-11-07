/// cobra_solution.mbt — Enhanced Solution (COBRApy core.Solution)
/// Solution with reduced costs, shadow prices, and DataFrame-like output.

///|
/// Enhanced FBA solution with dual information.
pub struct CobraSolution {
  objective_value : Double
  fluxes : Array[Double]
  reduced_costs : Array[Double]
  shadow_prices : Array[Double]
  status : String
}

///|
/// Compute enhanced solution from FBA result.
/// Reduced costs and shadow prices are approximated from sensitivity analysis.
pub fn cobra_enhanced_solution(
  model : FbaModel,
  sol : FbaSolution,
) -> CobraSolution raise FbaError {
  let n_r = model.reactions.length()
  let m = model.metabolites.length()
  // Approximate reduced costs: for each reaction, perturb its bound
  // and measure objective change (finite difference)
  let reduced_costs : Array[Double] = Array::make(n_r, 0.0)
  let epsilon = 1.0e-6
  for j = 0; j < n_r; j = j + 1 {
    let r = model.reactions[j]
    if r.objective != 0.0 {
      reduced_costs[j] = 0.0 // objective reaction
      continue
    }
    // Perturb upper bound
    let perturbed = fba_set_reaction_bound(
      model,
      r.id,
      r.lower_bound,
      r.upper_bound + epsilon,
    )
    let pert_sol = fba_solve(perturbed)
    if pert_sol.status == "optimal" && sol.status == "optimal" {
      reduced_costs[j] = (pert_sol.objective_value - sol.objective_value) /
        epsilon
    }
  }
  // Shadow prices: dual of steady-state constraints
  // Approximate by perturbing each metabolite's mass balance
  let shadow_prices : Array[Double] = Array::make(m, 0.0)
  for i = 0; i < m; i = i + 1 {
    // Perturb stoichiometry row slightly
    let pert_stoich : Array[Array[Double]] = Array(capacity=0)
    for row = 0; row < m; row = row + 1 {
      let row_arr : Array[Double] = Array(capacity=0)
      for col = 0; col < n_r; col = col + 1 {
        row_arr.push(model.stoichiometry[row][col])
      }
      pert_stoich.push(row_arr)
    }
    // Add small epsilon to first nonzero entry in row i
    for col = 0; col < n_r; col = col + 1 {
      if pert_stoich[i][col].abs() > 1.0e-9 {
        pert_stoich[i][col] = pert_stoich[i][col] + epsilon
        break
      }
    }
    let pert_model : FbaModel = {
      metabolites: model.metabolites,
      reactions: model.reactions,
      stoichiometry: pert_stoich,
    }
    let pert_sol = fba_solve(pert_model)
    if pert_sol.status == "optimal" && sol.status == "optimal" {
      shadow_prices[i] = (pert_sol.objective_value - sol.objective_value) /
        epsilon
    }
  }
  {
    objective_value: sol.objective_value,
    fluxes: sol.fluxes,
    reduced_costs,
    shadow_prices,
    status: sol.status,
  }
}

///|
/// Convert solution to tabular string (like pandas DataFrame).
pub fn CobraSolution::to_frame(
  self : CobraSolution,
  model : FbaModel,
) -> String {
  let lines : Array[String] = Array(capacity=0)
  lines.push("reaction\tflux\treduced_cost")
  for j = 0; j < model.reactions.length(); j = j + 1 {
    let r = model.reactions[j]
    let flux = if j < self.fluxes.length() { self.fluxes[j] } else { 0.0 }
    let rc = if j < self.reduced_costs.length() {
      self.reduced_costs[j]
    } else {
      0.0
    }
    lines.push("\{r.id}\t\{flux}\t\{rc}")
  }
  lines.join("\n")
}

///|
/// Get flux by reaction ID.
pub fn CobraSolution::get_flux(
  self : CobraSolution,
  model : FbaModel,
  reaction_id : String,
) -> Double {
  for j = 0; j < model.reactions.length(); j = j + 1 {
    if model.reactions[j].id == reaction_id {
      if j < self.fluxes.length() {
        return self.fluxes[j]
      }
    }
  }
  0.0
}
