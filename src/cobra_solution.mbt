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
/// 使用 cobra3 求解器返回的精确对偶解（shadow price / reduced cost）。
pub fn cobra_enhanced_solution(
  model : FbaModel,
  sol : FbaSolution,
) -> CobraSolution {
  let _ = sol
  let m3 = cobra3_from_fba_model(model)
  let opt = cobra3_model_optimize(m3)
  {
    objective_value: opt.objective_value,
    fluxes: opt.fluxes,
    reduced_costs: opt.reduced_costs,
    shadow_prices: opt.shadow_prices,
    status: if opt.status is Optimal {
      "optimal"
    } else {
      "not_optimal"
    },
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
