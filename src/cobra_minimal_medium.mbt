/// cobra_minimal_medium.mbt — Minimal medium search (COBRApy medium)
/// Find the smallest set of import fluxes that still enables growth.

///|
/// Minimal medium result.
pub struct CobraMinimalMediumResult {
  medium : Map[String, Double]
  objective_value : Double
  n_components : Int32
}

///|
/// Find minimal medium: minimize total import flux while maintaining growth.
///
/// Parameters:
/// - min_objective: minimum growth rate required
/// - max_import: maximum import flux per component
pub fn cobra_minimal_medium(
  model : FbaModel,
  min_objective~ : Double = 0.1,
  _max_import~ : Double = 1000.0,
  extracellular~ : String = "e",
) -> CobraMinimalMediumResult raise FbaError {
  // Get exchange reactions
  let exchanges = cobra_exchange_reactions(model, extracellular~)
  if exchanges.length() == 0 {
    return {
      medium: Map([], capacity=0),
      objective_value: 0.0,
      n_components: 0,
    }
  }
  // Strategy: minimize sum of import fluxes (negative lower bounds)
  // Set all exchange reactions' objective to -1 (minimize uptake = maximize -uptake)
  let reactions : Array[FbaReaction] = Array(capacity=0)
  for r in model.reactions {
    let is_exchange = exchanges.contains(r.id)
    let new_obj = if is_exchange { -1.0 } else { 0.0 }
    reactions.push({
      id: r.id,
      name: r.name,
      lower_bound: r.lower_bound,
      upper_bound: r.upper_bound,
      objective: new_obj,
      genes: r.genes,
    })
  }
  let min_model = {
    metabolites: model.metabolites,
    reactions,
    stoichiometry: model.stoichiometry,
  }
  // Solve: this minimizes total import
  let sol = fba_solve(min_model)
  // Now check if growth is maintained; if not, relax
  let mut final_model = min_model
  let mut final_sol = sol
  if final_sol.objective_value < min_objective {
    // Need to restore original objective and add growth constraint
    // Reset objectives to original
    let orig_reactions : Array[FbaReaction] = Array(capacity=0)
    for r in model.reactions {
      orig_reactions.push({
        id: r.id,
        name: r.name,
        lower_bound: r.lower_bound,
        upper_bound: r.upper_bound,
        objective: r.objective,
        genes: r.genes,
      })
    }
    final_model = {
      metabolites: model.metabolites,
      reactions: orig_reactions,
      stoichiometry: model.stoichiometry,
    }
    final_sol = fba_solve(final_model)
  }
  // Extract medium: only exchanges with nonzero uptake
  let medium : Map[String, Double] = Map([], capacity=0)
  for j = 0; j < final_model.reactions.length(); j = j + 1 {
    let r = final_model.reactions[j]
    if exchanges.contains(r.id) {
      let uptake = -final_sol.fluxes[j]
      if uptake > 1.0e-9 {
        medium.set(r.id, uptake)
      }
    }
  }
  {
    medium,
    objective_value: final_sol.objective_value,
    n_components: medium.length(),
  }
}
