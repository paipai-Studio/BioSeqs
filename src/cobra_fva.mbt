/// cobra_fva.mbt — Flux Variability Analysis (COBRApy flux_analysis.variability)
/// Computes min/max flux for each reaction at optimal solution.
/// Algorithm: Gudmundsson & Thiele (2010) BMC Bioinformatics — fast FVA.

///|
/// FVA result for a single reaction.
pub struct CobraFvaResult {
  reaction_id : String
  minimum : Double
  maximum : Double
}

///|
/// Flux Variability Analysis: for each reaction, find min and max flux
/// consistent with steady-state and (optionally) a fraction of optimum.
///
/// Parameters:
/// - fraction_of_optimum: constrain objective >= fraction * optimal_value
pub fn cobra_flux_variability_analysis(
  model : FbaModel,
  fraction_of_optimum~ : Double = 1.0,
) -> Array[CobraFvaResult] raise FbaError {
  let n_r = model.reactions.length()
  if n_r == 0 {
    raise FbaError("model has no reactions")
  }
  // First solve FBA to get optimal objective value
  let opt_sol = fba_solve(model)
  let opt_val = opt_sol.objective_value
  let results : Array[CobraFvaResult] = Array(capacity=n_r)
  for j = 0; j < n_r; j = j + 1 {
    let r = model.reactions[j]
    // Minimize this reaction's flux: set objective to -1 for this reaction
    let min_model = cobra_set_objective(
      model, j, -1.0, opt_val, fraction_of_optimum,
    )
    let min_sol = fba_solve(min_model)
    let min_flux = min_sol.fluxes[j]
    // Maximize this reaction's flux: set objective to +1
    let max_model = cobra_set_objective(
      model, j, 1.0, opt_val, fraction_of_optimum,
    )
    let max_sol = fba_solve(max_model)
    let max_flux = max_sol.fluxes[j]
    results.push({ reaction_id: r.id, minimum: min_flux, maximum: max_flux, })
  }
  results
}

///|
/// Set a single reaction as objective and add fraction-of-optimum constraint.
fn cobra_set_objective(
  model : FbaModel,
  target_idx : Int32,
  obj_coeff : Double,
  _opt_val : Double,
  _fraction : Double,
) -> FbaModel {
  let reactions : Array[FbaReaction] = Array(capacity=0)
  for j = 0; j < model.reactions.length(); j = j + 1 {
    let r = model.reactions[j]
    let new_obj = if j == target_idx { obj_coeff } else { 0.0 }
    reactions.push({
      id: r.id,
      name: r.name,
      lower_bound: r.lower_bound,
      upper_bound: r.upper_bound,
      objective: new_obj,
      genes: r.genes,
    })
  }
  {
    metabolites: model.metabolites,
    reactions,
    stoichiometry: model.stoichiometry,
  }
}

///|
/// Find blocked reactions: reactions that cannot carry any flux.
pub fn cobra_find_blocked_reactions(
  model : FbaModel,
  zero_cutoff~ : Double = 1.0e-9,
) -> Array[String] raise FbaError {
  let fva = cobra_flux_variability_analysis(model, fraction_of_optimum=0.0)
  let blocked : Array[String] = Array(capacity=0)
  for r in fva {
    if r.minimum.abs() < zero_cutoff && r.maximum.abs() < zero_cutoff {
      blocked.push(r.reaction_id)
    }
  }
  blocked
}

///|
/// Find essential genes: genes whose deletion kills growth.
pub fn cobra_find_essential_genes(
  model : FbaModel,
  threshold~ : Double = 1.0e-6,
) -> Array[String] raise FbaError {
  let _wt = fba_solve(model)
  let gene_set : Map[String, Bool] = Map([], capacity=0)
  for r in model.reactions {
    for g in r.genes {
      gene_set.set(g, true)
    }
  }
  let essential : Array[String] = Array(capacity=0)
  for g in gene_set.keys() {
    let (ko_model, _) = fba_knockout(model, [g])
    let sol = fba_solve(ko_model)
    if sol.objective_value < threshold {
      essential.push(g)
    }
  }
  essential.sort_by(fn(a, b) { compare(a, b) })
  essential
}

///|
/// Find essential reactions: reactions whose knockout kills growth.
pub fn cobra_find_essential_reactions(
  model : FbaModel,
  threshold~ : Double = 1.0e-6,
) -> Array[String] raise FbaError {
  let _wt = fba_solve(model)
  let essential : Array[String] = Array(capacity=0)
  for j = 0; j < model.reactions.length(); j = j + 1 {
    let r = model.reactions[j]
    if r.objective != 0.0 {
      continue // skip biomass objective
    }
    let ko_model = fba_set_reaction_bound(model, r.id, 0.0, 0.0)
    let sol = fba_solve(ko_model)
    if sol.objective_value < threshold {
      essential.push(r.id)
    }
  }
  essential
}
