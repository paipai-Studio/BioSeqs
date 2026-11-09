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
) -> Array[CobraFvaResult] {
  let m3 = cobra3_from_fba_model(model)
  let values = cobra3_fva(m3, fraction_of_optimum~)
  let results : Array[CobraFvaResult] = Array(capacity=values.length())
  for v in values {
    let item : CobraFvaResult = {
      reaction_id: v.reaction_id,
      minimum: v.minimum,
      maximum: v.maximum,
    }
    results.push(item)
  }
  results
}

///|
/// Find blocked reactions: reactions that cannot carry any flux.
pub fn cobra_find_blocked_reactions(
  model : FbaModel,
  zero_cutoff~ : Double = 1.0e-9,
) -> Array[String] {
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
