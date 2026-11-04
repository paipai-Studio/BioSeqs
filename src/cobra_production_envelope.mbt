/// cobra_production_envelope.mbt — Production envelope (COBRApy flux_analysis.phenotype_phase_plane)
/// Computes min/max objective flux conditioned on a target reaction flux.

///|
/// Production envelope point.
pub struct CobraEnvelopePoint {
  target_flux : Double
  objective_min : Double
  objective_max : Double
  carbon_yield_min : Double
  carbon_yield_max : Double
}

///|
/// Production envelope: vary a target reaction's flux and compute
/// the min/max objective at each point.
///
/// Parameters:
/// - target_reaction: reaction to vary
/// - points: number of grid points
/// - carbon_source: reaction for carbon uptake (for yield calculation)
pub fn cobra_production_envelope(
  model : FbaModel,
  target_reaction : String,
  points~ : Int32 = 20,
  carbon_source~ : String = "",
) -> Array[CobraEnvelopePoint] raise FbaError {
  // Find target reaction index
  let mut target_idx = -1
  for j = 0; j < model.reactions.length(); j = j + 1 {
    if model.reactions[j].id == target_reaction {
      target_idx = j
      break
    }
  }
  if target_idx < 0 {
    raise FbaError("target reaction not found: \{target_reaction}")
  }
  // Find carbon source index
  let mut carbon_idx = -1
  if carbon_source.length() > 0 {
    for j = 0; j < model.reactions.length(); j = j + 1 {
      if model.reactions[j].id == carbon_source {
        carbon_idx = j
        break
      }
    }
  }
  // Determine flux range via FVA
  let fva = cobra_flux_variability_analysis(model, fraction_of_optimum=0.0)
  let target_fva = fva[target_idx]
  let flux_min = target_fva.minimum
  let flux_max = target_fva.maximum
  if flux_max <= flux_min {
    return []
  }
  let results : Array[CobraEnvelopePoint] = Array(capacity=points)
  for p = 0; p < points; p = p + 1 {
    let t = if points > 1 {
      p.to_double() / (points - 1).to_double()
    } else {
      0.0
    }
    let target_flux = flux_min + t * (flux_max - flux_min)
    // Fix target reaction flux and optimize
    let fixed_model = fba_set_reaction_bound(
      model, target_reaction, target_flux, target_flux,
    )
    // Maximize objective
    let sol_max = fba_solve(fixed_model)
    let obj_max = sol_max.objective_value
    // Minimize objective (negate)
    let neg_model = cobra_negate_objective(fixed_model)
    let sol_min = fba_solve(neg_model)
    let obj_min = sol_min.objective_value
    // Carbon yield = objective / carbon_uptake
    let carbon_flux = if carbon_idx >= 0 {
      sol_max.fluxes[carbon_idx].abs()
    } else {
      1.0
    }
    let cy_min = if carbon_flux > 1.0e-9 { obj_min / carbon_flux } else { 0.0 }
    let cy_max = if carbon_flux > 1.0e-9 { obj_max / carbon_flux } else { 0.0 }
    results.push({
      target_flux,
      objective_min: obj_min,
      objective_max: obj_max,
      carbon_yield_min: cy_min,
      carbon_yield_max: cy_max,
    })
  }
  results
}

///|
/// Negate all objective coefficients (to minimize instead of maximize).
fn cobra_negate_objective(model : FbaModel) -> FbaModel {
  let reactions : Array[FbaReaction] = Array(capacity=0)
  for r in model.reactions {
    reactions.push({
      id: r.id,
      name: r.name,
      lower_bound: r.lower_bound,
      upper_bound: r.upper_bound,
      objective: -r.objective,
      genes: r.genes,
    })
  }
  {
    metabolites: model.metabolites,
    reactions,
    stoichiometry: model.stoichiometry,
  }
}
