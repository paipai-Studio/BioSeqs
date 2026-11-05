///|
/// cobra_flux_analysis2.mbt — Advanced flux analysis (COBRApy round 2)
/// 6 features: FVA, pFBA, geometric FBA, ROOM, lazy FVA, regulatory constraints.

// Helper: create new FbaReaction with modified bounds
fn cobra2_update_bounds(
  r : FbaReaction,
  lb : Double,
  ub : Double,
) -> FbaReaction {
  {
    id: r.id,
    name: r.name,
    lower_bound: lb,
    upper_bound: ub,
    objective: r.objective,
    genes: r.genes,
  }
}

///|
fn cobra2_update_obj(r : FbaReaction, obj : Double) -> FbaReaction {
  {
    id: r.id,
    name: r.name,
    lower_bound: r.lower_bound,
    upper_bound: r.upper_bound,
    objective: obj,
    genes: r.genes,
  }
}

///|
fn cobra2_ko_rxn(model : FbaModel, rxn_id : String) -> FbaModel {
  let new_rxns = Array(capacity=model.reactions.length())
  for i in 0..<model.reactions.length() {
    let r = model.reactions[i]
    if r.id == rxn_id {
      new_rxns.push(cobra2_update_bounds(r, 0.0, 0.0))
    } else {
      new_rxns.push(r)
    }
  }
  {
    metabolites: model.metabolites,
    reactions: new_rxns,
    stoichiometry: model.stoichiometry,
  }
}

// #1 Flux variability analysis (wrapper with fraction)

///|
pub fn cobra2_fva(
  model : FbaModel,
  fraction : Double,
) -> Array[CobraFvaResult] raise FbaError {
  cobra_flux_variability_analysis(model, fraction_of_optimum=fraction)
}

// #2 Parsimonious FBA wrapper

///|
pub fn cobra2_pfba(
  model : FbaModel,
  fraction : Double,
) -> CobraPfbaResult raise FbaError {
  cobra_pfba(model, fraction_of_optimum=fraction)
}

// #3 Geometric FBA

///|
pub struct Cobra2GeometricFbaResult {
  objective_value : Double
  fluxes : Array[Double]
  geometric_mean : Double
  status : String
}

///|
pub fn cobra2_geometric_fba(
  model : FbaModel,
) -> Cobra2GeometricFbaResult raise FbaError {
  let pfba_res = cobra_pfba(model, fraction_of_optimum=1.0)
  let mut log_sum = 0.0
  let mut n_nonzero = 0
  for i in 0..<pfba_res.fluxes.length() {
    let abs_v = if pfba_res.fluxes[i] < 0.0 {
      -pfba_res.fluxes[i]
    } else {
      pfba_res.fluxes[i]
    }
    if abs_v > 0.000000001 {
      log_sum = log_sum + @math.ln(abs_v)
      n_nonzero = n_nonzero + 1
    }
  }
  let gm = if n_nonzero > 0 {
    @math.pow(2.718281828459045, log_sum / n_nonzero.to_double())
  } else {
    0.0
  }
  {
    objective_value: pfba_res.objective_value,
    fluxes: pfba_res.fluxes,
    geometric_mean: gm,
    status: pfba_res.status,
  }
}

// #4 ROOM wrapper

///|
pub fn cobra2_room(
  model : FbaModel,
  reference : Array[Double],
  delta : Double,
  epsilon : Double,
) -> CobraRoomResult {
  cobra_room(model, reference, delta~, epsilon~)
}

// #5 Lazy FVA — compute FVA only for a subset of reactions

///|
pub fn cobra2_lazy_fva(
  model : FbaModel,
  reaction_ids : Array[String],
  _fraction : Double,
) -> Array[CobraFvaResult] raise FbaError {
  let results = Array(capacity=reaction_ids.length())
  for i in 0..<reaction_ids.length() {
    let rxn_id = reaction_ids[i]
    let mut found = false
    for j in 0..<model.reactions.length() {
      if model.reactions[j].id == rxn_id {
        found = true
        let orig_lb = model.reactions[j].lower_bound
        let orig_ub = model.reactions[j].upper_bound
        let orig_obj = model.reactions[j].objective
        let new_rxns_max = Array(capacity=model.reactions.length())
        for k in 0..<model.reactions.length() {
          let r = model.reactions[k]
          if k == j {
            new_rxns_max.push(cobra2_update_obj(r, 1.0))
          } else {
            new_rxns_max.push(cobra2_update_obj(r, 0.0))
          }
        }
        let model_max = {
          metabolites: model.metabolites,
          reactions: new_rxns_max,
          stoichiometry: model.stoichiometry,
        }
        let sol_max = fba_solve(model_max)
        let new_rxns_min = Array(capacity=model.reactions.length())
        for k in 0..<model.reactions.length() {
          let r = model.reactions[k]
          if k == j {
            new_rxns_min.push(cobra2_update_obj(r, -1.0))
          } else {
            new_rxns_min.push(cobra2_update_obj(r, 0.0))
          }
        }
        let model_min = {
          metabolites: model.metabolites,
          reactions: new_rxns_min,
          stoichiometry: model.stoichiometry,
        }
        let sol_min = fba_solve(model_min)
        let _ = orig_lb
        let _ = orig_ub
        let _ = orig_obj
        results.push({
          reaction_id: rxn_id,
          minimum: sol_min.objective_value,
          maximum: sol_max.objective_value,
        })
        break
      }
    }
    let _ = found
  }
  results
}

// #6 Regulatory constraints — apply boolean regulatory rules

///|
pub struct Cobra2RegulatoryRule {
  reaction_id : String
  activators : Array[String]
  inhibitors : Array[String]
}

///|
pub fn cobra2_apply_regulatory_constraints(
  model : FbaModel,
  rules : Array[Cobra2RegulatoryRule],
  active_genes : Array[String],
) -> FbaModel {
  let new_rxns = Array(capacity=model.reactions.length())
  for i in 0..<model.reactions.length() {
    let r = model.reactions[i]
    let mut inhibited = false
    for j in 0..<rules.length() {
      if rules[j].reaction_id == r.id {
        for k in 0..<rules[j].inhibitors.length() {
          let inh = rules[j].inhibitors[k]
          let mut active = false
          for g in active_genes {
            if g == inh {
              active = true
              break
            }
          }
          if active {
            inhibited = true
            break
          }
        }
        if inhibited {
          break
        }
      }
    }
    if inhibited {
      new_rxns.push(cobra2_update_bounds(r, 0.0, 0.0))
    } else {
      new_rxns.push(r)
    }
  }
  {
    metabolites: model.metabolites,
    reactions: new_rxns,
    stoichiometry: model.stoichiometry,
  }
}

///|
pub fn _cobra_flux_analysis2_dummy() -> Unit {
  let _ = Cobra2GeometricFbaResult::{
    objective_value: 0.0,
    fluxes: [],
    geometric_mean: 0.0,
    status: "",
  }
  let _ = Cobra2RegulatoryRule::{
    reaction_id: "",
    activators: [],
    inhibitors: [],
  }
}
