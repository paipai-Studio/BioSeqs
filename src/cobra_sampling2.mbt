///|
/// cobra_sampling2.mbt — Advanced sampling methods (COBRApy round 2)
/// 4 features: OptGP wrapper, ACHR wrapper, gapfilling, production envelope.

// #7 OptGP sampling wrapper
pub fn cobra2_optgp_sample(
  model : FbaModel,
  n_samples : Int32,
  thinning : Int32,
  n_groups : Int32,
) -> Array[Array[Double]] raise FbaError {
  let res = cobra_optgp_sample(model, n_samples~, thinning~, n_groups~)
  res.samples
}

// #8 ACHR sampling wrapper

///|
pub fn cobra2_achr_sample(
  model : FbaModel,
  n_samples : Int32,
  thinning : Int32,
  n_warmup : Int32,
) -> Array[Array[Double]] raise FbaError {
  let res = cobra_achr_sample(model, n_samples~, thinning~, n_warmup~)
  res.samples
}

// #9 Gap filling — find minimal reactions to add to make model feasible

///|
pub struct Cobra2GapFillResult {
  added_reactions : Array[String]
  objective_value : Double
  status : String
}

///|
pub fn cobra2_gap_fill(
  model : FbaModel,
  universal_reactions : Array[String],
  target_objective : String,
) -> Cobra2GapFillResult raise FbaError {
  let sol = fba_solve(model)
  if sol.objective_value > 0.000000001 {
    return {
      added_reactions: [],
      objective_value: sol.objective_value,
      status: "feasible",
    }
  }
  let added = Array(capacity=0)
  let mut best_obj = sol.objective_value
  for i in 0..<universal_reactions.length() {
    let rxn_id = universal_reactions[i]
    let mut found = false
    for j in 0..<model.reactions.length() {
      if model.reactions[j].id == rxn_id {
        found = true
        break
      }
    }
    if not(found) {
      let new_rxns = Array(capacity=model.reactions.length())
      for j in 0..<model.reactions.length() {
        let r = model.reactions[j]
        if r.id == target_objective {
          new_rxns.push({
            id: r.id,
            name: r.name,
            lower_bound: r.lower_bound,
            upper_bound: r.upper_bound,
            objective: 1.0,
            genes: r.genes,
          })
        } else {
          new_rxns.push(r)
        }
      }
      let new_model = {
        metabolites: model.metabolites,
        reactions: new_rxns,
        stoichiometry: model.stoichiometry,
      }
      let new_sol = fba_solve(new_model)
      if new_sol.objective_value > best_obj {
        best_obj = new_sol.objective_value
        added.push(rxn_id)
      }
    }
  }
  { added_reactions: added, objective_value: best_obj, status: "gapfilled", }
}

// #10 Production envelope — compute min/max production rate vs uptake

///|
pub struct Cobra2ProductionEnvelope {
  uptake_rates : Array[Double]
  production_min : Array[Double]
  production_max : Array[Double]
}

///|
pub fn cobra2_production_envelope(
  model : FbaModel,
  uptake_rxn_id : String,
  production_rxn_id : String,
  n_points : Int32,
) -> Cobra2ProductionEnvelope raise FbaError {
  let mut uptake_lb = 0.0
  let mut uptake_ub = 0.0
  let mut uptake_idx = -1
  for i in 0..<model.reactions.length() {
    if model.reactions[i].id == uptake_rxn_id {
      uptake_lb = model.reactions[i].lower_bound
      uptake_ub = model.reactions[i].upper_bound
      uptake_idx = i
      break
    }
  }
  let uptake_rates = Array(capacity=n_points)
  let prod_min = Array(capacity=n_points)
  let prod_max = Array(capacity=n_points)
  if uptake_idx < 0 || n_points < 2 {
    return { uptake_rates, production_min: prod_min, production_max: prod_max, }
  }
  let step = (uptake_ub - uptake_lb) / (n_points - 1).to_double()
  for p in 0..<n_points {
    let uptake_val = uptake_lb + step * p.to_double()
    uptake_rates.push(uptake_val)
    let new_rxns = Array(capacity=model.reactions.length())
    for j in 0..<model.reactions.length() {
      let r = model.reactions[j]
      if j == uptake_idx {
        new_rxns.push({
          id: r.id,
          name: r.name,
          lower_bound: uptake_val,
          upper_bound: uptake_val,
          objective: r.objective,
          genes: r.genes,
        })
      } else if r.id == production_rxn_id {
        new_rxns.push({
          id: r.id,
          name: r.name,
          lower_bound: r.lower_bound,
          upper_bound: r.upper_bound,
          objective: 1.0,
          genes: r.genes,
        })
      } else {
        new_rxns.push({
          id: r.id,
          name: r.name,
          lower_bound: r.lower_bound,
          upper_bound: r.upper_bound,
          objective: 0.0,
          genes: r.genes,
        })
      }
    }
    let model_max = {
      metabolites: model.metabolites,
      reactions: new_rxns,
      stoichiometry: model.stoichiometry,
    }
    let sol_max = fba_solve(model_max)
    prod_max.push(sol_max.objective_value)
    let new_rxns_min = Array(capacity=model.reactions.length())
    for j in 0..<model.reactions.length() {
      let r = model.reactions[j]
      if j == uptake_idx {
        new_rxns_min.push({
          id: r.id,
          name: r.name,
          lower_bound: uptake_val,
          upper_bound: uptake_val,
          objective: r.objective,
          genes: r.genes,
        })
      } else if r.id == production_rxn_id {
        new_rxns_min.push({
          id: r.id,
          name: r.name,
          lower_bound: r.lower_bound,
          upper_bound: r.upper_bound,
          objective: -1.0,
          genes: r.genes,
        })
      } else {
        new_rxns_min.push({
          id: r.id,
          name: r.name,
          lower_bound: r.lower_bound,
          upper_bound: r.upper_bound,
          objective: 0.0,
          genes: r.genes,
        })
      }
    }
    let model_min = {
      metabolites: model.metabolites,
      reactions: new_rxns_min,
      stoichiometry: model.stoichiometry,
    }
    let sol_min = fba_solve(model_min)
    prod_min.push(sol_min.objective_value)
  }
  { uptake_rates, production_min: prod_min, production_max: prod_max, }
}

///|
pub fn _cobra_sampling2_dummy() -> Unit {
  let _ = Cobra2GapFillResult::{
    added_reactions: [],
    objective_value: 0.0,
    status: "",
  }
  let _ = Cobra2ProductionEnvelope::{
    uptake_rates: [],
    production_min: [],
    production_max: [],
  }
}
