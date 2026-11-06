///|
/// cobra_model2.mbt — Model manipulation utilities (COBRApy round 2)
/// 5 features: remove metabolites, convert to irreversible, add boundary,
/// prune zero-flux, merge models.

// #27 Remove metabolites not involved in any reaction
pub fn cobra2_prune_metabolites(model : FbaModel) -> FbaModel {
  let used : Map[String, Bool] = Map([])
  for i in 0..<model.stoichiometry.length() {
    for j in 0..<model.stoichiometry[i].length() {
      if model.stoichiometry[i][j] != 0.0 {
        if i < model.metabolites.length() {
          used[model.metabolites[i]] = true
        }
        break
      }
    }
  }
  let new_mets = Array(capacity=0)
  let keep_idx = Array(capacity=0)
  for i in 0..<model.metabolites.length() {
    if used.contains(model.metabolites[i]) {
      new_mets.push(model.metabolites[i])
      keep_idx.push(i)
    }
  }
  let new_stoich = Array(capacity=new_mets.length())
  for i in 0..<keep_idx.length() {
    let old_row = model.stoichiometry[keep_idx[i]]
    new_stoich.push(old_row)
  }
  {
    metabolites: new_mets,
    reactions: model.reactions,
    stoichiometry: new_stoich,
  }
}

// #28 Convert to irreversible — split reversible reactions into forward + reverse

///|
pub fn cobra2_convert_to_irreversible(model : FbaModel) -> FbaModel {
  let new_rxns = Array(capacity=0)
  let new_stoich = Array(capacity=model.stoichiometry.length())
  for _i in 0..<model.stoichiometry.length() {
    new_stoich.push(Array(capacity=0))
  }
  for j in 0..<model.reactions.length() {
    let r = model.reactions[j]
    if r.lower_bound < 0.0 {
      new_rxns.push({
        id: r.id + "_forward",
        name: r.name + " (forward)",
        lower_bound: 0.0,
        upper_bound: r.upper_bound,
        objective: r.objective,
        genes: r.genes,
      })
      for i in 0..<new_stoich.length() {
        if j < model.stoichiometry[i].length() {
          new_stoich[i].push(model.stoichiometry[i][j])
        } else {
          new_stoich[i].push(0.0)
        }
      }
      new_rxns.push({
        id: r.id + "_reverse",
        name: r.name + " (reverse)",
        lower_bound: 0.0,
        upper_bound: -r.lower_bound,
        objective: -r.objective,
        genes: r.genes,
      })
      for i in 0..<new_stoich.length() {
        if j < model.stoichiometry[i].length() {
          new_stoich[i].push(-model.stoichiometry[i][j])
        } else {
          new_stoich[i].push(0.0)
        }
      }
    } else {
      new_rxns.push({
        id: r.id,
        name: r.name,
        lower_bound: r.lower_bound,
        upper_bound: r.upper_bound,
        objective: r.objective,
        genes: r.genes,
      })
      for i in 0..<new_stoich.length() {
        if j < model.stoichiometry[i].length() {
          new_stoich[i].push(model.stoichiometry[i][j])
        } else {
          new_stoich[i].push(0.0)
        }
      }
    }
  }
  {
    metabolites: model.metabolites,
    reactions: new_rxns,
    stoichiometry: new_stoich,
  }
}

// #29 Add boundary/exchange reaction for a metabolite

///|
pub fn cobra2_add_exchange_reaction(
  model : FbaModel,
  met_id : String,
  lb : Double,
  ub : Double,
) -> FbaModel raise FbaError {
  let rxn_id = "EX_" + met_id
  let new_rxn = {
    id: rxn_id,
    name: "Exchange for " + met_id,
    lower_bound: lb,
    upper_bound: ub,
    objective: 0.0,
    genes: [],
  }
  model.add_reaction(new_rxn, [(met_id, -1.0)])
}

// #30 Prune zero-flux reactions

///|
pub fn cobra2_prune_zero_flux(
  model : FbaModel,
  solution : FbaSolution,
  cutoff : Double,
) -> FbaModel {
  let keep_j = Array(capacity=0)
  let new_rxns = Array(capacity=0)
  for j in 0..<model.reactions.length() {
    let flux = if j < solution.fluxes.length() {
      solution.fluxes[j]
    } else {
      0.0
    }
    let abs_flux = if flux < 0.0 { -flux } else { flux }
    if abs_flux > cutoff {
      keep_j.push(j)
      new_rxns.push(model.reactions[j])
    }
  }
  let new_stoich = Array(capacity=model.stoichiometry.length())
  for i in 0..<model.stoichiometry.length() {
    let row = Array(capacity=keep_j.length())
    for k in 0..<keep_j.length() {
      let j = keep_j[k]
      if j < model.stoichiometry[i].length() {
        row.push(model.stoichiometry[i][j])
      } else {
        row.push(0.0)
      }
    }
    new_stoich.push(row)
  }
  {
    metabolites: model.metabolites,
    reactions: new_rxns,
    stoichiometry: new_stoich,
  }
}

// #31 Merge two models

///|
pub fn cobra2_merge_models(model_a : FbaModel, model_b : FbaModel) -> FbaModel {
  let met_set : Map[String, Int32] = Map([])
  let new_mets = Array(capacity=0)
  for i in 0..<model_a.metabolites.length() {
    if not(met_set.contains(model_a.metabolites[i])) {
      met_set[model_a.metabolites[i]] = new_mets.length()
      new_mets.push(model_a.metabolites[i])
    }
  }
  for i in 0..<model_b.metabolites.length() {
    if not(met_set.contains(model_b.metabolites[i])) {
      met_set[model_b.metabolites[i]] = new_mets.length()
      new_mets.push(model_b.metabolites[i])
    }
  }
  let new_rxns = Array(
    capacity=model_a.reactions.length() + model_b.reactions.length(),
  )
  for r in model_a.reactions {
    new_rxns.push(r)
  }
  for r in model_b.reactions {
    new_rxns.push(r)
  }
  let n_rxns_a = model_a.reactions.length()
  let new_stoich = Array(capacity=new_mets.length())
  for i in 0..<new_mets.length() {
    let row = Array(capacity=new_rxns.length())
    for j in 0..<new_rxns.length() {
      let mut val = 0.0
      if j < n_rxns_a {
        let a_idx_opt = met_set.get(model_a.metabolites[i])
        match a_idx_opt {
          Some(a_idx) =>
            if a_idx < model_a.stoichiometry.length() &&
              j < model_a.stoichiometry[a_idx].length() {
              val = model_a.stoichiometry[a_idx][j]
            }
          None => ()
        }
      } else {
        let b_j = j - n_rxns_a
        let b_idx_opt = met_set.get(model_b.metabolites[i])
        match b_idx_opt {
          Some(b_idx) =>
            if b_idx < model_b.stoichiometry.length() &&
              b_j < model_b.stoichiometry[b_idx].length() {
              val = model_b.stoichiometry[b_idx][b_j]
            }
          None => ()
        }
      }
      row.push(val)
    }
    new_stoich.push(row)
  }
  { metabolites: new_mets, reactions: new_rxns, stoichiometry: new_stoich, }
}

///|
pub fn _cobra_model2_dummy() -> Unit {
  ignore(0)
}
