/// cobra_mass_balance.mbt — Mass balance checking (COBRApy core.Reaction.check_mass_balance)
/// Check if reactions are element-balanced; convert to irreversible.

///|
/// Mass balance check result for a reaction.
pub struct CobraMassBalanceResult {
  reaction_id : String
  balanced : Bool
  unbalanced_elements : Array[String]
}

///|
/// Check mass balance for a single reaction.
/// Requires metabolite formulas to be provided.
pub fn cobra_check_reaction_mass_balance(
  model : FbaModel,
  reaction_idx : Int32,
  formulas : Map[String, String],
) -> CobraMassBalanceResult {
  let r = model.reactions[reaction_idx]
  let element_counts : Map[String, Double] = Map([], capacity=0)
  for i = 0; i < model.metabolites.length(); i = i + 1 {
    let s = model.stoichiometry[i][reaction_idx]
    if s.abs() < 1.0e-9 {
      continue
    }
    let met_id = model.metabolites[i]
    let formula = formulas.get(met_id).unwrap_or("")
    if formula.length() == 0 {
      continue
    }
    let elements = cobra_parse_formula(formula)
    for elem in elements.keys() {
      let count = elements.get(elem).unwrap_or(0)
      let prev = element_counts.get(elem).unwrap_or(0.0)
      element_counts.set(elem, prev + s * count.to_double())
    }
  }
  let unbalanced : Array[String] = Array(capacity=0)
  for elem in element_counts.keys() {
    if element_counts.get(elem).unwrap_or(0.0).abs() > 1.0e-6 {
      unbalanced.push(elem)
    }
  }
  {
    reaction_id: r.id,
    balanced: unbalanced.length() == 0,
    unbalanced_elements: unbalanced,
  }
}

///|
/// Check mass balance for all reactions.
pub fn cobra_check_all_mass_balance(
  model : FbaModel,
  formulas : Map[String, String],
) -> Array[CobraMassBalanceResult] {
  let results : Array[CobraMassBalanceResult] = Array(capacity=0)
  for j = 0; j < model.reactions.length(); j = j + 1 {
    results.push(cobra_check_reaction_mass_balance(model, j, formulas))
  }
  results
}

///|
/// Convert model to irreversible: split each reversible reaction into
/// forward and reverse reactions.
pub fn cobra_convert_to_irreversible(model : FbaModel) -> FbaModel {
  let reactions : Array[FbaReaction] = Array(capacity=0)
  let m = model.metabolites.length()
  let n_orig = model.reactions.length()
  // First pass: count new reactions
  let mut n_new = 0
  for j = 0; j < n_orig; j = j + 1 {
    if model.reactions[j].lower_bound < 0.0 {
      n_new = n_new + 2
    } else {
      n_new = n_new + 1
    }
  }
  // Build new stoichiometry
  let stoich : Array[Array[Double]] = Array(capacity=0)
  for i = 0; i < m; i = i + 1 {
    stoich.push(Array::make(n_new, 0.0))
  }
  let mut col = 0
  for j = 0; j < n_orig; j = j + 1 {
    let r = model.reactions[j]
    if r.lower_bound < 0.0 {
      // Forward reaction
      reactions.push({
        id: r.id + "_fwd",
        name: r.name + " (forward)",
        lower_bound: 0.0,
        upper_bound: r.upper_bound,
        objective: r.objective,
        genes: r.genes,
      })
      for i = 0; i < m; i = i + 1 {
        stoich[i][col] = model.stoichiometry[i][j]
      }
      col = col + 1
      // Reverse reaction
      reactions.push({
        id: r.id + "_rev",
        name: r.name + " (reverse)",
        lower_bound: 0.0,
        upper_bound: -r.lower_bound,
        objective: -r.objective,
        genes: r.genes,
      })
      for i = 0; i < m; i = i + 1 {
        stoich[i][col] = -model.stoichiometry[i][j]
      }
      col = col + 1
    } else {
      reactions.push({
        id: r.id,
        name: r.name,
        lower_bound: r.lower_bound,
        upper_bound: r.upper_bound,
        objective: r.objective,
        genes: r.genes,
      })
      for i = 0; i < m; i = i + 1 {
        stoich[i][col] = model.stoichiometry[i][j]
      }
      col = col + 1
    }
  }
  { metabolites: model.metabolites, reactions, stoichiometry: stoich, }
}

///|
/// Check charge balance for a reaction.
pub fn cobra_check_charge_balance(
  model : FbaModel,
  reaction_idx : Int32,
  charges : Map[String, Double],
) -> Bool {
  let mut total_charge = 0.0
  for i = 0; i < model.metabolites.length(); i = i + 1 {
    let s = model.stoichiometry[i][reaction_idx]
    if s.abs() < 1.0e-9 {
      continue
    }
    let met_id = model.metabolites[i]
    let charge = charges.get(met_id).unwrap_or(0.0)
    total_charge = total_charge + s * charge
  }
  total_charge.abs() < 1.0e-6
}
