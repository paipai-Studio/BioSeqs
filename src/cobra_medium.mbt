/// cobra_medium.mbt — Medium / boundary management (COBRApy core.Model.medium)
/// Detect and classify boundary reactions (exchange, demand, sink).

///|
/// Boundary reaction type.
pub enum CobraBoundaryType {
  CobraExchange
  CobraDemand
  CobraSink
  CobraInternal
}

///|
/// Classify a reaction based on its stoichiometry.
/// Exchange: extracellular metabolite <-> (one metabolite, extracellular compartment)
/// Demand: one metabolite --> (one metabolite, internal, irreversible out)
/// Sink: one metabolite <-> (one metabolite, internal, reversible)
pub fn cobra_classify_reaction(
  model : FbaModel,
  reaction_idx : Int32,
  extracellular~ : String = "e",
) -> CobraBoundaryType {
  let m = model.metabolites.length()
  // Count metabolites with nonzero stoichiometry
  let mut n_participants = 0
  let mut met_idx = -1
  let mut coeff = 0.0
  for i = 0; i < m; i = i + 1 {
    let s = model.stoichiometry[i][reaction_idx]
    if s.abs() > 1.0e-9 {
      n_participants = n_participants + 1
      met_idx = i
      coeff = s
    }
  }
  if n_participants != 1 {
    return CobraInternal
  }
  if met_idx < 0 {
    return CobraInternal
  }
  let met_id = model.metabolites[met_idx]
  // Check if metabolite is in extracellular compartment
  // Convention: metabolite IDs end with compartment suffix, e.g. "glc__e"
  let is_ext = met_id.has_suffix(extracellular) ||
    met_id.contains("__" + extracellular)
  let r = model.reactions[reaction_idx]
  let is_reversible = r.lower_bound < 0.0
  if is_ext {
    return CobraExchange
  }
  if !is_reversible && coeff < 0.0 {
    return CobraDemand
  }
  if is_reversible {
    return CobraSink
  }
  CobraDemand
}

///|
/// Get all boundary reactions (exchange + demand + sink).
pub fn cobra_boundary_reactions(
  model : FbaModel,
  extracellular~ : String = "e",
) -> Array[String] {
  let result : Array[String] = Array(capacity=0)
  for j = 0; j < model.reactions.length(); j = j + 1 {
    let btype = cobra_classify_reaction(model, j, extracellular~)
    match btype {
      CobraInternal => ()
      _ => result.push(model.reactions[j].id)
    }
  }
  result
}

///|
/// Get exchange reactions only.
pub fn cobra_exchange_reactions(
  model : FbaModel,
  extracellular~ : String = "e",
) -> Array[String] {
  let result : Array[String] = Array(capacity=0)
  for j = 0; j < model.reactions.length(); j = j + 1 {
    if cobra_classify_reaction(model, j, extracellular~) is CobraExchange {
      result.push(model.reactions[j].id)
    }
  }
  result
}

///|
/// Get demand reactions only.
pub fn cobra_demand_reactions(model : FbaModel) -> Array[String] {
  let result : Array[String] = Array(capacity=0)
  for j = 0; j < model.reactions.length(); j = j + 1 {
    if cobra_classify_reaction(model, j) is CobraDemand {
      result.push(model.reactions[j].id)
    }
  }
  result
}

///|
/// Get sink reactions only.
pub fn cobra_sink_reactions(model : FbaModel) -> Array[String] {
  let result : Array[String] = Array(capacity=0)
  for j = 0; j < model.reactions.length(); j = j + 1 {
    if cobra_classify_reaction(model, j) is CobraSink {
      result.push(model.reactions[j].id)
    }
  }
  result
}

///|
/// Get current medium: map of exchange reaction ID -> uptake bound.
pub fn cobra_get_medium(
  model : FbaModel,
  extracellular~ : String = "e",
) -> Map[String, Double] {
  let medium : Map[String, Double] = Map([], capacity=0)
  let exchanges = cobra_exchange_reactions(model, extracellular~)
  for ex_id in exchanges {
    for j = 0; j < model.reactions.length(); j = j + 1 {
      if model.reactions[j].id == ex_id {
        // Medium is the negative lower bound (uptake limit)
        medium.set(ex_id, -model.reactions[j].lower_bound)
      }
    }
  }
  medium
}

///|
/// Set medium: update exchange reaction lower bounds.
pub fn cobra_set_medium(
  model : FbaModel,
  medium : Map[String, Double],
) -> FbaModel raise FbaError {
  let mut result = model
  for ex_id in medium.keys() {
    let uptake = medium.get(ex_id).unwrap_or(0.0)
    // Find the exchange reaction
    let mut found = false
    for j = 0; j < result.reactions.length(); j = j + 1 {
      if result.reactions[j].id == ex_id {
        found = true
        break
      }
    }
    if found {
      result = fba_set_reaction_bound(
        result,
        ex_id,
        -uptake,
        result.reactions[0].upper_bound,
      )
    }
  }
  result
}

///|
/// Add a boundary (exchange) reaction for a metabolite.
pub fn cobra_add_boundary(
  model : FbaModel,
  metabolite_id : String,
  reaction_id~ : String = "",
  lower_bound~ : Double = -10.0,
  upper_bound~ : Double = 1000.0,
) -> FbaModel raise FbaError {
  if !model.metabolites.contains(metabolite_id) {
    raise FbaError("metabolite not found: \{metabolite_id}")
  }
  let rxn_id = if reaction_id.length() > 0 {
    reaction_id
  } else {
    "EX_" + metabolite_id
  }
  let reaction : FbaReaction = {
    id: rxn_id,
    name: "Exchange for " + metabolite_id,
    lower_bound,
    upper_bound,
    objective: 0.0,
    genes: [],
  }
  FbaModel::add_reaction(model, reaction, [(metabolite_id, -1.0)])
}
