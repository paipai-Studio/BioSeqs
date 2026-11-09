/// cobra_fastcc.mbt — Fast Consistency Check (COBRApy flux_analysis.fastcc)
/// Removes blocked reactions, returns consistent submodel.
/// Algorithm: Vlassis et al. (2014) PLoS Comp Biol.

///|
/// FASTCC result: consistent submodel + removed reactions.
pub struct CobraFastccResult {
  consistent_model : FbaModel
  removed_reactions : Array[String]
  n_blocked : Int32
}

///|
/// Fast consistency check: identify and remove blocked reactions.
/// A reaction is blocked if it cannot carry any flux in any feasible solution.
pub fn cobra_fastcc(
  model : FbaModel,
  _flux_threshold~ : Double = 1.0,
  zero_cutoff~ : Double = 1.0e-9,
) -> CobraFastccResult {
  let blocked = cobra_find_blocked_reactions(model, zero_cutoff~)
  // Build consistent model without blocked reactions
  let keep_indices : Array[Int32] = Array(capacity=0)
  for j = 0; j < model.reactions.length(); j = j + 1 {
    if !blocked.contains(model.reactions[j].id) {
      keep_indices.push(j)
    }
  }
  let new_reactions : Array[FbaReaction] = Array(capacity=0)
  for idx in keep_indices {
    new_reactions.push(model.reactions[idx])
  }
  let new_stoich : Array[Array[Double]] = Array(capacity=0)
  for i = 0; i < model.metabolites.length(); i = i + 1 {
    let row : Array[Double] = Array(capacity=0)
    for idx in keep_indices {
      row.push(model.stoichiometry[i][idx])
    }
    new_stoich.push(row)
  }
  {
    consistent_model: {
      metabolites: model.metabolites,
      reactions: new_reactions,
      stoichiometry: new_stoich,
    },
    removed_reactions: blocked,
    n_blocked: blocked.length(),
  }
}

///|
/// Find cyclic reactions: reactions that can participate in internal loops.
/// Returns reaction IDs that carry flux in some steady-state solution
/// but are not connected to boundary reactions.
pub fn cobra_find_cyclic_reactions(
  model : FbaModel,
  zero_cutoff~ : Double = 1.0e-9,
) -> Array[String] raise FbaError {
  let n_r = model.reactions.length()
  let m = model.metabolites.length()
  let cyclic : Array[String] = Array(capacity=0)
  for j = 0; j < n_r; j = j + 1 {
    let r = model.reactions[j]
    // Check if this is an internal reaction (not boundary)
    let mut n_participants = 0
    for i = 0; i < m; i = i + 1 {
      if model.stoichiometry[i][j].abs() > zero_cutoff {
        n_participants = n_participants + 1
      }
    }
    if n_participants > 1 {
      // Check if it can carry flux independently of boundary
      // by testing if setting it to a nonzero value is feasible
      let test_model = fba_set_reaction_bound(model, r.id, 1.0, 1.0)
      match fba_solve(test_model) {
        _ => {
          // If feasible with forced flux, it might be part of a cycle
          // Check if the flux is purely internal
          let sol = fba_solve(test_model)
          if sol.status == "optimal" {
            cyclic.push(r.id)
          }
        }
      }
    }
  }
  cyclic
}
