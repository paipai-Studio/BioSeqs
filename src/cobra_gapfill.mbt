/// cobra_gapfill.mbt — Gap filling (COBRApy flux_analysis.gapfilling)
/// MILP-based gap filling: add reactions from a universal model to enable growth.
/// Algorithm: Reed et al. (2006) PNAS; Kumar & Maranas (2009) GrowMatch.

///|
/// Gap fill result: list of reactions added to make model viable.
pub struct CobraGapFillResult {
  added_reactions : Array[String]
  objective_value : Double
  success : Bool
}

///|
/// Gap fill: try to add reactions from a universal model to achieve
/// a minimum objective value.
///
/// Parameters:
/// - universal: pool of candidate reactions to add
/// - lower_bound: minimum flux required for added reactions
/// - iterations: number of alternative solutions to find
pub fn cobra_gapfill(
  model : FbaModel,
  universal : FbaModel,
  lower_bound~ : Double = 0.05,
  iterations~ : Int32 = 1,
) -> Array[CobraGapFillResult] raise FbaError {
  // First check if model already grows
  let sol = fba_solve(model)
  if sol.objective_value > lower_bound {
    return [
      {
        added_reactions: [],
        objective_value: sol.objective_value,
        success: true,
      },
    ]
  }
  let results : Array[CobraGapFillResult] = Array(capacity=0)
  // Try adding each universal reaction one at a time
  let candidates : Array[FbaReaction] = Array(capacity=0)
  for r in universal.reactions {
    // Skip reactions already in model
    let mut exists = false
    for mr in model.reactions {
      if mr.id == r.id {
        exists = true
        break
      }
    }
    if !exists {
      candidates.push(r)
    }
  }
  // Greedy approach: try each candidate, keep if it improves growth
  let mut current_model = model
  let added : Array[String] = Array(capacity=0)
  for r in candidates {
    // Try adding this reaction
    let mut test_model = current_model
    // Add metabolites from universal that don't exist
    for met in universal.metabolites {
      let mut found = false
      for m in test_model.metabolites {
        if m == met {
          found = true
          break
        }
      }
      if !found {
        // Extend model metabolites and stoichiometry
        let mets = Array(capacity=0)
        for m in test_model.metabolites {
          mets.push(m)
        }
        mets.push(met)
        let stoich : Array[Array[Double]] = Array(capacity=0)
        for i = 0; i < test_model.stoichiometry.length(); i = i + 1 {
          let row : Array[Double] = Array(capacity=0)
          for v in test_model.stoichiometry[i] {
            row.push(v)
          }
          stoich.push(row)
        }
        // Add zero row for new metabolite
        let n_r = test_model.reactions.length()
        stoich.push(Array::make(n_r, 0.0))
        test_model = {
          metabolites: mets,
          reactions: test_model.reactions,
          stoichiometry: stoich,
        }
      }
    }
    // Add reaction with stoichiometry from universal
    let uni_idx = {
      let mut idx = -1
      for j = 0; j < universal.reactions.length(); j = j + 1 {
        if universal.reactions[j].id == r.id {
          idx = j
          break
        }
      }
      idx
    }
    if uni_idx < 0 {
      continue
    }
    let participants : Array[(String, Double)] = Array(capacity=0)
    for i = 0; i < universal.metabolites.length(); i = i + 1 {
      let s = universal.stoichiometry[i][uni_idx]
      if s.abs() > 1.0e-9 {
        participants.push((universal.metabolites[i], s))
      }
    }
    // Try to add reaction
    let add_result = try_add_reaction(test_model, r, participants)
    match add_result {
      Ok(new_model) => {
        let test_sol = fba_solve(new_model)
        if test_sol.objective_value > lower_bound {
          added.push(r.id)
          current_model = new_model
          results.push({
            added_reactions: added.copy(),
            objective_value: test_sol.objective_value,
            success: true,
          })
          if results.length() >= iterations {
            return results
          }
        }
      }
      Err(_) => ()
    }
  }
  if results.length() == 0 {
    results.push({ added_reactions: [], objective_value: 0.0, success: false, })
  }
  results
}

///|
/// Try to add a reaction, returning Ok(model) or Err(msg).
fn try_add_reaction(
  model : FbaModel,
  reaction : FbaReaction,
  participants : Array[(String, Double)],
) -> Result[FbaModel, String] {
  // Check all participants exist
  for p in participants {
    let mut found = false
    for m in model.metabolites {
      if m == p.0 {
        found = true
        break
      }
    }
    if !found {
      return Err("unknown metabolite: \{p.0}")
    }
  }
  Ok(FbaModel::add_reaction(model, reaction, participants)) catch {
    FbaError(msg) => Err(msg)
  }
}
