///|
/// cobra_deletion2.mbt — Gene/reaction deletion analysis (COBRApy round 2)
/// 5 features: single/double gene deletion, single/double reaction deletion, MOMA knockout.
fn cobra2_ko_gene(model : FbaModel, gene_id : String) -> FbaModel {
  let new_rxns = Array(capacity=model.reactions.length())
  for i in 0..<model.reactions.length() {
    let r = model.reactions[i]
    let mut hit = false
    for g in r.genes {
      if g == gene_id {
        hit = true
        break
      }
    }
    if hit {
      new_rxns.push({
        id: r.id,
        name: r.name,
        lower_bound: 0.0,
        upper_bound: 0.0,
        objective: r.objective,
        genes: r.genes,
      })
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

// #11 Single gene deletion

///|
pub struct Cobra2GeneDeletionResult {
  gene_id : String
  growth : Double
  status : String
}

///|
pub fn cobra2_single_gene_deletion(
  model : FbaModel,
) -> Array[Cobra2GeneDeletionResult] raise FbaError {
  let gene_set : Map[String, Bool] = Map([])
  for i in 0..<model.reactions.length() {
    for g in model.reactions[i].genes {
      gene_set[g] = true
    }
  }
  let results = Array(capacity=0)
  for gene in gene_set.keys() {
    let ko_model = cobra2_ko_gene(model, gene)
    let sol = fba_solve(ko_model)
    results.push({
      gene_id: gene,
      growth: sol.objective_value,
      status: sol.status,
    })
  }
  results
}

// #12 Double gene deletion

///|
pub struct Cobra2DoubleGeneDeletionResult {
  gene_pair : (String, String)
  growth : Double
  status : String
}

///|
pub fn cobra2_double_gene_deletion(
  model : FbaModel,
) -> Array[Cobra2DoubleGeneDeletionResult] raise FbaError {
  let gene_list = Array(capacity=0)
  let gene_set : Map[String, Bool] = Map([])
  for i in 0..<model.reactions.length() {
    for g in model.reactions[i].genes {
      if not(gene_set.contains(g)) {
        gene_set[g] = true
        gene_list.push(g)
      }
    }
  }
  let results = Array(capacity=0)
  for i in 0..<gene_list.length() {
    let start_j = i + 1
    for j in start_j..<gene_list.length() {
      let ko_model = cobra2_ko_gene(
        cobra2_ko_gene(model, gene_list[i]),
        gene_list[j],
      )
      let sol = fba_solve(ko_model)
      results.push({
        gene_pair: (gene_list[i], gene_list[j]),
        growth: sol.objective_value,
        status: sol.status,
      })
    }
  }
  results
}

// #13 Single reaction deletion

///|
pub struct Cobra2ReactionDeletionResult {
  reaction_id : String
  growth : Double
  status : String
}

///|
pub fn cobra2_single_reaction_deletion(
  model : FbaModel,
) -> Array[Cobra2ReactionDeletionResult] raise FbaError {
  let results = Array(capacity=model.reactions.length())
  for i in 0..<model.reactions.length() {
    let ko_model = cobra2_ko_rxn(model, model.reactions[i].id)
    let sol = fba_solve(ko_model)
    results.push({
      reaction_id: model.reactions[i].id,
      growth: sol.objective_value,
      status: sol.status,
    })
  }
  results
}

// #14 Double reaction deletion

///|
pub struct Cobra2DoubleReactionDeletionResult {
  reaction_pair : (String, String)
  growth : Double
  status : String
}

///|
pub fn cobra2_double_reaction_deletion(
  model : FbaModel,
) -> Array[Cobra2DoubleReactionDeletionResult] raise FbaError {
  let results = Array(capacity=0)
  for i in 0..<model.reactions.length() {
    let start_j = i + 1
    for j in start_j..<model.reactions.length() {
      let ko_model = cobra2_ko_rxn(
        cobra2_ko_rxn(model, model.reactions[i].id),
        model.reactions[j].id,
      )
      let sol = fba_solve(ko_model)
      results.push({
        reaction_pair: (model.reactions[i].id, model.reactions[j].id),
        growth: sol.objective_value,
        status: sol.status,
      })
    }
  }
  results
}

// #15 MOMA knockout

///|
pub fn cobra2_moma_knockout(
  model : FbaModel,
  gene_id : String,
) -> CobraMomaResult raise FbaError {
  let ref_sol = fba_solve(model)
  let ko_model = cobra2_ko_gene(model, gene_id)
  cobra_moma(ko_model, ref_sol.fluxes)
}

///|
pub fn _cobra_deletion2_dummy() -> Unit {
  let _ = Cobra2GeneDeletionResult::{ gene_id: "", growth: 0.0, status: "", }
  let _ = Cobra2DoubleGeneDeletionResult::{
    gene_pair: ("", ""),
    growth: 0.0,
    status: "",
  }
  let _ = Cobra2ReactionDeletionResult::{
    reaction_id: "",
    growth: 0.0,
    status: "",
  }
  let _ = Cobra2DoubleReactionDeletionResult::{
    reaction_pair: ("", ""),
    growth: 0.0,
    status: "",
  }
}
