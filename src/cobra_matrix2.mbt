///|
/// cobra_matrix2.mbt — Matrix / solver tools (COBRApy round 2)
/// 4 features: create_stoichiometric_matrix, constraint_matrices,
/// fix_objective_as_constraint, lexicographic FBA.

// #35 create_stoichiometric_matrix (from model.stoichiometry)
pub fn cobra2_create_stoichiometric_matrix(
  model : FbaModel,
) -> Array[Array[Double]] {
  let n_met = model.metabolites.length()
  let n_rxn = model.reactions.length()
  let s = Array(capacity=n_met)
  for i in 0..<n_met {
    let row = Array(capacity=n_rxn)
    for j in 0..<n_rxn {
      let mut val = 0.0
      if i < model.stoichiometry.length() && j < model.stoichiometry[i].length() {
        val = model.stoichiometry[i][j]
      }
      row.push(val)
    }
    s.push(row)
  }
  s
}

///|
pub fn cobra2_create_stoichiometric_matrix_sparse(
  model : FbaModel,
) -> Array[(Int32, Int32, Double)] {
  let entries = Array(capacity=0)
  for i in 0..<model.stoichiometry.length() {
    for j in 0..<model.stoichiometry[i].length() {
      let val = model.stoichiometry[i][j]
      if val != 0.0 {
        entries.push((i, j, val))
      }
    }
  }
  entries
}

// #36 constraint_matrices — full LP matrix form

///|
pub struct Cobra2ConstraintMatrices {
  equalities : Array[Array[Double]]
  b : Array[Double]
  inequalities : Array[Array[Double]]
  bounds : Array[(Double, Double)]
  variable_fixed : Array[Bool]
  variable_bounds : Array[(Double, Double)]
}

///|
pub fn cobra2_constraint_matrices(model : FbaModel) -> Cobra2ConstraintMatrices {
  let n_met = model.metabolites.length()
  let n_rxn = model.reactions.length()
  let equalities = cobra2_create_stoichiometric_matrix(model)
  let b = Array(capacity=n_met)
  for _i in 0..<n_met {
    b.push(0.0)
  }
  let bounds = Array(capacity=n_rxn)
  let var_fixed = Array(capacity=n_rxn)
  let var_bounds = Array(capacity=n_rxn)
  for i in 0..<n_rxn {
    bounds.push(
      (model.reactions[i].lower_bound, model.reactions[i].upper_bound),
    )
    var_bounds.push(
      (model.reactions[i].lower_bound, model.reactions[i].upper_bound),
    )
    var_fixed.push(
      model.reactions[i].lower_bound == model.reactions[i].upper_bound,
    )
  }
  {
    equalities,
    b,
    inequalities: [],
    bounds,
    variable_fixed: var_fixed,
    variable_bounds: var_bounds,
  }
}

// #37 fix_objective_as_constraint

///|
pub fn cobra2_fix_objective_as_constraint(
  model : FbaModel,
  fraction : Double,
) -> FbaModel raise FbaError {
  let sol = fba_solve(model)
  let target = sol.objective_value * fraction
  let new_rxns = Array(capacity=model.reactions.length())
  for i in 0..<model.reactions.length() {
    let r = model.reactions[i]
    if r.objective != 0.0 {
      let new_lb = if r.lower_bound < target / r.objective {
        r.lower_bound
      } else {
        target / r.objective
      }
      new_rxns.push({
        id: r.id,
        name: r.name,
        lower_bound: new_lb,
        upper_bound: r.upper_bound,
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

// #38 lexicographic FBA

///|
pub fn cobra2_lexicographic_fba(
  model : FbaModel,
  objectives : Array[Array[Double]],
  senses : Array[String],
) -> FbaSolution raise FbaError {
  let mut m = model
  let mut final_sol = fba_solve(m)
  for k in 0..<objectives.length() {
    let obj = objectives[k]
    let new_rxns = Array(capacity=m.reactions.length())
    for i in 0..<m.reactions.length() {
      let r = m.reactions[i]
      let new_obj = if i < obj.length() { obj[i] } else { 0.0 }
      let sense = if k < senses.length() { senses[k] } else { "max" }
      let actual_obj = if sense == "min" { -new_obj } else { new_obj }
      new_rxns.push({
        id: r.id,
        name: r.name,
        lower_bound: r.lower_bound,
        upper_bound: r.upper_bound,
        objective: actual_obj,
        genes: r.genes,
      })
    }
    m = {
      metabolites: m.metabolites,
      reactions: new_rxns,
      stoichiometry: m.stoichiometry,
    }
    final_sol = fba_solve(m)
    if k < objectives.length() - 1 {
      let achieved = final_sol.objective_value
      let constraint_rxns = Array(capacity=m.reactions.length())
      for i in 0..<m.reactions.length() {
        let r = m.reactions[i]
        if i < obj.length() && obj[i] != 0.0 {
          let new_lb = achieved / obj[i]
          constraint_rxns.push({
            id: r.id,
            name: r.name,
            lower_bound: new_lb,
            upper_bound: r.upper_bound,
            objective: r.objective,
            genes: r.genes,
          })
        } else {
          constraint_rxns.push(r)
        }
      }
      m = {
        metabolites: m.metabolites,
        reactions: constraint_rxns,
        stoichiometry: m.stoichiometry,
      }
    }
  }
  final_sol
}

///|
pub fn _cobra_matrix2_dummy() -> Unit {
  let _ = Cobra2ConstraintMatrices::{
    equalities: [],
    b: [],
    inequalities: [],
    bounds: [],
    variable_fixed: [],
    variable_bounds: [],
  }
}
