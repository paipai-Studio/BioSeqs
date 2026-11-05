///|
/// cobra_summary2.mbt — Summary / reporting (COBRApy round 2)
/// 3 features: Model.summary, Reaction.summary, Metabolite.summary.

// #32 Model.summary
pub struct Cobra2ModelSummary {
  objective_value : Double
  objective_direction : String
  producing : Array[(String, Double)]
  consuming : Array[(String, Double)]
  n_reactions : Int32
  n_metabolites : Int32
}

///|
pub fn cobra2_model_summary(
  model : FbaModel,
  solution : FbaSolution,
) -> Cobra2ModelSummary {
  let producing = Array(capacity=0)
  let consuming = Array(capacity=0)
  for i in 0..<solution.fluxes.length() {
    if solution.fluxes[i] > 0.0000000001 {
      let rxn_id = if i < model.reactions.length() {
        model.reactions[i].id
      } else {
        "rxn_" + i.to_string()
      }
      producing.push((rxn_id, solution.fluxes[i]))
    } else if solution.fluxes[i] < -0.0000000001 {
      let rxn_id = if i < model.reactions.length() {
        model.reactions[i].id
      } else {
        "rxn_" + i.to_string()
      }
      consuming.push((rxn_id, -solution.fluxes[i]))
    }
  }
  {
    objective_value: solution.objective_value,
    objective_direction: "max",
    producing,
    consuming,
    n_reactions: model.reactions.length(),
    n_metabolites: model.metabolites.length(),
  }
}

///|
pub fn cobra2_model_summary_to_string(summary : Cobra2ModelSummary) -> String {
  let sb = StringBuilder()
  let _ = sb.write(
    "Objective: " +
    summary.objective_direction +
    " " +
    summary.objective_value.to_string() +
    "\n",
  )
  let _ = sb.write("Reactions: " + summary.n_reactions.to_string() + "\n")
  let _ = sb.write("Metabolites: " + summary.n_metabolites.to_string() + "\n")
  let _ = sb.write("Producing:\n")
  for i in 0..<summary.producing.length() {
    let (rxn, flux) = summary.producing[i]
    let _ = sb.write("  " + rxn + ": " + flux.to_string() + "\n")
  }
  let _ = sb.write("Consuming:\n")
  for i in 0..<summary.consuming.length() {
    let (rxn, flux) = summary.consuming[i]
    let _ = sb.write("  " + rxn + ": " + flux.to_string() + "\n")
  }
  sb.to_string()
}

// #33 Reaction.summary

///|
pub struct Cobra2ReactionSummary {
  reaction_id : String
  flux : Double
  reactants : Array[(String, Double, Double)]
  products : Array[(String, Double, Double)]
  genes : Array[String]
}

///|
pub fn cobra2_reaction_summary(
  model : FbaModel,
  rxn_idx : Int32,
  solution : FbaSolution,
) -> Cobra2ReactionSummary {
  let r = model.reactions[rxn_idx]
  let flux = if rxn_idx < solution.fluxes.length() {
    solution.fluxes[rxn_idx]
  } else {
    0.0
  }
  let reactants = Array(capacity=0)
  let products = Array(capacity=0)
  for i in 0..<model.metabolites.length() {
    let coeff = if i < model.stoichiometry.length() &&
      rxn_idx < model.stoichiometry[i].length() {
      model.stoichiometry[i][rxn_idx]
    } else {
      0.0
    }
    let met_flux = coeff * flux
    if coeff < 0.0 {
      reactants.push((model.metabolites[i], -coeff, -met_flux))
    } else if coeff > 0.0 {
      products.push((model.metabolites[i], coeff, met_flux))
    }
  }
  { reaction_id: r.id, flux, reactants, products, genes: r.genes, }
}

///|
pub fn cobra2_reaction_summary_to_string(
  summary : Cobra2ReactionSummary,
) -> String {
  let sb = StringBuilder()
  let _ = sb.write("Reaction: " + summary.reaction_id + "\n")
  let _ = sb.write("Flux: " + summary.flux.to_string() + "\n")
  let _ = sb.write("Genes: ")
  for i in 0..<summary.genes.length() {
    let _ = sb.write(summary.genes[i])
    if i < summary.genes.length() - 1 {
      let _ = sb.write(", ")
    }
  }
  let _ = sb.write("\nReactants:\n")
  for i in 0..<summary.reactants.length() {
    let (met, coeff, flux) = summary.reactants[i]
    let _ = sb.write(
      "  " + met + " (" + coeff.to_string() + "): " + flux.to_string() + "\n",
    )
  }
  let _ = sb.write("Products:\n")
  for i in 0..<summary.products.length() {
    let (met, coeff, flux) = summary.products[i]
    let _ = sb.write(
      "  " + met + " (" + coeff.to_string() + "): " + flux.to_string() + "\n",
    )
  }
  sb.to_string()
}

// #34 Metabolite.summary

///|
pub struct Cobra2MetaboliteSummary {
  metabolite_id : String
  producing_flux : Double
  consuming_flux : Double
  net_flux : Double
  producing_reactions : Array[(String, Double)]
  consuming_reactions : Array[(String, Double)]
}

///|
pub fn cobra2_metabolite_summary(
  model : FbaModel,
  met_idx : Int32,
  solution : FbaSolution,
) -> Cobra2MetaboliteSummary {
  let met_id = if met_idx < model.metabolites.length() {
    model.metabolites[met_idx]
  } else {
    "met_" + met_idx.to_string()
  }
  let producing_rxns = Array(capacity=0)
  let consuming_rxns = Array(capacity=0)
  let mut total_producing = 0.0
  let mut total_consuming = 0.0
  for j in 0..<model.reactions.length() {
    let coeff = if met_idx < model.stoichiometry.length() &&
      j < model.stoichiometry[met_idx].length() {
      model.stoichiometry[met_idx][j]
    } else {
      0.0
    }
    let flux = if j < solution.fluxes.length() {
      solution.fluxes[j]
    } else {
      0.0
    }
    let contribution = coeff * flux
    if contribution > 0.0000000001 {
      producing_rxns.push((model.reactions[j].id, contribution))
      total_producing = total_producing + contribution
    } else if contribution < -0.0000000001 {
      consuming_rxns.push((model.reactions[j].id, -contribution))
      total_consuming = total_consuming + -contribution
    }
  }
  {
    metabolite_id: met_id,
    producing_flux: total_producing,
    consuming_flux: total_consuming,
    net_flux: total_producing - total_consuming,
    producing_reactions: producing_rxns,
    consuming_reactions: consuming_rxns,
  }
}

///|
pub fn _cobra_summary2_dummy() -> Unit {
  let _ = Cobra2ModelSummary::{
    objective_value: 0.0,
    objective_direction: "",
    producing: [],
    consuming: [],
    n_reactions: 0,
    n_metabolites: 0,
  }
  let _ = Cobra2ReactionSummary::{
    reaction_id: "",
    flux: 0.0,
    reactants: [],
    products: [],
    genes: [],
  }
  let _ = Cobra2MetaboliteSummary::{
    metabolite_id: "",
    producing_flux: 0.0,
    consuming_flux: 0.0,
    net_flux: 0.0,
    producing_reactions: [],
    consuming_reactions: [],
  }
}
