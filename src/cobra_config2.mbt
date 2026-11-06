///|
/// cobra_config2.mbt — Configuration / other (COBRApy round 2)
/// 3 features: Configuration singleton, find_external_compartment,
/// add_SBO annotations.

// #39 Configuration singleton
pub struct Cobra2Configuration {
  solver : String
  tolerance : Double
  lower_bound : Double
  upper_bound : Double
  processes : Int32
  cache_directory : String
  max_cache_size : Int32
  cache_expiration : Int32
}

///|
pub fn Cobra2Configuration::default() -> Cobra2Configuration {
  {
    solver: "glpk",
    tolerance: 0.000000001,
    lower_bound: -1000.0,
    upper_bound: 1000.0,
    processes: 1,
    cache_directory: "",
    max_cache_size: 100,
    cache_expiration: 3600,
  }
}

///|
pub fn cobra2_config_set_solver(
  config : Cobra2Configuration,
  solver : String,
) -> Cobra2Configuration {
  { ..config, solver, }
}

///|
pub fn cobra2_config_set_tolerance(
  config : Cobra2Configuration,
  tol : Double,
) -> Cobra2Configuration {
  { ..config, tolerance: tol, }
}

///|
pub fn cobra2_config_set_bounds(
  config : Cobra2Configuration,
  lb : Double,
  ub : Double,
) -> Cobra2Configuration {
  { ..config, lower_bound: lb, upper_bound: ub, }
}

///|
pub fn cobra2_config_set_processes(
  config : Cobra2Configuration,
  n : Int32,
) -> Cobra2Configuration {
  { ..config, processes: n, }
}

///|
pub fn cobra2_show_versions() -> String {
  let sb = StringBuilder()
  let _ = sb.write("System Information\n")
  let _ = sb.write("==================\n")
  let _ = sb.write("OS: linux\n")
  let _ = sb.write("Architecture: arm64\n")
  let _ = sb.write("BioSeqs: 0.1.13\n")
  let _ = sb.write("COBRA compatibility: COBRApy v0.32.1\n")
  sb.to_string()
}

// #40 find_external_compartment — heuristic identification

///|
fn cobra2_extract_compartment(met_id : String) -> String {
  let mut bracket_start = -1
  let mut bracket_end = -1
  for i in 0..<met_id.length() {
    let c = met_id[i].to_int()
    if c == 91 {
      bracket_start = i
    }
    if c == 93 {
      bracket_end = i
    }
  }
  if bracket_start >= 0 && bracket_end > bracket_start + 1 {
    met_id[bracket_start + 1:bracket_end].to_owned()
  } else {
    "c"
  }
}

///|
pub fn cobra2_find_external_compartment(model : FbaModel) -> String {
  let compartments = Array(capacity=0)
  let compartment_counts : Map[String, Int32] = Map([])
  for i in 0..<model.metabolites.length() {
    let comp = cobra2_extract_compartment(model.metabolites[i])
    let mut found = false
    for j in 0..<compartments.length() {
      if compartments[j] == comp {
        found = true
        break
      }
    }
    if not(found) {
      compartments.push(comp)
    }
    compartment_counts[comp] = compartment_counts.get(comp).unwrap_or(0) + 1
  }
  let external_shortlist = [
    "e", "ext", "extra", "extracellular", "out", "outside", "periplasm",
  ]
  for i in 0..<external_shortlist.length() {
    let candidate = external_shortlist[i]
    for j in 0..<compartments.length() {
      if compartments[j] == candidate {
        return candidate
      }
    }
  }
  let mut best_comp = "c"
  let mut best_count = -1
  for comp in compartment_counts.keys() {
    let count = compartment_counts[comp]
    if count > best_count {
      best_count = count
      best_comp = comp
    }
  }
  best_comp
}

///|
pub fn cobra2_is_boundary_type(
  reaction_id : String,
  boundary_type : String,
) -> Bool {
  let lower = reaction_id.to_lower()
  match boundary_type {
    "exchange" =>
      lower.has_prefix("ex_") ||
      lower.has_prefix("r_ex_") ||
      lower.has_prefix("exch_")
    "demand" =>
      lower.has_prefix("dm_") ||
      lower.has_prefix("demand_") ||
      (lower.has_prefix("dm") && not(lower.has_prefix("sk_")))
    "sink" => lower.has_prefix("sk_") || lower.has_prefix("sink_")
    _ => false
  }
}

///|
pub fn cobra2_find_boundary_types(
  model : FbaModel,
  boundary_type : String,
) -> Array[String] {
  let result = Array(capacity=0)
  for i in 0..<model.reactions.length() {
    if cobra2_is_boundary_type(model.reactions[i].id, boundary_type) {
      result.push(model.reactions[i].id)
    }
  }
  result
}

// #41 add_SBO annotations

///|
pub struct Cobra2SboTerm {
  term_id : String
  name : String
  description : String
}

///|
pub fn cobra2_sbo_exchange() -> Cobra2SboTerm {
  {
    term_id: "SBO:0000613",
    name: "exchange",
    description: "Exchange reaction",
  }
}

///|
pub fn cobra2_sbo_demand() -> Cobra2SboTerm {
  { term_id: "SBO:0000627", name: "demand", description: "Demand reaction", }
}

///|
pub fn cobra2_sbo_sink() -> Cobra2SboTerm {
  { term_id: "SBO:0000632", name: "sink", description: "Sink reaction", }
}

///|
pub fn cobra2_sbo_biomass() -> Cobra2SboTerm {
  {
    term_id: "SBO:0000629",
    name: "biomass",
    description: "Biomass production",
  }
}

///|
pub fn cobra2_sbo_pseudo() -> Cobra2SboTerm {
  { term_id: "SBO:0000642", name: "pseudo", description: "Pseudo reaction", }
}

///|
pub struct Cobra2AnnotatedReaction {
  reaction_id : String
  sbo_term : Cobra2SboTerm?
}

///|
pub fn cobra2_add_sbo_annotations(
  model : FbaModel,
) -> Array[Cobra2AnnotatedReaction] {
  let result = Array(capacity=model.reactions.length())
  for i in 0..<model.reactions.length() {
    let rxn_id = model.reactions[i].id
    let sbo = if cobra2_is_boundary_type(rxn_id, "exchange") {
      Some(cobra2_sbo_exchange())
    } else if cobra2_is_boundary_type(rxn_id, "demand") {
      Some(cobra2_sbo_demand())
    } else if cobra2_is_boundary_type(rxn_id, "sink") {
      Some(cobra2_sbo_sink())
    } else if rxn_id.to_lower().contains("biomass") {
      Some(cobra2_sbo_biomass())
    } else if rxn_id.to_lower().contains("atpm") ||
      rxn_id.to_lower().contains("maintenance") {
      Some(cobra2_sbo_pseudo())
    } else {
      None
    }
    result.push({ reaction_id: rxn_id, sbo_term: sbo, })
  }
  result
}

///|
pub fn cobra2_check_metabolite_compartment_formula(
  model : FbaModel,
) -> Array[(String, String)] {
  let issues = Array(capacity=0)
  for i in 0..<model.metabolites.length() {
    let met = model.metabolites[i]
    let mut has_non_alnum = false
    for j in 0..<met.length() {
      let c = met[j].to_int()
      if not(
          (c >= 65 && c <= 90) ||
          (c >= 97 && c <= 122) ||
          (c >= 48 && c <= 57) ||
          c == 95 ||
          c == 91 ||
          c == 93,
        ) {
        has_non_alnum = true
        break
      }
    }
    if has_non_alnum {
      issues.push((met, "Non-alphanumeric metabolite ID"))
    }
  }
  issues
}

///|
pub fn _cobra_config2_dummy() -> Unit {
  let _ = Cobra2Configuration::{
    solver: "",
    tolerance: 0.0,
    lower_bound: 0.0,
    upper_bound: 0.0,
    processes: 0,
    cache_directory: "",
    max_cache_size: 0,
    cache_expiration: 0,
  }
  let _ = Cobra2SboTerm::{ term_id: "", name: "", description: "", }
  let _ = Cobra2AnnotatedReaction::{ reaction_id: "", sbo_term: None, }
}
