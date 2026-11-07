/// biotite_rcsb_search.mbt — RCSB PDB advanced search query builder
/// Inspired by biotite.database.rcsb

///|
/// RCSB search condition types.
pub enum BiotiteRcsbSearchCondition {
  BiotiteRcsbExperimentalMethod(String)
  BiotiteRcsbResolution(Double, Double) // min, max
  BiotiteRcsbMolecularWeight(Double, Double) // min, max
  BiotiteRcsbOrganism(String)
  BiotiteRcsbLigand(String)
  BiotiteRcsbHasLigand(Bool)
  BiotiteRcsbDepositionDate(String, String) // from, to
  BiotiteRcsbReleaseDate(String, String) // from, to
  BiotiteRcsbNumChains(Int32, Int32) // min, max
  BiotiteRcsbNumResidues(Int32, Int32) // min, max
  BiotiteRcsbSpaceGroup(String)
  BiotiteRcsbPdbId(String)
}

///|
/// RCSB search query builder.
pub struct BiotiteRcsbSearchQuery {
  conditions : Array[BiotiteRcsbSearchCondition]
  return_type : String // "entry", "polymer_entity", etc.
  sort_by : String
  sort_dir : String
}

///|
/// Create a new RCSB search query.
pub fn BiotiteRcsbSearchQuery::new() -> BiotiteRcsbSearchQuery {
  { conditions: [], return_type: "entry", sort_by: "", sort_dir: "asc", }
}

///|
/// Add a search condition.
pub fn BiotiteRcsbSearchQuery::add_condition(
  self : BiotiteRcsbSearchQuery,
  condition : BiotiteRcsbSearchCondition,
) -> BiotiteRcsbSearchQuery {
  let new_conditions : Array[BiotiteRcsbSearchCondition] = []
  for i = 0; i < self.conditions.length(); i = i + 1 {
    new_conditions.push(self.conditions[i])
  }
  new_conditions.push(condition)
  {
    conditions: new_conditions,
    return_type: self.return_type,
    sort_by: self.sort_by,
    sort_dir: self.sort_dir,
  }
}

///|
/// Set return type.
pub fn BiotiteRcsbSearchQuery::set_return_type(
  self : BiotiteRcsbSearchQuery,
  return_type : String,
) -> BiotiteRcsbSearchQuery {
  {
    conditions: self.conditions,
    return_type,
    sort_by: self.sort_by,
    sort_dir: self.sort_dir,
  }
}

///|
/// Set sort.
pub fn BiotiteRcsbSearchQuery::set_sort(
  self : BiotiteRcsbSearchQuery,
  sort_by : String,
  sort_dir~ : String = "asc",
) -> BiotiteRcsbSearchQuery {
  {
    conditions: self.conditions,
    return_type: self.return_type,
    sort_by,
    sort_dir,
  }
}

///|
/// Build the RCSB search API URL.
pub fn BiotiteRcsbSearchQuery::build_url(
  self : BiotiteRcsbSearchQuery,
) -> String {
  let base = "https://search.rcsb.org/rcsbsearch/v2/query"
  let json = self.build_json()
  base + "?json=" + json
}

///|
/// Build the JSON query body for RCSB search API.
pub fn BiotiteRcsbSearchQuery::build_json(
  self : BiotiteRcsbSearchQuery,
) -> String {
  let sb = StringBuilder()
  let _ = sb.write(
    "{\"query\":{\"type\":\"group\",\"logical_operator\":\"and\",\"nodes\":[",
  )
  for i = 0; i < self.conditions.length(); i = i + 1 {
    if i > 0 {
      let _ = sb.write(",")
    }
    let _ = sb.write(self.condition_to_json(self.conditions[i]))
  }
  let _ = sb.write("]},\"return_type\":\"" + self.return_type + "\"")
  if self.sort_by != "" {
    let _ = sb.write(
      ",\"sort\":[{\"sort_by\":\"" +
      self.sort_by +
      "\",\"direction\":\"" +
      self.sort_dir +
      "\"}]",
    )
  }
  let _ = sb.write("}")
  sb.to_string()
}

///|
/// Convert a single condition to JSON.
fn BiotiteRcsbSearchQuery::condition_to_json(
  _self : BiotiteRcsbSearchQuery,
  cond : BiotiteRcsbSearchCondition,
) -> String {
  match cond {
    BiotiteRcsbExperimentalMethod(method_) =>
      "{\"type\":\"terminal\",\"service\":\"text\",\"parameters\":{\"attribute\":\"exptl.method\",\"operator\":\"equals\",\"value\":\"" +
      method_ +
      "\"}}"
    BiotiteRcsbResolution(min, max) =>
      "{\"type\":\"terminal\",\"service\":\"text\",\"parameters\":{\"attribute\":\"rcsb_entry_info.resolution_combined\",\"operator\":\"range\",\"from\":" +
      min.to_string() +
      ",\"to\":" +
      max.to_string() +
      "}}"
    BiotiteRcsbMolecularWeight(min, max) =>
      "{\"type\":\"terminal\",\"service\":\"text\",\"parameters\":{\"attribute\":\"rcsb_entry_info.molecular_weight\",\"operator\":\"range\",\"from\":" +
      min.to_string() +
      ",\"to\":" +
      max.to_string() +
      "}}"
    BiotiteRcsbOrganism(name) =>
      "{\"type\":\"terminal\",\"service\":\"text\",\"parameters\":{\"attribute\":\"rcsb_entity_source_organism.ncbi_scientific_name\",\"operator\":\"contains\",\"value\":\"" +
      name +
      "\"}}"
    BiotiteRcsbLigand(name) =>
      "{\"type\":\"terminal\",\"service\":\"text\",\"parameters\":{\"attribute\":\"rcsb_nonpolymer_entity_container_identifiers.comp_id\",\"operator\":\"equals\",\"value\":\"" +
      name +
      "\"}}"
    BiotiteRcsbHasLigand(has) =>
      "{\"type\":\"terminal\",\"service\":\"text\",\"parameters\":{\"attribute\":\"rcsb_entry_info.has_nonpolymer_entities\",\"operator\":\"equals\",\"value\":" +
      has.to_string() +
      "}}"
    BiotiteRcsbDepositionDate(from, to) =>
      "{\"type\":\"terminal\",\"service\":\"text\",\"parameters\":{\"attribute\":\"rcsb_accession_info.deposit_date\",\"operator\":\"range\",\"from\":\"" +
      from +
      "\",\"to\":\"" +
      to +
      "\"}}"
    BiotiteRcsbReleaseDate(from, to) =>
      "{\"type\":\"terminal\",\"service\":\"text\",\"parameters\":{\"attribute\":\"rcsb_accession_info.initial_release_date\",\"operator\":\"range\",\"from\":\"" +
      from +
      "\",\"to\":\"" +
      to +
      "\"}}"
    BiotiteRcsbNumChains(min, max) =>
      "{\"type\":\"terminal\",\"service\":\"text\",\"parameters\":{\"attribute\":\"rcsb_entry_info.chain_count\",\"operator\":\"range\",\"from\":" +
      min.to_string() +
      ",\"to\":" +
      max.to_string() +
      "}}"
    BiotiteRcsbNumResidues(min, max) =>
      "{\"type\":\"terminal\",\"service\":\"text\",\"parameters\":{\"attribute\":\"rcsb_entry_info.residue_count\",\"operator\":\"range\",\"from\":" +
      min.to_string() +
      ",\"to\":" +
      max.to_string() +
      "}}"
    BiotiteRcsbSpaceGroup(sg) =>
      "{\"type\":\"terminal\",\"service\":\"text\",\"parameters\":{\"attribute\":\"rcsb_entry_info.space_group\",\"operator\":\"equals\",\"value\":\"" +
      sg +
      "\"}}"
    BiotiteRcsbPdbId(id) =>
      "{\"type\":\"terminal\",\"service\":\"text\",\"parameters\":{\"attribute\":\"rcsb_entry_container_identifiers.entry_id\",\"operator\":\"equals\",\"value\":\"" +
      id +
      "\"}}"
  }
}

///|
/// Convenience: search by experimental method.
pub fn biotite_rcsb_search_by_method(method_ : String) -> String {
  let search = BiotiteRcsbSearchQuery::new()
  let _ = search.add_condition(BiotiteRcsbExperimentalMethod(method_))
  search.build_url()
}

///|
/// Convenience: search by resolution range.
pub fn biotite_rcsb_search_by_resolution(
  min_res : Double,
  max_res : Double,
) -> String {
  let search = BiotiteRcsbSearchQuery::new()
  let _ = search.add_condition(BiotiteRcsbResolution(min_res, max_res))
  search.build_url()
}

///|
/// Use unused constructors.
fn _biotite_rcsb_use_constructors() -> Unit {
  ignore(BiotiteRcsbMolecularWeight(0.0, 0.0))
  ignore(BiotiteRcsbOrganism(""))
  ignore(BiotiteRcsbLigand(""))
  ignore(BiotiteRcsbHasLigand(true))
  ignore(BiotiteRcsbDepositionDate("", ""))
  ignore(BiotiteRcsbReleaseDate("", ""))
  ignore(BiotiteRcsbNumChains(0, 0))
  ignore(BiotiteRcsbNumResidues(0, 0))
  ignore(BiotiteRcsbSpaceGroup(""))
  ignore(BiotiteRcsbPdbId(""))
}
