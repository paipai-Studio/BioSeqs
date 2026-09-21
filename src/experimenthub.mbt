///|
/// ExperimentHub (Bioconductor data hub) core: resource metadata plus the
/// query/filter layer over the resource table.

///|
/// Metadata for one ExperimentHub resource.
pub struct EhResource {
  title : String
  species : String
  taxonomy_id : Int
  genome : String
  description : String
  dataprovider : String
}

///|
pub fn EhResource::new(
  title : String,
  species : String,
  taxonomy_id : Int,
  genome : String,
  description : String,
  dataprovider : String,
) -> EhResource {
  EhResource::{ title, species, taxonomy_id, genome, description, dataprovider }
}

///|
/// Indices of resources whose species field contains `species`
/// (case-insensitive substring, matching ExperimentHub's `query` semantics).
pub fn eh_filter_species(resources : Array[EhResource], species : String) -> Array[Int] {
  let key = to_lowercase(species)
  let out : Array[Int] = Array::new()
  for i = 0; i < resources.length(); i = i + 1 {
    if contains(to_lowercase(resources[i].species), key) {
      out.push(i)
    }
  }
  out
}

///|
/// Indices of resources with an exact matching `taxonomy_id`.
pub fn eh_filter_taxonomy(resources : Array[EhResource], taxid : Int) -> Array[Int] {
  let out : Array[Int] = Array::new()
  for i = 0; i < resources.length(); i = i + 1 {
    if resources[i].taxonomy_id == taxid {
      out.push(i)
    }
  }
  out
}

///|
/// Indices of resources matching `keyword` (case-insensitive substring) in any
/// of title / species / genome / description / dataprovider.
pub fn eh_search(resources : Array[EhResource], keyword : String) -> Array[Int] {
  let key = to_lowercase(keyword)
  let out : Array[Int] = Array::new()
  for i = 0; i < resources.length(); i = i + 1 {
    let r = resources[i]
    let hay = to_lowercase(
      r.title + " " + r.species + " " + r.genome + " " + r.description + " " + r.dataprovider,
    )
    if contains(hay, key) {
      out.push(i)
    }
  }
  out
}

///|
/// Default ExperimentHub cache root (`<home>/.cache/ExperimentHub`), matching
/// the BiocFileCache default used by ExperimentHub.
pub fn eh_cache_root(home : String) -> String {
  home + "/.cache/ExperimentHub"
}