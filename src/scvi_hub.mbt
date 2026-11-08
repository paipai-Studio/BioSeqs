// scvi_hub.mbt - Model hub metadata mirroring scvi.hub
//
// Provides HubMetadata / HubModel-style model cards, validation, JSON
// serialisation and local cache paths. Network transfer is intentionally kept
// out of the library; callers can feed the generated URL to their own client.

///|
pub struct ScviHubMetadata {
  model_name : String
  version : String
  description : String
  url : String
  license : String
  scvi_version : String
  anndata_version : String
  n_obs : Int32
  n_vars : Int32
}

///|
pub fn ScviHubMetadata::new(
  model_name : String,
  version : String,
  description : String,
  url : String,
  license : String,
  scvi_version : String,
  anndata_version : String,
  n_obs : Int32,
  n_vars : Int32,
) -> ScviHubMetadata {
  {
    model_name,
    version,
    description,
    url,
    license,
    scvi_version,
    anndata_version,
    n_obs,
    n_vars,
  }
}

///|
pub fn ScviHubMetadata::empty() -> ScviHubMetadata {
  ScviHubMetadata::new(
    "unnamed", "0.0.0", "", "", "unknown", "0.0.0", "0.0.0", 0, 0,
  )
}

///|
pub fn scvi_hub_validate(meta : ScviHubMetadata) -> Bool {
  if meta.model_name == "" || meta.version == "" {
    return false
  }
  if meta.url == "" && meta.description == "" {
    return false
  }
  if meta.n_obs < 0 || meta.n_vars < 0 {
    return false
  }
  true
}

///|
/// Flat key/value model card, suitable for hf-hub style metadata.
pub fn scvi_hub_model_card(meta : ScviHubMetadata) -> Array[(String, String)] {
  [
    ("model_name", meta.model_name),
    ("version", meta.version),
    ("description", meta.description),
    ("url", meta.url),
    ("license", meta.license),
    ("scvi_version", meta.scvi_version),
    ("anndata_version", meta.anndata_version),
    ("n_obs", meta.n_obs.to_string()),
    ("n_vars", meta.n_vars.to_string()),
  ]
}

///|
pub fn scvi_hub_to_json(meta : ScviHubMetadata) -> String {
  "{\"model_name\":\"" +
  meta.model_name +
  "\",\"version\":\"" +
  meta.version +
  "\",\"description\":\"" +
  meta.description +
  "\",\"url\":\"" +
  meta.url +
  "\",\"license\":\"" +
  meta.license +
  "\",\"scvi_version\":\"" +
  meta.scvi_version +
  "\",\"anndata_version\":\"" +
  meta.anndata_version +
  "\",\"n_obs\":" +
  meta.n_obs.to_string() +
  ",\"n_vars\":" +
  meta.n_vars.to_string() +
  "}"
}

///|
fn scvi_hub_parse_int(text : String) -> Int32 {
  let mut value = 0
  for i in 0..<text.length() {
    let c = text.unsafe_get(i)
    let d = c.to_int() - '0'.to_int()
    if d >= 0 && d <= 9 {
      value = value * 10 + d
    }
  }
  value
}

///|
/// Rebuild metadata from a flat key/value card (missing fields keep defaults).
pub fn scvi_hub_from_pairs(pairs : Array[(String, String)]) -> ScviHubMetadata {
  let base = ScviHubMetadata::empty()
  let mut model_name = base.model_name
  let mut version = base.version
  let mut description = base.description
  let mut url = base.url
  let mut license = base.license
  let mut scvi_version = base.scvi_version
  let mut anndata_version = base.anndata_version
  let mut n_obs = base.n_obs
  let mut n_vars = base.n_vars
  for i in 0..<pairs.length() {
    let (k, v) = pairs[i]
    if k == "model_name" {
      model_name = v
    } else if k == "version" {
      version = v
    } else if k == "description" {
      description = v
    } else if k == "url" {
      url = v
    } else if k == "license" {
      license = v
    } else if k == "scvi_version" {
      scvi_version = v
    } else if k == "anndata_version" {
      anndata_version = v
    } else if k == "n_obs" {
      n_obs = scvi_hub_parse_int(v)
    } else if k == "n_vars" {
      n_vars = scvi_hub_parse_int(v)
    }
  }
  ScviHubMetadata::new(
    model_name, version, description, url, license, scvi_version, anndata_version,
    n_obs, n_vars,
  )
}

///|
/// Local cache path for a hub model.
pub fn scvi_hub_cache_path(
  cache_dir : String,
  model_name : String,
  version : String,
) -> String {
  let dir = if cache_dir.has_suffix("/") { cache_dir } else { cache_dir + "/" }
  dir + model_name + "/" + version + "/model.pt"
}

///|
/// Download URL for a hub model.
pub fn scvi_hub_download_url(
  base_url : String,
  model_name : String,
  version : String,
) -> String {
  let base = if base_url.has_suffix("/") { base_url } else { base_url + "/" }
  base + model_name + "/resolve/" + version + "/model.pt"
}

///|
pub fn scvi_hub_matches(meta : ScviHubMetadata, query : String) -> Bool {
  meta.model_name.contains(query) || meta.description.contains(query)
}

///|
pub fn scvi_hub_summary(meta : ScviHubMetadata) -> String {
  meta.model_name +
  "@" +
  meta.version +
  " (" +
  meta.n_obs.to_string() +
  "x" +
  meta.n_vars.to_string() +
  ", scvi " +
  meta.scvi_version +
  ", anndata " +
  meta.anndata_version +
  ")"
}
