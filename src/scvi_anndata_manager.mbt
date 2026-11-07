// scvi_anndata_manager.mbt - AnnData field registry and validation
// Inspired by scvi-tools AnnDataManager and field system
// Manages registration of AnnData fields for model setup

///|
pub struct AnnDataField {
  registry_key : String
  attr_key : String
  attr_name : String
  is_count_data : Bool
  is_categorical : Bool
  is_required : Bool
  n_cats : Int32
}

///|
pub struct AnnDataManager {
  fields : Array[AnnDataField]
  setup_args : Map[String, String]
  mut validated : Bool
  n_vars : Int32
  n_obs : Int32
  mut n_batch : Int32
  mut n_labels : Int32
}

///|
pub fn AnnDataManager::new() -> AnnDataManager {
  {
    fields: Array(capacity=0),
    setup_args: Map([], capacity=0),
    validated: false,
    n_vars: 0,
    n_obs: 0,
    n_batch: 1,
    n_labels: 1,
  }
}

///|
/// Register a layer field (count data)
pub fn anndata_register_layer(
  manager : AnnDataManager,
  layer_key : String,
  is_count_data : Bool,
) -> AnnDataManager {
  let m = manager
  m.fields.push({
    registry_key: "X",
    attr_key: layer_key,
    attr_name: "layer",
    is_count_data,
    is_categorical: false,
    is_required: true,
    n_cats: 0,
  })
  m
}

///|
/// Register a categorical obs field (e.g. batch)
pub fn anndata_register_categorical(
  manager : AnnDataManager,
  registry_key : String,
  obs_key : String,
  n_cats : Int32,
) -> AnnDataManager {
  let m = manager
  m.fields.push({
    registry_key,
    attr_key: obs_key,
    attr_name: "obs",
    is_count_data: false,
    is_categorical: true,
    is_required: true,
    n_cats,
  })
  if registry_key == "batch" {
    m.n_batch = n_cats
  }
  if registry_key == "labels" {
    m.n_labels = n_cats
  }
  m
}

///|
/// Register a numerical obs field (e.g. size factor)
pub fn anndata_register_numerical(
  manager : AnnDataManager,
  registry_key : String,
  obs_key : String,
  is_required : Bool,
) -> AnnDataManager {
  let m = manager
  m.fields.push({
    registry_key,
    attr_key: obs_key,
    attr_name: "obs",
    is_count_data: false,
    is_categorical: false,
    is_required,
    n_cats: 0,
  })
  m
}

///|
/// Register continuous covariate keys
pub fn anndata_register_continuous_covs(
  manager : AnnDataManager,
  cov_keys : Array[String],
) -> AnnDataManager {
  let m = manager
  for i in 0..<cov_keys.length() {
    m.fields.push({
      registry_key: "continuous_covs",
      attr_key: cov_keys[i],
      attr_name: "obs",
      is_count_data: false,
      is_categorical: false,
      is_required: true,
      n_cats: 0,
    })
  }
  m
}

///|
/// Validate that all required fields are present
pub fn scvi_anndata_validate_field(
  manager : AnnDataManager,
  available_keys : Array[String],
) -> (AnnDataManager, Array[String]) {
  let m = manager
  let errors = Array(capacity=0)
  for i in 0..<m.fields.length() {
    let field = m.fields[i]
    if field.is_required {
      let mut found = false
      for j in 0..<available_keys.length() {
        if available_keys[j] == field.attr_key {
          found = true
          break
        }
      }
      if !found {
        errors.push(
          "Missing required field: {field.attr_key} (registry: {field.registry_key})",
        )
      }
    }
  }
  if errors.length() == 0 {
    m.validated = true
  }
  (m, errors)
}

///|
/// Get field by registry key
pub fn anndata_get_field(
  manager : AnnDataManager,
  registry_key : String,
) -> AnnDataField? {
  for i in 0..<manager.fields.length() {
    if manager.fields[i].registry_key == registry_key {
      return Some(manager.fields[i])
    }
  }
  None
}

///|
/// Get summary stats
pub fn anndata_summary_stats(manager : AnnDataManager) -> String {
  let _n_fields = manager.fields.length()
  "AnnDataManager Summary:\n  Fields: {n_fields}\n  n_obs: {manager.n_obs}\n  n_vars: {manager.n_vars}\n  n_batch: {manager.n_batch}\n  n_labels: {manager.n_labels}\n  Validated: {manager.validated}"
}

///|
/// Create registry from setup args
pub fn anndata_setup_from_args(
  batch_key : String?,
  labels_key : String?,
  layer_key : String?,
  size_factor_key : String?,
  categorical_covariate_keys : Array[String]?,
  continuous_covariate_keys : Array[String]?,
) -> AnnDataManager {
  let mut manager = AnnDataManager::new()
  // Layer field
  match layer_key {
    Some(k) => manager = anndata_register_layer(manager, k, true)
    None => manager = anndata_register_layer(manager, "X", true)
  }
  // Batch field
  match batch_key {
    Some(k) => manager = anndata_register_categorical(manager, "batch", k, 1)
    None => ()
  }
  // Labels field
  match labels_key {
    Some(k) => manager = anndata_register_categorical(manager, "labels", k, 1)
    None => ()
  }
  // Size factor
  match size_factor_key {
    Some(k) =>
      manager = anndata_register_numerical(manager, "size_factor", k, false)
    None => ()
  }
  // Categorical covariates
  match categorical_covariate_keys {
    Some(keys) =>
      for i in 0..<keys.length() {
        manager = anndata_register_categorical(
          manager,
          "cat_cov_{i}",
          keys[i],
          1,
        )
      }
    None => ()
  }
  // Continuous covariates
  match continuous_covariate_keys {
    Some(keys) => manager = anndata_register_continuous_covs(manager, keys)
    None => ()
  }
  manager
}
