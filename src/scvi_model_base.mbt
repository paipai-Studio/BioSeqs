// scvi_model_base.mbt - Shared model bookkeeping mirroring scvi.model.base
//
// MoonBit has no class inheritance, so the BaseModelClass / mixin surface is
// provided as a small state object plus free functions that every model can
// embed or delegate to.

///|
pub struct ScviModelBase {
  n_obs : Int32
  n_vars : Int32
  n_batch : Int32
  mut is_trained : Bool
  mut train_indices : Array[Int32]
  mut validation_indices : Array[Int32]
  mut test_indices : Array[Int32]
}

///|
pub fn scvi_model_base_init(
  n_obs : Int32,
  n_vars : Int32,
  n_batch : Int32,
) -> ScviModelBase {
  {
    n_obs,
    n_vars,
    n_batch,
    is_trained: false,
    train_indices: Array(capacity=0),
    validation_indices: Array(capacity=0),
    test_indices: Array(capacity=0),
  }
}

///|
pub fn scvi_model_base_from_registry(
  registry : ScviDataRegistry,
) -> ScviModelBase {
  scvi_model_base_init(
    registry.n_obs,
    registry.n_vars,
    registry.batch_categories.length(),
  )
}

///|
pub fn scvi_model_base_set_split(
  obj : ScviModelBase,
  split : ScviDataSplit,
) -> Unit {
  obj.train_indices = split.train_indices
  obj.validation_indices = split.validation_indices
  obj.test_indices = split.test_indices
}

///|
pub fn scvi_model_base_split(
  obj : ScviModelBase,
  config : ScviDataSplitterConfig,
) -> Unit {
  let split = scvi_data_splitter_split(obj.n_obs, config)
  scvi_model_base_set_split(obj, split)
}

///|
pub fn scvi_model_base_mark_trained(obj : ScviModelBase) -> Unit {
  obj.is_trained = true
}

///|
pub fn scvi_model_base_is_trained(obj : ScviModelBase) -> Bool {
  obj.is_trained
}

///|
pub fn scvi_model_base_n_train(obj : ScviModelBase) -> Int32 {
  if obj.train_indices.length() > 0 {
    obj.train_indices.length()
  } else {
    obj.n_obs
  }
}

///|
pub fn scvi_model_base_n_validation(obj : ScviModelBase) -> Int32 {
  obj.validation_indices.length()
}

///|
pub fn scvi_model_base_n_test(obj : ScviModelBase) -> Int32 {
  obj.test_indices.length()
}

///|
/// Validate shape and split consistency.
pub fn scvi_model_base_validate(obj : ScviModelBase) -> Bool {
  if obj.n_obs <= 0 || obj.n_vars <= 0 || obj.n_batch <= 0 {
    return false
  }
  let seen = Array::make(obj.n_obs, 0)
  for i in 0..<obj.train_indices.length() {
    let idx = obj.train_indices[i]
    if idx < 0 || idx >= obj.n_obs {
      return false
    }
    seen[idx] = seen[idx] + 1
  }
  for i in 0..<obj.validation_indices.length() {
    let idx = obj.validation_indices[i]
    if idx < 0 || idx >= obj.n_obs {
      return false
    }
    seen[idx] = seen[idx] + 1
  }
  for i in 0..<obj.test_indices.length() {
    let idx = obj.test_indices[i]
    if idx < 0 || idx >= obj.n_obs {
      return false
    }
    seen[idx] = seen[idx] + 1
  }
  for i in 0..<seen.length() {
    if seen[i] > 1 {
      return false
    }
  }
  true
}

///|
pub fn scvi_model_base_summary(obj : ScviModelBase) -> String {
  "ScviModelBase(n_obs=" +
  obj.n_obs.to_string() +
  ", n_vars=" +
  obj.n_vars.to_string() +
  ", n_batch=" +
  obj.n_batch.to_string() +
  ", trained=" +
  obj.is_trained.to_string() +
  ")"
}

///|
/// Check a latent matrix against the declared observation count.
pub fn scvi_model_base_check_latent(
  obj : ScviModelBase,
  latent : Array[Array[Double]],
) -> Bool {
  latent.length() == obj.n_obs
}

///|
/// scvi-tools `is_trained` guard used before inference.
pub fn scvi_model_base_require_trained(obj : ScviModelBase) -> Bool {
  obj.is_trained
}
