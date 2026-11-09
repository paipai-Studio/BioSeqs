// scvi_settings.mbt - Global runtime settings and small utilities
// Mirrors scvi.settings and scvi.utils (attrdict, DEFAULT_KEYS, subsampling).

///|
pub struct ScviRuntimeSettings {
  mut seed : Int32
  mut verbosity : Int32
  mut dl_num_workers : Int32
  mut progress_bar : Bool
  mut num_threads : Int32
  mut save_dir : String
}

///|
pub fn scvi_runtime_settings_default() -> ScviRuntimeSettings {
  {
    seed: 0,
    verbosity: 1,
    dl_num_workers: 0,
    progress_bar: true,
    num_threads: 1,
    save_dir: "./scvi_model",
  }
}

///|
pub fn ScviRuntimeSettings::set_seed(
  self : ScviRuntimeSettings,
  seed : Int32,
) -> Unit {
  self.seed = seed
}

///|
pub fn ScviRuntimeSettings::set_verbosity(
  self : ScviRuntimeSettings,
  verbosity : Int32,
) -> Unit {
  self.verbosity = verbosity
}

///|
pub fn ScviRuntimeSettings::set_dl_num_workers(
  self : ScviRuntimeSettings,
  workers : Int32,
) -> Unit {
  self.dl_num_workers = workers
}

///|
pub fn ScviRuntimeSettings::set_progress_bar(
  self : ScviRuntimeSettings,
  enabled : Bool,
) -> Unit {
  self.progress_bar = enabled
}

///|
pub fn ScviRuntimeSettings::set_num_threads(
  self : ScviRuntimeSettings,
  threads : Int32,
) -> Unit {
  self.num_threads = threads
}

///|
pub fn ScviRuntimeSettings::set_save_dir(
  self : ScviRuntimeSettings,
  dir : String,
) -> Unit {
  self.save_dir = dir
}

///|
pub fn ScviRuntimeSettings::summary(self : ScviRuntimeSettings) -> String {
  "seed=" +
  self.seed.to_string() +
  ";verbosity=" +
  self.verbosity.to_string() +
  ";dl_num_workers=" +
  self.dl_num_workers.to_string() +
  ";num_threads=" +
  self.num_threads.to_string() +
  ";save_dir=" +
  self.save_dir
}

///|
/// Keys registered by scvi-tools `setup_anndata` implementations.
pub fn scvi_utils_default_keys() -> Array[String] {
  [
    "X", "batch", "labels", "size_factor", "protein_expression", "cat_covs", "cont_covs",
  ]
}

///|
/// Lightweight attrdict replacement: ordered (key, value) pairs.
pub fn scvi_utils_attrdict(
  pairs : Array[(String, String)],
) -> Array[(String, String)] {
  pairs
}

///|
pub fn scvi_utils_get(pairs : Array[(String, String)], key : String) -> String? {
  for i in 0..<pairs.length() {
    let (k, v) = pairs[i]
    if k == key {
      return Some(v)
    }
  }
  None
}

///|
/// Uniform subsample of `0..<n` keeping approximately `fraction` of rows.
pub fn scvi_utils_subsample_indices(
  n : Int32,
  fraction : Double,
  seed : Int32,
) -> Array[Int32] {
  let frac = if fraction < 0.0 {
    0.0
  } else if fraction > 1.0 {
    1.0
  } else {
    fraction
  }
  let out = Array(capacity=n)
  for i in 0..<n {
    if scvi_dist_uniform(seed + i * 7919) < frac {
      out.push(i)
    }
  }
  out
}

///|
/// Deterministic per-label subsample used by semi-supervised training.
pub fn scvi_utils_stratified_indices(
  labels : Array[Int32],
  per_label : Int32,
  seed : Int32,
) -> Array[Int32] {
  let mut min_label = 0
  let mut max_label = 0
  for i in 0..<labels.length() {
    if labels[i] < min_label {
      min_label = labels[i]
    }
    if labels[i] > max_label {
      max_label = labels[i]
    }
  }
  let out = Array(capacity=labels.length())
  for lab in min_label..<=max_label {
    let mut taken = 0
    for i in 0..<labels.length() {
      if labels[i] == lab && taken < per_label {
        out.push(i)
        taken = taken + 1
      }
    }
  }
  ignore(seed)
  out
}

///|
/// Format seconds as `HH:MM:SS`.
pub fn scvi_utils_format_duration(seconds : Int32) -> String {
  let total = if seconds < 0 { 0 } else { seconds }
  let h = total / 3600
  let m = total % 3600 / 60
  let s = total % 60
  let hh = if h < 10 { "0" + h.to_string() } else { h.to_string() }
  let mm = if m < 10 { "0" + m.to_string() } else { m.to_string() }
  let ss = if s < 10 { "0" + s.to_string() } else { s.to_string() }
  hh + ":" + mm + ":" + ss
}

///|
/// Classification accuracy between predictions and ground truth.
pub fn scvi_utils_accuracy(
  predictions : Array[Int32],
  truth : Array[Int32],
) -> Double {
  let n = if predictions.length() < truth.length() {
    predictions.length()
  } else {
    truth.length()
  }
  if n == 0 {
    return 0.0
  }
  let mut correct = 0
  for i in 0..<n {
    if predictions[i] == truth[i] {
      correct = correct + 1
    }
  }
  correct.to_double() / n.to_double()
}
