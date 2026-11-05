///|
/// Other tools: AnnTorchDataset, ModelHub, Reticulate, Settings, Logging.
pub struct AnnTorchDataset {
  n_obs : Int32
  n_vars : Int32
  data : Array[Array[Double]]
  obs_indices : Array[Int32]
  var_indices : Array[Int32]
  is_sparse : Bool
}

///|
pub fn ann_torch_dataset_init(data : Array[Array[Double]]) -> AnnTorchDataset {
  let n_obs = data.length()
  let n_vars = if n_obs > 0 { data[0].length() } else { 0 }
  let obs_idx = Array(capacity=n_obs)
  for i in 0..<n_obs {
    obs_idx.push(i)
  }
  let var_idx = Array(capacity=n_vars)
  for i in 0..<n_vars {
    var_idx.push(i)
  }
  {
    n_obs,
    n_vars,
    data,
    obs_indices: obs_idx,
    var_indices: var_idx,
    is_sparse: false,
  }
}

///|
pub fn ann_torch_dataset_getitem(
  dataset : AnnTorchDataset,
  idx : Int32,
) -> Array[Double] {
  if idx >= 0 && idx < dataset.data.length() {
    dataset.data[idx]
  } else {
    []
  }
}

///|
pub fn ann_torch_dataset_subset(
  dataset : AnnTorchDataset,
  obs_mask : Array[Int32],
) -> AnnTorchDataset {
  let subset_data = Array(capacity=obs_mask.length())
  for i in 0..<obs_mask.length() {
    let idx = obs_mask[i]
    if idx >= 0 && idx < dataset.data.length() {
      subset_data.push(dataset.data[idx])
    }
  }
  {
    n_obs: obs_mask.length(),
    n_vars: dataset.n_vars,
    data: subset_data,
    obs_indices: obs_mask,
    var_indices: dataset.var_indices,
    is_sparse: dataset.is_sparse,
  }
}

///|
pub fn ann_torch_dataset_batch(
  dataset : AnnTorchDataset,
  indices : Array[Int32],
) -> Array[Array[Double]] {
  let batch = Array(capacity=indices.length())
  for i in 0..<indices.length() {
    batch.push(ann_torch_dataset_getitem(dataset, indices[i]))
  }
  batch
}

///|
pub fn ann_torch_dataset_to_sparse(
  dataset : AnnTorchDataset,
) -> Array[(Int32, Int32, Double)] {
  let sparse = Array(capacity=0)
  for i in 0..<dataset.n_obs {
    for j in 0..<dataset.n_vars {
      if i < dataset.data.length() && j < dataset.data[i].length() {
        if dataset.data[i][j] != 0.0 {
          sparse.push((i, j, dataset.data[i][j]))
        }
      }
    }
  }
  sparse
}

///|
pub struct ModelHubEntry {
  name : String
  author : String
  url : String
  n_genes : Int32
  description : String
  tags : Array[String]
}

///|
pub fn model_hub_search(
  entries : Array[ModelHubEntry],
  query : String,
) -> Array[ModelHubEntry] {
  let results = Array(capacity=0)
  for i in 0..<entries.length() {
    let entry = entries[i]
    if entry.name.contains(query) || entry.description.contains(query) {
      results.push(entry)
    } else {
      for tag in entry.tags {
        if tag.contains(query) {
          results.push(entry)
          break
        }
      }
    }
  }
  results
}

///|
pub fn model_hub_format_entry(entry : ModelHubEntry) -> String {
  let sb = StringBuilder()
  let _ = sb.write("Name: " + entry.name + "\n")
  let _ = sb.write("Author: " + entry.author + "\n")
  let _ = sb.write("URL: " + entry.url + "\n")
  let _ = sb.write("Genes: " + entry.n_genes.to_string() + "\n")
  let _ = sb.write("Description: " + entry.description + "\n")
  let _ = sb.write("Tags: ")
  for i in 0..<entry.tags.length() {
    let _ = sb.write(entry.tags[i])
    if i < entry.tags.length() - 1 {
      let _ = sb.write(", ")
    }
  }
  let _ = sb.write("\n")
  sb.to_string()
}

///|
pub struct ReticulateConfig {
  r_home : String
  python_path : String
  is_connected : Bool
}

///|
pub fn reticulate_init(
  r_home : String,
  python_path : String,
) -> ReticulateConfig {
  { r_home, python_path, is_connected: true, }
}

///|
pub fn reticulate_call_r(
  config : ReticulateConfig,
  function_name : String,
  args : Array[String],
) -> String {
  if not(config.is_connected) {
    return "Error: R not connected"
  }
  let sb = StringBuilder()
  let _ = sb.write(function_name + "(")
  for i in 0..<args.length() {
    let _ = sb.write(args[i])
    if i < args.length() - 1 {
      let _ = sb.write(", ")
    }
  }
  let _ = sb.write(")")
  sb.to_string()
}

///|
pub fn reticulate_convert_anndata_to_sce(
  config : ReticulateConfig,
  n_obs : Int32,
  n_vars : Int32,
) -> String {
  if not(config.is_connected) {
    return "Error: R not connected"
  }
  "SingleCellExperiment(assays=list(counts=matrix(0, " +
  n_vars.to_string() +
  ", " +
  n_obs.to_string() +
  ")))"
}

///|
pub struct ScviSettings {
  seed : Int32
  logging_level : String
  n_threads : Int32
  progress_bar : Bool
  fig_format : String
  fig_dpi : Int32
}

///|
pub fn ScviSettings::default() -> ScviSettings {
  {
    seed: 0,
    logging_level: "INFO",
    n_threads: 1,
    progress_bar: true,
    fig_format: "png",
    fig_dpi: 100,
  }
}

///|
pub fn scvi_settings_set_seed(
  settings : ScviSettings,
  seed : Int32,
) -> ScviSettings {
  { ..settings, seed, }
}

///|
pub fn scvi_settings_set_logging(
  settings : ScviSettings,
  level : String,
) -> ScviSettings {
  { ..settings, logging_level: level, }
}

///|
pub fn scvi_settings_set_threads(
  settings : ScviSettings,
  n : Int32,
) -> ScviSettings {
  { ..settings, n_threads: n, }
}

///|
pub struct ScviLogger {
  level : String
  messages : Array[(String, String)]
}

///|
pub fn ScviLogger::new(level : String) -> ScviLogger {
  { level, messages: [], }
}

///|
pub fn scvi_logger_log(
  logger : ScviLogger,
  level : String,
  msg : String,
) -> ScviLogger {
  let l = logger
  let level_priority = match level {
    "DEBUG" => 0
    "INFO" => 1
    "WARNING" => 2
    "ERROR" => 3
    _ => 1
  }
  let threshold = match l.level {
    "DEBUG" => 0
    "INFO" => 1
    "WARNING" => 2
    "ERROR" => 3
    _ => 1
  }
  if level_priority >= threshold {
    l.messages.push((level, msg))
  }
  l
}

///|
pub fn scvi_logger_info(logger : ScviLogger, msg : String) -> ScviLogger {
  scvi_logger_log(logger, "INFO", msg)
}

///|
pub fn scvi_logger_warning(logger : ScviLogger, msg : String) -> ScviLogger {
  scvi_logger_log(logger, "WARNING", msg)
}

///|
pub fn scvi_logger_error(logger : ScviLogger, msg : String) -> ScviLogger {
  scvi_logger_log(logger, "ERROR", msg)
}

///|
pub fn scvi_logger_get_messages(logger : ScviLogger) -> Array[(String, String)] {
  logger.messages
}

///|
pub fn scvi_logger_format(logger : ScviLogger) -> String {
  let sb = StringBuilder()
  for i in 0..<logger.messages.length() {
    let (level, msg) = logger.messages[i]
    let _ = sb.write("[" + level + "] " + msg + "\n")
  }
  sb.to_string()
}

///|
pub fn _scvi_tools_dummy() -> Unit {
  let _ = AnnTorchDataset::{
    n_obs: 0,
    n_vars: 0,
    data: [],
    obs_indices: [],
    var_indices: [],
    is_sparse: false,
  }
  let _ = ModelHubEntry::{
    name: "",
    author: "",
    url: "",
    n_genes: 0,
    description: "",
    tags: [],
  }
  let _ = ReticulateConfig::{
    r_home: "",
    python_path: "",
    is_connected: false,
  }
  let _ = ScviSettings::{
    seed: 0,
    logging_level: "",
    n_threads: 0,
    progress_bar: false,
    fig_format: "",
    fig_dpi: 0,
  }
  let _ = ScviLogger::{ level: "", messages: [], }
}
