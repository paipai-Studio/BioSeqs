///|
/// Model lifecycle: save, load, convert_legacy_save, minify_adata,
/// setup_annbatch, deregister_manager, transfer_fields.
pub struct ScviModelMetadata {
  model_name : String
  version : String
  n_genes : Int32
  n_latent : Int32
  n_batches : Int32
  is_trained : Bool
  init_params : Map[String, String]
  registry : Map[String, String]
}

///|
pub fn scvi_model_save(
  path : String,
  model_name : String,
  encoder_weights : Array[Array[Double]],
  decoder_weights : Array[Array[Double]],
  metadata : ScviModelMetadata,
) -> String {
  let sb = StringBuilder()
  let _ = sb.write("--- scvi-tools model save ---\n")
  let _ = sb.write("path: " + path + "\n")
  let _ = sb.write("model: " + model_name + "\n")
  let _ = sb.write("version: " + metadata.version + "\n")
  let _ = sb.write("n_genes: " + metadata.n_genes.to_string() + "\n")
  let _ = sb.write("n_latent: " + metadata.n_latent.to_string() + "\n")
  let _ = sb.write("n_batches: " + metadata.n_batches.to_string() + "\n")
  let _ = sb.write("is_trained: " + metadata.is_trained.to_string() + "\n")
  let _ = sb.write("--- encoder weights ---\n")
  for i in 0..<encoder_weights.length() {
    for j in 0..<encoder_weights[i].length() {
      let _ = sb.write(encoder_weights[i][j].to_string())
      if j < encoder_weights[i].length() - 1 {
        let _ = sb.write(",")
      }
    }
    let _ = sb.write("\n")
  }
  let _ = sb.write("--- decoder weights ---\n")
  for i in 0..<decoder_weights.length() {
    for j in 0..<decoder_weights[i].length() {
      let _ = sb.write(decoder_weights[i][j].to_string())
      if j < decoder_weights[i].length() - 1 {
        let _ = sb.write(",")
      }
    }
    let _ = sb.write("\n")
  }
  sb.to_string()
}

///|
pub fn scvi_model_load(text : String) -> (String, ScviModelMetadata) {
  let lines = text.split("\n").to_array()
  let mut model_name = ""
  let mut version = ""
  let mut n_genes = 0
  let mut n_latent = 0
  let mut n_batches = 0
  let mut is_trained = false
  for i in 0..<lines.length() {
    let line = lines[i].to_owned()
    if line.has_prefix("model: ") {
      model_name = line[7:].to_owned()
    } else if line.has_prefix("version: ") {
      version = line[9:].to_owned()
    } else if line.has_prefix("n_genes: ") {
      n_genes = parse_int(line[9:].to_owned())
    } else if line.has_prefix("n_latent: ") {
      n_latent = parse_int(line[10:].to_owned())
    } else if line.has_prefix("n_batches: ") {
      n_batches = parse_int(line[11:].to_owned())
    } else if line.has_prefix("is_trained: ") {
      is_trained = line[12:].to_owned() == "true"
    }
  }
  let meta = {
    model_name,
    version,
    n_genes,
    n_latent,
    n_batches,
    is_trained,
    init_params: Map([]),
    registry: Map([]),
  }
  (model_name, meta)
}

///|
pub fn scvi_convert_legacy_save(legacy_text : String) -> String {
  let lines = legacy_text.split("\n").to_array()
  let sb = StringBuilder()
  let _ = sb.write("--- scvi-tools model save ---\n")
  for i in 0..<lines.length() {
    let line = lines[i].to_owned()
    if line.has_prefix("model_name:") {
      let _ = sb.write("model: " + line[11:].trim().to_owned() + "\n")
    } else if line.has_prefix("trained:") {
      let _ = sb.write("is_trained: " + line[8:].trim().to_owned() + "\n")
    } else if !line.has_prefix("#") && line.length() > 0 {
      let _ = sb.write(line + "\n")
    }
  }
  sb.to_string()
}

///|
pub fn scvi_minify_adata(
  latent : Array[Array[Double]],
  library_sizes : Array[Double],
  obs_names : Array[String],
) -> Map[String, Array[Double]] {
  let minified = Map([])
  for i in 0..<latent.length() {
    let key = if i < obs_names.length() {
      obs_names[i]
    } else {
      "cell_" + i.to_string()
    }
    let values = Array(capacity=latent[i].length() + 1)
    for j in 0..<latent[i].length() {
      values.push(latent[i][j])
    }
    let lib = if i < library_sizes.length() { library_sizes[i] } else { 0.0 }
    values.push(lib)
    minified[key] = values
  }
  minified
}

///|
pub struct AnnBatchConfig {
  batch_size : Int32
  n_obs : Int32
  shuffle : Bool
  drop_last : Bool
  mut current_idx : Int32
  indices : Array[Int32]
}

///|
pub fn scvi_setup_annbatch(
  n_obs : Int32,
  batch_size : Int32,
  shuffle : Bool,
) -> AnnBatchConfig {
  let indices = Array(capacity=n_obs)
  for i in 0..<n_obs {
    indices.push(i)
  }
  if shuffle {
    let mut rng = 42
    for i in 0..<indices.length() {
      rng = (rng * 1103515245 + 12345) % 2147483647
      let j = rng % indices.length()
      if j >= 0 {
        let tmp = indices[i]
        indices[i] = indices[j]
        indices[j] = tmp
      }
    }
  }
  { batch_size, n_obs, shuffle, drop_last: false, current_idx: 0, indices, }
}

///|
pub fn annbatch_next(config : AnnBatchConfig) -> (AnnBatchConfig, Array[Int32]) {
  let c = config
  if c.current_idx >= c.n_obs {
    return (c, [])
  }
  let end_idx = if c.current_idx + c.batch_size < c.n_obs {
    c.current_idx + c.batch_size
  } else if c.drop_last {
    c.current_idx = c.n_obs
    return (c, [])
  } else {
    c.n_obs
  }
  let batch = Array(capacity=end_idx - c.current_idx)
  for i in c.current_idx..<end_idx {
    if i < c.indices.length() {
      batch.push(c.indices[i])
    }
  }
  c.current_idx = end_idx
  (c, batch)
}

///|
pub fn annbatch_reset(config : AnnBatchConfig) -> AnnBatchConfig {
  let c = config
  c.current_idx = 0
  if c.shuffle {
    let mut rng = 42
    for i in 0..<c.indices.length() {
      rng = (rng * 1103515245 + 12345) % 2147483647
      let j = rng % c.indices.length()
      if j >= 0 {
        let tmp = c.indices[i]
        c.indices[i] = c.indices[j]
        c.indices[j] = tmp
      }
    }
  }
  c
}

///|
pub fn scvi_deregister_manager(
  registry : Map[String, String],
  keys_to_remove : Array[String],
) -> Map[String, String] {
  let r = registry
  for i in 0..<keys_to_remove.length() {
    let key = keys_to_remove[i]
    r.remove(key)
  }
  r
}

///|
pub fn scvi_transfer_fields(
  source_registry : Map[String, String],
  target_keys : Array[String],
) -> Map[String, String] {
  let transferred = Map([])
  for i in 0..<target_keys.length() {
    let key = target_keys[i]
    match source_registry.get(key) {
      Some(val) => transferred[key] = val
      None => ()
    }
  }
  transferred
}

///|
pub fn _scvi_model_io_dummy() -> Unit {
  let _ = ScviModelMetadata::{
    model_name: "",
    version: "",
    n_genes: 0,
    n_latent: 0,
    n_batches: 0,
    is_trained: false,
    init_params: Map([]),
    registry: Map([]),
  }
  let _ = AnnBatchConfig::{
    batch_size: 0,
    n_obs: 0,
    shuffle: false,
    drop_last: false,
    current_idx: 0,
    indices: [],
  }
}
