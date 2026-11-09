// scvi_dataloaders.mbt - Mini-batch loaders and splits mirroring scvi.dataloaders
//
// Provides AnnDataLoader-style batching, train/validation/test splitting,
// semi-supervised labelled/unlabelled partitioning, masked-gene batches and
// concatenated multi-dataset loading.

// ===========================================================================
// AnnDataLoader
// ===========================================================================

///|
pub struct ScviDataLoader {
  data : Array[Array[Double]]
  batch_size : Int32
  shuffle : Bool
  seed : Int32
  order : Array[Int32]
  mut epoch : Int32
}

///|
pub fn scvi_dataloader_init(
  data : Array[Array[Double]],
  batch_size : Int32,
  shuffle : Bool,
  seed : Int32,
) -> ScviDataLoader {
  let bs = if batch_size < 1 { 1 } else { batch_size }
  let order = Array(capacity=data.length())
  for i in 0..<data.length() {
    order.push(i)
  }
  let loader : ScviDataLoader = {
    data,
    batch_size: bs,
    shuffle,
    seed,
    order,
    epoch: 0,
  }
  if shuffle {
    scvi_dataloader_shuffle(loader)
  }
  loader
}

///|
/// Fisher-Yates shuffle of the sample order (deterministic given `seed`).
pub fn scvi_dataloader_shuffle(loader : ScviDataLoader) -> Unit {
  loader.epoch = loader.epoch + 1
  let base = loader.seed + loader.epoch * 100003
  let n = loader.order.length()
  for i in 0..<n {
    let j = (scvi_dist_uniform(base + i * 131) * (n - i).to_double()).to_int()
    let jj = if j < 0 { 0 } else if j >= n - i { n - i - 1 } else { j }
    let a = loader.order[i]
    loader.order[i] = loader.order[i + jj]
    loader.order[i + jj] = a
  }
}

///|
pub fn scvi_dataloader_n_obs(obj : ScviDataLoader) -> Int32 {
  obj.data.length()
}

///|
pub fn scvi_dataloader_n_batches(obj : ScviDataLoader) -> Int32 {
  (obj.data.length() + obj.batch_size - 1) / obj.batch_size
}

///|
pub fn scvi_dataloader_batch_indices(
  obj : ScviDataLoader,
  batch_idx : Int32,
) -> Array[Int32] {
  let start = batch_idx * obj.batch_size
  let out = Array(capacity=obj.batch_size)
  let mut i = start
  while i < start + obj.batch_size && i < obj.order.length() {
    out.push(obj.order[i])
    i = i + 1
  }
  out
}

///|
pub fn scvi_dataloader_batch(
  obj : ScviDataLoader,
  batch_idx : Int32,
) -> Array[Array[Double]] {
  let idx = scvi_dataloader_batch_indices(obj, batch_idx)
  let out = Array(capacity=idx.length())
  for i in 0..<idx.length() {
    out.push(obj.data[idx[i]])
  }
  out
}

///|
/// Materialise all mini-batches for one epoch.
pub fn scvi_dataloader_all_batches(
  obj : ScviDataLoader,
) -> Array[Array[Array[Double]]] {
  let n = scvi_dataloader_n_batches(obj)
  let out = Array(capacity=n)
  for b in 0..<n {
    out.push(scvi_dataloader_batch(obj, b))
  }
  out
}

// ===========================================================================
// DataSplitter
// ===========================================================================

///|
pub struct ScviDataSplitterConfig {
  train_size : Double
  validation_size : Double
  test_size : Double
  shuffle : Bool
  seed : Int32
}

///|
pub fn ScviDataSplitterConfig::new(
  train_size : Double,
  validation_size : Double,
  test_size : Double,
  shuffle : Bool,
  seed : Int32,
) -> ScviDataSplitterConfig {
  { train_size, validation_size, test_size, shuffle, seed, }
}

///|
pub fn ScviDataSplitterConfig::default() -> ScviDataSplitterConfig {
  {
    train_size: 0.9,
    validation_size: 0.1,
    test_size: 0.0,
    shuffle: true,
    seed: 42,
  }
}

///|
pub struct ScviDataSplit {
  train_indices : Array[Int32]
  validation_indices : Array[Int32]
  test_indices : Array[Int32]
}

///|
pub fn scvi_data_splitter_split(
  n_obs : Int32,
  config : ScviDataSplitterConfig,
) -> ScviDataSplit {
  let order = Array(capacity=n_obs)
  for i in 0..<n_obs {
    order.push(i)
  }
  if config.shuffle {
    for i in 0..<n_obs {
      let j = (scvi_dist_uniform(config.seed + i * 97) * (n_obs - i).to_double()).to_int()
      let jj = if j < 0 {
        0
      } else if j >= n_obs - i {
        n_obs - i - 1
      } else {
        j
      }
      let a = order[i]
      order[i] = order[i + jj]
      order[i + jj] = a
    }
  }
  let n_train = (n_obs.to_double() * config.train_size).to_int()
  let n_val = (n_obs.to_double() * config.validation_size).to_int()
  let train = Array(capacity=n_train)
  let validation = Array(capacity=n_val)
  let test_split = Array(capacity=n_obs - n_train - n_val)
  for i in 0..<n_obs {
    if i < n_train {
      train.push(order[i])
    } else if i < n_train + n_val {
      validation.push(order[i])
    } else {
      test_split.push(order[i])
    }
  }
  {
    train_indices: train,
    validation_indices: validation,
    test_indices: test_split,
  }
}

///|
/// Split each label group independently so that every split keeps the class
/// composition of the full data set.
pub fn scvi_data_splitter_stratified(
  labels : Array[Int32],
  config : ScviDataSplitterConfig,
) -> ScviDataSplit {
  let n_obs = labels.length()
  let mut min_label = 0
  let mut max_label = 0
  for i in 0..<n_obs {
    if labels[i] < min_label {
      min_label = labels[i]
    }
    if labels[i] > max_label {
      max_label = labels[i]
    }
  }
  let train = Array(capacity=n_obs)
  let validation = Array(capacity=0)
  let test_split = Array(capacity=0)
  for lab in min_label..<=max_label {
    let group = Array(capacity=n_obs)
    for i in 0..<n_obs {
      if labels[i] == lab {
        group.push(i)
      }
    }
    let n_g = group.length()
    let n_tr = (n_g.to_double() * config.train_size).to_int()
    let n_va = (n_g.to_double() * config.validation_size).to_int()
    for i in 0..<n_g {
      if i < n_tr {
        train.push(group[i])
      } else if i < n_tr + n_va {
        validation.push(group[i])
      } else {
        test_split.push(group[i])
      }
    }
  }
  {
    train_indices: train,
    validation_indices: validation,
    test_indices: test_split,
  }
}

// ===========================================================================
// Semi-supervised loader
// ===========================================================================

///|
pub struct ScviSemiSupervisedSplit {
  labeled_indices : Array[Int32]
  unlabeled_indices : Array[Int32]
}

///|
/// Partition cells by label; `unlabeled_label` (usually -1) marks unlabelled
/// observations.
pub fn scvi_semisupervised_split(
  labels : Array[Int32],
  unlabeled_label : Int32,
) -> ScviSemiSupervisedSplit {
  let labeled = Array(capacity=labels.length())
  let unlabeled = Array(capacity=labels.length())
  for i in 0..<labels.length() {
    if labels[i] == unlabeled_label {
      unlabeled.push(i)
    } else {
      labeled.push(i)
    }
  }
  { labeled_indices: labeled, unlabeled_indices: unlabeled, }
}

///|
pub fn scvi_semisupervised_batch(
  data : Array[Array[Double]],
  indices : Array[Int32],
  batch_size : Int32,
  batch_idx : Int32,
) -> Array[Array[Double]] {
  let start = batch_idx * batch_size
  let out = Array(capacity=batch_size)
  let mut i = start
  while i < start + batch_size && i < indices.length() {
    out.push(data[indices[i]])
    i = i + 1
  }
  out
}

// ===========================================================================
// Masked-gene loader
// ===========================================================================

///|
pub struct ScviMaskedBatch {
  values : Array[Array[Double]]
  mask : Array[Array[Double]]
}

///|
/// Randomly mask a fraction of entries. `mask=1` keeps a value, `mask=0`
/// marks a masked entry (the value itself is zeroed).
pub fn scvi_masked_loader_mask(
  data : Array[Array[Double]],
  mask_prob : Double,
  seed : Int32,
) -> ScviMaskedBatch {
  let values = Array(capacity=data.length())
  let mask = Array(capacity=data.length())
  for i in 0..<data.length() {
    let row_v = Array(capacity=data[i].length())
    let row_m = Array(capacity=data[i].length())
    for j in 0..<data[i].length() {
      let u = scvi_dist_uniform(seed + i * 7919 + j * 131)
      if u < mask_prob {
        row_v.push(0.0)
        row_m.push(0.0)
      } else {
        row_v.push(data[i][j])
        row_m.push(1.0)
      }
    }
    values.push(row_v)
    mask.push(row_m)
  }
  { values, mask, }
}

///|
/// Apply a fixed mask to a subset of cells.
pub fn scvi_masked_loader_apply(
  data : Array[Array[Double]],
  indices : Array[Int32],
  mask_prob : Double,
  seed : Int32,
) -> ScviMaskedBatch {
  let subset = Array(capacity=indices.length())
  for i in 0..<indices.length() {
    subset.push(data[indices[i]])
  }
  scvi_masked_loader_mask(subset, mask_prob, seed)
}

// ===========================================================================
// ConcatenateLoader
// ===========================================================================

///|
pub struct ScviConcatenateLoader {
  datasets : Array[Array[Array[Double]]]
  batch_size : Int32
}

///|
pub fn scvi_concatenate_loader_init(
  datasets : Array[Array[Array[Double]]],
  batch_size : Int32,
) -> ScviConcatenateLoader {
  { datasets, batch_size: if batch_size < 1 { 1 } else { batch_size }, }
}

///|
pub fn scvi_concatenate_loader_n_obs(obj : ScviConcatenateLoader) -> Int32 {
  let mut n = 0
  for i in 0..<obj.datasets.length() {
    n = n + obj.datasets[i].length()
  }
  n
}

///|
pub fn scvi_concatenate_loader_n_batches(obj : ScviConcatenateLoader) -> Int32 {
  (scvi_concatenate_loader_n_obs(obj) + obj.batch_size - 1) / obj.batch_size
}

///|
pub fn scvi_concatenate_loader_flat(
  obj : ScviConcatenateLoader,
) -> Array[Array[Double]] {
  let out = Array(capacity=scvi_concatenate_loader_n_obs(obj))
  for d in 0..<obj.datasets.length() {
    for i in 0..<obj.datasets[d].length() {
      out.push(obj.datasets[d][i])
    }
  }
  out
}

///|
pub fn scvi_concatenate_loader_batch(
  obj : ScviConcatenateLoader,
  batch_idx : Int32,
) -> Array[Array[Double]] {
  let flat = scvi_concatenate_loader_flat(obj)
  let start = batch_idx * obj.batch_size
  let out = Array(capacity=obj.batch_size)
  let mut i = start
  while i < start + obj.batch_size && i < flat.length() {
    out.push(flat[i])
    i = i + 1
  }
  out
}
