/// AnnData chunked iteration.
///
/// `for chunk in adata.chunked(0, 1000):` equivalent.
///
/// Iterates over batches of obs or var for memory-efficient processing
/// of large datasets.

///|
/// A chunk of AnnData data along one axis.
pub struct AnnDataChunk {
  axis : Int32 // 0 = obs, 1 = var
  start : Int32 // start index
  end_ : Int32 // end index (exclusive)
  data : AnnData // subset AnnData for this chunk
}

///|
pub fn AnnDataChunk::start(self : AnnDataChunk) -> Int32 {
  self.start
}

///|
pub fn AnnDataChunk::end_(self : AnnDataChunk) -> Int32 {
  self.end_
}

///|
pub fn AnnDataChunk::data(self : AnnDataChunk) -> AnnData {
  self.data
}

///|
/// Create chunk indices for a given size and chunk_size.
fn anndata_chunk_indices(
  n : Int32,
  chunk_size : Int32,
) -> Array[(Int32, Int32)] {
  if chunk_size <= 0 || n <= 0 {
    return []
  }
  let out : Array[(Int32, Int32)] = Array(capacity=0)
  let mut i = 0
  while i < n {
    let end_ = if i + chunk_size < n { i + chunk_size } else { n }
    out.push((i, end_))
    i = i + chunk_size
  }
  out
}

///|
/// Iterate over obs chunks. Returns all chunks as an array.
/// Each chunk is a subset AnnData containing chunk_size observations.
pub fn anndata_chunked_obs(
  adata : AnnData,
  chunk_size : Int32,
) -> Array[AnnDataChunk] {
  let n_obs = adata.n_obs()
  let n_vars = adata.n_vars()
  let all_vars : Array[Int32] = Array(capacity=n_vars)
  for j = 0; j < n_vars; j = j + 1 {
    all_vars.push(j)
  }
  let chunks : Array[AnnDataChunk] = Array(capacity=0)
  for (start, end_) in anndata_chunk_indices(n_obs, chunk_size) {
    let obs_idx : Array[Int32] = Array(capacity=end_ - start)
    for i = start; i < end_; i = i + 1 {
      obs_idx.push(i)
    }
    let subset = adata.subset(obs_idx, all_vars)
    chunks.push({ axis: 0, start, end_, data: subset, })
  }
  chunks
}

///|
/// Iterate over var chunks. Returns all chunks as an array.
pub fn anndata_chunked_var(
  adata : AnnData,
  chunk_size : Int32,
) -> Array[AnnDataChunk] {
  let n_obs = adata.n_obs()
  let n_vars = adata.n_vars()
  let all_obs : Array[Int32] = Array(capacity=n_obs)
  for i = 0; i < n_obs; i = i + 1 {
    all_obs.push(i)
  }
  let chunks : Array[AnnDataChunk] = Array(capacity=0)
  for (start, end_) in anndata_chunk_indices(n_vars, chunk_size) {
    let var_idx : Array[Int32] = Array(capacity=end_ - start)
    for j = start; j < end_; j = j + 1 {
      var_idx.push(j)
    }
    let subset = adata.subset(all_obs, var_idx)
    chunks.push({ axis: 1, start, end_, data: subset, })
  }
  chunks
}

///|
/// Apply a function to each obs chunk and collect results.
/// fn_id selects a built-in operation: 0=row_sums, 1=row_means, 2=nnz
pub fn anndata_chunked_apply_obs(
  adata : AnnData,
  chunk_size : Int32,
  fn_id : Int32,
) -> Array[Double] {
  let chunks = anndata_chunked_obs(adata, chunk_size)
  let out : Array[Double] = Array(capacity=0)
  for chunk in chunks {
    let sub = chunk.data
    match fn_id {
      0 => {
        let sums = sub.row_sums()
        for s in sums {
          out.push(s)
        }
      }
      1 => {
        let means = sub.row_means()
        for m in means {
          out.push(m)
        }
      }
      2 => out.push(sub.nnz().to_double())
      _ => ()
    }
  }
  out
}

///|
/// Number of chunks for given size and chunk_size.
pub fn anndata_n_chunks(n : Int32, chunk_size : Int32) -> Int32 {
  if chunk_size <= 0 || n <= 0 {
    return 0
  }
  (n + chunk_size - 1) / chunk_size
}
