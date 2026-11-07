/// anndata_chunked2.mbt — Chunked iteration over X and layers
/// Item 21: chunked_X / chunked_Layer with axis selection

///|
/// Iterate over X in chunks along obs axis (rows).
/// Returns (chunk_data, start_index, end_index).
pub fn anndata_chunked_x(
  adata : AnnData,
  chunk_size : Int32,
) -> Array[(Array[Array[Double]], Int32, Int32)] {
  let x = adata.x_matrix()
  let n = x.length()
  let chunks : Array[(Array[Array[Double]], Int32, Int32)] = Array(
    capacity=n / chunk_size + 1,
  )
  let mut start = 0
  while start < n {
    let end = if start + chunk_size < n { start + chunk_size } else { n }
    let chunk : Array[Array[Double]] = Array(capacity=end - start)
    for i = start; i < end; i = i + 1 {
      chunk.push(x[i])
    }
    chunks.push((chunk, start, end))
    start = end
  }
  chunks
}

///|
/// Iterate over X in chunks along var axis (columns).
pub fn anndata_chunked_x_cols(
  adata : AnnData,
  chunk_size : Int32,
) -> Array[(Array[Array[Double]], Int32, Int32)] {
  let x = adata.x_matrix()
  let n_obs = x.length()
  let n_vars = if n_obs > 0 { x[0].length() } else { 0 }
  let chunks : Array[(Array[Array[Double]], Int32, Int32)] = Array(
    capacity=n_vars / chunk_size + 1,
  )
  let mut start = 0
  while start < n_vars {
    let end = if start + chunk_size < n_vars {
      start + chunk_size
    } else {
      n_vars
    }
    // Each chunk is n_obs × (end-start)
    let chunk : Array[Array[Double]] = Array(capacity=n_obs)
    for i = 0; i < n_obs; i = i + 1 {
      let row : Array[Double] = Array(capacity=end - start)
      for j = start; j < end; j = j + 1 {
        row.push(x[i][j])
      }
      chunk.push(row)
    }
    chunks.push((chunk, start, end))
    start = end
  }
  chunks
}

///|
/// Iterate over a layer in chunks.
pub fn anndata_chunked_layer(
  adata : AnnData,
  layer_key : String,
  chunk_size : Int32,
) -> Array[(Array[Array[Double]], Int32, Int32)]? {
  if not(adata.layers.contains(layer_key)) {
    return None
  }
  let m = adata.layers[layer_key]
  let n = m.length()
  let chunks : Array[(Array[Array[Double]], Int32, Int32)] = Array(
    capacity=n / chunk_size + 1,
  )
  let mut start = 0
  while start < n {
    let end = if start + chunk_size < n { start + chunk_size } else { n }
    let chunk : Array[Array[Double]] = Array(capacity=end - start)
    for i = start; i < end; i = i + 1 {
      chunk.push(m[i])
    }
    chunks.push((chunk, start, end))
    start = end
  }
  Some(chunks)
}

///|
/// Random chunk sampling (chunk_X with replacement).
pub fn anndata_chunk_x_random(
  adata : AnnData,
  select : Int32,
  rng~ : (Double, Double, Double) = (0.3, 0.7, 0.5),
) -> Array[Array[Double]] {
  let x = adata.x_matrix()
  let n = x.length()
  if n == 0 {
    return Array(capacity=0)
  }
  let result : Array[Array[Double]] = Array(capacity=select)
  // Simple deterministic "random" sampling using provided rng values
  let rngs = [rng.0, rng.1, rng.2]
  for i = 0; i < select; i = i + 1 {
    let r = rngs[i % 3]
    let idx = (r * n.to_double()).floor().to_int()
    let clamped = if idx >= n { n - 1 } else { idx }
    result.push(x[clamped])
  }
  result
}

///|
/// Chunked iteration with a callback (simplified: returns all chunk sums).
pub fn anndata_chunked_x_sum(
  adata : AnnData,
  chunk_size : Int32,
) -> Array[Double] {
  let chunks = anndata_chunked_x(adata, chunk_size)
  let sums : Array[Double] = Array(capacity=chunks.length())
  for chunk in chunks {
    let mut s = 0.0
    for row in chunk.0 {
      for v in row {
        s = s + v
      }
    }
    sums.push(s)
  }
  sums
}
