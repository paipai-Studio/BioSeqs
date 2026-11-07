/// anndata_delayed.mbt — chunked lazy evaluation over AnnData X.
///
/// MoonBit analogue of `anndata.experimental.read_elem_as_dask` / dask
/// arrays: an `AnnDataDelayed` wraps an AnnData with a chunk size and a
/// chain of elementwise ops. Chunks are materialized on demand
/// (`anndata_delayed_chunk`), streaming reductions never build the full
/// matrix (`anndata_delayed_row_sums` / `col_sums`), and
/// `anndata_delayed_collect` materializes everything.

///|
/// Elementwise ops applicable lazily per chunk.
pub(all) enum AnnDataDelayedOp {
  DIdentity
  DLog1p
  DAddScalar(Double)
  DMulScalar(Double)
  DSquare
}

///|
/// A lazy, chunked view of an AnnData's X with a pending op chain.
pub struct AnnDataDelayed {
  source : AnnData
  chunk_size : Int32
  ops : Array[AnnDataDelayedOp]
}

///|
/// Create a delayed view of X (no data copied).
pub fn anndata_delayed_new(
  adata : AnnData,
  chunk_size~ : Int32 = 1000,
) -> AnnDataDelayed {
  AnnDataDelayed::{
    source: adata,
    chunk_size: if chunk_size <= 0 {
      1000
    } else {
      chunk_size
    },
    ops: Array(capacity=0),
  }
}

///|
/// Append an elementwise op (returns a new delayed handle; the source is
/// shared, ops accumulate).
pub fn anndata_delayed_map(
  d : AnnDataDelayed,
  op : AnnDataDelayedOp,
) -> AnnDataDelayed {
  let ops = d.ops.copy()
  ops.push(op)
  AnnDataDelayed::{ source: d.source, chunk_size: d.chunk_size, ops, }
}

///|
/// Number of row chunks.
pub fn anndata_delayed_n_chunks(d : AnnDataDelayed) -> Int32 {
  anndata_zarr_chunk_count(d.source.n_obs(), d.chunk_size)
}

///|
/// Apply the op chain to a single value.
fn delayed_apply_one(ops : Array[AnnDataDelayedOp], v : Double) -> Double {
  let mut x = v
  for op in ops {
    match op {
      DIdentity => ()
      DLog1p => x = @math.ln(1.0 + x)
      DAddScalar(delta) => x = x + delta
      DMulScalar(f) => x = x * f
      DSquare => x = x * x
    }
  }
  x
}

///|
/// Materialize chunk `i` (rows i*chunk .. (i+1)*chunk) with all pending
/// ops applied. Out-of-range chunk indices yield an empty matrix.
pub fn anndata_delayed_chunk(
  d : AnnDataDelayed,
  i : Int32,
) -> Array[Array[Double]] {
  if i < 0 || i >= anndata_delayed_n_chunks(d) {
    return []
  }
  let start = i * d.chunk_size
  let end_ = if start + d.chunk_size > d.source.n_obs() {
    d.source.n_obs()
  } else {
    start + d.chunk_size
  }
  let out : Array[Array[Double]] = Array(capacity=0)
  for r = start; r < end_; r = r + 1 {
    let row : Array[Double] = Array(capacity=d.source.n_vars())
    for c = 0; c < d.source.n_vars(); c = c + 1 {
      row.push(delayed_apply_one(d.ops, d.source.get_x(r, c)))
    }
    out.push(row)
  }
  out
}

///|
/// Materialize the full X with all pending ops applied.
pub fn anndata_delayed_collect(d : AnnDataDelayed) -> Array[Array[Double]] {
  let out : Array[Array[Double]] = Array(capacity=d.source.n_obs())
  for i = 0; i < d.source.n_obs(); i = i + 1 {
    let row : Array[Double] = Array(capacity=d.source.n_vars())
    for c = 0; c < d.source.n_vars(); c = c + 1 {
      row.push(delayed_apply_one(d.ops, d.source.get_x(i, c)))
    }
    out.push(row)
  }
  out
}

///|
/// Materialize the delayed object into a full AnnData.
pub fn anndata_delayed_to_adata(d : AnnDataDelayed) -> AnnData {
  AnnData::new(
    anndata_delayed_collect(d),
    d.source.obs_names.copy(),
    d.source.var_names.copy(),
  )
}

///|
/// Streaming row sums: one chunk at a time, never materializing the
/// full matrix.
pub fn anndata_delayed_row_sums(d : AnnDataDelayed) -> Array[Double] {
  let out : Array[Double] = Array(capacity=d.source.n_obs())
  let n_chunks = anndata_delayed_n_chunks(d)
  for i = 0; i < n_chunks; i = i + 1 {
    let chunk = anndata_delayed_chunk(d, i)
    for row in chunk {
      let mut sum = 0.0
      for v in row {
        sum = sum + v
      }
      out.push(sum)
    }
  }
  out
}

///|
/// Streaming column sums (chunk-at-a-time).
pub fn anndata_delayed_col_sums(d : AnnDataDelayed) -> Array[Double] {
  let sums : Array[Double] = Array::make(d.source.n_vars(), 0.0)
  let n_chunks = anndata_delayed_n_chunks(d)
  for i = 0; i < n_chunks; i = i + 1 {
    let chunk = anndata_delayed_chunk(d, i)
    for r = 0; r < chunk.length(); r = r + 1 {
      for c = 0; c < chunk[r].length(); c = c + 1 {
        sums[c] = sums[c] + chunk[r][c]
      }
    }
  }
  sums
}

///|
/// Single element access through the op chain (no chunk materialization).
pub fn anndata_delayed_elem(d : AnnDataDelayed, i : Int32, j : Int32) -> Double {
  delayed_apply_one(d.ops, d.source.get_x(i, j))
}

///|
/// `anndata.experimental.read_elem_as_dask` analogue: wrap an X element
/// (JSON dense array) as a delayed, chunked object. Returns `None` when
/// the element is not a dense numeric array.
pub fn anndata_read_elem_as_delayed(
  x_elem : JsonValue,
  chunk_size~ : Int32 = 1000,
) -> AnnDataDelayed? {
  match x_elem {
    JsonArr(rows) => {
      let mat : Array[Array[Double]] = Array(capacity=rows.length())
      for row in rows {
        match row {
          JsonArr(cols) => {
            let r : Array[Double] = Array(capacity=cols.length())
            for c in cols {
              match c {
                JsonNum(v) => r.push(v)
                _ => r.push(0.0)
              }
            }
            mat.push(r)
          }
          _ => return None
        }
      }
      let n_obs = mat.length()
      let n_vars = if n_obs > 0 { mat[0].length() } else { 0 }
      let obs_names : Array[String] = Array(capacity=n_obs)
      for i = 0; i < n_obs; i = i + 1 {
        obs_names.push("obs_" + i.to_string())
      }
      let var_names : Array[String] = Array(capacity=n_vars)
      for j = 0; j < n_vars; j = j + 1 {
        var_names.push("var_" + j.to_string())
      }
      Some(
        anndata_delayed_new(
          AnnData::new(mat, obs_names, var_names),
          chunk_size~,
        ),
      )
    }
    _ => None
  }
}
