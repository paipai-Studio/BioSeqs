/// anndata_inplace.mbt — Inplace operation semantics
/// Item 13: inplace vs copy semantics

///|
/// Inplace operation wrapper: applies a function to X in-place or returns a copy.
pub fn anndata_inplace_log(adata : AnnData, inplace~ : Bool = false) -> AnnData {
  if inplace {
    let x = adata.x_matrix()
    for i = 0; i < x.length(); i = i + 1 {
      for j = 0; j < x[i].length(); j = j + 1 {
        x[i][j] = @math.ln(1.0 + x[i][j])
      }
    }
    adata
  } else {
    let result = adata.copy()
    let x = result.x_matrix()
    for i = 0; i < x.length(); i = i + 1 {
      for j = 0; j < x[i].length(); j = j + 1 {
        x[i][j] = @math.ln(1.0 + x[i][j])
      }
    }
    result
  }
}

///|
/// Inplace scalar multiplication.
pub fn anndata_inplace_mul_scalar(
  adata : AnnData,
  scalar : Double,
  inplace~ : Bool = false,
) -> AnnData {
  if inplace {
    let x = adata.x_matrix()
    for i = 0; i < x.length(); i = i + 1 {
      for j = 0; j < x[i].length(); j = j + 1 {
        x[i][j] = x[i][j] * scalar
      }
    }
    adata
  } else {
    let result = adata.copy()
    let x = result.x_matrix()
    for i = 0; i < x.length(); i = i + 1 {
      for j = 0; j < x[i].length(); j = j + 1 {
        x[i][j] = x[i][j] * scalar
      }
    }
    result
  }
}

///|
/// Inplace scalar addition.
pub fn anndata_inplace_add_scalar(
  adata : AnnData,
  scalar : Double,
  inplace~ : Bool = false,
) -> AnnData {
  if inplace {
    let x = adata.x_matrix()
    for i = 0; i < x.length(); i = i + 1 {
      for j = 0; j < x[i].length(); j = j + 1 {
        x[i][j] = x[i][j] + scalar
      }
    }
    adata
  } else {
    let result = adata.copy()
    let x = result.x_matrix()
    for i = 0; i < x.length(); i = i + 1 {
      for j = 0; j < x[i].length(); j = j + 1 {
        x[i][j] = x[i][j] + scalar
      }
    }
    result
  }
}

///|
/// Inplace normalization (divide by row sum).
pub fn anndata_inplace_normalize(
  adata : AnnData,
  target_sum~ : Double = 10000.0,
  inplace~ : Bool = false,
) -> AnnData {
  let result = if inplace { adata } else { adata.copy() }
  let x = result.x_matrix()
  for i = 0; i < x.length(); i = i + 1 {
    let mut row_sum = 0.0
    for j = 0; j < x[i].length(); j = j + 1 {
      row_sum = row_sum + x[i][j]
    }
    if row_sum > 0.0 {
      let scale = target_sum / row_sum
      for j = 0; j < x[i].length(); j = j + 1 {
        x[i][j] = x[i][j] * scale
      }
    }
  }
  result
}

///|
/// Inplace filter genes by minimum count.
pub fn anndata_inplace_filter_genes(
  adata : AnnData,
  min_count~ : Double = 0.0,
  inplace~ : Bool = false,
) -> AnnData {
  let x = adata.x_matrix()
  let n_vars = adata.n_vars()
  // Compute column sums
  let col_sums : Array[Double] = Array(capacity=n_vars)
  for j = 0; j < n_vars; j = j + 1 {
    let mut s = 0.0
    for i = 0; i < x.length(); i = i + 1 {
      s = s + x[i][j]
    }
    col_sums.push(s)
  }
  // Find kept columns
  let keep : Array[Int32] = Array(capacity=n_vars)
  for j = 0; j < n_vars; j = j + 1 {
    if col_sums[j] > min_count {
      keep.push(j)
    }
  }
  if keep.length() == n_vars {
    return if inplace { adata } else { adata.copy() }
  }
  // Subset
  let keep_i32 : Array[Int32] = Array(capacity=keep.length())
  for k in keep {
    keep_i32.push(k)
  }
  adata.subset(Array(capacity=0), keep_i32)
}
