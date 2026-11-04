/// AnnData element-wise X operations and layer arithmetic.
///
/// `adata.X = np.log1p(adata.X)` and
/// `adata.layers["diff"] = layers["a"] - layers["b"]`

///|
/// Apply log (natural) to every element of X. Returns new matrix.
pub fn anndata_x_log(adata : AnnData) -> Array[Array[Double]] {
  let n_obs = adata.n_obs()
  let n_vars = adata.n_vars()
  let out : Array[Array[Double]] = Array::make(n_obs, Array::make(n_vars, 0.0))
  for i = 0; i < n_obs; i = i + 1 {
    for j = 0; j < n_vars; j = j + 1 {
      let v = adata.get_x(i, j)
      out[i][j] = if v > 0.0 { @math.ln(v) } else { 0.0 }
    }
  }
  out
}

///|
/// Apply log1p (log(1+x)) to every element of X.
pub fn anndata_x_log1p(adata : AnnData) -> Array[Array[Double]] {
  let n_obs = adata.n_obs()
  let n_vars = adata.n_vars()
  let out : Array[Array[Double]] = Array::make(n_obs, Array::make(n_vars, 0.0))
  for i = 0; i < n_obs; i = i + 1 {
    for j = 0; j < n_vars; j = j + 1 {
      let v = adata.get_x(i, j)
      if v > -1.0 {
        out[i][j] = @math.ln(1.0 + v)
      } else {
        out[i][j] = 0.0
      }
    }
  }
  out
}

///|
/// Apply log2 to every element of X.
pub fn anndata_x_log2(adata : AnnData) -> Array[Array[Double]] {
  let n_obs = adata.n_obs()
  let n_vars = adata.n_vars()
  let out : Array[Array[Double]] = Array::make(n_obs, Array::make(n_vars, 0.0))
  for i = 0; i < n_obs; i = i + 1 {
    for j = 0; j < n_vars; j = j + 1 {
      let v = adata.get_x(i, j)
      out[i][j] = if v > 0.0 { @math.log2(v) } else { 0.0 }
    }
  }
  out
}

///|
/// Apply sqrt to every element of X.
pub fn anndata_x_sqrt(adata : AnnData) -> Array[Array[Double]] {
  let n_obs = adata.n_obs()
  let n_vars = adata.n_vars()
  let out : Array[Array[Double]] = Array::make(n_obs, Array::make(n_vars, 0.0))
  for i = 0; i < n_obs; i = i + 1 {
    for j = 0; j < n_vars; j = j + 1 {
      let v = adata.get_x(i, j)
      out[i][j] = if v > 0.0 { v.sqrt() } else { 0.0 }
    }
  }
  out
}

///|
/// Add scalar to every element of X.
pub fn anndata_x_add_scalar(
  adata : AnnData,
  scalar : Double,
) -> Array[Array[Double]] {
  let n_obs = adata.n_obs()
  let n_vars = adata.n_vars()
  let out : Array[Array[Double]] = Array::make(n_obs, Array::make(n_vars, 0.0))
  for i = 0; i < n_obs; i = i + 1 {
    for j = 0; j < n_vars; j = j + 1 {
      out[i][j] = adata.get_x(i, j) + scalar
    }
  }
  out
}

///|
/// Multiply every element of X by scalar.
pub fn anndata_x_mul_scalar(
  adata : AnnData,
  scalar : Double,
) -> Array[Array[Double]] {
  let n_obs = adata.n_obs()
  let n_vars = adata.n_vars()
  let out : Array[Array[Double]] = Array::make(n_obs, Array::make(n_vars, 0.0))
  for i = 0; i < n_obs; i = i + 1 {
    for j = 0; j < n_vars; j = j + 1 {
      out[i][j] = adata.get_x(i, j) * scalar
    }
  }
  out
}

///|
/// Apply arbitrary element-wise function to X.
/// fn_id: 0=log, 1=log1p, 2=log2, 3=sqrt, 4=exp, 5=square, 6=abs
pub fn anndata_x_apply(adata : AnnData, fn_id : Int32) -> Array[Array[Double]] {
  let n_obs = adata.n_obs()
  let n_vars = adata.n_vars()
  let out : Array[Array[Double]] = Array::make(n_obs, Array::make(n_vars, 0.0))
  for i = 0; i < n_obs; i = i + 1 {
    for j = 0; j < n_vars; j = j + 1 {
      let v = adata.get_x(i, j)
      out[i][j] = match fn_id {
        0 => if v > 0.0 { @math.ln(v) } else { 0.0 }
        1 => if v > -1.0 { @math.ln(1.0 + v) } else { 0.0 }
        2 => if v > 0.0 { @math.log2(v) } else { 0.0 }
        3 => if v > 0.0 { v.sqrt() } else { 0.0 }
        4 => @math.exp(v)
        5 => v * v
        6 => if v < 0.0 { -v } else { v }
        _ => v
      }
    }
  }
  out
}

// === Layer arithmetic ===

///|
/// Add two matrices element-wise.
pub fn anndata_mat_add(
  a : Array[Array[Double]],
  b : Array[Array[Double]],
) -> Array[Array[Double]] {
  let n = a.length()
  if n == 0 {
    return []
  }
  let m = a[0].length()
  let out : Array[Array[Double]] = Array::make(n, Array::make(m, 0.0))
  for i = 0; i < n; i = i + 1 {
    for j = 0; j < m; j = j + 1 {
      let bv = if i < b.length() && j < b[i].length() { b[i][j] } else { 0.0 }
      out[i][j] = a[i][j] + bv
    }
  }
  out
}

///|
/// Subtract matrix b from a element-wise.
pub fn anndata_mat_sub(
  a : Array[Array[Double]],
  b : Array[Array[Double]],
) -> Array[Array[Double]] {
  let n = a.length()
  if n == 0 {
    return []
  }
  let m = a[0].length()
  let out : Array[Array[Double]] = Array::make(n, Array::make(m, 0.0))
  for i = 0; i < n; i = i + 1 {
    for j = 0; j < m; j = j + 1 {
      let bv = if i < b.length() && j < b[i].length() { b[i][j] } else { 0.0 }
      out[i][j] = a[i][j] - bv
    }
  }
  out
}

///|
/// Multiply two matrices element-wise (Hadamard product).
pub fn anndata_mat_mul(
  a : Array[Array[Double]],
  b : Array[Array[Double]],
) -> Array[Array[Double]] {
  let n = a.length()
  if n == 0 {
    return []
  }
  let m = a[0].length()
  let out : Array[Array[Double]] = Array::make(n, Array::make(m, 0.0))
  for i = 0; i < n; i = i + 1 {
    for j = 0; j < m; j = j + 1 {
      let bv = if i < b.length() && j < b[i].length() { b[i][j] } else { 0.0 }
      out[i][j] = a[i][j] * bv
    }
  }
  out
}

///|
/// Divide matrix a by b element-wise.
pub fn anndata_mat_div(
  a : Array[Array[Double]],
  b : Array[Array[Double]],
) -> Array[Array[Double]] {
  let n = a.length()
  if n == 0 {
    return []
  }
  let m = a[0].length()
  let out : Array[Array[Double]] = Array::make(n, Array::make(m, 0.0))
  for i = 0; i < n; i = i + 1 {
    for j = 0; j < m; j = j + 1 {
      let bv = if i < b.length() && j < b[i].length() { b[i][j] } else { 1.0 }
      out[i][j] = if bv != 0.0 { a[i][j] / bv } else { 0.0 }
    }
  }
  out
}

///|
/// Compute layer from X and store: result = X op layer_b.
/// op: 0=add, 1=sub, 2=mul, 3=div
pub fn anndata_layer_from_x(
  adata : AnnData,
  layer_b_name : String,
  result_name : String,
  op : Int32,
) -> Unit {
  match adata.get_layer(layer_b_name) {
    Some(b) => {
      let x = adata.x_matrix()
      let result = match op {
        0 => anndata_mat_add(x, b)
        1 => anndata_mat_sub(x, b)
        2 => anndata_mat_mul(x, b)
        3 => anndata_mat_div(x, b)
        _ => x
      }
      ignore(adata.set_layer(result_name, result))
    }
    None => ()
  }
}
