/// AnnData numeric obs/var columns.
///
/// Python AnnData stores obs/var as pandas DataFrames with mixed types
/// (string, float, int, bool). BioSeqs' AnnData only stores String.
/// This module adds numeric annotation support alongside string columns.

///|
/// Numeric annotation storage for obs and var.
/// Keys are column names, values are Double arrays.
pub struct AnnDataNumericAnno {
  obs_num : Map[String, Array[Double]]
  var_num : Map[String, Array[Double]]
}

///|
pub fn AnnDataNumericAnno::new() -> AnnDataNumericAnno {
  { obs_num: Map([], capacity=4), var_num: Map([], capacity=4), }
}

///|
/// Set a numeric obs column (e.g. n_genes, total_counts).
pub fn AnnDataNumericAnno::set_obs_num(
  self : AnnDataNumericAnno,
  key : String,
  values : Array[Double],
) -> Unit {
  self.obs_num[key] = values
}

///|
/// Set a numeric var column (e.g. mean_expression, dispersion).
pub fn AnnDataNumericAnno::set_var_num(
  self : AnnDataNumericAnno,
  key : String,
  values : Array[Double],
) -> Unit {
  self.var_num[key] = values
}

///|
/// Get a numeric obs column.
pub fn AnnDataNumericAnno::get_obs_num(
  self : AnnDataNumericAnno,
  key : String,
) -> Array[Double]? {
  match self.obs_num.get(key) {
    Some(v) => Some(v)
    None => None
  }
}

///|
/// Get a numeric var column.
pub fn AnnDataNumericAnno::get_var_num(
  self : AnnDataNumericAnno,
  key : String,
) -> Array[Double]? {
  match self.var_num.get(key) {
    Some(v) => Some(v)
    None => None
  }
}

///|
/// Compute QC metrics and store as numeric obs columns.
/// n_genes: number of non-zero genes per cell.
/// total_counts: sum of expression per cell.
pub fn anndata_compute_qc_obs(adata : AnnData) -> AnnDataNumericAnno {
  let anno = AnnDataNumericAnno::new()
  let n_obs = adata.n_obs()
  let n_genes : Array[Double] = Array::make(n_obs, 0.0)
  let total_counts : Array[Double] = Array::make(n_obs, 0.0)
  for i = 0; i < n_obs; i = i + 1 {
    let mut ng = 0.0
    let mut tc = 0.0
    for j = 0; j < adata.n_vars(); j = j + 1 {
      let v = adata.get_x(i, j)
      if v != 0.0 {
        ng = ng + 1.0
        tc = tc + v
      }
    }
    n_genes[i] = ng
    total_counts[i] = tc
  }
  anno.set_obs_num("n_genes", n_genes)
  anno.set_obs_num("total_counts", total_counts)
  // Per-gene stats
  let n_vars = adata.n_vars()
  let n_cells : Array[Double] = Array::make(n_vars, 0.0)
  let mean_expr : Array[Double] = Array::make(n_vars, 0.0)
  for j = 0; j < n_vars; j = j + 1 {
    let mut nc = 0.0
    let mut s = 0.0
    for i = 0; i < n_obs; i = i + 1 {
      let v = adata.get_x(i, j)
      if v != 0.0 {
        nc = nc + 1.0
      }
      s = s + v
    }
    n_cells[j] = nc
    mean_expr[j] = if n_obs > 0 { s / n_obs.to_double() } else { 0.0 }
  }
  anno.set_var_num("n_cells", n_cells)
  anno.set_var_num("mean_expression", mean_expr)
  anno
}

///|
/// Filter obs by numeric column threshold.
/// Returns boolean mask: true if value >= min and value <= max.
pub fn anndata_filter_numeric(
  anno : AnnDataNumericAnno,
  key : String,
  min_val~ : Double = -999999999999999.0,
  max_val~ : Double = 999999999999999.0,
) -> Array[Bool] {
  match anno.get_obs_num(key) {
    Some(col) => {
      let mask : Array[Bool] = Array::make(col.length(), false)
      for i = 0; i < col.length(); i = i + 1 {
        mask[i] = col[i] >= min_val && col[i] <= max_val
      }
      mask
    }
    None => []
  }
}

///|
/// Subset numeric annotations by index.
pub fn AnnDataNumericAnno::subset(
  self : AnnDataNumericAnno,
  obs_idx : Array[Int32],
  var_idx : Array[Int32],
) -> AnnDataNumericAnno {
  let result = AnnDataNumericAnno::new()
  let obs_keys = self.obs_num.keys().collect()
  for k in obs_keys {
    match self.get_obs_num(k) {
      Some(col) => {
        let sub : Array[Double] = Array(capacity=obs_idx.length())
        for i in obs_idx {
          if i >= 0 && i < col.length() {
            sub.push(col[i])
          } else {
            sub.push(0.0)
          }
        }
        result.set_obs_num(k, sub)
      }
      None => ()
    }
  }
  let var_keys = self.var_num.keys().collect()
  for k in var_keys {
    match self.get_var_num(k) {
      Some(col) => {
        let sub : Array[Double] = Array(capacity=var_idx.length())
        for i in var_idx {
          if i >= 0 && i < col.length() {
            sub.push(col[i])
          } else {
            sub.push(0.0)
          }
        }
        result.set_var_num(k, sub)
      }
      None => ()
    }
  }
  result
}

///|
/// Summary statistics for a numeric obs column.
/// Returns (mean, std, min, max, median).
pub fn anndata_numeric_summary(
  anno : AnnDataNumericAnno,
  key : String,
) -> (Double, Double, Double, Double, Double) {
  match anno.get_obs_num(key) {
    Some(col) => {
      if col.length() == 0 {
        return (0.0, 0.0, 0.0, 0.0, 0.0)
      }
      let mut sum = 0.0
      let mut min_v = col[0]
      let mut max_v = col[0]
      for v in col {
        sum = sum + v
        if v < min_v {
          min_v = v
        }
        if v > max_v {
          max_v = v
        }
      }
      let mean = sum / col.length().to_double()
      let mut var_sum = 0.0
      for v in col {
        let d = v - mean
        var_sum = var_sum + d * d
      }
      let std = (var_sum / col.length().to_double()).sqrt()
      // Simple median (sort copy)
      let sorted = col.copy()
      sorted.sort()
      let mid = sorted.length() / 2
      let median = if sorted.length() % 2 == 0 {
        (sorted[mid - 1] + sorted[mid]) * 0.5
      } else {
        sorted[mid]
      }
      (mean, std, min_v, max_v, median)
    }
    None => (0.0, 0.0, 0.0, 0.0, 0.0)
  }
}
