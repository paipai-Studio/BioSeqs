/// AnnData groupby aggregation.
///
/// `adata.obs.groupby("leiden").mean()` equivalent.
///
/// Groups observations by a categorical column and aggregates X
/// (mean, sum, median, std) per group.

///|
/// Aggregation result: group names + aggregated matrix
pub struct AnnDataGroupbyResult {
  group_names : Array[String] // unique group labels
  group_sizes : Array[Int32] // number of obs per group
  result : Array[Array[Double]] // n_groups × n_vars matrix
}

///|
pub fn AnnDataGroupbyResult::group_names(
  self : AnnDataGroupbyResult,
) -> Array[String] {
  self.group_names
}

///|
pub fn AnnDataGroupbyResult::group_sizes(
  self : AnnDataGroupbyResult,
) -> Array[Int32] {
  self.group_sizes
}

///|
pub fn AnnDataGroupbyResult::result(
  self : AnnDataGroupbyResult,
) -> Array[Array[Double]] {
  self.result
}

///|
/// Group observations by an obs column and aggregate X.
/// agg: 0=mean, 1=sum, 2=median, 3=std
pub fn anndata_groupby(
  adata : AnnData,
  key : String,
  agg : Int32,
) -> AnnDataGroupbyResult? {
  match adata.obs_column(key) {
    None => None
    Some(col) => {
      let cat = anndata_encode_categorical(col)
      let n_groups = cat.categories.length()
      let n_vars = adata.n_vars()
      let n_obs = adata.n_obs()
      // Compute group sizes
      let sizes = cat.counts()
      // For each group, accumulate X values
      let group_sums : Array[Array[Double]] = Array::make(
        n_groups,
        Array::make(n_vars, 0.0),
      )
      for i = 0; i < n_obs; i = i + 1 {
        let g = cat.codes[i]
        if g >= 0 && g < n_groups {
          for j = 0; j < n_vars; j = j + 1 {
            group_sums[g][j] = group_sums[g][j] + adata.get_x(i, j)
          }
        }
      }
      let result : Array[Array[Double]] = Array::make(
        n_groups,
        Array::make(n_vars, 0.0),
      )
      match agg {
        1 =>
          // Sum
          for g = 0; g < n_groups; g = g + 1 {
            for j = 0; j < n_vars; j = j + 1 {
              result[g][j] = group_sums[g][j]
            }
          }
        0 =>
          // Mean
          for g = 0; g < n_groups; g = g + 1 {
            let sz = if sizes[g] > 0 { sizes[g].to_double() } else { 1.0 }
            for j = 0; j < n_vars; j = j + 1 {
              result[g][j] = group_sums[g][j] / sz
            }
          }
        2 =>
          // Median: collect all values per group per gene, then median
          for g = 0; g < n_groups; g = g + 1 {
            for j = 0; j < n_vars; j = j + 1 {
              let vals : Array[Double] = Array(capacity=0)
              for i = 0; i < n_obs; i = i + 1 {
                if cat.codes[i] == g {
                  vals.push(adata.get_x(i, j))
                }
              }
              if vals.length() > 0 {
                vals.sort()
                let mid = vals.length() / 2
                result[g][j] = if vals.length() % 2 == 0 {
                  (vals[mid - 1] + vals[mid]) * 0.5
                } else {
                  vals[mid]
                }
              }
            }
          }
        3 =>
          // Std: compute mean first, then std
          for g = 0; g < n_groups; g = g + 1 {
            let sz = if sizes[g] > 0 { sizes[g].to_double() } else { 1.0 }
            for j = 0; j < n_vars; j = j + 1 {
              let mean = group_sums[g][j] / sz
              let mut var_sum = 0.0
              for i = 0; i < n_obs; i = i + 1 {
                if cat.codes[i] == g {
                  let d = adata.get_x(i, j) - mean
                  var_sum = var_sum + d * d
                }
              }
              result[g][j] = (var_sum / sz).sqrt()
            }
          }
        _ =>
          // Default: mean
          for g = 0; g < n_groups; g = g + 1 {
            let sz = if sizes[g] > 0 { sizes[g].to_double() } else { 1.0 }
            for j = 0; j < n_vars; j = j + 1 {
              result[g][j] = group_sums[g][j] / sz
            }
          }
      }
      Some({ group_names: cat.categories, group_sizes: sizes, result, })
    }
  }
}

///|
/// Groupby mean (convenience wrapper).
pub fn anndata_groupby_mean(
  adata : AnnData,
  key : String,
) -> AnnDataGroupbyResult? {
  anndata_groupby(adata, key, 0)
}

///|
/// Groupby sum (convenience wrapper).
pub fn anndata_groupby_sum(
  adata : AnnData,
  key : String,
) -> AnnDataGroupbyResult? {
  anndata_groupby(adata, key, 1)
}

///|
/// Groupby median (convenience wrapper).
pub fn anndata_groupby_median(
  adata : AnnData,
  key : String,
) -> AnnDataGroupbyResult? {
  anndata_groupby(adata, key, 2)
}

///|
/// Groupby std (convenience wrapper).
pub fn anndata_groupby_std(
  adata : AnnData,
  key : String,
) -> AnnDataGroupbyResult? {
  anndata_groupby(adata, key, 3)
}

///|
/// Aggregate an obsm matrix by group (e.g. mean PCA per cluster).
pub fn anndata_groupby_obsm(
  adata : AnnData,
  obsm_key : String,
  group_key : String,
  agg : Int32,
) -> AnnDataGroupbyResult? {
  match (adata.get_obsm(obsm_key), adata.obs_column(group_key)) {
    (Some(mat), Some(col)) => {
      let cat = anndata_encode_categorical(col)
      let n_groups = cat.categories.length()
      let n_cols = if mat.length() > 0 { mat[0].length() } else { 0 }
      let sizes = cat.counts()
      let sums : Array[Array[Double]] = Array::make(
        n_groups,
        Array::make(n_cols, 0.0),
      )
      for i = 0; i < mat.length(); i = i + 1 {
        let g = if i < cat.codes.length() { cat.codes[i] } else { -1 }
        if g >= 0 && g < n_groups {
          for j = 0; j < n_cols && j < mat[i].length(); j = j + 1 {
            sums[g][j] = sums[g][j] + mat[i][j]
          }
        }
      }
      let result : Array[Array[Double]] = Array::make(
        n_groups,
        Array::make(n_cols, 0.0),
      )
      for g = 0; g < n_groups; g = g + 1 {
        let sz = if sizes[g] > 0 { sizes[g].to_double() } else { 1.0 }
        for j = 0; j < n_cols; j = j + 1 {
          if agg == 1 {
            result[g][j] = sums[g][j]
          } else {
            result[g][j] = sums[g][j] / sz
          }
        }
      }
      Some({ group_names: cat.categories, group_sizes: sizes, result, })
    }
    _ => None
  }
}
