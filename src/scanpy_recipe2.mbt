///|
/// scanpy `pp.recipe_weinreb17` — the Weinreb et al. (2018) preprocessing
/// recipe for haematopoiesis single-cell data.  Pipeline: filter genes
/// (`min_cells`), normalize to the median count depth, log1p, select highly
/// variable genes by the coefficient of variation (cv^2), scale to zero mean
/// / unit variance, then PCA and a diffusion map on the kNN graph of the
/// scaled HVG matrix (scanpy stores both in `obsm` for trajectory analysis).

///|
/// Median of a non-empty slice (average of the two middle values for even
/// lengths).
fn scanpy_recipe_median(a : Array[Double]) -> Double {
  let n = a.length()
  if n == 0 {
    return 0.0
  }
  let s = a.copy()
  s.sort()
  if n % 2 == 1 {
    s[n / 2]
  } else {
    (s[n / 2 - 1] + s[n / 2]) / 2.0
  }
}

///|
/// Result of the Weinreb 2017 recipe: the scaled HVG matrix (the PCA input),
/// the selected HVG indices (into the gene-filtered space), their PCA, the
/// diffusion map of the kNN graph over the PCA coordinates and the median
/// count depth used for normalization.
pub struct ScanpyWeinreb17 {
  scaled : Array[Array[Double]]
  hvg_indices : Array[Int32]
  pca : ScanpyPca
  diffusion : ScanpyDiffmap
  counts_median : Double
}

///|
/// Run the Weinreb 2017 recipe on a raw count matrix (cells x genes):
///
/// 1. keep genes expressed in at least `min_cells` cells,
/// 2. normalize each cell to `target_sum` (default: the median total count
///    over all cells, Weinreb's choice) and apply log1p,
/// 3. select the `n_top_genes` with the highest coefficient of variation
///    squared (cv^2 = var / mean^2) among kept genes,
/// 4. z-score the HVG columns (zero mean, unit variance),
/// 5. PCA to `n_comps` components,
/// 6. build the kNN connectivities (`n_neighbors`) on the PCA coordinates and
///    compute a diffusion map with `n_comps` components.
pub fn scanpy_recipe_weinreb17(
  x : Array[Array[Double]],
  min_cells~ : Int32 = 3,
  target_sum~ : Double = 0.0,
  n_top_genes~ : Int32 = 2000,
  n_comps~ : Int32 = 50,
  n_neighbors~ : Int32 = 15,
) -> ScanpyWeinreb17 {
  let n_cells = x.length()
  let n_genes = if n_cells > 0 { x[0].length() } else { 0 }
  if n_cells == 0 || n_genes == 0 {
    return ScanpyWeinreb17::{
      scaled: [],
      hvg_indices: [],
      pca: ScanpyPca::{ coords: [], loadings: [], variance_ratio: [], },
      diffusion: ScanpyDiffmap::{ eigen_basis: [], eigen_values: [], },
      counts_median: 0.0,
    }
  }
  // step 1: gene filter
  let kept : Array[Int32] = Array(capacity=n_genes)
  for j = 0; j < n_genes; j = j + 1 {
    let mut cnt = 0
    for i = 0; i < n_cells; i = i + 1 {
      if x[i][j] > 0.0 {
        cnt = cnt + 1
      }
    }
    if cnt >= min_cells {
      kept.push(j)
    }
  }
  let n_kept = kept.length()
  if n_kept == 0 {
    return ScanpyWeinreb17::{
      scaled: [],
      hvg_indices: [],
      pca: ScanpyPca::{ coords: [], loadings: [], variance_ratio: [], },
      diffusion: ScanpyDiffmap::{ eigen_basis: [], eigen_values: [], },
      counts_median: 0.0,
    }
  }
  // step 2: median count depth + normalize_total + log1p
  let totals : Array[Double] = Array(capacity=n_cells)
  for i = 0; i < n_cells; i = i + 1 {
    let mut s = 0.0
    for j = 0; j < n_genes; j = j + 1 {
      s = s + x[i][j]
    }
    totals.push(s)
  }
  let median = if target_sum > 0.0 {
    target_sum
  } else {
    let m = scanpy_recipe_median(totals)
    if m > 0.0 {
      m
    } else {
      1.0
    }
  }
  let norm : Array[Array[Double]] = Array(capacity=n_cells)
  for i = 0; i < n_cells; i = i + 1 {
    let row : Array[Double] = Array(capacity=n_kept)
    let scale = if totals[i] > 0.0 { median / totals[i] } else { 0.0 }
    for t = 0; t < n_kept; t = t + 1 {
      row.push(@math.ln(1.0 + x[i][kept[t]] * scale))
    }
    norm.push(row)
  }
  // step 3: HVG by cv^2 = var / mean^2 (skip zero-mean genes)
  let cv2 : Array[(Double, Int32)] = Array(capacity=n_kept)
  for t = 0; t < n_kept; t = t + 1 {
    let mut s = 0.0
    for i = 0; i < n_cells; i = i + 1 {
      s = s + norm[i][t]
    }
    let mean = s / n_cells.to_double()
    let mut v = 0.0
    for i = 0; i < n_cells; i = i + 1 {
      let d = norm[i][t] - mean
      v = v + d * d
    }
    let variance = v / n_cells.to_double()
    let score = if mean > 1.0e-12 { variance / (mean * mean) } else { 0.0 }
    cv2.push((score, t))
  }
  cv2.sort_by(fn(a, b) {
    if a.0 > b.0 {
      -1
    } else if a.0 < b.0 {
      1
    } else if a.1 < b.1 {
      -1
    } else if a.1 > b.1 {
      1
    } else {
      0
    }
  })
  let n_select = if n_top_genes < n_kept { n_top_genes } else { n_kept }
  let hvg_indices : Array[Int32] = Array(capacity=n_select)
  let hvg_mask = Array::make(n_kept, false)
  for t = 0; t < n_select; t = t + 1 {
    hvg_indices.push(cv2[t].1)
    hvg_mask[cv2[t].1] = true
  }
  // step 4: scale HVG columns to zero mean / unit variance
  let hvg_means : Array[Double] = Array(capacity=n_select)
  let hvg_stds : Array[Double] = Array(capacity=n_select)
  for t = 0; t < n_select; t = t + 1 {
    let j = hvg_indices[t]
    let mut s = 0.0
    for i = 0; i < n_cells; i = i + 1 {
      s = s + norm[i][j]
    }
    let mean = s / n_cells.to_double()
    let mut v = 0.0
    for i = 0; i < n_cells; i = i + 1 {
      let d = norm[i][j] - mean
      v = v + d * d
    }
    hvg_means.push(mean)
    hvg_stds.push((v / n_cells.to_double()).sqrt())
  }
  let scaled : Array[Array[Double]] = Array(capacity=n_cells)
  for i = 0; i < n_cells; i = i + 1 {
    let row : Array[Double] = Array(capacity=n_select)
    for t = 0; t < n_select; t = t + 1 {
      let j = hvg_indices[t]
      let std = hvg_stds[t]
      if std > 1.0e-12 {
        row.push((norm[i][j] - hvg_means[t]) / std)
      } else {
        row.push(0.0)
      }
    }
    scaled.push(row)
  }
  ignore(hvg_mask)
  // step 5: PCA on the scaled HVG matrix
  let comps = if n_comps < n_select { n_comps } else { n_select }
  let pca = scanpy_pca(scaled, comps)
  // step 6: neighbors + diffusion map on the PCA coordinates
  let conn = scanpy_neighbors(pca.coords, n_neighbors)
  let diff = scanpy_diffmap(conn, comps)
  ScanpyWeinreb17::{
    scaled,
    hvg_indices,
    pca,
    diffusion: diff,
    counts_median: median,
  }
}
