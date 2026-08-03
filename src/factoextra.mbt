///|
/// factoextra - PCA/factor analysis utilities inspired by Bioconductor factoextra.
///
/// Provides:
/// - get_eigenvalue: Extract eigenvalues and explained variance
/// - get_pca_ind: Extract individual coordinates
/// - get_pca_var: Extract variable coordinates
/// - contrib: Contribution of variables to dimensions
/// - cos: Cosine quality of representation
/// - Dimdesc: Dimension description (quantitative)
/// - PCA computation via eigenvalue decomposition

///|
pub struct EigenValue {
  dimension : Int
  eigenvalue : Double
  variance_explained : Double
  cumulative_variance : Double
}

///|
pub struct PcaIndResult {
  coord : Array[Array[Double]]
  cos2 : Array[Array[Double]]
  contrib : Array[Array[Double]]
}

///|
pub struct PcaVarResult {
  coord : Array[Array[Double]]
  cos2 : Array[Array[Double]]
  contrib : Array[Array[Double]]
}

///|
pub struct PcaResult {
  eigenvalues : Array[EigenValue]
  ind : PcaIndResult
  var : PcaVarResult
  n_ind : Int
  n_var : Int
  n_dims : Int
}

///|
/// Perform PCA on a data matrix.
/// Returns eigenvalues, individual coordinates, and variable loadings.
pub fn facto_pca(
  data : Array[Array[Double]],
  center? : Bool = true,
  scale? : Bool = true,
  ncp? : Int = 0,
) -> PcaResult {
  let n = data.length()
  if n == 0 {
    return PcaResult::{
      eigenvalues: [],
      ind: PcaIndResult::{ coord: [], cos2: [], contrib: [] },
      var: PcaVarResult::{ coord: [], cos2: [], contrib: [] },
      n_ind: 0,
      n_var: 0,
      n_dims: 0,
    }
  }
  let p = data[0].length()
  // Center and scale
  let centered = facto_center_scale(data, center, scale)
  // Compute covariance matrix (p x p)
  let cov : Array[Array[Double]] = Array::make(p, Array::make(p, 0.0))
  let mut i = 0
  while i < n {
    let mut j = 0
    while j < p {
      let mut k = 0
      while k < p {
        cov[j][k] = cov[j][k] + centered[i][j] * centered[i][k]
        k = k + 1
      }
      j = j + 1
    }
    i = i + 1
  }
  let scale_factor = 1.0 / (n.to_double() - 1.0)
  let mut ii = 0
  while ii < p {
    let mut jj = 0
    while jj < p { cov[ii][jj] = cov[ii][jj] * scale_factor; jj = jj + 1 }
    ii = ii + 1
  }
  // Compute eigenvalues via power iteration
  let n_components = if ncp > 0 {
    if ncp > p { p } else { ncp }
  } else {
    p
  }
  let (eigenvalues, eigenvectors) = facto_power_iteration(cov, n_components)
  // Sort by descending eigenvalue
  let sorted = facto_sort_eigen(eigenvalues, eigenvectors)
  let eig_vals = sorted.0
  let eig_vecs = sorted.1
  // Compute explained variance
  let mut tv = 0.0
  for ev in eig_vals { tv = tv + ev }
  // Compute individual scores: centered %*% eigenvectors
  let ind_coord : Array[Array[Double]] = Array::make(n, Array::make(n_components, 0.0))
  let mut si = 0
  while si < n {
    let mut sj = 0
    while sj < n_components {
      let mut s = 0.0
      let mut sk = 0
      while sk < p { s = s + centered[si][sk] * eig_vecs[sk][sj]; sk = sk + 1 }
      ind_coord[si][sj] = s
      sj = sj + 1
    }
    si = si + 1
  }
  // Variable coordinates: eigenvectors * sqrt(eigenvalues)
  let var_coord : Array[Array[Double]] = Array::make(p, Array::make(n_components, 0.0))
  let mut vi = 0
  while vi < p {
    let mut vj = 0
    while vj < n_components {
      let eig = eig_vals[vj]
      let sqrt_eig = if eig > 0.0 { eig.sqrt() } else { 0.0 }
      var_coord[vi][vj] = eig_vecs[vi][vj] * sqrt_eig
      vj = vj + 1
    }
    vi = vi + 1
  }
  // Build eigenvalue list
  let eigenvalues_list : Array[EigenValue] = Array::new()
  let mut cv : Double = 0.0
  let mut ei = 0
  while ei < n_components {
    let ev = eig_vals[ei]
    let var_explained = if tv > 0.0 { ev / tv * 100.0 } else { 0.0 }
    cv = cv + var_explained
    eigenvalues_list.push(EigenValue::{
      dimension: ei + 1,
      eigenvalue: ev,
      variance_explained: var_explained,
      cumulative_variance: cv,
    })
    ei = ei + 1
  }
  // Cosine quality for individuals
  let ind_cos2 = facto_cos2(ind_coord)
  // Cosine quality for variables
  let var_cos2 = facto_cos2(var_coord)
  // Contribution of individuals
  let ind_contrib = facto_contrib_ind(ind_coord)
  // Contribution of variables
  let var_contrib = facto_contrib_var(var_coord, eig_vals)
  PcaResult::{
    eigenvalues: eigenvalues_list,
    ind: PcaIndResult::{ coord: ind_coord, cos2: ind_cos2, contrib: ind_contrib },
    var: PcaVarResult::{ coord: var_coord, cos2: var_cos2, contrib: var_contrib },
    n_ind: n,
    n_var: p,
    n_dims: n_components,
  }
}

fn facto_center_scale(
  data : Array[Array[Double]],
  center : Bool,
  scale : Bool,
) -> Array[Array[Double]] {
  let n = data.length()
  if n == 0 { return [] }
  let p = data[0].length()
  let means : Array[Double] = Array::make(p, 0.0)
  let sds : Array[Double] = Array::make(p, 1.0)
  let mut j = 0
  while j < p {
    let mut s = 0.0
    let mut i = 0
    while i < n { s = s + data[i][j]; i = i + 1 }
    means[j] = if center { s / n.to_double() } else { 0.0 }
    if scale {
      let mut s2 = 0.0
      let mut i2 = 0
      while i2 < n {
        let d = data[i2][j] - means[j]
        s2 = s2 + d * d
        i2 = i2 + 1
      }
      let variance = if n < 2 { 1.0 } else { s2 / (n.to_double() - 1.0) }
      sds[j] = if variance > 0.0 { variance.sqrt() } else { 1.0 }
    }
    j = j + 1
  }
  let result : Array[Array[Double]] = Array::make(n, Array::make(p, 0.0))
  let mut i = 0
  while i < n {
    let mut j2 = 0
    while j2 < p {
      result[i][j2] = (data[i][j2] - means[j2]) / sds[j2]
      j2 = j2 + 1
    }
    i = i + 1
  }
  result
}

fn facto_power_iteration(
  a : Array[Array[Double]],
  n_components : Int,
) -> (Array[Double], Array[Array[Double]]) {
  let n = a.length()
  if n == 0 { return ([], []) }
  let max_iter = 200
  let tol = 1.0e-10
  let eigenvalues : Array[Double] = Array::make(n_components, 0.0)
  let eigenvectors : Array[Array[Double]] = Array::make(n, Array::make(n_components, 0.0))
  let mut i = 0
  while i < n_components {
    // Initialize vector
    let v : Array[Double] = Array::make(n, 0.0)
    v[0] = 1.0
    // Power iteration
    let mut lambda = 0.0
    let mut iter = 0
    while iter < max_iter {
      // v_new = A * v
      let v_new : Array[Double] = Array::make(n, 0.0)
      let mut j = 0
      while j < n {
        let mut s = 0.0
        let mut k = 0
        while k < n { s = s + a[j][k] * v[k]; k = k + 1 }
        v_new[j] = s
        j = j + 1
      }
      // Find infinity norm
      let mut max_val = 0.0
      let mut mj = 0
      while mj < n {
        let av = if v_new[mj] < 0.0 { -v_new[mj] } else { v_new[mj] }
        if av > max_val { max_val = av }
        mj = mj + 1
      }
      lambda = max_val
      // Normalize
      if max_val > 0.0 {
        let mut nj = 0
        while nj < n { v[nj] = v_new[nj] / max_val; nj = nj + 1 }
      }
      // Check convergence
      if iter > 0 && max_val > 0.0 {
        let mut diff = 0.0
        let mut dj = 0
        while dj < n {
          let d = v_new[dj] / max_val - eigenvectors[dj][i]
          let ad = if d < 0.0 { -d } else { d }
          if ad > diff { diff = ad }
          dj = dj + 1
        }
        if diff < tol { break }
      }
      iter = iter + 1
    }
    eigenvalues[i] = lambda
    let mut ei = 0
    while ei < n { eigenvectors[ei][i] = v[ei]; ei = ei + 1 }
    // Deflation: subtract contribution of found eigenvector
    if lambda > 0.0 {
      let mut di = 0
      while di < n {
        let mut dj = 0
        while dj < n {
          let contrib = v[di] * v[dj] * lambda
          a[di][dj] = a[di][dj] - contrib
          dj = dj + 1
        }
        di = di + 1
      }
    }
    i = i + 1
  }
  (eigenvalues, eigenvectors)
}

fn facto_sort_eigen(
  eigenvalues : Array[Double],
  eigenvectors : Array[Array[Double]],
) -> (Array[Double], Array[Array[Double]]) {
  let n = eigenvalues.length()
  let indices : Array[Int] = Array::make(n, 0)
  let mut i = 0
  while i < n { indices[i] = i; i = i + 1 }
  // Bubble sort by eigenvalue descending
  let mut iter = 0
  while iter < n - 1 {
    let mut swapped = false
    let mut j = 0
    while j < n - iter - 1 {
      if eigenvalues[indices[j]] < eigenvalues[indices[j + 1]] {
        let tmp = indices[j]
        indices[j] = indices[j + 1]
        indices[j + 1] = tmp
        swapped = true
      }
      j = j + 1
    }
    if !swapped { break }
    iter = iter + 1
  }
  let sorted_eigenvalues : Array[Double] = Array::make(n, 0.0)
  let sorted_eigenvectors : Array[Array[Double]] = Array::make(eigenvectors.length(), Array::make(n, 0.0))
  let mut k = 0
  while k < n {
    sorted_eigenvalues[k] = eigenvalues[indices[k]]
    let mut col = 0
    while col < eigenvectors.length() {
      sorted_eigenvectors[col][k] = eigenvectors[col][indices[k]]
      col = col + 1
    }
    k = k + 1
  }
  (sorted_eigenvalues, sorted_eigenvectors)
}

fn facto_cos2(coords : Array[Array[Double]]) -> Array[Array[Double]] {
  let n = coords.length()
  if n == 0 { return [] }
  let ncp = coords[0].length()
  let cos2 : Array[Array[Double]] = Array::make(n, Array::make(ncp, 0.0))
  let mut i = 0
  while i < n {
    let mut total = 0.0
    let mut j = 0
    while j < ncp { total = total + coords[i][j] * coords[i][j]; j = j + 1 }
    let mut k = 0
    while k < ncp {
      cos2[i][k] = if total > 0.0 { coords[i][k] * coords[i][k] / total } else { 0.0 }
      k = k + 1
    }
    i = i + 1
  }
  cos2
}

fn facto_contrib_ind(
  coords : Array[Array[Double]],
) -> Array[Array[Double]] {
  let n = coords.length()
  if n == 0 { return [] }
  let ncp = coords[0].length()
  let contrib : Array[Array[Double]] = Array::make(n, Array::make(ncp, 0.0))
  let total_contrib : Array[Double] = Array::make(ncp, 0.0)
  let mut i = 0
  while i < n {
    let mut j = 0
    while j < ncp {
      let c = coords[i][j] * coords[i][j]
      contrib[i][j] = c
      total_contrib[j] = total_contrib[j] + c
      j = j + 1
    }
    i = i + 1
  }
  let mut k = 0
  while k < ncp {
    if total_contrib[k] > 0.0 {
      let mut l = 0
      while l < n { contrib[l][k] = contrib[l][k] / total_contrib[k] * 100.0; l = l + 1 }
    }
    k = k + 1
  }
  contrib
}

fn facto_contrib_var(
  coords : Array[Array[Double]],
  eigenvalues : Array[Double],
) -> Array[Array[Double]] {
  let n = coords.length()
  if n == 0 { return [] }
  let ncp = coords[0].length()
  let contrib : Array[Array[Double]] = Array::make(n, Array::make(ncp, 0.0))
  let mut i = 0
  while i < n {
    let mut j = 0
    while j < ncp {
      let eig = eigenvalues[j]
      let c = if eig > 0.0 { coords[i][j] * coords[i][j] / eig * 100.0 } else { 0.0 }
      contrib[i][j] = c
      j = j + 1
    }
    i = i + 1
  }
  contrib
}

///|
/// Get eigenvalues and explained variance from PCA result.
pub fn facto_get_eigenvalue(result : PcaResult) -> Array[EigenValue] {
  result.eigenvalues
}

///|
/// Get individual coordinates.
pub fn facto_get_pca_ind(result : PcaResult) -> PcaIndResult {
  result.ind
}

///|
/// Get variable coordinates.
pub fn facto_get_pca_var(result : PcaResult) -> PcaVarResult {
  result.var
}

///|
/// Get contribution of variables to a specific dimension.
pub fn facto_var_contrib_dim(result : PcaResult, dim : Int) -> Array[(String, Double)] {
  let d = dim - 1
  if d < 0 || d >= result.n_dims { return [] }
  let contribs : Array[(String, Double)] = Array::new()
  let mut i = 0
  while i < result.n_var {
    contribs.push(("var_\{i + 1}", result.var.contrib[i][d]))
    i = i + 1
  }
  // Sort by contribution descending
  let n = contribs.length()
  let mut iter = 0
  while iter < n - 1 {
    let mut swapped = false
    let mut j = 0
    while j < n - iter - 1 {
      if contribs[j].1 < contribs[j + 1].1 {
        let tmp = contribs[j]
        contribs[j] = contribs[j + 1]
        contribs[j + 1] = tmp
        swapped = true
      }
      j = j + 1
    }
    if !swapped { break }
    iter = iter + 1
  }
  contribs
}

///|
/// Get contribution of individuals to a specific dimension.
pub fn facto_ind_contrib_dim(result : PcaResult, dim : Int) -> Array[(String, Double)] {
  let d = dim - 1
  if d < 0 || d >= result.n_dims { return [] }
  let contribs : Array[(String, Double)] = Array::new()
  let mut i = 0
  while i < result.n_ind {
    contribs.push(("ind_\{i + 1}", result.ind.contrib[i][d]))
    i = i + 1
  }
  let n = contribs.length()
  let mut iter = 0
  while iter < n - 1 {
    let mut swapped = false
    let mut j = 0
    while j < n - iter - 1 {
      if contribs[j].1 < contribs[j + 1].1 {
        let tmp = contribs[j]
        contribs[j] = contribs[j + 1]
        contribs[j + 1] = tmp
        swapped = true
      }
      j = j + 1
    }
    if !swapped { break }
    iter = iter + 1
  }
  contribs
}

///|
/// Get cosine quality of representation for variables on a specific dimension.
pub fn facto_var_cos2_dim(result : PcaResult, dim : Int) -> Array[(String, Double)] {
  let d = dim - 1
  if d < 0 || d >= result.n_dims { return [] }
  let cos2s : Array[(String, Double)] = Array::new()
  let mut i = 0
  while i < result.n_var {
    cos2s.push(("var_\{i + 1}", result.var.cos2[i][d]))
    i = i + 1
  }
  cos2s
}

///|
/// Get cosine quality of representation for individuals on a specific dimension.
pub fn facto_ind_cos2_dim(result : PcaResult, dim : Int) -> Array[(String, Double)] {
  let d = dim - 1
  if d < 0 || d >= result.n_dims { return [] }
  let cos2s : Array[(String, Double)] = Array::new()
  let mut i = 0
  while i < result.n_ind {
    cos2s.push(("ind_\{i + 1}", result.ind.cos2[i][d]))
    i = i + 1
  }
  cos2s
}

///|
pub struct DimDesc {
  dimension : Int
  var_name : String
  correlation : Double
  p_value : Double
}

///|
pub fn facto_dimdesc(
  data : Array[Array[Double]],
  result : PcaResult,
  dim : Int,
) -> Array[DimDesc] {
  let d = dim - 1
  if d < 0 || d >= result.n_dims { return [] }
  let n = data.length()
  let p = data[0].length()
  let desc : Array[DimDesc] = Array::new()
  let pc_scores : Array[Double] = Array::make(n, 0.0)
  let mut i = 0
  while i < n { pc_scores[i] = result.ind.coord[i][d]; i = i + 1 }
  let mut j = 0
  while j < p {
    let var_vals : Array[Double] = Array::make(n, 0.0)
    let mut k = 0
    while k < n { var_vals[k] = data[k][j]; k = k + 1 }
    let corr = facto_pearson_correlation(pc_scores, var_vals)
    let denom_base = 1.0 - corr * corr
    let t_stat = if denom_base > 0.0 { corr * ((n.to_double() - 2.0) / denom_base).sqrt() } else { 0.0 }
    let p_val = if t_stat >= 0.0 { @math.exp(-t_stat) } else { 1.0 }
    desc.push(DimDesc::{
      dimension: dim,
      var_name: "var_\{j + 1}",
      correlation: corr,
      p_value: p_val,
    })
    j = j + 1
  }
  let n_desc = desc.length()
  let mut iter = 0
  while iter < n_desc - 1 {
    let mut swapped = false
    let mut l = 0
    while l < n_desc - iter - 1 {
      let abs_a = if desc[l].correlation < 0.0 { -desc[l].correlation } else { desc[l].correlation }
      let abs_b = if desc[l + 1].correlation < 0.0 { -desc[l + 1].correlation } else { desc[l + 1].correlation }
      if abs_a < abs_b {
        let tmp = desc[l]
        desc[l] = desc[l + 1]
        desc[l + 1] = tmp
        swapped = true
      }
      l = l + 1
    }
    if !swapped { break }
    iter = iter + 1
  }
  desc
}

fn facto_pearson_correlation(x : Array[Double], y : Array[Double]) -> Double {
  let n = x.length()
  if n < 3 { return 0.0 }
  let mut mx = 0.0
  let mut my = 0.0
  let mut i = 0
  while i < n { mx = mx + x[i]; my = my + y[i]; i = i + 1 }
  mx = mx / n.to_double()
  my = my / n.to_double()
  let mut num = 0.0
  let mut dx = 0.0
  let mut dy = 0.0
  let mut j = 0
  while j < n {
    let vx = x[j] - mx
    let vy = y[j] - my
    num = num + vx * vy
    dx = dx + vx * vx
    dy = dy + vy * vy
    j = j + 1
  }
  let denom = (dx * dy).sqrt()
  if denom == 0.0 { return 0.0 }
  num / denom
}

///|
/// Compute total inertia (sum of eigenvalues).
pub fn facto_total_inertia(result : PcaResult) -> Double {
  let mut total = 0.0
  for ev in result.eigenvalues {
    total = total + ev.eigenvalue
  }
  total
}

///|
/// Get the number of dimensions to retain based on threshold.
pub fn facto_nb_dim(result : PcaResult, threshold : Double) -> Int {
  let mut cum = 0.0
  let mut i = 0
  while i < result.eigenvalues.length() {
    cum = cum + result.eigenvalues[i].variance_explained
    if cum >= threshold { return i + 1 }
    i = i + 1
  }
  result.eigenvalues.length()
}

///|
pub struct PcaSummary {
  n_individuals : Int
  n_variables : Int
  n_dimensions : Int
  total_inertia : Double
  first_dim_variance : Double
}

///|
pub fn facto_summary(result : PcaResult) -> PcaSummary {
  let total = facto_total_inertia(result)
  let first_var = if result.eigenvalues.length() > 0 { result.eigenvalues[0].variance_explained } else { 0.0 }
  PcaSummary::{
    n_individuals: result.n_ind,
    n_variables: result.n_var,
    n_dimensions: result.n_dims,
    total_inertia: total,
    first_dim_variance: first_var,
  }
}

///|
/// Compute loadings (rotation matrix) from variable coordinates.
pub fn facto_loadings(result : PcaResult) -> Array[Array[Double]] {
  result.var.coord
}

///|
/// Print a summary of the PCA result.
pub fn facto_print_summary(result : PcaResult) -> String {
  let s = facto_summary(result)
  let lines : Array[String] = Array::new()
  lines.push("PCA Summary:")
  lines.push("  Individuals: \{s.n_individuals}")
  lines.push("  Variables: \{s.n_variables}")
  lines.push("  Dimensions: \{s.n_dimensions}")
  lines.push("  Total inertia: \{s.total_inertia}")
  lines.push("  First dimension variance: \{s.first_dim_variance}%")
  lines.push("")
  lines.push("Eigenvalues:")
  for ev in result.eigenvalues {
    lines.push("  Dim \{ev.dimension}: eigenvalue=\{ev.eigenvalue}, var=\{ev.variance_explained}%, cum=\{ev.cumulative_variance}%")
  }
  lines.join("\n")
}
