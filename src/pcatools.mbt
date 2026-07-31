///|
/// Advanced PCA analysis & visualization - inspired by Bioconductor {PCAtools}.
///
/// References:
///   - https://bioconductor.org/packages/release/bioc/html/PCAtools.html
///   - Blighe K, Lun A (2024). PCAtools: Everything Principal Components Analysis.
///     R package version 2.16.0.

// ---------- helper linear algebra / stats primitives ----------

fn pca_is_na(x : Double) -> Bool { x.is_nan() }

fn pca_transpose(a : Array[Array[Double]]) -> Array[Array[Double]] {
  let n = a.length()
  if n == 0 { return [] }
  let m = a[0].length()
  let out : Array[Array[Double]] = Array::new()
  let mut j = 0
  while j < m {
    let row : Array[Double] = Array::new()
    let mut i = 0
    while i < n {
      row.push(a[i][j])
      i = i + 1
    }
    out.push(row)
    j = j + 1
  }
  out
}

fn pca_copy(a : Array[Array[Double]]) -> Array[Array[Double]] {
  let out : Array[Array[Double]] = Array::new()
  for row in a {
    let nr : Array[Double] = Array::new()
    for x in row { nr.push(x) }
    out.push(nr)
  }
  out
}

fn center_and_scale_cols(
  data : Array[Array[Double]],
  scale : Bool,
) -> (Array[Array[Double]], Array[Double], Array[Double]) {
  let n = data.length()
  if n == 0 { return ([], [], []) }
  let p = data[0].length()
  let means : Array[Double] = Array::new()
  let sds : Array[Double] = Array::new()
  let mut j = 0
  while j < p {
    let mut s = 0.0
    let mut c = 0
    let mut i = 0
    while i < n {
      let v = data[i][j]
      if !pca_is_na(v) {
        s = s + v
        c = c + 1
      }
      i = i + 1
    }
    let mu = if c == 0 { 0.0 } else { s / c.to_double() }
    means.push(mu)
    if scale {
      let mut s2 = 0.0
      let mut c2 = 0
      let mut i2 = 0
      while i2 < n {
        let v = data[i2][j]
        if !pca_is_na(v) {
          let d = v - mu
          s2 = s2 + d * d
          c2 = c2 + 1
        }
        i2 = i2 + 1
      }
      let variance = if c2 < 2 {
        1.0
      } else {
        s2 / (c2 - 1).to_double()
      }
      let sd = if variance <= 0.0 {
        1.0
      } else {
        variance.sqrt()
      }
      sds.push(sd)
    } else {
      sds.push(1.0)
    }
    j = j + 1
  }
  let out = pca_copy(data)
  let mut ii = 0
  while ii < n {
    let mut jj = 0
    while jj < p {
      let v = out[ii][jj]
      if pca_is_na(v) {
        out[ii][jj] = 0.0
      } else {
        out[ii][jj] = (v - means[jj]) / sds[jj]
      }
      jj = jj + 1
    }
    ii = ii + 1
  }
  (out, means, sds)
}

fn compute_covariance_matrix(
  x : Array[Array[Double]],
) -> Array[Array[Double]] {
  let n = x.length()
  if n == 0 { return [] }
  let p = x[0].length()
  let out : Array[Array[Double]] = Array::new()
  let denom = (if n - 1 < 1 { 1 } else { n - 1 }).to_double()
  let mut i = 0
  while i < p {
    let row : Array[Double] = Array::new()
    let mut j = 0
    while j < p {
      let mut s = 0.0
      let mut k = 0
      while k < n {
        s = s + x[k][i] * x[k][j]
        k = k + 1
      }
      row.push(s / denom)
      j = j + 1
    }
    out.push(row)
    i = i + 1
  }
  out
}

fn power_iteration_pca(
  cov : Array[Array[Double]],
  n_components : Int,
) -> (Array[Array[Double]], Array[Double]) {
  let p = cov.length()
  if p == 0 { return ([], []) }
  let k = if n_components < 1 {
    1
  } else if n_components > p {
    p
  } else {
    n_components
  }
  let eigvecs_rows : Array[Array[Double]] = Array::new()
  let eigvals : Array[Double] = Array::new()
  let c = pca_copy(cov)
  let zero_eps = 0.000000000001
  let mut comp = 0
  while comp < k {
    let mut v : Array[Double] = Array::new()
    let mut ii = 0
    while ii < p {
      let vv = if ii == comp {
        1.0
      } else {
        0.5.pow(comp.to_double())
      }
      v.push(vv)
      ii = ii + 1
    }
    let mut nrm_v = 0.0
    let mut t = 0
    while t < p {
      nrm_v = nrm_v + v[t] * v[t]
      t = t + 1
    }
    nrm_v = if nrm_v == 0.0 { 1.0 } else { nrm_v.sqrt() }
    let mut tt = 0
    while tt < p {
      v[tt] = v[tt] / nrm_v
      tt = tt + 1
    }
    let mut eig = 0.0
    let mut iter_count = 0
    while iter_count < 200 {
      let v_new : Array[Double] = Array::new()
      let mut iii = 0
      while iii < p {
        let mut dot = 0.0
        let mut jjj = 0
        while jjj < p {
          dot = dot + c[iii][jjj] * v[jjj]
          jjj = jjj + 1
        }
        v_new.push(dot)
        iii = iii + 1
      }
      let mut nrm2 = 0.0
      let mut kk = 0
      while kk < p {
        nrm2 = nrm2 + v_new[kk] * v_new[kk]
        kk = kk + 1
      }
      let nrm = if nrm2 == 0.0 { 0.0 } else { nrm2.sqrt() }
      if nrm < zero_eps { break }
      let mut kk2 = 0
      while kk2 < p {
        v_new[kk2] = v_new[kk2] / nrm
        kk2 = kk2 + 1
      }
      let mut rq = 0.0
      let mut iii2 = 0
      while iii2 < p {
        let mut cv = 0.0
        let mut jjj2 = 0
        while jjj2 < p {
          cv = cv + c[iii2][jjj2] * v_new[jjj2]
          jjj2 = jjj2 + 1
        }
        rq = rq + v_new[iii2] * cv
        iii2 = iii2 + 1
      }
      let mut diff = 0.0
      let mut i3 = 0
      while i3 < p {
        let d = v_new[i3] - v[i3]
        diff = diff + d * d
        i3 = i3 + 1
      }
      v = v_new
      eig = rq
      if diff.sqrt() < 0.00000001 { break }
      iter_count = iter_count + 1
    }
    if eig < 0.0 { eig = 0.0 }
    eigvecs_rows.push(v)
    eigvals.push(eig)
    // deflate
    let mut id = 0
    while id < p {
      let mut jd = 0
      while jd < p {
        c[id][jd] = c[id][jd] - eig * v[id] * v[jd]
        jd = jd + 1
      }
      id = id + 1
    }
    comp = comp + 1
  }
  let v_cols = pca_transpose(eigvecs_rows)
  (v_cols, eigvals)
}

fn pca_matmul(
  a : Array[Array[Double]],
  b : Array[Array[Double]],
) -> Array[Array[Double]] {
  let n = a.length()
  if n == 0 { return [] }
  let ka = a[0].length()
  let kb = b.length()
  let p = if kb == 0 { 0 } else { b[0].length() }
  let k = if ka < kb { ka } else { kb }
  let out : Array[Array[Double]] = Array::new()
  let mut i = 0
  while i < n {
    let row : Array[Double] = Array::new()
    let mut j = 0
    while j < p {
      let mut s = 0.0
      let mut t = 0
      while t < k {
        s = s + a[i][t] * b[t][j]
        t = t + 1
      }
      row.push(s)
      j = j + 1
    }
    out.push(row)
    i = i + 1
  }
  out
}

fn compute_variance_explained(eigvals : Array[Double]) -> Array[Double] {
  let mut total = 0.0
  for e in eigvals { total = total + e }
  let out : Array[Double] = Array::new()
  if total == 0.0 {
    for e in eigvals {
      let _e2 = e
      out.push(0.0)
    }
  } else {
    for e in eigvals { out.push(e / total) }
  }
  out
}

fn compute_cumulative(x : Array[Double]) -> Array[Double] {
  let out : Array[Double] = Array::new()
  let mut s = 0.0
  for v in x {
    s = s + v
    out.push(s)
  }
  out
}

// ---------- FullPCAResult ----------

pub struct FullPCAResult {
  pub scores : Array[Array[Double]]
  pub loadings : Array[Array[Double]]
  pub eigenvalues : Array[Double]
  pub variance_explained : Array[Double]
  pub cumulative_variance : Array[Double]
  pub n_samples : Int
  pub n_variables : Int
  pub n_components : Int
  pub center : Array[Double]
  pub scale : Array[Double]
  pub used_scaling : Bool
}

pub fn pcatools_run_pca(
  data : Array[Array[Double]],
  n_components : Int,
  scale : Bool,
) -> FullPCAResult {
  let (centered, means, sds) = center_and_scale_cols(data, scale)
  if centered.length() == 0 {
    return FullPCAResult::{
      scores: [],
      loadings: [],
      eigenvalues: [],
      variance_explained: [],
      cumulative_variance: [],
      n_samples: 0,
      n_variables: 0,
      n_components: 0,
      center: [],
      scale: [],
      used_scaling: scale,
    }
  }
  let cov = compute_covariance_matrix(centered)
  let (eigvecs, eigvals) = power_iteration_pca(cov, n_components)
  let scores = pca_matmul(centered, eigvecs)
  let variance_explained = compute_variance_explained(eigvals)
  let cumulative_variance = compute_cumulative(variance_explained)
  FullPCAResult::{
    scores: scores,
    loadings: eigvecs,
    eigenvalues: eigvals,
    variance_explained: variance_explained,
    cumulative_variance: cumulative_variance,
    n_samples: data.length(),
    n_variables: if data.length() == 0 { 0 } else { data[0].length() },
    n_components: eigvals.length(),
    center: means,
    scale: sds,
    used_scaling: scale,
  }
}

pub fn pcatools_summary(res : FullPCAResult) -> String {
  let buf = StringBuilder::new()
  buf.write_string("PCAResult: \{res.n_samples} samples, \{res.n_variables} variables\n")
  buf.write_string("  components = \{res.n_components}, scaling = \{res.used_scaling}\n")
  let k = if res.n_components < 8 { res.n_components } else { 8 }
  let mut i = 0
  while i < k {
    let ev = if i < res.eigenvalues.length() { res.eigenvalues[i] } else { 0.0 }
    let ve = if i < res.variance_explained.length() { res.variance_explained[i] } else { 0.0 }
    let cv = if i < res.cumulative_variance.length() { res.cumulative_variance[i] } else { 0.0 }
    buf.write_string("  PC\{i+1}: eig=\{ev}  var=\{ve}  cum=\{cv}\n")
    i = i + 1
  }
  if res.n_components > k {
    buf.write_string("  ... (\{res.n_components - k} more PCs)\n")
  }
  buf.to_string()
}

// ---------- scree_plot_ascii ----------

pub fn scree_plot_ascii(
  res : FullPCAResult,
  height? : Int = 8,
  width_per_pc? : Int = 3,
) -> String {
  let k = res.variance_explained.length()
  if k == 0 { return "(no PCs available for scree plot)".to_string() }
  let buf = StringBuilder::new()
  let mut max_ve = 0.0
  if res.variance_explained.length() > 0 {
    max_ve = res.variance_explained[0]
    for x in res.variance_explained {
      if x > max_ve { max_ve = x }
    }
  }
  let mm = if max_ve == 0.0 { 1.0 } else { max_ve }
  let h_use = if height < 2 { 2 } else { height }
  let mut r = 0
  while r < h_use {
    let frac = (h_use - 1 - r).to_double() / (h_use - 1).to_double()
    let threshold = frac * mm
    let mut i = 0
    while i < k {
      let cell = if res.variance_explained[i] >= threshold {
        "#"
      } else {
        " "
      }
      let mut w = 0
      while w < width_per_pc {
        buf.write_string(cell)
        w = w + 1
      }
      if i + 1 < k { buf.write_string(" ") }
      i = i + 1
    }
    buf.write_string("  \{threshold}\n")
    r = r + 1
  }
  buf.write_string("--PCs--")
  let mut i2 = 0
  while i2 < k {
    let label = "\{i2+1}"
    buf.write_string(" ")
    let pad = width_per_pc - label.length()
    let mut p = 0
    while p < pad {
      buf.write_string(" ")
      p = p + 1
    }
    buf.write_string(label)
    i2 = i2 + 1
  }
  buf.write_char('\n')
  buf.to_string()
}

// ---------- pca_biplot_ascii ----------

pub struct BiplotOptions {
  pub pc_x : Int
  pub pc_y : Int
  pub width : Int
  pub height : Int
  pub max_loadings : Int
}

pub fn BiplotOptions::new() -> BiplotOptions {
  BiplotOptions::{
    pc_x: 1,
    pc_y: 2,
    width: 60,
    height: 20,
    max_loadings: 8,
  }
}

pub fn pca_biplot_ascii(
  res : FullPCAResult,
  opts? : BiplotOptions = BiplotOptions::new(),
) -> String {
  let buf = StringBuilder::new()
  if res.scores.length() == 0 || res.n_components < 2 {
    return "(insufficient data for PCA biplot)".to_string()
  }
  let ix = if opts.pc_x < 1 {
    0
  } else if opts.pc_x > res.n_components {
    res.n_components - 1
  } else {
    opts.pc_x - 1
  }
  let iy = if opts.pc_y < 1 {
    1
  } else if opts.pc_y > res.n_components {
    res.n_components - 1
  } else {
    opts.pc_y - 1
  }
  let mut xmin = res.scores[0][ix]
  let mut xmax = xmin
  let mut ymin = res.scores[0][iy]
  let mut ymax = ymin
  for s in res.scores {
    if s[ix] < xmin { xmin = s[ix] }
    if s[ix] > xmax { xmax = s[ix] }
    if s[iy] < ymin { ymin = s[iy] }
    if s[iy] > ymax { ymax = s[iy] }
  }
  let xr = if xmax - xmin == 0.0 { 1.0 } else { xmax - xmin }
  let yr = if ymax - ymin == 0.0 { 1.0 } else { ymax - ymin }
  let w = if opts.width < 10 { 10 } else { opts.width }
  let h = if opts.height < 5 { 5 } else { opts.height }
  let grid : Array[Array[Int]] = Array::new()
  let mut _r = 0
  while _r < h {
    let row : Array[Int] = Array::new()
    let mut _cc = 0
    while _cc < w {
      row.push(' '.to_int())
      _cc = _cc + 1
    }
    grid.push(row)
    _r = _r + 1
  }
  // scores = 'o'
  for s in res.scores {
    let gx = ((s[ix] - xmin) / xr * (w - 1).to_double()).to_int()
    let gy = ((s[iy] - ymin) / yr * (h - 1).to_double()).to_int()
    let xc = if gx < 0 { 0 } else if gx >= w { w - 1 } else { gx }
    let gyy = if gy < 0 { 0 } else if gy >= h { h - 1 } else { gy }
    let gr = h - 1 - gyy
    grid[gr][xc] = 'o'.to_int()
  }
  // axes origin
  let ox = ((0.0 - xmin) / xr * (w - 1).to_double()).to_int()
  let oy = ((0.0 - ymin) / yr * (h - 1).to_double()).to_int()
  let oxc = if ox < 0 { 0 } else if ox >= w { w - 1 } else { ox }
  let oyc = if oy < 0 { 0 } else if oy >= h { h - 1 } else { oy }
  let oyr = h - 1 - oyc
  let mut c_axis = 0
  while c_axis < w {
    let ch = grid[oyr][c_axis]
    if ch == ' '.to_int() { grid[oyr][c_axis] = '-'.to_int() }
    c_axis = c_axis + 1
  }
  let mut r_axis = 0
  while r_axis < h {
    let ch = grid[r_axis][oxc]
    if ch == ' '.to_int() { grid[r_axis][oxc] = '|'.to_int() }
    r_axis = r_axis + 1
  }
  // loadings
  let loadings_nrow = res.loadings.length()
  if loadings_nrow > 0 && iy < res.loadings[0].length() {
    let ranked : Array[(Int, Double)] = Array::new()
    let mut idx = 0
    while idx < loadings_nrow {
      let lx = res.loadings[idx][ix]
      let ly = res.loadings[idx][iy]
      ranked.push((idx, lx * lx + ly * ly))
      idx = idx + 1
    }
    ranked.sort_by(fn(p : (Int, Double), q : (Int, Double)) -> Int {
      if p.1 > q.1 { -1 } else if p.1 < q.1 { 1 } else { 0 }
    })
    let lm = if opts.max_loadings < 1 { 1 } else { opts.max_loadings }
    let ntop = if lm < ranked.length() { lm } else { ranked.length() }
    let mut lmax = 0.0
    let mut t = 0
    while t < ntop {
      let idxi = ranked[t].0
      let lx = res.loadings[idxi][ix]
      let ly = res.loadings[idxi][iy]
      let m_val = (lx * lx + ly * ly).sqrt()
      if m_val > lmax { lmax = m_val }
      t = t + 1
    }
    if lmax == 0.0 { lmax = 1.0 }
    let mut t2 = 0
    while t2 < ntop {
      let idxi = ranked[t2].0
      let lx = res.loadings[idxi][ix] / lmax
      let ly = res.loadings[idxi][iy] / lmax
      let steps = 30
      let mut s_step = 0
      while s_step < steps {
        let frac = (s_step + 1).to_double() / steps.to_double()
        let lx_s = lx * frac
        let ly_s = ly * frac
        let gx = ((lx_s - xmin) / xr * (w - 1).to_double()).to_int()
        let gy = ((ly_s - ymin) / yr * (h - 1).to_double()).to_int()
        let xcc = if gx < 0 { 0 } else if gx >= w { w - 1 } else { gx }
        let gyy2 = if gy < 0 { 0 } else if gy >= h { h - 1 } else { gy }
        let grc = h - 1 - gyy2
        let cur = grid[grc][xcc]
        if cur == ' '.to_int() || cur == '-'.to_int() || cur == '|'.to_int() {
          grid[grc][xcc] = '*'.to_int()
        }
        s_step = s_step + 1
      }
      t2 = t2 + 1
    }
  }
  buf.write_string("PCA biplot - PC\{opts.pc_x} vs PC\{opts.pc_y}\n")
  buf.write_string("  PC\{opts.pc_x} range: [\{xmin}, \{xmax}]\n")
  buf.write_string("  PC\{opts.pc_y} range: [\{ymin}, \{ymax}]\n")
  let mut rr = 0
  while rr < h {
    let mut cc = 0
    while cc < w {
      let ch_int = grid[rr][cc]
      match ch_int.to_char() {
        Some(cv) => buf.write_char(cv)
        None => buf.write_char(' ')
      }
      cc = cc + 1
    }
    buf.write_char('\n')
    rr = rr + 1
  }
  buf.to_string()
}

// ---------- find_pca_outliers ----------

pub struct PCAOutlierResult {
  pub indices : Array[Int]
  pub scores : Array[Double]
  pub n_tested : Int
  pub cutoff : Double
}

fn quantile_simple(xs_in : Array[Double], q : Double) -> Double {
  if xs_in.length() == 0 { return 0.0 }
  let xs : Array[Double] = Array::new()
  for x in xs_in { xs.push(x) }
  xs.sort_by(fn(a : Double, b : Double) -> Int {
    if a < b { -1 } else if a > b { 1 } else { 0 }
  })
  let n = xs.length()
  let qq = if q < 0.0 { 0.0 } else if q > 1.0 { 1.0 } else { q }
  let pos = qq * (n - 1).to_double()
  let lo = pos.to_int()
  let hi = if lo + 1 < n { lo + 1 } else { n - 1 }
  let frac = pos - lo.to_double()
  xs[lo] * (1.0 - frac) + xs[hi] * frac
}

pub fn find_pca_outliers(
  res : FullPCAResult,
  n_pcs? : Int = 5,
  quantile_cutoff? : Double = 0.975,
) -> PCAOutlierResult {
  let n = res.scores.length()
  let k_total = res.n_components
  let k = if n_pcs < 1 {
    1
  } else if n_pcs > k_total {
    k_total
  } else {
    n_pcs
  }
  if n == 0 || k == 0 {
    return PCAOutlierResult::{
      indices: [],
      scores: [],
      n_tested: 0,
      cutoff: 0.0,
    }
  }
  let distances : Array[Double] = Array::new()
  let mut _tot = 0.0
  for s in res.scores {
    let mut d2 = 0.0
    let mut t = 0
    while t < k {
      let ev = if t < res.eigenvalues.length() && res.eigenvalues[t] > 0.0 {
        res.eigenvalues[t]
      } else {
        1.0
      }
      d2 = d2 + s[t] * s[t] / ev
      t = t + 1
    }
    distances.push(d2)
    _tot = _tot + d2
  }
  let cutoff = quantile_simple(distances, quantile_cutoff)
  let outliers : Array[Int] = Array::new()
  let out_scores : Array[Double] = Array::new()
  let mut i = 0
  while i < distances.length() {
    if distances[i] > cutoff {
      outliers.push(i)
      out_scores.push(distances[i])
    }
    i = i + 1
  }
  let _used = _tot
  PCAOutlierResult::{
    indices: outliers,
    scores: out_scores,
    n_tested: n,
    cutoff: cutoff,
  }
}

// ---------- variable_correlations ----------

pub fn variable_correlations(
  data : Array[Array[Double]],
  res : FullPCAResult,
  which_pcs? : Int = 5,
) -> Array[Array[Double]] {
  let n = data.length()
  if n == 0 { return [] }
  let p = data[0].length()
  let k = if which_pcs < 1 {
    1
  } else if which_pcs > res.n_components {
    res.n_components
  } else {
    which_pcs
  }
  let means : Array[Double] = Array::new()
  let sds : Array[Double] = Array::new()
  let mut j = 0
  while j < p {
    let mut s = 0.0
    let mut c = 0
    let mut i = 0
    while i < n {
      let v = data[i][j]
      if !pca_is_na(v) {
        s = s + v
        c = c + 1
      }
      i = i + 1
    }
    let m = if c == 0 { 0.0 } else { s / c.to_double() }
    means.push(m)
    let mut s2 = 0.0
    let mut c2 = 0
    let mut i2 = 0
    while i2 < n {
      let v = data[i2][j]
      if !pca_is_na(v) {
        let d = v - m
        s2 = s2 + d * d
        c2 = c2 + 1
      }
      i2 = i2 + 1
    }
    let variance = if c2 < 2 { 1.0 } else { s2 / (c2 - 1).to_double() }
    let sd = if variance <= 0.0 { 1.0 } else { variance.sqrt() }
    sds.push(sd)
    j = j + 1
  }
  let pc_means : Array[Double] = Array::new()
  let pc_sds : Array[Double] = Array::new()
  let mut t = 0
  while t < k {
    let mut s = 0.0
    let mut i = 0
    while i < n {
      s = s + res.scores[i][t]
      i = i + 1
    }
    let m = s / n.to_double()
    pc_means.push(m)
    let mut s2 = 0.0
    let mut i2 = 0
    while i2 < n {
      let d = res.scores[i2][t] - m
      s2 = s2 + d * d
      i2 = i2 + 1
    }
    let variance = if n < 2 { 1.0 } else { s2 / (n - 1).to_double() }
    let sd = if variance <= 0.0 { 1.0 } else { variance.sqrt() }
    pc_sds.push(sd)
    t = t + 1
  }
  let out : Array[Array[Double]] = Array::new()
  let mut jj = 0
  while jj < p {
    let row : Array[Double] = Array::new()
    let mut tt = 0
    while tt < k {
      let mut sum_prod = 0.0
      let mut c3 = 0
      let mut ii = 0
      while ii < n {
        let v = data[ii][jj]
        if !pca_is_na(v) {
          let zx = (v - means[jj]) / sds[jj]
          let zy = (res.scores[ii][tt] - pc_means[tt]) / pc_sds[tt]
          sum_prod = sum_prod + zx * zy
          c3 = c3 + 1
        }
        ii = ii + 1
      }
      let r = if c3 < 2 {
        0.0
      } else {
        sum_prod / (c3 - 1).to_double()
      }
      let rclamp = if r > 1.0 {
        1.0
      } else if r < -1.0 {
        -1.0
      } else {
        r
      }
      row.push(rclamp)
      tt = tt + 1
    }
    out.push(row)
    jj = jj + 1
  }
  out
}
