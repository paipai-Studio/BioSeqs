///|
/// Variance Stabilizing Normalization - inspired by the Bioconductor {vsn} package.
///
/// References:
///   - https://www.bioconductor.org/packages/release/bioc/html/vsn.html
///   - Huber W, von Heydebreck A, Sueltmann H, Poustka A, Vingron M. (2002)
///     Variance stabilization applied to microarray data calibration and to the
///     quantification of differential expression. Bioinformatics 18 Suppl 1:S96-S104.

// ---- helper ---------------------------------------------------------------

fn vsn_is_na(x : Double) -> Bool { x.is_nan() }

fn vsn_copy(mat : Array[Array[Double]]) -> Array[Array[Double]] {
  let n = mat.length()
  let out : Array[Array[Double]] = Array::new()
  let mut i = 0
  while i < n {
    let row = mat[i]
    let m = row.length()
    let nr : Array[Double] = Array::new()
    let mut j = 0
    while j < m {
      nr.push(row[j])
      j = j + 1
    }
    out.push(nr)
    i = i + 1
  }
  out
}

fn vsn_col_mean_sd(
  mat : Array[Array[Double]],
  j : Int,
) -> (Double, Double) {
  let mut s = 0.0
  let mut s2 = 0.0
  let mut c = 0
  let n = mat.length()
  let mut i = 0
  while i < n {
    let v = mat[i][j]
    if !vsn_is_na(v) {
      s = s + v
      s2 = s2 + v * v
      c = c + 1
    }
    i = i + 1
  }
  if c < 2 { return (0.0, 0.0) }
  let mu = s / c.to_double()
  let myvar = (s2 - s * s / c.to_double()) / (c - 1).to_double()
  let sd = if myvar > 0.0 { myvar.sqrt() } else { 0.0 }
  (mu, sd)
}

// ---- parameters & result --------------------------------------------------

pub struct VSNColParam {
  a : Double
  b : Double
}

pub struct VSNResult {
  params : Array[VSNColParam]
  reference_mean_sd : Array[(Double, Double)]
}

pub struct VSNControl {
  pub n_bins : Int
  pub min_bin_count : Int
  pub max_iter : Int
  pub lts_frac : Double
}

pub fn VSNControl::new() -> VSNControl {
  VSNControl::{
    n_bins: 20,
    min_bin_count: 5,
    max_iter: 20,
    lts_frac: 0.9,
  }
}

// ---- core transform -------------------------------------------------------

pub fn glog(x : Double, a : Double, b : Double) -> Double {
  let v = (a + b * x) / 2.0
  let inner = v + (v * v + 1.0).sqrt()
  if inner <= 0.0 { 0.0 } else { @math.ln(inner) }
}

pub fn glog_inv(y : Double, a : Double, b : Double) -> Double {
  let s = ((@math.exp(y) - @math.exp(-y)) / 2.0)
  if b == 0.0 { 0.0 } else { (2.0 * s - a) / b }
}

// ---- per-column parameter fitting ----------------------------------------

fn vsn_lm(
  x : Array[Double],
  y : Array[Double],
  lts_frac : Double,
  max_iter : Int,
) -> (Double, Double) {
  let n = x.length()
  if n == 0 { return (0.0, 1.0) }
  let mut sx = 0.0
  let mut sy = 0.0
  let mut sxx = 0.0
  let mut sxy = 0.0
  let mut i = 0
  while i < n {
    sx = sx + x[i]
    sy = sy + y[i]
    sxx = sxx + x[i] * x[i]
    sxy = sxy + x[i] * y[i]
    i = i + 1
  }
  let denom = (n.to_double() * sxx - sx * sx)
  let mut b = if denom == 0.0 {
    1.0
  } else {
    (n.to_double() * sxy - sx * sy) / denom
  }
  let mut a = (sy - b * sx) / n.to_double()
  let mut _it = 0
  while _it < max_iter {
    let resid : Array[(Int, Double)] = Array::new()
    let mut i2 = 0
    while i2 < n {
      let e = y[i2] - (a + b * x[i2])
      resid.push((i2, e * e))
      i2 = i2 + 1
    }
    resid.sort_by(fn(p : (Int, Double), q : (Int, Double)) -> Int {
      if p.1 < q.1 { -1 } else if p.1 > q.1 { 1 } else { 0 }
    })
    let keep = if n < 2 {
      n
    } else {
      let k = (lts_frac * n.to_double()).to_int()
      if k < 2 { 2 } else { k }
    }
    let mut sx2 = 0.0
    let mut sy2 = 0.0
    let mut sxx2 = 0.0
    let mut sxy2 = 0.0
    let k2 = if keep < n { keep } else { n }
    let mut t = 0
    while t < k2 {
      let i3 = resid[t].0
      sx2 = sx2 + x[i3]
      sy2 = sy2 + y[i3]
      sxx2 = sxx2 + x[i3] * x[i3]
      sxy2 = sxy2 + x[i3] * y[i3]
      t = t + 1
    }
    let kd = k2.to_double()
    let d2 = (kd * sxx2 - sx2 * sx2)
    let nb = if d2 == 0.0 { b } else { (kd * sxy2 - sx2 * sy2) / d2 }
    let na = if kd == 0.0 { a } else { (sy2 - nb * sx2) / kd }
    let da = (na - a).abs() + (nb - b).abs()
    a = na
    b = nb
    if da < 0.00000001 { break }
    _it = _it + 1
  }
  (a, b)
}

fn vsn_fit_col(
  mat : Array[Array[Double]],
  j : Int,
  ctrl : VSNControl,
) -> (VSNColParam, Array[(Double, Double)]) {
  let n = mat.length()
  let row_means : Array[Double] = Array::new()
  let mut i = 0
  while i < n {
    let mut s = 0.0
    let mut c = 0
    let row = mat[i]
    let m = row.length()
    let mut k = 0
    while k < m {
      let v = row[k]
      if !vsn_is_na(v) {
        s = s + v
        c = c + 1
      }
      k = k + 1
    }
    row_means.push(if c == 0 { 0.0 } else { s / c.to_double() })
    i = i + 1
  }
  let xs : Array[Double] = Array::new()
  let ys : Array[Double] = Array::new()
  let mut i2 = 0
  while i2 < n {
    let v = mat[i2][j]
    if !vsn_is_na(v) {
      let cv = v - row_means[i2]
      xs.push(row_means[i2])
      ys.push(cv.abs())
    }
    i2 = i2 + 1
  }
  let n_bins = if ctrl.n_bins < 1 { 1 } else { ctrl.n_bins }
  let pts = xs.length()
  let order : Array[(Double, Double)] = Array::new()
  let mut ii = 0
  while ii < pts {
    order.push((xs[ii], ys[ii]))
    ii = ii + 1
  }
  order.sort_by(fn(p : (Double, Double), q : (Double, Double)) -> Int {
    if p.0 < q.0 { -1 } else if p.0 > q.0 { 1 } else { 0 }
  })
  let bin_x : Array[Double] = Array::new()
  let bin_y : Array[Double] = Array::new()
  if pts > 0 {
    let bin_size = if pts / n_bins < ctrl.min_bin_count {
      if ctrl.min_bin_count < 1 { 1 } else { ctrl.min_bin_count }
    } else {
      pts / n_bins
    }
    let mut start = 0
    while start < pts {
      let end = if start + bin_size > pts { pts } else { start + bin_size }
      let cnt = end - start
      if cnt < ctrl.min_bin_count && bin_y.length() > 0 {
        let last_idx = bin_y.length() - 1
        let prev_cnt_d = bin_size.to_double()
        let mut mx = 0.0
        let mut my = 0.0
        let mut kk = 0
        let mut tt = start
        while tt < end {
          mx = mx + order[tt].0
          my = my + order[tt].1
          kk = kk + 1
          tt = tt + 1
        }
        if kk > 0 {
          let kd = kk.to_double()
          let new_len = (prev_cnt_d + kd)
          bin_x[last_idx] = (bin_x[last_idx] * prev_cnt_d + mx / kd * kd) / new_len
          bin_y[last_idx] = (bin_y[last_idx] * prev_cnt_d + my / kd * kd) / new_len
        }
        start = end
        continue
      }
      let mut mx = 0.0
      let mut my = 0.0
      let mut tt2 = start
      while tt2 < end {
        mx = mx + order[tt2].0
        my = my + order[tt2].1
        tt2 = tt2 + 1
      }
      let cnt_d = cnt.to_double()
      bin_x.push(mx / cnt_d)
      bin_y.push(my / cnt_d)
      start = end
    }
  }
  let (a, b_fit) = vsn_lm(bin_x, bin_y, ctrl.lts_frac, ctrl.max_iter)
  let bb = if b_fit <= 0.0 { 1.0 } else { b_fit }
  let ref_ms : Array[(Double, Double)] = Array::new()
  let mut t = 0
  while t < bin_x.length() {
    ref_ms.push((bin_x[t], bin_y[t]))
    t = t + 1
  }
  (VSNColParam::{ a: a, b: bb, }, ref_ms)
}

fn vsn_fit(mat : Array[Array[Double]], ctrl : VSNControl) -> VSNResult {
  if mat.length() == 0 {
    return VSNResult::{ params: [], reference_mean_sd: [], }
  }
  let cols = mat[0].length()
  let params : Array[VSNColParam] = Array::new()
  let mut ref_out : Array[(Double, Double)] = []
  let mut j = 0
  while j < cols {
    let (p, rms) = vsn_fit_col(mat, j, ctrl)
    params.push(p)
    if j == 0 { ref_out = rms }
    j = j + 1
  }
  VSNResult::{
    params: params,
    reference_mean_sd: ref_out,
  }
}

fn vsn_transform(
  mat : Array[Array[Double]],
  fit : VSNResult,
) -> Array[Array[Double]] {
  let out = vsn_copy(mat)
  let cols = fit.params.length()
  let n = out.length()
  let mut i = 0
  while i < n {
    let mut j = 0
    while j < cols {
      let v = out[i][j]
      let p = fit.params[j]
      if !vsn_is_na(v) {
        out[i][j] = glog(v, p.a, p.b)
      }
      j = j + 1
    }
    i = i + 1
  }
  out
}

pub fn vsn2(mat : Array[Array[Double]]) -> Array[Array[Double]] {
  let fit = vsn_fit(mat, VSNControl::new())
  vsn_transform(mat, fit)
}

pub fn vsn2_with_control(
  mat : Array[Array[Double]],
  ctrl : VSNControl,
) -> Array[Array[Double]] {
  let fit = vsn_fit(mat, ctrl)
  vsn_transform(mat, fit)
}

pub fn vsn_fit_and_report(mat : Array[Array[Double]]) -> VSNResult {
  vsn_fit(mat, VSNControl::new())
}

// ---- ASCII Mean-SD plot ---------------------------------------------------

pub struct MeanSDBin {
  mean_x : Double
  sd_y : Double
  count : Int
}

pub fn mean_sd_bins(
  mat : Array[Array[Double]],
  n_bins : Int,
) -> Array[MeanSDBin] {
  let n = mat.length()
  if n == 0 { return Array::new() }
  let cols = mat[0].length()
  let pts : Array[(Double, Double)] = Array::new()
  let mut i = 0
  while i < n {
    let mut j = 0
    while j < cols {
      let v = mat[i][j]
      if !vsn_is_na(v) {
        let (mu, sd) = vsn_col_mean_sd(mat, j)
        pts.push((mu, sd))
      }
      j = j + 1
    }
    i = i + 1
  }
  pts.sort_by(fn(p : (Double, Double), q : (Double, Double)) -> Int {
    if p.0 < q.0 { -1 } else if p.0 > q.0 { 1 } else { 0 }
  })
  let nb = if n_bins < 1 { 1 } else { n_bins }
  let total = pts.length()
  let out : Array[MeanSDBin] = Array::new()
  if total == 0 { return out }
  let step = total / nb
  let mut start = 0
  while start < total {
    let end = if start + step > total { total } else { start + step }
    let mut sx = 0.0
    let mut sy = 0.0
    let cnt = end - start
    let mut t = start
    while t < end {
      sx = sx + pts[t].0
      sy = sy + pts[t].1
      t = t + 1
    }
    let cnt_d = cnt.to_double()
    out.push(MeanSDBin::{
      mean_x: sx / cnt_d,
      sd_y: sy / cnt_d,
      count: cnt,
    })
    start = end
  }
  out
}

pub fn mean_sd_ascii(
  mat : Array[Array[Double]],
  n_bins? : Int = 40,
) -> String {
  let bins = mean_sd_bins(mat, n_bins)
  let width = 40
  let height = 10
  if bins.length() == 0 { return "(empty mean-SD plot)".to_string() }
  let mut xmin = bins[0].mean_x
  let mut xmax = bins[0].mean_x
  let mut ymin = bins[0].sd_y
  let mut ymax = bins[0].sd_y
  for b in bins {
    if b.mean_x < xmin { xmin = b.mean_x }
    if b.mean_x > xmax { xmax = b.mean_x }
    if b.sd_y < ymin { ymin = b.sd_y }
    if b.sd_y > ymax { ymax = b.sd_y }
  }
  let xr = if xmax - xmin == 0.0 { 1.0 } else { xmax - xmin }
  let yr = if ymax - ymin == 0.0 { 1.0 } else { ymax - ymin }
  let grid : Array[Array[Int]] = Array::make(
    height,
    Array::make(width, ' '.to_int())
  )
  let k = bins.length()
  let mut idx = 0
  while idx < k {
    let b = bins[idx]
    let xi = ((b.mean_x - xmin) / xr * (width - 1).to_double()).to_int()
    let yi = ((b.sd_y - ymin) / yr * (height - 1).to_double()).to_int()
    let xii = if xi < 0 { 0 } else { if xi >= width { width - 1 } else { xi } }
    let yii = if yi < 0 { 0 } else { if yi >= height { height - 1 } else { yi } }
    let row_idx = height - 1 - yii
    let c = grid[row_idx][xii]
    let nc = if c == ' '.to_int() {
      '*'.to_int()
    } else {
      if c == '*'.to_int() { '+'.to_int() } else { '#'.to_int() }
    }
    grid[row_idx][xii] = nc
    idx = idx + 1
  }
  let buf = StringBuilder::new()
  buf.write_string("mean-SD plot (mean=[\{xmin},\{xmax}] sd=[\{ymin},\{ymax}])\n")
  let mut r = 0
  while r < height {
    let row = grid[r]
    let mut c = 0
    while c < width {
      let ch = row[c]
      if ch == ' '.to_int() { buf.write_char(' ') }
      else if ch == '*'.to_int() { buf.write_char('*') }
      else if ch == '+'.to_int() { buf.write_char('+') }
      else { buf.write_char('#') }
      c = c + 1
    }
    buf.write_char('\n')
    r = r + 1
  }
  buf.to_string()
}

pub fn summarize_vsn_fit(fit : VSNResult) -> String {
  let buf = StringBuilder::new()
  let n = fit.params.length()
  buf.write_string("VSNResult: \{n} columns fitted\n")
  let mut i = 0
  while i < 5 && i < n {
    let p = fit.params[i]
    buf.write_string("  col \{i}: a=\{p.a}  b=\{p.b}\n")
    i = i + 1
  }
  if n > 5 { buf.write_string("  ... (\{n - 5} more)\n") }
  buf.to_string()
}

pub fn vsn_denoise(
  mat : Array[Array[Double]],
  frac_of_median_var? : Double = 0.05,
) -> Array[Array[Double]] {
  let out = vsn_copy(mat)
  if out.length() == 0 { return out }
  let cols = out[0].length()
  let vars : Array[Double] = Array::new()
  let mut j = 0
  while j < cols {
    let (_, sd) = vsn_col_mean_sd(out, j)
    vars.push(sd * sd)
    j = j + 1
  }
  let sorted_vars : Array[Double] = Array::new()
  for v in vars { sorted_vars.push(v) }
  sorted_vars.sort_by(fn(a : Double, b : Double) -> Int {
    if a < b { -1 } else if a > b { 1 } else { 0 }
  })
  let med = if sorted_vars.length() == 0 {
    0.0
  } else {
    sorted_vars[sorted_vars.length() / 2]
  }
  let threshold = med * frac_of_median_var
  let n = out.length()
  let mut i = 0
  while i < n {
    let mut jj = 0
    while jj < cols {
      if vars[jj] < threshold { out[i][jj] = 0.0 }
      jj = jj + 1
    }
    i = i + 1
  }
  out
}
