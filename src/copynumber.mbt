// Copy-number segmentation - port of Bioconductor `copynumber`.
//
// Implements Piecewise Constant Fitting (PCF) of Nilsen and Liestol et al.
// (2012), BMC Genomics: the exact Potts-filter dynamic program used by
// copynumber::pcf()/pcfPlain(), the multi-sample joint segmentation used by
// multipcf(), and the MAD/PCF winsorization of outliers. PCF complements the
// Circular Binary Segmentation implementation in dnacopy.mbt: CBS finds
// changepoints by permutation tests, while PCF minimizes a penalized least
// squares criterion (each break costs `gamma`).
//
// Every numeric result here was cross-checked against the unmodified
// Bioconductor copynumber R functions (exactPcf, doMultiPCF, madWins,
// pcfWins, medianFilter/getMad) to within floating point round-off.

///|
/// Scale factor of the MAD estimator, identical to R `mad()` default
/// (1.4826 makes the MAD a consistent estimator of the Gaussian SD).
pub let copy_mad_const : Double = 1.4826

///|
/// Median of three values, port of the `med3` helper of R `smoothEnds()`.
fn copy_med3(a : Double, b : Double, c : Double) -> Double {
  if a < b {
    if c < b {
      if a >= c {
        a
      } else {
        c
      }
    } else {
      b
    }
  } else if c > b {
    if a <= c {
      a
    } else {
      c
    }
  } else {
    b
  }
}

///|
/// Exact median of an odd-length sub-array, equivalent to R
/// `sort(x, partial = h)[h]`.
fn copy_med_odd_range(
  x : Array[Double],
  lo : Int,
  count : Int,
) -> Double {
  let window = x[lo:lo + count].to_owned()
  window.sort()
  window[count / 2]
}

///|
/// Running median smoother used by the copynumber package: a literal port of
/// `copynumber::medianFilter`, i.e. R `runmed(x, k = 2*k+1, endrule =
/// "median")`. Interior probes get the exact median of the full window; the
/// edges are first kept ("keep" rule) and then post-processed by
/// `stats::smoothEnds`, including Tukey's end-point rule with linear
/// extrapolation. `k` is the half-window size (default 25, window width 51).
pub fn copy_median_filter(x : Array[Double], k? : Int = 25) -> Array[Double] {
  let n = x.length()
  if n == 0 {
    return []
  }
  let mut width = 2 * k + 1
  if width > n {
    // runmed requires an odd width not exceeding n
    width = if n % 2 == 0 { n - 1 } else { n }
  }
  let h = width / 2
  // end_rule = "keep": edge observations are copied verbatim
  let y = x.copy()
  let mut i = h
  while i < n - h {
    y[i] = copy_med_odd_range(x, i - h, 2 * h + 1)
    i = i + 1
  }
  if h < 1 {
    return y
  }
  // smoothEnds(y, k = width)
  let sm = y.copy()
  if h >= 2 {
    sm[1] = copy_med3(y[0], y[1], y[2])
    sm[n - 2] = copy_med3(y[n - 1], y[n - 2], y[n - 3])
    if h >= 3 {
      let mut ii = 3 // 1-based index as in the R for-loop
      while ii <= h {
        let j = 2 * ii - 1
        sm[ii - 1] = copy_med_odd_range(y, 0, j)
        sm[n - ii] = copy_med_odd_range(y, n - j, j)
        ii = ii + 1
      }
    }
  }
  // Tukey end points: median of the observation, the next smoothed value and
  // its linear extrapolation
  sm[0] = copy_med3(y[0], sm[1], sm[1] - 2.0 * (sm[2] - sm[1]))
  sm[n - 1] = copy_med3(
    y[n - 1],
    sm[n - 2],
    sm[n - 2] - 2.0 * (sm[n - 3] - sm[n - 2]),
  )
  sm
}

///|
/// R `mad(x)`: 1.4826 times the median absolute deviation from the median.
pub fn copy_mad(x : Array[Double]) -> Double {
  let center = stat_median(x)
  let dev : Array[Double] = []
  for v in x {
    dev.push((v - center).abs())
  }
  copy_mad_const * stat_median(dev)
}

///|
/// Robust SD estimate `copynumber::getMad`: zero observations (typically
/// imputed missing values) are discarded, the remainder is smoothed with the
/// running median, and the MAD of the residuals is returned.
pub fn copy_get_mad(x : Array[Double], k? : Int = 25) -> Double {
  let nonzero : Array[Double] = []
  for v in x {
    if v != 0.0 {
      nonzero.push(v)
    }
  }
  if nonzero.length() == 0 {
    return 0.0
  }
  let runmed = copy_median_filter(nonzero, k=k)
  let dif : Array[Double] = []
  let mut i = 0
  while i < nonzero.length() {
    dif.push(nonzero[i] - runmed[i])
    i = i + 1
  }
  copy_mad(dif)
}

///|
/// One fitted constant region (segment) of a PCF fit.
pub struct PcfSegment {
  /// 1-based index of the first probe in the segment
  start : Int
  /// number of probes in the segment
  length : Int
  /// fitted value (segment mean)
  mean : Double
}

///|
/// Result of a single-sample piecewise constant fit, analogous to the list
/// returned by copynumber `exactPcf()`/`pcfPlain(return.est = TRUE)`.
pub struct PcfResult {
  segments : Array[PcfSegment]
  n_intervals : Int
  /// fitted value at every probe (yhat)
  yhat : Array[Double]
}

///|
/// Construct a one-segment result covering all probes - used for short
/// sequences and for zero variance.
fn pcf_single_segment(y : Array[Double]) -> PcfResult {
  let mut sum = 0.0
  for v in y {
    sum += v
  }
  let mean = sum / y.length().to_double()
  let yhat = Array::make(y.length(), mean)
  { segments: [{ start: 1, length: y.length(), mean }], n_intervals: 1, yhat }
}

///|
/// Exact Piecewise Constant Fitting, a literal port of copynumber
/// `exactPcf(y, kmin, gamma, yest = TRUE)`. The sequence is partitioned into
/// intervals of at least `kmin` probes minimizing
/// `sum SSE + gamma * n_breaks`; the Potts-filter dynamic program evaluates
/// every admissible last breakpoint. `gamma` is the penalty per break
/// (default 40, as in copynumber).
pub fn exact_pcf(
  y : Array[Double],
  kmin? : Int = 5,
  gamma? : Double = 40.0,
) -> PcfResult raise SeqError {
  let n = y.length()
  if n == 0 {
    raise SeqError("exact_pcf requires at least one observation")
  }
  if kmin < 1 {
    raise SeqError("kmin must be positive")
  }
  // R: fewer than 2*kmin probes -> one segment, the overall mean
  if n < 2 * kmin {
    return pcf_single_segment(y)
  }
  // DP vectors use 1-based indexing, mirroring the R source verbatim.
  let mut init_sum = 0.0
  let mut i = 0
  while i < kmin {
    init_sum += y[i]
    i = i + 1
  }
  let init_ave = init_sum / kmin.to_double()
  let best_cost = Array::make(n + 1, 0.0)
  let best_split = Array::make(n + 1, 0)
  let best_aver = Array::make(n + 1, 0.0)
  best_cost[kmin] = -init_sum * init_ave
  best_aver[kmin] = init_ave
  let sums = Array::make(n + 1, 0.0)
  let aver = Array::make(n + 1, 0.0)
  let cand = Array::make(n + 1, 0.0)
  let kmin1 = kmin + 1
  let mut kk = kmin1
  while kk < 2 * kmin {
    let mut q = kmin1
    while q <= kk {
      sums[q] += y[kk - 1]
      aver[q] = sums[q] / (kk - q + 1).to_double()
      q = q + 1
    }
    best_aver[kk] = (init_sum + sums[kmin1]) / kk.to_double()
    best_cost[kk] = -kk.to_double() * best_aver[kk] * best_aver[kk]
    kk = kk + 1
  }
  let mut nn = 2 * kmin
  while nn <= n {
    let yn = y[nn - 1]
    let mut q = kmin1
    while q <= nn {
      sums[q] += yn
      aver[q] = sums[q] / (nn - q + 1).to_double()
      q = q + 1
    }
    let upper = nn - kmin + 1
    let mut t = kmin1
    while t <= upper {
      cand[t] = best_cost[t - 1] - sums[t] * aver[t] + gamma
      t = t + 1
    }
    // which.min(): first position attaining the minimum
    let mut pos = kmin1
    t = kmin1 + 1
    while t <= upper {
      if cand[t] < cand[pos] {
        pos = t
      }
      t = t + 1
    }
    let mut cost = cand[pos]
    let mut mu = aver[pos]
    // compare with taking the whole prefix as a single interval
    let tot_ave = (sums[kmin1] + init_sum) / nn.to_double()
    let tot_cost = -nn.to_double() * tot_ave * tot_ave
    if tot_cost < cost {
      pos = 1
      cost = tot_cost
      mu = tot_ave
    }
    best_cost[nn] = cost
    best_aver[nn] = mu
    best_split[nn] = pos - 1
    nn = nn + 1
  }
  // Fill fitted values by backtracking
  let yhat = Array::make(n + 1, 0.0)
  let mut cur = n
  let mut nint = 0
  while cur > 0 {
    let val = best_aver[cur]
    let mut ii = best_split[cur] + 1
    while ii <= cur {
      yhat[ii] = val
      ii = ii + 1
    }
    cur = best_split[cur]
    nint = nint + 1
  }
  // Recover segment starts/lengths/means (backtracking yields them reversed)
  let starts = Array::make(nint, 0)
  let lengths = Array::make(nint, 0)
  let means = Array::make(nint, 0.0)
  cur = n
  let mut old = cur
  let mut cnt = nint
  while cur > 0 {
    starts[cnt - 1] = best_split[cur] + 1
    lengths[cnt - 1] = old - best_split[cur]
    means[cnt - 1] = best_aver[cur]
    cur = best_split[cur]
    old = cur
    cnt = cnt - 1
  }
  let segments : Array[PcfSegment] = []
  let mut s = 0
  while s < nint {
    segments.push({ start: starts[s], length: lengths[s], mean: means[s] })
    s = s + 1
  }
  { segments, n_intervals: nint, yhat: yhat[1:n + 1].to_owned() }
}

///|
/// User-level single-sample PCF, port of the per-sample logic in
/// `copynumber::pcfPlain`. With `normalize = TRUE` (the copynumber default)
/// the penalty is scaled by the robust variance, `gamma * getMad(y)^2`; with
/// zero variance the sample is returned as one constant segment. The exact
/// dynamic program is always used (copynumber switches to a fast
/// approximation above 400 probes; the exact solution is the reference
/// optimum).
pub fn pcf_plain(
  y : Array[Double],
  kmin? : Int = 5,
  gamma? : Double = 40.0,
  normalize? : Bool = true,
  k? : Int = 25,
) -> PcfResult raise SeqError {
  let mut use_gamma = gamma
  if normalize {
    let sd = copy_get_mad(y, k=k)
    use_gamma = gamma * sd * sd
    if use_gamma == 0.0 || use_gamma.is_nan() {
      return pcf_single_segment(y)
    }
  }
  exact_pcf(y, kmin=kmin, gamma=use_gamma)
}

///|
/// Round to `digits` decimal places using round-half-to-even, matching R
/// `round()` for the outlier comparisons.
fn copy_round(x : Double, digits : Int) -> Double {
  let scale = @math.pow(10.0, digits.to_double())
  let q = x * scale
  let floor_q = q.floor()
  let rem = q - floor_q
  let nearest = if rem < 0.5 {
    floor_q
  } else if rem > 0.5 {
    floor_q + 1.0
  } else if floor_q.to_int() % 2 == 0 {
    floor_q
  } else {
    floor_q + 1.0
  }
  nearest / scale
}

///|
/// Winsorization result: clipped values, the robust SD and the per-probe
/// outlier flag (+1 high outlier, -1 low outlier, 0 ordinary probe).
pub struct WinsResult {
  ywin : Array[Double]
  sdev : Double
  outliers : Array[Int]
}

///|
/// MAD winsorization, a literal port of `copynumber::madWins`: each probe is
/// replaced by its running median plus its residual clipped to
/// `+- tau * MAD(residuals)` (copynumber default `tau = 2.5`).
pub fn winsorize_mad(
  x : Array[Double],
  tau? : Double = 2.5,
  k? : Int = 25,
  digits? : Int = 4,
) -> WinsResult {
  let n = x.length()
  let xhat = copy_median_filter(x, k=k)
  let d : Array[Double] = []
  let mut i = 0
  while i < n {
    d.push(x[i] - xhat[i])
    i = i + 1
  }
  let sdev = copy_mad(d)
  let z = tau * sdev
  let ywin : Array[Double] = []
  while ywin.length() < n {
    let idx = ywin.length()
    let clipped = if d[idx] < -z {
      -z
    } else if d[idx] > z {
      z
    } else {
      d[idx]
    }
    ywin.push(xhat[idx] + clipped)
  }
  let outliers = Array::make(n, 0)
  i = 0
  while i < n {
    let rx = copy_round(x[i], digits)
    let rw = copy_round(ywin[i], digits)
    if rx > rw {
      outliers[i] = 1
    } else if rx < rw {
      outliers[i] = -1
    }
    i = i + 1
  }
  { ywin, sdev, outliers }
}

///|
/// PCF winsorization, a literal port of `copynumber::pcfWins`: residuals are
/// clipped against the running median and then a PCF fit is used as the
/// baseline for `iter` iterations (copynumber default one), after which final
/// outliers are flagged against the piecewise constant fit.
pub fn winsorize_pcf(
  x : Array[Double],
  tau? : Double = 2.5,
  k? : Int = 25,
  gamma? : Double = 40.0,
  iter? : Int = 1,
  digits? : Int = 4,
) -> WinsResult raise SeqError {
  let n = x.length()
  let mut xhat = copy_median_filter(x, k=k)
  let mut xwin = x.copy()
  let mut sdev = 0.0
  let mut j = 0
  while j < iter {
    let d : Array[Double] = []
    let mut i = 0
    while i < n {
      d.push(x[i] - xhat[i])
      i = i + 1
    }
    sdev = copy_mad(d)
    let z = tau * sdev
    xwin = []
    i = 0
    while i < n {
      let clipped = if d[i] < -z {
        -z
      } else if d[i] > z {
        z
      } else {
        d[i]
      }
      xwin.push(xhat[i] + clipped)
      i = i + 1
    }
    xhat = exact_pcf(xwin, kmin=5, gamma=gamma * sdev * sdev).yhat
    j = j + 1
  }
  // final clipping against the PCF baseline
  let d : Array[Double] = []
  let mut i = 0
  while i < n {
    d.push(xwin[i] - xhat[i])
    i = i + 1
  }
  sdev = copy_mad(d)
  let z = tau * sdev
  let ywin : Array[Double] = []
  i = 0
  while i < n {
    let residual = x[i] - xhat[i]
    let clipped = if residual < -z {
      -z
    } else if residual > z {
      z
    } else {
      residual
    }
    ywin.push(xhat[i] + clipped)
    i = i + 1
  }
  let outliers = Array::make(n, 0)
  i = 0
  while i < n {
    let rx = copy_round(x[i], digits)
    let rw = copy_round(ywin[i], digits)
    if rx > rw {
      outliers[i] = 1
    } else if rx < rw {
      outliers[i] = -1
    }
    i = i + 1
  }
  { ywin, sdev, outliers }
}

///|
/// One joint segment of a multi-sample PCF fit; `means` holds one fitted
/// value per sample.
pub struct MultiPcfSegment {
  start : Int
  length : Int
  means : Array[Double]
}

///|
/// Result of a multi-sample joint piecewise constant fit
/// (`copynumber::multipcf`): all samples share the same breakpoints.
pub struct MultiPcfResult {
  segments : Array[MultiPcfSegment]
  n_intervals : Int
  /// fitted values, indexed `yhat[sample][probe]`
  yhat : Array[Array[Double]]
  /// per-sample robust SD used for normalization (1 when not normalized)
  sds : Array[Double]
}

///|
/// Multi-sample joint PCF, a port of `copynumber::multipcf`/`doMultiPCF`.
/// `samples[s][p]` is the copy-number value of sample `s` at probe `p`; all
/// samples must contain the same number of complete (non-NaN) probes. The
/// samples share breakpoints; the penalty is scaled by the number of samples
/// (`gamma * n_samples`, as in copynumber). With `normalize = TRUE` each
/// sample is divided by its MAD-SD before fitting (and optionally multiplied
/// by `weights`); fitted means are converted back afterwards. A zero
/// variance sample yields a single constant segment, mirroring copynumber.
pub fn multi_pcf(
  samples : Array[Array[Double]],
  gamma? : Double = 40.0,
  normalize? : Bool = true,
  k? : Int = 25,
  weights? : Array[Double] = [],
) -> MultiPcfResult raise SeqError {
  let ns = samples.length()
  if ns == 0 {
    raise SeqError("multi_pcf requires at least one sample")
  }
  let n = samples[0].length()
  if n == 0 {
    raise SeqError("multi_pcf requires at least one probe")
  }
  for sample in samples {
    if sample.length() != n {
      raise SeqError("all samples must have the same number of probes")
    }
    for v in sample {
      if v.is_nan() {
        raise SeqError(
          "multi_pcf cannot be run with missing values; impute first",
        )
      }
    }
  }
  let w = if weights.length() == 0 {
    Array::make(ns, 1.0)
  } else if weights.length() != ns {
    raise SeqError("weights must have one entry per sample")
  } else {
    weights
  }
  let sd = Array::make(ns, 1.0)
  if normalize {
    let mut j = 0
    while j < ns {
      sd[j] = copy_get_mad(samples[j], k=k)
      j = j + 1
    }
  }
  let total_mean : Array[Double] = []
  let mut zero_sd = false
  for sdv in sd {
    if sdv == 0.0 || sdv.is_nan() {
      zero_sd = true
    }
  }
  // degenerate case: copynumber returns the per-sample overall mean
  if zero_sd {
    let mut j = 0
    while j < ns {
      let mut sum = 0.0
      for v in samples[j] {
        sum += v
      }
      total_mean.push(sum / n.to_double())
      j = j + 1
    }
    let yhat0 : Array[Array[Double]] = []
    j = 0
    while j < ns {
      yhat0.push(Array::make(n, total_mean[j]))
      j = j + 1
    }
    let segments0 : Array[MultiPcfSegment] = []
    segments0.push({ start: 1, length: n, means: total_mean })
    return { segments: segments0, n_intervals: 1, yhat: yhat0, sds: sd }
  }
  // scale once: z[j][p] = sample / sd * weight
  let z : Array[Array[Double]] = []
  let mut j = 0
  while j < ns {
    let row = Array::make(n, 0.0)
    let mut p = 0
    while p < n {
      row[p] = samples[j][p] / sd[j] * w[j]
      p = p + 1
    }
    z.push(row)
    j = j + 1
  }
  let scaled_gamma = gamma * ns.to_double()
  // doMultiPCF: samples in rows, probes in columns.
  let best_cost = Array::make(n + 1, 0.0)
  let best_split = Array::make(n + 1, 0)
  let best_aver : Array[Array[Double]] = []
  let sums : Array[Array[Double]] = []
  let aver : Array[Array[Double]] = []
  let mut jj = 0
  while jj < ns {
    best_aver.push(Array::make(n + 1, 0.0))
    sums.push(Array::make(n + 1, 0.0))
    aver.push(Array::make(n + 1, 0.0))
    jj = jj + 1
  }
  let mut nn = 1
  while nn <= n {
    let mut q = 1
    while q <= nn {
      let mut s = 0
      while s < ns {
        sums[s][q] += z[s][nn - 1]
        aver[s][q] = sums[s][q] / (nn - q + 1).to_double()
        s = s + 1
      }
      q = q + 1
    }
    // candidate costs: q == 1 is the whole prefix (no break, no gamma)
    let mut best_q = 1
    let mut best = 0.0
    let mut s = 0
    while s < ns {
      best -= sums[s][1] * aver[s][1]
      s = s + 1
    }
    q = 2
    while q <= nn {
      let mut c = best_cost[q - 1] + scaled_gamma
      s = 0
      while s < ns {
        c -= sums[s][q] * aver[s][q]
        s = s + 1
      }
      if c < best {
        best = c
        best_q = q
      }
      q = q + 1
    }
    best_cost[nn] = best
    best_split[nn] = best_q - 1
    s = 0
    while s < ns {
      best_aver[s][nn] = aver[s][best_q]
      s = s + 1
    }
    nn = nn + 1
  }
  // backtrack fitted values (each sample needs its own mutable row)
  let yhat_z : Array[Array[Double]] = []
  let mut row_i = 0
  while row_i < ns {
    yhat_z.push(Array::make(n + 1, 0.0))
    row_i = row_i + 1
  }
  let mut nint = 0
  let mut cur = n
  while cur > 0 {
    let mut s = 0
    while s < ns {
      let val = best_aver[s][cur]
      let mut p = best_split[cur] + 1
      while p <= cur {
        yhat_z[s][p] = val
        p = p + 1
      }
      s = s + 1
    }
    cur = best_split[cur]
    nint = nint + 1
  }
  // recover segments in probe order
  let starts = Array::make(nint, 0)
  let lengths = Array::make(nint, 0)
  let seg_means : Array[Array[Double]] = Array::make(nint, [])
  cur = n
  let mut old = cur
  let mut cnt = nint
  while cur > 0 {
    starts[cnt - 1] = best_split[cur] + 1
    lengths[cnt - 1] = old - best_split[cur]
    let row : Array[Double] = []
    let mut s = 0
    while s < ns {
      // un-weight and un-normalize: mean / w * sd
      row.push(best_aver[s][cur] / w[s] * sd[s])
      s = s + 1
    }
    seg_means[cnt - 1] = row
    cur = best_split[cur]
    old = cur
    cnt = cnt - 1
  }
  let yhat1 : Array[Array[Double]] = []
  let mut s = 0
  while s < ns {
    let row : Array[Double] = []
    let mut p = 1
    while p <= n {
      row.push(yhat_z[s][p] / w[s] * sd[s])
      p = p + 1
    }
    yhat1.push(row)
    s = s + 1
  }
  let segments1 : Array[MultiPcfSegment] = []
  let mut si = 0
  while si < nint {
    segments1.push({
      start: starts[si],
      length: lengths[si],
      means: seg_means[si],
    })
    si = si + 1
  }
  { segments: segments1, n_intervals: nint, yhat: yhat1, sds: sd }
}
