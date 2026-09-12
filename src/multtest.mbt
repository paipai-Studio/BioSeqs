///|
/// # multtest: Resampling-based multiple testing procedures
///
/// MoonBit port of the Bioconductor `multtest` package (Pollard, Gilbert, Ge,
/// Taylor & Dudoit), which implements the multiple testing procedures of
/// Dudoit, Shaffer & Boldrick (2003). Two core families are provided:
///
/// 1. Marginal raw-p adjustments, `mt.rawp2adjp`, covering all nine
///    procedures shipped with the package:
///    Bonferroni, Holm, Hochberg, Sidak single-step (`SidakSS`),
///    Sidak step-down (`SidakSD`), Benjamini-Hochberg (`BH`),
///    Benjamini-Yekutieli (`BY`), adaptive BH (`ABH`) and the two-stage BH
///    (`TSBH`, one adjusted-p column per nominal level `alpha`).
///    `mt.reject` counts and locates rejected hypotheses.
/// 2. Westfall-Young permutation maxT (`mt_max_t`) and minP (`mt_min_p`)
///    step-down adjusted p-values, which take observed test statistics and a
///    matrix of null statistics computed under label permutations/bootstrap
///    (the caller generates the matrix; the two-sample statistics live in
///    other modules).
///
/// Numerical behaviour (column order, sort order, monotone enforcement,
/// adaptive `h0` estimators, NaN propagation) is verified bitwise against the
/// reference R implementation on fixed input grids.

///|
/// Names of the nine adjustment procedures understood by `mt_rawp2adjp`,
/// in the reference order used by `mt.rawp2adjp`.
pub let mt_procs : Array[String] = [
  "Bonferroni", "Holm", "Hochberg", "SidakSS", "SidakSD", "BH", "BY", "ABH", "TSBH",
]

///|
/// Result of `mt_rawp2adjp`.
///
/// - `rawp`: the input p-values sorted ascending, NaN entries last
///   (reference rows are always returned in sorted order).
/// - `adjp`: one column per requested procedure (keyed by procedure name;
///   TSBH columns are keyed `"TSBH_<alpha>"`), in the same sorted order.
/// - `index`: zero-based indices into the original input giving the sort
///   order.
/// - `h0_abh`: ABH estimate of the number of true null hypotheses (NaN when
///   the estimator is not defined).
/// - `h0_tsbh`: `(alpha, h0)` pairs, one per nominal TSBH level.
pub struct MTPAdjustResult {
  rawp : Array[Double]
  adjp : Map[String, Array[Double]]
  index : Array[Int]
  h0_abh : Double
  h0_tsbh : Array[(Double, Double)]
}

///|
/// Adjust marginal raw p-values. Port of `multtest::mt.rawp2adjp`.
///
/// `procs` selects procedures (default: all nine); unknown names are ignored,
/// as in R. `alpha` is the nominal level (or levels) used by TSBH, producing
/// one `"TSBH_<alpha>"` column per level. `na_rm` only changes the effective
/// number of hypotheses (`mgood`) used in the multipliers; NaN p-values are
/// always sorted last and propagated as NaN.
pub fn mt_rawp2adjp(
  rawp~ : Array[Double],
  procs? : Array[String] = [],
  alpha? : Array[Double] = [0.05],
  na_rm? : Bool = false,
) -> MTPAdjustResult raise SeqError {
  let procs = if procs.is_empty() { mt_procs } else { procs }
  let m = rawp.length()
  if m == 0 {
    raise SeqError("mt_rawp2adjp: rawp must have length > 0")
  }
  for p in rawp {
    if !p.is_nan() && (p < 0.0 || p > 1.0) {
      raise SeqError(
        "mt_rawp2adjp: raw p-values must lie in [0,1] (got " +
        p.to_string() +
        ")",
      )
    }
  }
  let mgood = if na_rm {
    rawp.fold(init=0, fn(acc : Int, p : Double) -> Int {
      if p.is_nan() {
        acc
      } else {
        acc + 1
      }
    })
  } else {
    m
  }
  // order(rawp): ascending finite values first (stable), NaN last.
  let index = Array::make(m, 0)
  let mut k = 0
  for i = 0; i < m; i = i + 1 {
    if !rawp[i].is_nan() {
      index[k] = i
      k = k + 1
    }
  }
  let nfinite = k
  for i = 0; i < m; i = i + 1 {
    if rawp[i].is_nan() {
      index[k] = i
      k = k + 1
    }
  }
  // insertion sort the finite prefix in ascending order (stable)
  for i = 1; i < nfinite; i = i + 1 {
    let key = index[i]
    let mut j = i - 1
    while j >= 0 && rawp[index[j]] > rawp[key] {
      index[j + 1] = index[j]
      j = j - 1
    }
    index[j + 1] = key
  }
  let spval = index.map(fn(i) { rawp[i] })
  let cols : Map[String, Array[Double]] = Map([], capacity=16)
  cols.set("rawp", spval.copy())

  let want = fn(name : String) -> Bool { procs.contains(name) }

  // First-pass BH statistic shared by ABH and TSBH.
  let bh_first = fn() -> Array[Double] {
    let tmp = spval.copy()
    let mut i = m - 2
    while i >= 0 {
      if spval[i].is_nan() {
        tmp[i] = nan_value()
      } else {
        let x = mgood.to_double() / (i + 1).to_double() * spval[i]
        tmp[i] = rmin_ignore_nan(tmp[i + 1], rmin_ignore_nan(x, 1.0))
      }
      i = i - 1
    }
    tmp
  }

  // TSBH must be prepared first in the reference so its columns slot in
  // before ABH; map order is irrelevant for MoonBit callers.
  let h0_tsbh : Array[(Double, Double)] = []
  if want("TSBH") {
    let tmp = bh_first()
    for a in alpha {
      let mut h0 = mgood.to_double()
      for x in tmp {
        if x < a / (1.0 + a) {
          h0 = h0 - 1.0
        }
      }
      let col = tmp.map(fn(x) { x * h0 / mgood.to_double() })
      cols.set("TSBH_" + a.to_string(), col)
      h0_tsbh.push((a, h0))
    }
  }

  if want("Bonferroni") {
    let col = spval.map(fn(p) {
      let v = mgood.to_double() * p
      if v > 1.0 {
        1.0
      } else {
        v
      }
    })
    cols.set("Bonferroni", col)
  }

  if want("Holm") {
    let tmp = spval.copy()
    tmp[0] = rmin1(mgood.to_double() * spval[0])
    for i = 1; i < m; i = i + 1 {
      let w = (mgood - i).to_double()
      tmp[i] = rmax_propagate_nan(tmp[i - 1], rmin1(w * spval[i]))
    }
    cols.set("Holm", tmp)
  }

  if want("Hochberg") {
    let tmp = spval.copy()
    let mut i = m - 2
    while i >= 0 {
      if spval[i].is_nan() {
        tmp[i] = nan_value()
      } else {
        let w = (mgood - i).to_double()
        tmp[i] = rmin_ignore_nan(tmp[i + 1], rmin_ignore_nan(w * spval[i], 1.0))
      }
      i = i - 1
    }
    cols.set("Hochberg", tmp)
  }

  if want("SidakSS") {
    let col = spval.map(fn(p) { 1.0 - sidak_pow(1.0 - p, mgood) })
    cols.set("SidakSS", col)
  }

  if want("SidakSD") {
    let tmp = spval.copy()
    tmp[0] = 1.0 - sidak_pow(1.0 - spval[0], mgood)
    for i = 1; i < m; i = i + 1 {
      let q = 1.0 - sidak_pow(1.0 - spval[i], mgood - i)
      tmp[i] = rmax_propagate_nan(tmp[i - 1], q)
    }
    cols.set("SidakSD", tmp)
  }

  if want("BH") {
    cols.set("BH", bh_first())
  }

  if want("BY") {
    let tmp = spval.copy()
    let a = harmonic_sum(mgood)
    tmp[m - 1] = rmin1(a * spval[m - 1])
    let mut i = m - 2
    while i >= 0 {
      if spval[i].is_nan() {
        tmp[i] = nan_value()
      } else {
        let x = mgood.to_double() * a / (i + 1).to_double() * spval[i]
        tmp[i] = rmin_ignore_nan(tmp[i + 1], rmin_ignore_nan(x, 1.0))
      }
      i = i - 1
    }
    cols.set("BY", tmp)
  }

  let mut h0_abh = nan_value()
  if want("ABH") {
    let tmp = bh_first()
    // h0.m[k] = (mgood+1-k)/(1-p[k]), k = 1..mgood (finite, ascending p).
    if mgood > 0 {
      let h0m : Array[Double] = Array::make(mgood, 0.0)
      for k = 0; k < mgood; k = k + 1 {
        h0m[k] = (mgood - k).to_double() / (1.0 - spval[k])
      }
      // first strict increase of h0.m; R's grab indexes the earlier point
      let mut grab = -1
      for k = 0; k < mgood - 1; k = k + 1 {
        if h0m[k + 1] > h0m[k] && grab < 0 {
          grab = k
        }
      }
      if grab >= 0 {
        h0_abh = (if h0m[grab] < mgood.to_double() {
          h0m[grab]
        } else {
          mgood.to_double()
        }).ceil()
      }
    }
    let col = tmp.map(fn(x) { x * h0_abh / mgood.to_double() })
    cols.set("ABH", col)
  }

  { rawp: spval, adjp: cols, index, h0_abh, h0_tsbh, }
}

///|
/// Count rejections and locate rejected hypotheses. Port of `mt.reject`.
///
/// For an adjusted-p vector and a single `alpha`, returns `(r, which)` where
/// `r` is the number of values `<= alpha` and `which[i]` is true when the
/// i-th hypothesis is rejected (NaN when the adjusted p is NaN).
pub fn mt_reject(
  adjp~ : Array[Double],
  alpha? : Double = 0.05,
) -> (Int, Array[Bool?]) raise SeqError {
  if adjp.length() == 0 {
    raise SeqError("mt_reject: adjp must have length > 0")
  }
  let mut r = 0
  let mut seen_nan = false
  let which = adjp.map(fn(p) -> Bool? {
    if p.is_nan() {
      seen_nan = true
      None
    } else if p <= alpha {
      r = r + 1
      Some(true)
    } else {
      Some(false)
    }
  })
  if seen_nan {
    (0, which) // R returns NA for the total when any comparison is NA
  } else {
    (r, which)
  }
}

///|
/// Westfall-Young single-step/maxT permutation adjusted p-values
/// (`mt.maxT` engine).
///
/// `observed` holds the m observed test statistics; `null_stats` is the m x B
/// matrix of test statistics under B label permutations (rows hypotheses).
/// Two-sided statistics are compared by absolute value. The adjusted
/// p-value of hypothesis i is the fraction of permutations whose largest
/// null statistic (over all hypotheses) is at least as extreme as the
/// observed statistic of i, enforced monotone by descending observed
/// statistic. No +1 correction is applied, matching the reference C engine
/// (`get_maxT`).
pub fn mt_max_t(
  observed~ : Array[Double],
  null_stats~ : Array[Array[Double]],
) -> Array[Double] raise SeqError {
  let (m, b) = check_null_matrix(observed~, null_stats~)
  let t = observed.map(fn(x) { x.abs() })
  // per-permutation maximum absolute null statistic
  let mmax : Array[Double] = Array::make(b, 0.0)
  for j = 0; j < b; j = j + 1 {
    let mut v = 0.0
    for i = 0; i < m; i = i + 1 {
      let z = null_stats[i][j].abs()
      if z > v {
        v = z
      }
    }
    mmax[j] = v
  }
  let rawp : Array[Double] = Array::make(m, 0.0)
  for i = 0; i < m; i = i + 1 {
    let mut c = 0
    for j = 0; j < b; j = j + 1 {
      if mmax[j] >= t[i] {
        c = c + 1
      }
    }
    rawp[i] = c.to_double() / b.to_double()
  }
  // order hypotheses by descending observed statistic, enforce monotonicity
  let order : Array[Int] = []
  for i = 0; i < m; i = i + 1 {
    order.push(i)
  }
  order.sort_by(fn(a, b2) { t[b2].compare(t[a]) })
  let adj : Array[Double] = Array::make(m, 0.0)
  let mut run = 0.0
  for i in order {
    if rawp[i] > run {
      run = rawp[i]
    }
    adj[i] = run
  }
  adj
}

///|
/// Westfall-Young single-step/minP permutation adjusted p-values
/// (`mt.minP` engine).
///
/// For each permutation, every null statistic is converted to a marginal
/// tail probability against that hypothesis' own empirical null distribution,
/// and the smallest probability over hypotheses forms the permutation
/// minP statistic. Adjusted p-values are the fraction of permutations whose
/// minP statistic is no larger than the observed marginal p-value, enforced
/// monotone by ascending marginal p-value.
pub fn mt_min_p(
  observed~ : Array[Double],
  null_stats~ : Array[Array[Double]],
) -> Array[Double] raise SeqError {
  let (m, b) = check_null_matrix(observed~, null_stats~)
  let t = observed.map(fn(x) { x.abs() })
  let v = null_stats.map(fn(row) { row.map(fn(x) { x.abs() }) })
  // marginal empirical p-values: p[i][j] = #{j2: V[i][j2] >= V[i][j]} / B
  let pp : Array[Array[Double]] = Array::make(m, [])
  for i = 0; i < m; i = i + 1 {
    pp[i] = Array::make(b, 0.0)
    for j = 0; j < b; j = j + 1 {
      let mut c = 0
      for j2 = 0; j2 < b; j2 = j2 + 1 {
        if v[i][j2] >= v[i][j] {
          c = c + 1
        }
      }
      pp[i][j] = c.to_double() / b.to_double()
    }
  }
  // per-permutation minimum
  let mmin : Array[Double] = Array::make(b, 1.0)
  for j = 0; j < b; j = j + 1 {
    let mut x = 1.0
    for i = 0; i < m; i = i + 1 {
      if pp[i][j] < x {
        x = pp[i][j]
      }
    }
    mmin[j] = x
  }
  let rawp : Array[Double] = Array::make(m, 0.0)
  for i = 0; i < m; i = i + 1 {
    let mut c = 0
    for j2 = 0; j2 < b; j2 = j2 + 1 {
      if v[i][j2] >= t[i] {
        c = c + 1
      }
    }
    let pobs = c.to_double() / b.to_double()
    let mut cc = 0
    for j = 0; j < b; j = j + 1 {
      if mmin[j] <= pobs {
        cc = cc + 1
      }
    }
    rawp[i] = cc.to_double() / b.to_double()
  }
  let order : Array[Int] = []
  for i = 0; i < m; i = i + 1 {
    order.push(i)
  }
  order.sort_by(fn(a, b2) { rawp[a].compare(rawp[b2]) })
  let adj : Array[Double] = Array::make(m, 0.0)
  let mut run = 0.0
  for i in order {
    if rawp[i] > run {
      run = rawp[i]
    }
    adj[i] = run
  }
  adj
}

///|
fn check_null_matrix(
  observed~ : Array[Double],
  null_stats~ : Array[Array[Double]],
) -> (Int, Int) raise SeqError {
  let m = observed.length()
  if m == 0 {
    raise SeqError("mt_max_t/mt_min_p: observed must be non-empty")
  }
  if null_stats.length() != m {
    raise SeqError(
      "mt_max_t/mt_min_p: null_stats must have one row per hypothesis",
    )
  }
  let b = null_stats[0].length()
  if b == 0 {
    raise SeqError("mt_max_t/mt_min_p: at least one permutation is required")
  }
  for row in null_stats {
    if row.length() != b {
      raise SeqError("mt_max_t/mt_min_p: ragged null_stats matrix")
    }
  }
  (m, b)
}

///|
/// R `min(x, 1)`: NaN propagates.
fn rmin1(x : Double) -> Double {
  if x.is_nan() {
    x
  } else if x < 1.0 {
    x
  } else {
    1.0
  }
}

///|
/// R `min(..., na.rm = TRUE)`: NaN arguments ignored.
fn rmin_ignore_nan(x : Double, y : Double) -> Double {
  if x.is_nan() {
    y
  } else if y.is_nan() {
    x
  } else if x < y {
    x
  } else {
    y
  }
}

///|
/// R `max(...)` without na.rm: NaN propagates.
fn rmax_propagate_nan(x : Double, y : Double) -> Double {
  if x.is_nan() || y.is_nan() {
    nan_value()
  } else if x > y {
    x
  } else {
    y
  }
}

///|
/// R `(1-p)^n`: note `NA^0 == 1` in R, so exponent zero is exactly 1.
fn sidak_pow(base : Double, n : Int) -> Double {
  if n == 0 {
    1.0
  } else {
    @math.pow(base, n.to_double())
  }
}

///|
/// Harmonic number H_n = sum(1/k, k=1..n); H_0 = 0 (as in R `sum(1/(1:0))`).
fn harmonic_sum(n : Int) -> Double {
  let mut s = 0.0
  for k = 1; k <= n; k = k + 1 {
    s = s + 1.0 / k.to_double()
  }
  s
}

///|
fn nan_value() -> Double {
  0.0 / 0.0
}
