///|
/// Bioconductor apcluster - Affinity Propagation Clustering.
///
/// Port of the Bioconductor `apcluster` package (Bodenhofer et al.,
/// Bioinformatics 2011), implementing the affinity propagation algorithm
/// of Frey & Dueck (Science 2007). Affinity propagation identifies
/// representative data points ("exemplars") by passing two kinds of real-
/// valued messages between pairs of points:
///   - responsibility r(i,k): how well-suited k is to serve as exemplar
///     for i, compared to other candidate exemplars;
///   - availability a(i,k): how appropriate it is for i to choose k as
///     exemplar, taking into account support from other points.
/// Messages are updated iteratively with a damping factor until the set of
/// exemplars remains stable for `convits` iterations.
///
/// Provided:
///   - Similarity matrices: ap_neg_dist_mat (negative Minkowski distance^r,
///     default negative squared Euclidean), ap_exp_sim_mat (Gaussian/RBF),
///     ap_lin_sim_mat (linear), ap_cor_sim_mat (Pearson / Spearman).
///   - apcluster: run affinity propagation with a scalar or per-point
///     preference (default: median similarity).
///   - apcluster_k: bisect the preference to obtain exactly K clusters,
///     with deterministic pruning merges as a fallback.
///   - Reproducible noise runs via apcluster_seed.

///|
/// Error raised by affinity propagation for invalid arguments, constant
/// input rows, or when an exact-K solution is unreachable.
pub suberror ApError {
  ApError(String)
}

///|
fn[T] ap_fail(msg : String) -> T raise ApError {
  raise ApError(msg)
}

// ---------------------------------------------------------------------------
// Random helpers (Park-Miller), used for the optional noise runs
// ---------------------------------------------------------------------------

///|
let ap_rng_state : Ref[Int] = Ref(777)

///|
/// Seed the internal random number generator used for noise runs. The
/// generator is deterministic (Park-Miller MINSTD).
pub fn apcluster_seed(seed : Int) -> Unit {
  ap_rng_state.val = if seed <= 0 { 777 } else { seed % 2147483647 }
}

///|
fn ap_uniform_signed() -> Double {
  let a = 48271
  let m = 2147483647
  let q = 44488
  let r = 3399
  let s = ap_rng_state.val
  let hi = s / q
  let lo = s % q
  let mut t = a * lo - r * hi
  if t <= 0 {
    t = t + m
  }
  ap_rng_state.val = t
  // uniform in (-1, 1)
  2.0 * (t.to_double() / m.to_double()) - 1.0
}

// ---------------------------------------------------------------------------
// Small numeric helpers
// ---------------------------------------------------------------------------

///|
fn ap_is_finite(v : Double) -> Bool {
  v == v && v.abs() <= 1.0e300
}

///|
/// R's default quantile (type 7) of a sorted array: h = (n-1)*q + 1 with
/// linear interpolation between order statistics.
fn ap_quantile_sorted(
  sorted : Array[Double],
  q : Double,
) -> Double raise ApError {
  let n = sorted.length()
  if n == 0 {
    ap_fail("cannot compute quantile of an empty vector")
  }
  if n == 1 {
    return sorted[0]
  }
  let h = (n - 1).to_double() * q + 1.0
  let kf = h.floor().to_int() // 1-based lower order statistic
  let k = kf - 1 // 0-based
  if k >= n - 1 {
    return sorted[n - 1]
  }
  let frac = h - kf.to_double()
  sorted[k] + frac * (sorted[k + 1] - sorted[k])
}

///|
/// Median of an array (R's median == type-7 50% quantile).
fn ap_median(v : Array[Double]) -> Double raise ApError {
  let s = v.copy()
  s.sort()
  ap_quantile_sorted(s, 0.5)
}

///|
/// Check that x is a non-empty rectangular matrix; return (n, p).
fn ap_check_matrix(x : Array[Array[Double]]) -> (Int, Int) raise ApError {
  let n = x.length()
  if n == 0 {
    ap_fail("empty input matrix")
  }
  let p = x[0].length()
  if p == 0 {
    ap_fail("input matrix has zero columns")
  }
  for row in x {
    if row.length() != p {
      ap_fail("ragged input matrix: all rows must have the same length")
    }
  }
  (n, p)
}

///|
/// Allocate an n x n matrix filled with 0.0 (independent rows).
fn ap_zero_matrix(n : Int) -> Array[Array[Double]] {
  let m : Array[Array[Double]] = Array(capacity=n)
  for _i in 0..<n {
    m.push(Array::make(n, 0.0))
  }
  m
}

///|
/// Deep copy of a square matrix.
fn ap_copy_matrix(m : Array[Array[Double]]) -> Array[Array[Double]] {
  let n = m.length()
  let out : Array[Array[Double]] = Array(capacity=n)
  for i = 0; i < n; i = i + 1 {
    out.push(m[i].copy())
  }
  out
}

// ---------------------------------------------------------------------------
// Similarity matrix constructors
// ---------------------------------------------------------------------------

///|
/// Pairwise Minkowski-type dissimilarity d(i,j)^r = sum_k |x_ik - x_jk|^r.
/// With r = 2 (default) this is the squared Euclidean distance, as used by
/// the default argument negDistMat(r=2) of the R package.
pub fn ap_minkowski_power(
  x : Array[Array[Double]],
  r? : Double = 2.0,
) -> Array[Array[Double]] raise ApError {
  let (n, p) = ap_check_matrix(x)
  let d = ap_zero_matrix(n)
  for i = 0; i < n; i = i + 1 {
    for j = i + 1; j < n; j = j + 1 {
      let mut sum = 0.0
      for k in 0..<p {
        let diff = x[i][k] - x[j][k]
        let ad = if diff < 0.0 { -diff } else { diff }
        sum = sum + @math.pow(ad, r)
      }
      d[i][j] = sum
      d[j][i] = sum
    }
  }
  d
}

///|
/// Negative distance similarity matrix (R: negDistMat). S(i,j) = -d(i,j)^r,
/// diagonal = 0 (self-distance zero). The diagonal is later overwritten by
/// the preference when running apcluster.
pub fn ap_neg_dist_mat(
  x : Array[Array[Double]],
  r? : Double = 2.0,
) -> Array[Array[Double]] raise ApError {
  let d = ap_minkowski_power(x, r~)
  let n = d.length()
  let s = ap_zero_matrix(n)
  for i in 0..<n {
    for j in 0..<n {
      if i != j {
        s[i][j] = -d[i][j]
      }
    }
  }
  s
}

///|
/// Exponential (Gaussian/RBF) similarity matrix (R: expSimMat).
/// S(i,j) = exp(-d(i,j)^r / w). If w is not supplied it is set to the
/// median of the off-diagonal d(i,j)^r values. Diagonal = 1.
pub fn ap_exp_sim_mat(
  x : Array[Array[Double]],
  r? : Double = 2.0,
  w? : Double? = None,
) -> Array[Array[Double]] raise ApError {
  let d = ap_minkowski_power(x, r~)
  let n = d.length()
  let width = match w {
    Some(v) => if v <= 0.0 { ap_fail("w must be positive") } else { v }
    None => {
      let off : Array[Double] = []
      for i = 0; i < n; i = i + 1 {
        for j = i + 1; j < n; j = j + 1 {
          off.push(d[i][j])
        }
      }
      let med = ap_median(off)
      if med <= 0.0 {
        ap_fail("cannot default w: all points are identical")
      }
      med
    }
  }
  let s = ap_zero_matrix(n)
  for i = 0; i < n; i = i + 1 {
    s[i][i] = 1.0
    for j = i + 1; j < n; j = j + 1 {
      let v = @math.exp(-d[i][j] / width)
      s[i][j] = v
      s[j][i] = v
    }
  }
  s
}

///|
/// Linear similarity matrix (R: linSimMat).
/// S(i,j) = 1 - d(i,j)^r / w. If w is not supplied it is set to the
/// maximum off-diagonal d(i,j)^r, so similarities lie in [0, 1].
/// Diagonal = 1.
pub fn ap_lin_sim_mat(
  x : Array[Array[Double]],
  r? : Double = 2.0,
  w? : Double? = None,
) -> Array[Array[Double]] raise ApError {
  let d = ap_minkowski_power(x, r~)
  let n = d.length()
  let width = match w {
    Some(v) => if v <= 0.0 { ap_fail("w must be positive") } else { v }
    None => {
      let mut mx = 0.0
      for i = 0; i < n; i = i + 1 {
        for j = i + 1; j < n; j = j + 1 {
          if d[i][j] > mx {
            mx = d[i][j]
          }
        }
      }
      if mx <= 0.0 {
        ap_fail("cannot default w: all points are identical")
      }
      mx
    }
  }
  let s = ap_zero_matrix(n)
  for i = 0; i < n; i = i + 1 {
    s[i][i] = 1.0
    for j = i + 1; j < n; j = j + 1 {
      let v = 1.0 - d[i][j] / width
      s[i][j] = v
      s[j][i] = v
    }
  }
  s
}

///|
/// Average ranks (R rank(..., ties.method="average")) of a vector.
fn ap_average_ranks(v : Array[Double]) -> Array[Double] {
  let p = v.length()
  // index-sort by value (stable; ties get adjacent positions)
  let idx : Array[Int] = Array(capacity=p)
  for i = 0; i < p; i = i + 1 {
    idx.push(i)
  }
  idx.sort_by(fn(a : Int, b : Int) -> Int {
    if v[a] < v[b] {
      -1
    } else if v[a] > v[b] {
      1
    } else {
      0
    }
  })
  let ranks = Array::make(p, 0.0)
  let mut pos = 0
  while pos < p {
    let mut end = pos + 1
    while end < p && v[idx[end]] == v[idx[pos]] {
      end = end + 1
    }
    // positions pos..end-1 are tied; R ranks are 1-based, average the ranks
    let avg_rank = (pos + 1 + end).to_double() / 2.0
    let mut k = pos
    while k < end {
      ranks[idx[k]] = avg_rank
      k = k + 1
    }
    pos = end
  }
  ranks
}

///|
/// Pearson correlation between two vectors (means precomputed). Zero
/// variance raises ApError, as cor() would produce NA.
fn ap_pearson(
  a : Array[Double],
  b : Array[Double],
  mean_a : Double,
  mean_b : Double,
) -> Double raise ApError {
  let p = a.length()
  let mut saa = 0.0
  let mut sbb = 0.0
  let mut sab = 0.0
  for k = 0; k < p; k = k + 1 {
    let da = a[k] - mean_a
    let db = b[k] - mean_b
    saa = saa + da * da
    sbb = sbb + db * db
    sab = sab + da * db
  }
  if saa == 0.0 || sbb == 0.0 {
    ap_fail("cannot compute correlation of a constant vector")
  }
  sab / (saa * sbb).sqrt()
}

///|
/// Correlation-based similarity matrix (R: corSimMat). Similarity between
/// two data points (rows) is their correlation across features (columns).
/// method: "pearson" (default) or "spearman" (Pearson correlation of the
/// average ranks of features within each row). Diagonal = 1.
pub fn ap_cor_sim_mat(
  x : Array[Array[Double]],
  method? : String = "pearson",
) -> Array[Array[Double]] raise ApError {
  let (n, p) = ap_check_matrix(x)
  // precompute transformed rows and their means
  let rows : Array[Array[Double]] = Array(capacity=n)
  let means : Array[Double] = Array(capacity=n)
  for i = 0; i < n; i = i + 1 {
    let row = match method {
      "pearson" => x[i].copy()
      "spearman" => ap_average_ranks(x[i])
      _ => ap_fail("unknown correlation method: \{method}")
    }
    let mut mean = 0.0
    for v in row {
      mean = mean + v
    }
    mean = mean / p.to_double()
    rows.push(row)
    means.push(mean)
  }
  let s = ap_zero_matrix(n)
  for i = 0; i < n; i = i + 1 {
    s[i][i] = 1.0
    for j = i + 1; j < n; j = j + 1 {
      let c = ap_pearson(rows[i], rows[j], means[i], means[j])
      s[i][j] = c
      s[j][i] = c
    }
  }
  s
}

// ---------------------------------------------------------------------------
// Result object
// ---------------------------------------------------------------------------

///|
/// Result of an affinity propagation run.
///
/// - exemplars: indices (0-based) of the chosen exemplar points, ascending.
/// - clusters: member indices per exemplar (same order as `exemplars`),
///   each sorted ascending.
/// - labels: exemplar index assigned to every input point (labels[i] is the
///   exemplar of point i).
/// - preferences: per-point preferences actually used.
/// - iters / converged: iterations performed and whether the exemplar set
///   was stable for `convits` iterations.
/// - dpsim: sum of similarities of non-exemplar points to their exemplars.
/// - expref: sum of preferences of the chosen exemplars.
/// - netsim: total net similarity (objective), dpsim + expref.
/// - sim: the similarity matrix used (diagonal overwritten with the
///   preferences).
pub struct ApResult {
  exemplars : Array[Int]
  clusters : Array[Array[Int]]
  labels : Array[Int]
  preferences : Array[Double]
  iters : Int
  converged : Bool
  dpsim : Double
  expref : Double
  netsim : Double
  sim : Array[Array[Double]]
}

///|
/// Accessor: exemplar indices (0-based), ascending.
pub fn ApResult::exemplars(self : ApResult) -> Array[Int] {
  self.exemplars
}

///|
/// Accessor: cluster member lists (same order as exemplars()), each sorted
/// ascending.
pub fn ApResult::clusters(self : ApResult) -> Array[Array[Int]] {
  self.clusters
}

///|
/// Accessor: label (exemplar index) of every input point.
pub fn ApResult::labels(self : ApResult) -> Array[Int] {
  self.labels
}

///|
/// Accessor: number of clusters found.
pub fn ApResult::num_clusters(self : ApResult) -> Int {
  self.exemplars.length()
}

///|
/// Accessor: per-point preferences used.
pub fn ApResult::preferences(self : ApResult) -> Array[Double] {
  self.preferences
}

///|
/// Accessor: iterations performed.
pub fn ApResult::iters(self : ApResult) -> Int {
  self.iters
}

///|
/// Accessor: whether the run converged.
pub fn ApResult::converged(self : ApResult) -> Bool {
  self.converged
}

///|
/// Accessor: data-point part of the net similarity.
pub fn ApResult::dpsim(self : ApResult) -> Double {
  self.dpsim
}

///|
/// Accessor: sum of exemplar preferences.
pub fn ApResult::expref(self : ApResult) -> Double {
  self.expref
}

///|
/// Accessor: total net similarity (objective value).
pub fn ApResult::netsim(self : ApResult) -> Double {
  self.netsim
}

///|
/// Accessor: the similarity matrix used (diagonal holds the preferences).
pub fn ApResult::sim(self : ApResult) -> Array[Array[Double]] {
  self.sim
}

// ---------------------------------------------------------------------------
// Core affinity propagation (one run)
// ---------------------------------------------------------------------------

///|
/// Resolve the per-point preferences: explicit vector, explicit scalar, or
/// the type-7 quantile q (default median) of the finite similarities.
fn ap_resolve_preferences(
  sim : Array[Array[Double]],
  p~ : Double?,
  q~ : Double,
  preferences~ : Array[Double]?,
) -> Array[Double] raise ApError {
  let n = sim.length()
  match preferences {
    Some(v) => {
      if v.length() != n {
        ap_fail("preferences vector must have one entry per point")
      }
      v.copy()
    }
    None => {
      let pref = match p {
        Some(v) => v
        None => {
          let finite : Array[Double] = []
          for i = 0; i < n; i = i + 1 {
            for j = 0; j < n; j = j + 1 {
              if ap_is_finite(sim[i][j]) {
                finite.push(sim[i][j])
              }
            }
          }
          if finite.length() == 0 {
            ap_fail("similarity matrix contains no finite values")
          }
          finite.sort()
          ap_quantile_sorted(finite, q)
        }
      }
      Array::make(n, pref)
    }
  }
}

///|
/// One affinity propagation run on a similarity matrix whose diagonal is
/// set to the per-point preferences. Returns the working matrix, the
/// availability matrix, the iteration count, the convergence flag and the
/// per-point labels.
fn ap_run(
  sim0 : Array[Array[Double]],
  preferences : Array[Double],
  maxits : Int,
  convits : Int,
  lam : Double,
) -> (Array[Array[Double]], Array[Array[Double]], Int, Bool, Array[Int]) {
  let n = sim0.length()
  // working similarity matrix with diagonal = preferences
  let s = ap_copy_matrix(sim0)
  for i = 0; i < n; i = i + 1 {
    s[i][i] = preferences[i]
  }
  let a = ap_zero_matrix(n)
  let r = ap_zero_matrix(n)
  let mut prev_exemplars : Array[Int] = []
  let mut stable = 0
  let mut iters = 0
  let mut converged = false
  for it = 1; it <= maxits; it = it + 1 {
    iters = it
    // --- responsibility update (with top-two trick per row) ---
    let rnew = ap_zero_matrix(n)
    for i = 0; i < n; i = i + 1 {
      // largest and second largest of a(i,k)+s(i,k) and their columns
      let mut v1 = -1.7976931348623157e308
      let mut v2 = -1.7976931348623157e308
      let mut i1 = -1
      for k = 0; k < n; k = k + 1 {
        let val = a[i][k] + s[i][k]
        if val > v1 {
          v2 = v1
          v1 = val
          i1 = k
        } else if val > v2 {
          v2 = val
        }
      }
      for k = 0; k < n; k = k + 1 {
        let best_other = if k == i1 { v2 } else { v1 }
        rnew[i][k] = s[i][k] - best_other
      }
    }
    // damping of R
    for i = 0; i < n; i = i + 1 {
      for k = 0; k < n; k = k + 1 {
        r[i][k] = lam * r[i][k] + (1.0 - lam) * rnew[i][k]
      }
    }
    // --- availability update (per column) ---
    let anew = ap_zero_matrix(n)
    for k = 0; k < n; k = k + 1 {
      let rkk = r[k][k]
      let mut sumpos = 0.0
      let pos = Array::make(n, 0.0)
      for j = 0; j < n; j = j + 1 {
        if j != k && r[j][k] > 0.0 {
          pos[j] = r[j][k]
          sumpos = sumpos + r[j][k]
        }
      }
      anew[k][k] = sumpos
      for i = 0; i < n; i = i + 1 {
        if i != k {
          let v = rkk + sumpos - pos[i]
          anew[i][k] = if v < 0.0 { v } else { 0.0 }
        }
      }
    }
    // damping of A
    for i = 0; i < n; i = i + 1 {
      for k = 0; k < n; k = k + 1 {
        a[i][k] = lam * a[i][k] + (1.0 - lam) * anew[i][k]
      }
    }
    // --- current exemplar set ---
    let exemplars : Array[Int] = []
    for i = 0; i < n; i = i + 1 {
      if a[i][i] + r[i][i] > 0.0 {
        exemplars.push(i)
      }
    }
    if exemplars.length() == prev_exemplars.length() {
      let mut same = true
      for idx in 0..<exemplars.length() {
        if exemplars[idx] != prev_exemplars[idx] {
          same = false
        }
      }
      if same {
        stable = stable + 1
      } else {
        stable = 0
      }
    } else {
      stable = 0
    }
    prev_exemplars = exemplars
    if stable >= convits && exemplars.length() > 0 {
      converged = true
      break
    }
  }
  // --- assignment: each point chooses the exemplar maximizing a(i,e)+r(i,e)
  let mut exemplars = prev_exemplars
  if exemplars.length() == 0 {
    // degenerate fallback: the single best point becomes the exemplar
    let mut best = 0
    let mut bestv = -1.7976931348623157e308
    for i = 0; i < n; i = i + 1 {
      let v = a[i][i] + r[i][i]
      if v > bestv {
        bestv = v
        best = i
      }
    }
    exemplars = [best]
  }
  let labels = Array::make(n, -1)
  for i = 0; i < n; i = i + 1 {
    let mut best_e = exemplars[0]
    let mut best_v = a[i][best_e] + r[i][best_e]
    for ei in 1..<exemplars.length() {
      let e = exemplars[ei]
      let v = a[i][e] + r[i][e]
      if v > best_v {
        best_v = v
        best_e = e
      }
    }
    labels[i] = best_e
  }
  (s, a, iters, converged, labels)
}

///|
/// Build the ApResult from a completed run (exemplars, clusters, objective).
fn ap_assemble(
  s : Array[Array[Double]],
  labels : Array[Int],
  preferences : Array[Double],
  iters : Int,
  converged : Bool,
) -> ApResult {
  let n = s.length()
  // distinct exemplars, ascending
  let is_ex : Array[Bool] = Array::make(n, false)
  for i = 0; i < n; i = i + 1 {
    is_ex[labels[i]] = true
  }
  let exemplars : Array[Int] = []
  for i = 0; i < n; i = i + 1 {
    if is_ex[i] {
      exemplars.push(i)
    }
  }
  let clusters : Array[Array[Int]] = Array(capacity=exemplars.length())
  for e in exemplars {
    let members : Array[Int] = []
    for i = 0; i < n; i = i + 1 {
      if labels[i] == e {
        members.push(i)
      }
    }
    clusters.push(members)
  }
  // objective: sum over non-exemplar points of s(i, e_i) + sum over
  // exemplars of their preference
  let mut dpsim = 0.0
  let mut expref = 0.0
  for i = 0; i < n; i = i + 1 {
    if is_ex[i] {
      expref = expref + preferences[i]
    } else {
      dpsim = dpsim + s[i][labels[i]]
    }
  }
  {
    exemplars,
    clusters,
    labels,
    preferences,
    iters,
    converged,
    dpsim,
    expref,
    netsim: dpsim + expref,
    sim: s,
  }
}

///|
/// Run affinity propagation.
///
/// - sim: n x n similarity matrix (symmetric); its diagonal is ignored and
///   overwritten with the preferences.
/// - p: scalar preference (all points share it). If neither `p` nor
///   `preferences` is given, the preference defaults to the median
///   (quantile q, default 0.5) of the finite similarities.
/// - preferences: per-point preferences (overrides `p`).
/// - maxits (1000) / convits (100) / lam (0.9 damping): R package defaults.
/// - nonoise: number of extra runs with small random noise added to the
///   similarities (R: nonoise); the best net-similarity solution wins.
///   Use apcluster_seed() for reproducibility.
/// - noise_r: noise magnitude relative to the similarity range (R: r,
///   default 0.001).
pub fn apcluster(
  sim~ : Array[Array[Double]],
  p? : Double? = None,
  q? : Double = 0.5,
  preferences? : Array[Double]? = None,
  maxits? : Int = 1000,
  convits? : Int = 100,
  lam? : Double = 0.9,
  nonoise? : Int = 0,
  noise_r? : Double = 0.001,
) -> ApResult raise ApError {
  let n = sim.length()
  if n == 0 {
    ap_fail("empty similarity matrix")
  }
  for row in sim {
    if row.length() != n {
      ap_fail("similarity matrix must be square")
    }
  }
  if lam < 0.5 || lam >= 1.0 {
    ap_fail("lam must lie in [0.5, 1)")
  }
  if maxits <= 0 || convits <= 0 {
    ap_fail("maxits and convits must be positive")
  }
  let prefs = ap_resolve_preferences(sim, p~, q~, preferences~)
  // deterministic run first (matches R when nonoise = 0)
  let (s0, _a0, it0, c0, l0) = ap_run(sim, prefs, maxits, convits, lam)
  let mut best = ap_assemble(s0, l0, prefs, it0, c0)
  // optional noise runs
  if nonoise > 0 {
    // finite range over all entries
    let mut lo = 1.7976931348623157e308
    let mut hi = -1.7976931348623157e308
    for i = 0; i < n; i = i + 1 {
      for j = 0; j < n; j = j + 1 {
        if ap_is_finite(sim[i][j]) {
          if sim[i][j] < lo {
            lo = sim[i][j]
          }
          if sim[i][j] > hi {
            hi = sim[i][j]
          }
        }
      }
    }
    let range = hi - lo
    for _run in 0..<nonoise {
      let noisy = ap_copy_matrix(sim)
      for i = 0; i < n; i = i + 1 {
        for j = 0; j < n; j = j + 1 {
          if ap_is_finite(noisy[i][j]) {
            noisy[i][j] = noisy[i][j] + ap_uniform_signed() * noise_r * range
          }
        }
      }
      let (s, _a, iters, converged, labels) = ap_run(
        noisy, prefs, maxits, convits, lam,
      )
      let cand = ap_assemble(s, labels, prefs, iters, converged)
      if cand.netsim > best.netsim {
        best = cand
      }
    }
  }
  best
}

///|
/// Convenience wrapper: build the similarity matrix with `sim_fn` (e.g.
/// ap_neg_dist_mat, ap_exp_sim_mat, ap_cor_sim_mat) and run apcluster.
/// All other arguments are forwarded.
pub fn apcluster_on_data(
  x~ : Array[Array[Double]],
  sim_fn~ : (Array[Array[Double]]) -> Array[Array[Double]] raise ApError,
  p? : Double? = None,
  q? : Double = 0.5,
  preferences? : Array[Double]? = None,
  maxits? : Int = 1000,
  convits? : Int = 100,
  lam? : Double = 0.9,
  nonoise? : Int = 0,
  noise_r? : Double = 0.001,
) -> ApResult raise ApError {
  let sim = sim_fn(x)
  apcluster(
    sim~,
    p~,
    q~,
    preferences~,
    maxits~,
    convits~,
    lam~,
    nonoise~,
    noise_r~,
  )
}

// ---------------------------------------------------------------------------
// apclusterK: preference bisection for a requested number of clusters
// ---------------------------------------------------------------------------

///|
/// Reassign every point formerly labelled `from` to exemplar `to`
/// (used as a deterministic pruning fallback when bisection cannot produce
/// exactly K clusters).
fn ap_merge_labels(labels : Array[Int], from : Int, to : Int) -> Array[Int] {
  let n = labels.length()
  let out = labels.copy()
  for i in 0..<n {
    if out[i] == from {
      out[i] = to
    }
  }
  out
}

///|
/// Run affinity propagation while targeting exactly K clusters (R:
/// apclusterK). The preference is bisected between the minimum and maximum
/// finite similarity; a run yielding K clusters is returned. If too many
/// clusters remain, clusters are pruned by repeatedly merging the pair of
/// exemplars with the highest mutual similarity. If too few clusters remain
/// and exact is true, ApError is raised (R: "could not find K clusters").
pub fn apcluster_k(
  sim~ : Array[Array[Double]],
  k~ : Int,
  p? : Double? = None,
  bimaxit? : Int = 100,
  exact? : Bool = true,
  maxits? : Int = 1000,
  convits? : Int = 100,
  lam? : Double = 0.9,
) -> ApResult raise ApError {
  let n = sim.length()
  if k <= 0 {
    ap_fail("k must be positive")
  }
  if k > n {
    ap_fail("k cannot exceed the number of points (\{n})")
  }
  if n == 1 {
    return apcluster(sim~, p~, maxits~, convits~, lam~)
  }
  let finite : Array[Double] = []
  for i = 0; i < n; i = i + 1 {
    for j = 0; j < n; j = j + 1 {
      if ap_is_finite(sim[i][j]) {
        finite.push(sim[i][j])
      }
    }
  }
  finite.sort()
  let mut p_lo = finite[0]
  let mut p_hi = finite[finite.length() - 1]
  let mut p_cur = match p {
    Some(v) => v
    None => ap_quantile_sorted(finite, 0.5)
  }
  let mut result : ApResult? = None
  for _it in 0..<bimaxit {
    let run = apcluster(sim~, p=Some(p_cur), maxits~, convits~, lam~)
    let ne = run.num_clusters()
    if ne == k {
      return run
    }
    // remember the run with the closest cluster count
    result = match result {
      Some(prev) =>
        if (ne - k).abs() < (prev.num_clusters() - k).abs() {
          Some(run)
        } else {
          Some(prev)
        }
      None => Some(run)
    }
    if ne > k {
      // too many exemplars -> lower the preference
      p_hi = p_cur
    } else {
      // too few -> raise it
      p_lo = p_cur
    }
    let next = (p_lo + p_hi) / 2.0
    if next == p_cur {
      break
    }
    p_cur = next
  }
  // bisection exhausted without an exact-K solution
  let mut run = match result {
    Some(r) => r
    None => ap_fail("apcluster_k failed to produce any result")
  }
  if run.num_clusters() > k {
    // prune: merge closest exemplar pairs until K clusters remain
    let mut labels = run.labels().copy()
    let mut exemplars = run.exemplars().copy()
    let sim_mat = run.sim()
    while exemplars.length() > k {
      // pick the exemplar pair with the highest mutual similarity
      let mut best_a = exemplars[0]
      let mut best_b = exemplars[1]
      let mut best_s = -1.7976931348623157e308
      for ii in 0..<exemplars.length() {
        for jj = ii + 1; jj < exemplars.length(); jj = jj + 1 {
          let ea = exemplars[ii]
          let eb = exemplars[jj]
          if sim_mat[ea][eb] > best_s {
            best_s = sim_mat[ea][eb]
            best_a = ea
            best_b = eb
          }
        }
      }
      labels = ap_merge_labels(labels, best_b, best_a)
      let new_ex : Array[Int] = []
      for e in exemplars {
        if e != best_b {
          new_ex.push(e)
        }
      }
      exemplars = new_ex
    }
    // rebuild result with the pruned assignment
    let prefs = run.preferences().copy()
    let s = run.sim()
    run = ap_assemble(s, labels, prefs, run.iters(), run.converged())
  }
  if run.num_clusters() < k && exact {
    let got = run.num_clusters()
    ap_fail(
      "could not find a solution with \{k} clusters (got \{got}); use exact=false to accept it",
    )
  }
  run
}

///|
/// Convenience wrapper for apcluster_k: build similarities with `sim_fn`
/// and target K clusters.
pub fn apcluster_k_on_data(
  x~ : Array[Array[Double]],
  sim_fn~ : (Array[Array[Double]]) -> Array[Array[Double]] raise ApError,
  k~ : Int,
  bimaxit? : Int = 100,
  exact? : Bool = true,
  maxits? : Int = 1000,
  convits? : Int = 100,
  lam? : Double = 0.9,
) -> ApResult raise ApError {
  let sim = sim_fn(x)
  apcluster_k(sim~, k~, bimaxit~, exact~, maxits~, convits~, lam~)
}
