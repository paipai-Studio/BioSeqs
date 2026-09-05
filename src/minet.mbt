// Mutual-information gene-network inference - port of Bioconductor `minet`
// (Meyer, Lafitte & Bontempi, BMC Bioinformatics 2008).
//
// Pipeline: discretize the expression matrix -> build the mutual
// information matrix (MIM) with one of four entropy estimators
// (empirical plug-in, Miller-Madow bias correction, Hausser-Strimmer
// shrinkage, Schuermann-Grassberger Dirichlet) -> infer the network with
// ARACNE (data-processing-inequality pruning), CLR (background
// z-score correction) or MRNET (mRMR forward selection). A network can
// then be scored against a reference adjacency matrix with
// precision/recall curves, AUROC, AUPR and the maximal F-score.
//
// Matrix convention (consistent with the WGCNA port): rows are
// variables (genes), columns are samples. All information values are
// returned in bits (log base 2).

///|
/// Result of scoring an inferred network against a reference network,
/// analog of minet's `validate()` output.
pub struct MinetValidation {
  thresholds : Array[Double]
  precision : Array[Double]
  recall : Array[Double]
  fpr : Array[Double]
  fscore : Array[Double]
  auc_roc : Double
  auc_pr : Double
  max_fscore : Double
}

///|
fn mt_log2(x : Double) -> Double {
  @math.log2(x)
}

///|
/// Discretize one continuous variable (samples along the row) into
/// `nbin` factor levels. `"equalfreq"` puts (about) equally many samples
/// in each bin (rank/quantile based); `"equalwidth"` uses bins of equal
/// value range. Returns bin indices in 1..nbin. Ties straddling an
/// equal-frequency boundary share a single bin, as with R's cut().
fn minet_discretize_row(row : Array[Double], nbin : Int, disc_method : String) -> Array[Int] {
  let n = row.length()
  let values : Array[Double] = []
  let mut i = 0
  while i < n {
    values.push(row[i])
    i = i + 1
  }
  // insertion-sort indices by value (n is small, gene profile length)
  let order : Array[Int] = []
  i = 0
  while i < n {
    order.push(i)
    i = i + 1
  }
  i = 1
  while i < n {
    let key = order[i]
    let kv = values[key]
    let mut j = i - 1
    while j >= 0 && values[order[j]] > kv {
      order[j + 1] = order[j]
      j = j - 1
    }
    order[j + 1] = key
    i = i + 1
  }
  let bins = Array::make(n, 1)
  if disc_method == "equalwidth" {
    let mut lo = row[0]
    let mut hi = row[0]
    i = 1
    while i < n {
      if row[i] < lo {
        lo = row[i]
      }
      if row[i] > hi {
        hi = row[i]
      }
      i = i + 1
    }
    let width = (hi - lo) / nbin.to_double()
    i = 0
    while i < n {
      let mut b = if width > 0.0 {
        ((row[i] - lo) / width).floor().to_int() + 1
      } else {
        1
      }
      if b < 1 {
        b = 1
      }
      if b > nbin {
        b = nbin
      }
      bins[i] = b
      i = i + 1
    }
  } else {
    // equalfreq: rank based bins of ~equal size; merge tie blocks that
    // straddle a bin boundary into the lower bin.
    i = 0
    while i < n {
      let rank_pos = i // 0-based sorted position
      let mut b = (rank_pos * nbin) / n + 1
      if b > nbin {
        b = nbin
      }
      bins[order[i]] = b
      i = i + 1
    }
    // enforce ties get identical bins (walk sorted sequence, propagate
    // the lower bin across equal values)
    i = 1
    while i < n {
      if values[order[i]] == values[order[i - 1]] {
        let idx = order[i]
        let prev = order[i - 1]
        let b = if bins[idx] < bins[prev] { bins[idx] } else { bins[prev] }
        bins[idx] = b
        bins[prev] = b
      }
      i = i + 1
    }
  }
  bins
}

///|
/// Discretize a gene x sample expression matrix. `nbin` defaults to
/// floor(sqrt(n_samples)), the minet default; `disc_method` is
/// "equalfreq" (default) or "equalwidth".
pub fn minet_discretize(
  data : Array[Array[Double]],
  disc_method? : String = "equalfreq",
  nbin? : Int = 0,
) -> Array[Array[Int]] raise SeqError {
  if data.length() == 0 || data[0].length() == 0 {
    raise SeqError("minet_discretize: empty expression matrix")
  }
  let nsamples = data[0].length()
  let mut i = 0
  while i < data.length() {
    if data[i].length() != nsamples {
      raise SeqError("minet_discretize: ragged expression matrix (row lengths differ)")
    }
    i = i + 1
  }
  if disc_method != "equalfreq" && disc_method != "equalwidth" {
    raise SeqError(
      "minet_discretize: unknown method '\{disc_method}' (expected equalfreq|equalwidth)",
    )
  }
  let nb = if nbin == 0 {
    (nsamples.to_double().sqrt()).floor().to_int()
  } else {
    nbin
  }
  if nb < 2 {
    raise SeqError("minet_discretize: nbin must be >= 2")
  }
  if nb > nsamples {
    raise SeqError("minet_discretize: nbin (\{nb}) exceeds number of samples (\{nsamples})")
  }
  let result : Array[Array[Int]] = []
  i = 0
  while i < data.length() {
    result.push(minet_discretize_row(data[i], nb, disc_method))
    i = i + 1
  }
  result
}

///|
/// Shannon entropy (bits) of a discretized variable, with the four
/// estimators offered by minet:
///  - "empirical"  : plug-in maximum likelihood
///  - "mi.mm"      : Miller-Madow bias correction (default)
///  - "mi.shrink"  : Hausser-Strimmer shrinkage toward the uniform
///  - "mi.sg"      : Schuermann-Grassberger Dirichlet prior (beta=1/K)
///
/// `k_levels` is the total number of factor levels (empty levels count);
/// when 0 it is taken as the number of observed levels.
pub fn minet_entropy(
  x : Array[Int],
  estimator? : String = "mi.mm",
  k_levels? : Int = 0,
) -> Double raise SeqError {
  let n = x.length()
  if n == 0 {
    raise SeqError("minet_entropy: empty variable")
  }
  let mut max_level = 0
  let counts : Map[Int, Int] = Map([])
  let mut i = 0
  while i < n {
    let v = x[i]
    if v < 1 {
      raise SeqError("minet_entropy: bin indices must be >= 1")
    }
    if v > max_level {
      max_level = v
    }
    counts.set(v, counts.get_or_default(v, 0) + 1)
    i = i + 1
  }
  let k_tot = if k_levels > 0 { k_levels } else { max_level }
  let nd = n.to_double()
  match estimator {
    "empirical" => {
      let mut h = 0.0
      for kv in counts {
        let c = kv.1.to_double()
        let p = c / nd
        if p > 0.0 {
          h = h - p * mt_log2(p)
        }
      }
      h
    }
    "mi.mm" => {
      let mut h = 0.0
      let nonempty = counts.length()
      for kv in counts {
        let c = kv.1.to_double()
        let p = c / nd
        if p > 0.0 {
          h = h - p * mt_log2(p)
        }
      }
      // Miller-Madow: H_mm = H_emp + (K_hat - 1) / (2 n)
      h + (nonempty.to_double() - 1.0) / (2.0 * nd)
    }
    "mi.shrink" => {
      let u = 1.0 / k_tot.to_double()
      let mut sum_p2 = 0.0
      let mut sum_dev2 = 0.0
      let mut level = 1
      while level <= k_tot {
        let p = counts.get_or_default(level, 0).to_double() / nd
        sum_p2 = sum_p2 + p * p
        let d = p - u
        sum_dev2 = sum_dev2 + d * d
        level = level + 1
      }
      let mut lambda = if sum_dev2 > 0.0 {
        (1.0 - sum_p2) / ((nd - 1.0) * sum_dev2)
      } else {
        0.0
      }
      if lambda < 0.0 {
        lambda = 0.0
      }
      if lambda > 1.0 {
        lambda = 1.0
      }
      let mut h = 0.0
      level = 1
      while level <= k_tot {
        let p_ml = counts.get_or_default(level, 0).to_double() / nd
        let p = (1.0 - lambda) * p_ml + lambda * u
        if p > 0.0 {
          h = h - p * mt_log2(p)
        }
        level = level + 1
      }
      h
    }
    "mi.sg" => {
      // Schuermann-Grassberger: beta = 1/K, p = (c + beta)/(n + K*beta)
      let beta = 1.0 / k_tot.to_double()
      let denom = nd + k_tot.to_double() * beta
      let mut h = 0.0
      let mut level = 1
      while level <= k_tot {
        let p = (counts.get_or_default(level, 0).to_double() + beta) / denom
        if p > 0.0 {
          h = h - p * mt_log2(p)
        }
        level = level + 1
      }
      h
    }
    _ => raise SeqError(
      "minet_entropy: unknown estimator '\{estimator}' (expected empirical|mi.mm|mi.shrink|mi.sg)",
    )
  }
}

///|
/// Joint entropy H(X, Y) over the `kx * ky` contingency table.
fn minet_joint_entropy(
  x : Array[Int],
  y : Array[Int],
  estimator : String,
  kx : Int,
  ky : Int,
) -> Double raise SeqError {
  let n = x.length()
  // contingency counts keyed by (bx, by)
  let joint : Map[Int, Int] = Map([])
  let mut nonempty = 0
  let mut max_bx = 0
  let mut max_by = 0
  let mut i = 0
  while i < n {
    let bx = x[i]
    let by = y[i]
    if bx > max_bx {
      max_bx = bx
    }
    if by > max_by {
      max_by = by
    }
    let key = bx * 100003 + by
    let prev = joint.get_or_default(key, 0)
    if prev == 0 {
      nonempty = nonempty + 1
    }
    joint.set(key, prev + 1)
    i = i + 1
  }
  let ktot_x = if kx > 0 { kx } else { max_bx }
  let ktot_y = if ky > 0 { ky } else { max_by }
  let k_tot = ktot_x * ktot_y
  let nd = n.to_double()
  match estimator {
    "empirical" => {
      let mut h = 0.0
      for kv in joint {
        let c = kv.1.to_double()
        let p = c / nd
        if p > 0.0 {
          h = h - p * mt_log2(p)
        }
      }
      h
    }
    "mi.mm" => {
      let mut h = 0.0
      for kv in joint {
        let c = kv.1.to_double()
        let p = c / nd
        if p > 0.0 {
          h = h - p * mt_log2(p)
        }
      }
      h + (nonempty.to_double() - 1.0) / (2.0 * nd)
    }
    "mi.shrink" => {
      let u = 1.0 / k_tot.to_double()
      let mut sum_p2 = 0.0
      let mut sum_dev2 = 0.0
      // iterate over all cells including empty ones
      let mut bx = 1
      while bx <= ktot_x {
        let mut by = 1
        while by <= ktot_y {
          let p = joint.get_or_default(bx * 100003 + by, 0).to_double() / nd
          sum_p2 = sum_p2 + p * p
          let d = p - u
          sum_dev2 = sum_dev2 + d * d
          by = by + 1
        }
        bx = bx + 1
      }
      let mut lambda = if sum_dev2 > 0.0 {
        (1.0 - sum_p2) / ((nd - 1.0) * sum_dev2)
      } else {
        0.0
      }
      if lambda < 0.0 {
        lambda = 0.0
      }
      if lambda > 1.0 {
        lambda = 1.0
      }
      let mut h = 0.0
      bx = 1
      while bx <= ktot_x {
        let mut by = 1
        while by <= ktot_y {
          let p_ml = joint.get_or_default(bx * 100003 + by, 0).to_double() / nd
          let p = (1.0 - lambda) * p_ml + lambda * u
          if p > 0.0 {
            h = h - p * mt_log2(p)
          }
          by = by + 1
        }
        bx = bx + 1
      }
      h
    }
    "mi.sg" => {
      let beta = 1.0 / k_tot.to_double()
      let denom = nd + k_tot.to_double() * beta
      let mut h = 0.0
      let mut bx = 1
      while bx <= ktot_x {
        let mut by = 1
        while by <= ktot_y {
          let p = (joint.get_or_default(bx * 100003 + by, 0).to_double() + beta) / denom
          if p > 0.0 {
            h = h - p * mt_log2(p)
          }
          by = by + 1
        }
        bx = bx + 1
      }
      h
    }
    _ => raise SeqError(
      "minet_joint_entropy: unknown estimator '\{estimator}' (expected empirical|mi.mm|mi.shrink|mi.sg)",
    )
  }
}

///|
/// Mutual information I(X;Y) = H(X) + H(Y) - H(X,Y) in bits, using the
/// given estimator consistently for marginals and the joint table.
pub fn minet_mutual_information(
  x : Array[Int],
  y : Array[Int],
  estimator? : String = "mi.mm",
) -> Double raise SeqError {
  if x.length() != y.length() || x.length() == 0 {
    raise SeqError("minet_mutual_information: variables must be non-empty and equal length")
  }
  let hx = minet_entropy(x, estimator=estimator)
  let hy = minet_entropy(y, estimator=estimator)
  let hxy = minet_joint_entropy(x, y, estimator, 0, 0)
  let mi = hx + hy - hxy
  // guard tiny negative values from floating error
  if mi < 0.0 && mi > -1.0e-9 {
    0.0
  } else {
    mi
  }
}

///|
/// Build the mutual information matrix (MIM) from a discretized
/// gene x sample matrix. The diagonal carries each variable's entropy,
/// matching minet's build.mim().
pub fn minet_build_mim(
  disc : Array[Array[Int]],
  estimator? : String = "mi.mm",
) -> Array[Array[Double]] raise SeqError {
  let n_genes = disc.length()
  if n_genes == 0 {
    raise SeqError("minet_build_mim: empty discretized matrix")
  }
  let nsamples = disc[0].length()
  let mut g = 0
  while g < n_genes {
    if disc[g].length() != nsamples {
      raise SeqError("minet_build_mim: ragged discretized matrix (row lengths differ)")
    }
    g = g + 1
  }
  let mim : Array[Array[Double]] = []
  let mut i = 0
  while i < n_genes {
    let row : Array[Double] = Array::make(n_genes, 0.0)
    mim.push(row)
    i = i + 1
  }
  i = 0
  while i < n_genes {
    mim[i][i] = minet_entropy(disc[i], estimator=estimator)
    let mut j = i + 1
    while j < n_genes {
      let mi = minet_mutual_information(disc[i], disc[j], estimator=estimator)
      mim[i][j] = mi
      mim[j][i] = mi
      j = j + 1
    }
    i = i + 1
  }
  mim
}

///|
/// ARACNE network: apply the Data Processing Inequality to every
/// triplet of variables and remove the weakest edge when it is
/// significantly smaller than the second-weakest (tolerance `eps`,
/// default 0.1). The MIM is not modified; a pruned copy is returned.
pub fn minet_aracne(mim : Array[Array[Double]], eps? : Double = 0.1) -> Array[Array[Double]] {
  let n = mim.length()
  let net : Array[Array[Double]] = []
  let mut i = 0
  while i < n {
    let row : Array[Double] = []
    let mut j = 0
    while j < mim[i].length() {
      row.push(mim[i][j])
      j = j + 1
    }
    net.push(row)
    i = i + 1
  }
  i = 0
  while i < n {
    let mut j = i + 1
    while j < n {
      let mut k = j + 1
      while k < n {
        let e_ij = net[i][j]
        let e_ik = net[i][k]
        let e_jk = net[j][k]
        // identify a unique weakest edge and prune it when it is weaker
        // than the second-weakest edge beyond the relative tolerance eps
        // (minet rule: wmin < w2 - eps * w2 = w2 * (1 - eps)).
        if e_ij < e_ik && e_ij < e_jk {
          if e_ij < e_ik * (1.0 - eps) && e_ij < e_jk * (1.0 - eps) {
            net[i][j] = 0.0
            net[j][i] = 0.0
          }
        } else if e_ik < e_ij && e_ik < e_jk {
          if e_ik < e_ij * (1.0 - eps) && e_ik < e_jk * (1.0 - eps) {
            net[i][k] = 0.0
            net[k][i] = 0.0
          }
        } else if e_jk < e_ij && e_jk < e_ik {
          if e_jk < e_ij * (1.0 - eps) && e_jk < e_ik * (1.0 - eps) {
            net[j][k] = 0.0
            net[k][j] = 0.0
          }
        }
        k = k + 1
      }
      j = j + 1
    }
    i = i + 1
  }
  net
}

///|
/// CLR network: normalize each MIM entry against the row and column
/// background distributions of MI values (mean and standard deviation
/// over each variable's n-1 off-diagonal entries; the standard deviation
/// follows R's sample-sd convention, dividing the n-1 deviations by
/// n-2), then combine the two positive z-scores quadratically.
pub fn minet_clr(mim : Array[Array[Double]]) -> Array[Array[Double]] {
  let n = mim.length()
  let means : Array[Double] = Array::make(n, 0.0)
  let sds : Array[Double] = Array::make(n, 0.0)
  let mut i = 0
  while i < n {
    let mut sum = 0.0
    let mut j = 0
    while j < n {
      if j != i {
        sum = sum + mim[i][j]
      }
      j = j + 1
    }
    let mean = sum / (n - 1).to_double()
    means[i] = mean
    let mut ss = 0.0
    j = 0
    while j < n {
      if j != i {
        let d = mim[i][j] - mean
        ss = ss + d * d
      }
      j = j + 1
    }
    // R sd() uses n-1 denominator; n-1 off-diagonal entries -> /(n-2)
    sds[i] = if n > 2 { (ss / (n - 2).to_double()).sqrt() } else { 0.0 }
    i = i + 1
  }
  let net : Array[Array[Double]] = []
  i = 0
  while i < n {
    let row : Array[Double] = Array::make(n, 0.0)
    net.push(row)
    i = i + 1
  }
  i = 0
  while i < n {
    let mut j = i + 1
    while j < n {
      let z_i = if sds[i] > 0.0 { (mim[i][j] - means[i]) / sds[i] } else { 0.0 }
      let z_j = if sds[j] > 0.0 { (mim[j][i] - means[j]) / sds[j] } else { 0.0 }
      let p_i = if z_i > 0.0 { z_i } else { 0.0 }
      let p_j = if z_j > 0.0 { z_j } else { 0.0 }
      let score = (p_i * p_i + p_j * p_j).sqrt()
      net[i][j] = score
      net[j][i] = score
      j = j + 1
    }
    i = i + 1
  }
  net
}

///|
/// MRNET network: for each target variable run the mRMR forward
/// selection (relevance = MI with target, redundancy = mean MI with
/// already selected variables; stop when the best score <= 0) and
/// symmetrize by taking the maximum of the two directed scores.
pub fn minet_mrnet(mim : Array[Array[Double]]) -> Array[Array[Double]] {
  let n = mim.length()
  let net : Array[Array[Double]] = []
  let mut i = 0
  while i < n {
    net.push(Array::make(n, 0.0))
    i = i + 1
  }
  let mut target = 0
  while target < n {
    // candidate list (gene indices except target)
    let cands : Array[Int] = []
    let mut g = 0
    while g < n {
      if g != target {
        cands.push(g)
      }
      g = g + 1
    }
    let selected : Array[Int] = []
    while cands.length() > 0 {
      let mut best = 0
      let mut best_score = -1.0e300
      let mut ci = 0
      while ci < cands.length() {
        let cand = cands[ci]
        let mut score = mim[target][cand]
        if selected.length() > 0 {
          let mut redundancy = 0.0
          let mut si = 0
          while si < selected.length() {
            redundancy = redundancy + mim[cand][selected[si]]
            si = si + 1
          }
          score = score - redundancy / selected.length().to_double()
        }
        if score > best_score {
          best_score = score
          best = ci
        }
        ci = ci + 1
      }
      if best_score <= 0.0 {
        break
      }
      let chosen = cands[best]
      net[target][chosen] = best_score
      selected.push(chosen)
      ignore(cands.remove(best))
    }
    target = target + 1
  }
  // symmetrize with max
  let mut i = 0
  while i < n {
    let mut j = i + 1
    while j < n {
      let v = if net[i][j] > net[j][i] { net[i][j] } else { net[j][i] }
      net[i][j] = v
      net[j][i] = v
      j = j + 1
    }
    i = i + 1
  }
  net
}

///|
/// Full minet pipeline: discretize a continuous expression matrix, build
/// the MIM and infer the network with `net_method` ("aracne", "clr" or
/// "mrnet", default "mrnet"). Rows are genes, columns samples.
pub fn minet(
  data : Array[Array[Double]],
  net_method? : String = "mrnet",
  estimator? : String = "mi.mm",
  disc? : String = "equalfreq",
  nbin? : Int = 0,
  eps? : Double = 0.1,
) -> Array[Array[Double]] raise SeqError {
  let d = minet_discretize(data, disc_method=disc, nbin=nbin)
  let mim = minet_build_mim(d, estimator=estimator)
  match net_method {
    "aracne" => minet_aracne(mim, eps=eps)
    "clr" => minet_clr(mim)
    "mrnet" => minet_mrnet(mim)
    _ => raise SeqError(
      "minet: unknown method '\{net_method}' (expected aracne|clr|mrnet)",
    )
  }
}

///|
/// Score an inferred weighted network against a binary reference
/// adjacency matrix (1 = true edge, 0 = no edge; diagonal ignored).
/// Returns precision/recall/FPR/F-score curves over all distinct edge
/// weights, plus AUROC (trapezoid), AUPR (trapezoid) and the maximal
/// F-score. Analogous to minet's validate().
pub fn minet_validate(
  net : Array[Array[Double]],
  reference : Array[Array[Int]],
) -> MinetValidation raise SeqError {
  let n = net.length()
  if reference.length() != n {
    raise SeqError("minet_validate: net and reference dimensions differ")
  }
  // collect upper-triangle edges
  let weights : Array[Double] = []
  let truths : Array[Int] = []
  let mut i = 0
  while i < n {
    if net[i].length() != n || reference[i].length() != n {
      raise SeqError("minet_validate: networks must be square matrices")
    }
    let mut j = i + 1
    while j < n {
      weights.push(net[i][j])
      truths.push(reference[i][j])
      j = j + 1
    }
    i = i + 1
  }
  let m = weights.length()
  // sort indices by descending weight
  let order : Array[Int] = []
  let mut k = 0
  while k < m {
    order.push(k)
    k = k + 1
  }
  k = 1
  while k < m {
    let key = order[k]
    let kv = weights[key]
    let mut j = k - 1
    while j >= 0 && weights[order[j]] < kv {
      order[j + 1] = order[j]
      j = j - 1
    }
    order[j + 1] = key
    k = k + 1
  }
  let total_pos = {
    let mut tp = 0
    let mut idx = 0
    while idx < m {
      if truths[idx] == 1 {
        tp = tp + 1
      }
      idx = idx + 1
    }
    tp
  }
  let total_neg = m - total_pos
  let thresholds : Array[Double] = []
  let precision : Array[Double] = []
  let recall : Array[Double] = []
  let fpr : Array[Double] = []
  let fscore : Array[Double] = []
  // initial point: threshold above the max weight -> empty prediction
  thresholds.push(if m > 0 { weights[order[0]] + 1.0 } else { 1.0 })
  precision.push(1.0)
  recall.push(0.0)
  fpr.push(0.0)
  fscore.push(0.0)
  let mut tp = 0
  let mut fp = 0
  k = 0
  while k < m {
    let idx = order[k]
    if truths[idx] == 1 {
      tp = tp + 1
    } else {
      fp = fp + 1
    }
    let prec = if (tp + fp) > 0 {
      tp.to_double() / (tp + fp).to_double()
    } else {
      1.0
    }
    let rec = if total_pos > 0 { tp.to_double() / total_pos.to_double() } else { 0.0 }
    let f_p = if total_neg > 0 { fp.to_double() / total_neg.to_double() } else { 0.0 }
    let fs = if prec + rec > 0.0 { 2.0 * prec * rec / (prec + rec) } else { 0.0 }
    thresholds.push(weights[idx])
    precision.push(prec)
    recall.push(rec)
    fpr.push(f_p)
    fscore.push(fs)
    k = k + 1
  }
  // AUROC: trapezoid over FPR sorted ascending (curve built descending
  // thresholds; reverse order), extend to (1, last_TPR).
  let mut auc_roc = 0.0
  let mut auc_pr = 0.0
  let mut max_f = 0.0
  k = 0
  while k < fscore.length() {
    if fscore[k] > max_f {
      max_f = fscore[k]
    }
    k = k + 1
  }
  // points are stored in order of decreasing threshold; iterate reverse
  // for ascending x (fpr / recall): x0 is the point with the larger x.
  let mut p = fpr.length() - 1
  while p > 0 {
    let x0 = fpr[p]
    let x1 = fpr[p - 1]
    let y0 = recall[p]
    let y1 = recall[p - 1]
    auc_roc = auc_roc + (x0 - x1) * (y0 + y1) / 2.0
    let rx0 = recall[p]
    let rx1 = recall[p - 1]
    let py0 = precision[p]
    let py1 = precision[p - 1]
    auc_pr = auc_pr + (rx0 - rx1) * (py0 + py1) / 2.0
    p = p - 1
  }
  // extend ROC to FPR = 1 at the last TPR reached
  let last_tpr = recall[recall.length() - 1]
  let last_fpr = fpr[fpr.length() - 1]
  auc_roc = auc_roc + (1.0 - last_fpr) * last_tpr
  MinetValidation::{
    thresholds,
    precision,
    recall,
    fpr,
    fscore,
    auc_roc,
    auc_pr,
    max_fscore: max_f,
  }
}
