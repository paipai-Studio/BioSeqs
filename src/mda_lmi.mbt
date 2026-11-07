///|
/// Linear mutual information (MDAnalysis `lmi`): mutual information between
/// two scalar signals estimated on equal-width histograms with natural
/// logarithms, plus the pairwise MI matrix and a normalized variant.

///|
/// Bin index of `v` on the equal-width histogram of `nb` bins over
/// `[lo, hi]`. The right endpoint `hi` is folded into the last bin.
fn mda_lmi_bin_index(v : Double, lo : Double, hi : Double, nb : Int32) -> Int32 {
  let mut b = ((v - lo) / (hi - lo) * nb.to_double()).to_int()
  if b >= nb {
    b = nb - 1
  }
  if b < 0 {
    b = 0
  }
  b
}

///|
/// Equal-width histogram counts of `x` over `[min(x), max(x)]`.
/// Returns `(counts, lo, hi)`; for a constant sequence `lo == hi` and all
/// counts collapse (the caller decides the semantics).
fn mda_lmi_hist_counts(
  x : Array[Double],
  nb : Int32,
) -> (Array[Int32], Double, Double) {
  let n = x.length()
  let counts = Array::make(nb, 0)
  if n == 0 || nb <= 0 {
    return (counts, 0.0, 0.0)
  }
  let mut lo = x[0]
  let mut hi = x[0]
  for v in x {
    if v < lo {
      lo = v
    }
    if v > hi {
      hi = v
    }
  }
  if hi <= lo {
    return (counts, lo, hi)
  }
  for v in x {
    let b = mda_lmi_bin_index(v, lo, hi, nb)
    counts[b] = counts[b] + 1
  }
  (counts, lo, hi)
}

///|
/// Shannon entropy (natural logarithm) of the equal-width histogram of
/// `x` with `nbins` bins: `H(X) = -sum_b p_b ln p_b`. An empty input, a
/// non-positive bin count, or a constant sequence yields `0.0`.
pub fn mda_entropy_hist(x : Array[Double], nbins~ : Int32 = 20) -> Double {
  let nb = if nbins < 1 { 1 } else { nbins }
  let n = x.length()
  if n == 0 {
    return 0.0
  }
  let (counts, lo, hi) = mda_lmi_hist_counts(x, nb)
  if hi <= lo {
    return 0.0
  }
  let mut h = 0.0
  for c in counts {
    if c > 0 {
      let p = c.to_double() / n.to_double()
      h = h - p * @math.ln(p)
    }
  }
  h
}

///|
/// Core mutual information on equal-width joint histograms (natural log),
/// assuming non-empty equal-length inputs. Constant signals yield `0.0`.
fn mda_lmi_core(x : Array[Double], y : Array[Double], nb : Int32) -> Double {
  let n = x.length()
  let (cx, lox, hix) = mda_lmi_hist_counts(x, nb)
  let (cy, loy, hiy) = mda_lmi_hist_counts(y, nb)
  if hix <= lox || hiy <= loy {
    return 0.0
  }
  let joint : Array[Array[Int32]] = Array(capacity=nb)
  for i = 0; i < nb; i = i + 1 {
    joint.push(Array::make(nb, 0))
  }
  for k = 0; k < n; k = k + 1 {
    let bx = mda_lmi_bin_index(x[k], lox, hix, nb)
    let by = mda_lmi_bin_index(y[k], loy, hiy, nb)
    joint[bx][by] = joint[bx][by] + 1
  }
  let mut mi = 0.0
  let nf = n.to_double()
  for i = 0; i < nb; i = i + 1 {
    if cx[i] > 0 {
      let px = cx[i].to_double() / nf
      for j = 0; j < nb; j = j + 1 {
        let c = joint[i][j]
        if c > 0 && cy[j] > 0 {
          let pxy = c.to_double() / nf
          let py = cy[j].to_double() / nf
          mi = mi + pxy * @math.ln(pxy / (px * py))
        }
      }
    }
  }
  mi
}

///|
/// Linear mutual information `I(X;Y)` between two equal-length signals
/// (MDAnalysis `lmi`): estimated on an `nbins x nbins` equal-width joint
/// histogram with natural logarithms,
/// `I = sum pxy ln(pxy / (px py))`.
///
/// When either marginal entropy vanishes (constant sequence) the result is
/// `0.0`. Raises if `x` and `y` have different lengths or are empty.
pub fn mda_lmi(
  x : Array[Double],
  y : Array[Double],
  nbins~ : Int32 = 20,
) -> Double raise {
  if x.length() != y.length() {
    raise SeqError("mda_lmi: x and y must have the same length")
  }
  if x.length() == 0 {
    raise SeqError("mda_lmi: input signals must be non-empty")
  }
  let nb = if nbins < 1 { 1 } else { nbins }
  mda_lmi_core(x, y, nb)
}

///|
/// Pairwise linear mutual information matrix (MDAnalysis `lmi` over all
/// signal pairs). The result is symmetric with `M[i][i] = H(X_i) = I(X_i;
/// X_i)`; the off-diagonal entries are `I(X_i; X_j)`. Rows or columns of
/// length zero yield an empty matrix.
pub fn mda_lmi_matrix(
  signals : Array[Array[Double]],
  nbins~ : Int32 = 20,
) -> Array[Array[Double]] {
  let n = signals.length()
  let nb = if nbins < 1 { 1 } else { nbins }
  let m : Array[Array[Double]] = Array(capacity=n)
  for i = 0; i < n; i = i + 1 {
    m.push(Array::make(n, 0.0))
  }
  for i = 0; i < n; i = i + 1 {
    let si = signals[i]
    if si.length() > 0 {
      m[i][i] = mda_entropy_hist(si, nbins~)
      for j = i + 1; j < n; j = j + 1 {
        let v = mda_lmi_core(si, signals[j], nb)
        m[i][j] = v
        m[j][i] = v
      }
    }
  }
  m
}

///|
/// Normalized linear mutual information `I / sqrt(H(X) H(Y))` in `[0, 1]`;
/// returns `0.0` when either entropy (or their product) is zero.
pub fn mda_lmi_normalized(i : Double, hx : Double, hy : Double) -> Double {
  let den = (hx * hy).sqrt()
  if den <= 0.0 {
    0.0
  } else {
    i / den
  }
}
