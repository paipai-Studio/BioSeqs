///|
/// swfdr demo:
///   1. Covariate-modulated pi0 estimation (lm_pi0) and q-values
///      (lm_qvalue) with unit.spline and smooth.spline smoothing.
///   2. Science-wise false discovery rate for rounded and truncated
///      p-values (swfdr_calculate, the EM model of the `science` table).
///
/// Run with: moon run examples/swfdr_demo
fn round4(x : Double) -> Double {
  (x * 10000.0).round() / 10000.0
}

///|
fn main raise {
  // ===== 1. Covariate-modulated pi0 and q-values =========================
  println("== covariate-modulated FDR (40 hypotheses) ==")
  // p-values and a single covariate on a grid from -1 to 2; small p-values
  // cluster at the high-covariate end, so estimated pi0 decreases with x.
  let n = 40
  let p = [
    0.02, 0.97, 0.163581, 0.872852, 0.553713, 0.910457, 0.668837, 0.385033, 0.623886,
    0.170607, 0.391942, 0.851453, 0.44818, 0.41567, 0.279233, 0.638824, 0.654353,
    0.605966, 0.650483, 0.061336, 0.74109, 0.452955, 0.348697, 0.529169, 0.075983,
    0.48949, 0.753233, 0.047598, 0.715859, 0.019905, 0.191923, 0.456433, 0.278204,
    0.690245, 0.074963, 0.426591, 0.372464, 0.327325, 0.63941, 0.306681,
  ]
  let cov = [
    -1.0, -0.923077, -0.846154, -0.769231, -0.692308, -0.615385, -0.538462, -0.461538,
    -0.384615, -0.307692, -0.230769, -0.153846, -0.076923, 0.0, 0.076923, 0.153846,
    0.230769, 0.307692, 0.384615, 0.461538, 0.538462, 0.615385, 0.692308, 0.769231,
    0.846154, 0.923077, 1.0, 1.076923, 1.153846, 1.230769, 1.307692, 1.384615, 1.461538,
    1.538462, 1.615385, 1.692308, 1.769231, 1.846154, 1.923077, 2.0,
  ]
  // design matrix: one row per hypothesis (x is n x p1, row oriented).
  let xrows : Array[Array[Double]] = []
  for v in cov {
    xrows.push([v])
  }
  let xm = Some(xrows)
  let mut i = 0

  for smoothing in ["unit.spline", "smooth.spline"] {
    let fit = @bio.swfdr_lm_pi0(p, x=xm, smoothing~, smooth_df=3)
    let qv = @bio.swfdr_lm_qvalue(p, x=xm, pi0=Some(fit))
    // average pi0 in the low- and high-covariate halves.
    let mut lo_sum = 0.0
    let mut hi_sum = 0.0
    let mut lo_n = 0
    let mut hi_n = 0
    let mut sig = 0
    i = 0
    while i < n {
      if cov[i] < 0.0 {
        lo_sum += fit.pi0()[i]
        lo_n += 1
      } else {
        hi_sum += fit.pi0()[i]
        hi_n += 1
      }
      if qv.qvalues()[i] < 0.05 {
        sig += 1
      }
      i = i + 1
    }
    println(
      "\{smoothing}: mean pi0 low-cov=\{round4(lo_sum / lo_n.to_double())} " +
      "high-cov=\{round4(hi_sum / hi_n.to_double())}, q < 0.05: \{sig} hypotheses",
    )
  }

  // ===== 2. Rounded / truncated p-values (science table EM) ===============
  println("")
  println("== science-wise FDR with rounded and truncated p-values ==")
  // 60 p-values; indices 41-50 and 56-57 (1-based) are reported rounded to
  // the coarse 0.001/0.005 grid, and 51-55 are truncated upper bounds.
  let pp = [
    0.58471185, 0.11378168, 0.68426474, 0.99250878, 0.53499358, 0.96661406, 0.67142756,
    0.29457771, 0.35836298, 0.17531475, 0.54881739, 0.5054517, 0.19383647, 0.63690411,
    0.68780009, 0.64019077, 0.35788536, 0.102585, 0.09779092, 0.18288626, 0.22790347,
    0.08052415, 0.82161842, 0.59111398, 0.77338902, 0.35008598, 0.00606132, 0.81450622,
    0.00117923, 0.20069018, 0.01314413, 0.00013078, 0.01661963, 0.05042985, 0.00544384,
    3.049e-05, 3.492e-05, 0.00217114, 0.02335548, 0.00793186, 0.003, 0.011, 0.022,
    0.031, 0.047, 0.001, 0.018, 0.027, 0.039, 0.008, 0.05, 0.04, 0.03, 0.02, 0.01,
    0.0, 0.05, 0.0051082, 0.01139009, 0.00418146,
  ]
  let rounded = Array::make(60, 0)
  i = 40
  while i < 50 {
    rounded[i] = 1
    i = i + 1
  }
  rounded[55] = 1
  rounded[56] = 1
  let truncated = Array::make(60, 0)
  i = 50
  while i < 55 {
    truncated[i] = 1
    i = i + 1
  }

  let science = @bio.swfdr_calculate(
    pp,
    truncated,
    rounded,
    num_em_iterations=100,
  )
  println(
    "EM fit (100 iterations): pi0=\{round4(science.pi0())} " +
    "alpha=\{round4(science.alpha())} beta=\{round4(science.beta())}",
  )
  // posterior null probability is only defined for precise p-values;
  // rounded/truncated positions report None.
  let z_small = science.z()[26].unwrap() // p = 0.00606132
  let z_big = science.z()[0].unwrap() // p = 0.58471185
  println(
    "posterior null probability: exact p=0.00606 -> \{round4(z_small)}, " +
    "exact p=0.58471 -> \{round4(z_big)}",
  )
  println("rounded p-values: posterior null indicator is None (by design)")
  println("coarse bins (0,.005] ... (.045,.05]:")
  let mut n0line = "  expected null n0:"
  let mut nline = "  observed count n :"
  let mut k = 0
  while k < science.n0().length() {
    n0line += " " + round4(science.n0()[k]).to_string()
    nline += " " + science.n()[k].to_string()
    k += 1
  }
  println(n0line)
  println(nline)

  let one = @bio.swfdr_calculate(pp, truncated, rounded, num_em_iterations=1)
  println(
    "single EM iteration:    pi0=\{round4(one.pi0())} alpha=\{round4(one.alpha())} " +
    "beta=\{round4(one.beta())}",
  )
}
