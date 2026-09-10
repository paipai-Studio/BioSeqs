///|
/// copynumber demo: Piecewise Constant Fitting (PCF) of copy-number data,
/// ported from the Bioconductor copynumber package (Nilsen & Liestoel et al.
/// 2012). Shows MAD winsorization of outliers, single-sample PCF and joint
/// multi-sample multiPCF with a shared segmentation.
///
/// Run: moon run examples/copynumber_demo
fn demo_rand(seed : Int) -> Double {
  let s = (1103515245 * seed + 12345) & 0x7FFFFFFF
  s.to_double() / 2147483647.0
}

///|
/// One noisy 120-probe sample with four underlying levels (edges 30/60/90).
fn make_profile(
  levels : Array[Double],
  seed : Int,
  spikes : Array[(Int, Double)],
) -> Array[Double] {
  let n = 120
  let edges = [0, 30, 60, 90, n]
  let y = Array::new()
  let mut i = 0
  while i < n {
    let mut level = levels[3]
    let mut b = 0
    while b < 4 {
      if i >= edges[b] && i < edges[b + 1] {
        level = levels[b]
      }
      b = b + 1
    }
    y.push(level + (demo_rand(seed + i) - 0.5) * 0.1)
    i = i + 1
  }
  let mut t = 0
  while t < spikes.length() {
    let (idx, delta) = spikes[t]
    y[idx] = y[idx] + delta
    t = t + 1
  }
  y
}

///|
/// One character per probe, binned by log-ratio magnitude.
fn ascii_track(y : Array[Double]) -> String {
  let sb = StringBuilder()
  for v in y {
    sb.write_char(if v > 0.25 { '+' } else if v < -0.25 { '-' } else { '.' })
  }
  sb.to_string()
}

///|
/// Right-align an integer in a fixed-width field.
fn field_int(i : Int, width : Int) -> String {
  let sb = StringBuilder()
  let s = i.to_string()
  let mut n = s.length()
  while n < width {
    sb.write_char(' ')
    n = n + 1
  }
  sb.write_string(s)
  sb.to_string()
}

///|
fn main raise {
  println("=== copynumber Demo (PCF copy-number segmentation) ===\n")

  // 1. A noisy four-level profile with two spikes (a gain and a loss)
  let y = make_profile([0.0, 0.45, -0.35, 0.1], 17, [(12, 1.2), (75, -1.0)])
  println("1. Synthetic profile: 120 probes, levels 0 / +0.45 / -0.35 / +0.1")
  println(
    "   raw log-ratios:        \{ascii_track(y)}  (+ gain, - loss, . neutral)",
  )

  // 2. MAD winsorization clips single-probe outliers against a median filter
  let w = @bio.winsorize_mad(y)
  let flagged = StringBuilder()
  let mut i = 0
  while i < w.outliers.length() {
    if w.outliers[i] != 0 {
      flagged.write_string(" \{i + 1}")
    }
    i = i + 1
  }
  println("\n2. winsorize_mad (tau=2.5, k=25): MAD = \{w.sdev}")
  println("   outlier probes (1-based):\{flagged}")
  println("   winsorized log-ratios: \{ascii_track(w.ywin)}")

  // 3. PCF on the winsorized profile; gamma penalizes segment count
  let fit = @bio.pcf_plain(w.ywin)
  println(
    "\n3. pcf_plain (kmin=5, gamma=40, MAD-normalized): \{fit.n_intervals} segments",
  )
  for seg in fit.segments {
    println(
      "   probes \{seg.start}-\{seg.start + seg.length - 1}  n=\{seg.length}  mean=\{seg.mean}",
    )
  }
  let fitted = Array::new()
  for seg in fit.segments {
    let mut j = 0
    while j < seg.length {
      fitted.push(seg.mean)
      j = j + 1
    }
  }
  println("   PCF fitted track:       \{ascii_track(fitted)}")

  // 4. multiPCF: three samples share one segmentation, per-sample levels
  let s1 = make_profile([0.0, 0.45, -0.35, 0.1], 17, [(12, 1.2), (75, -1.0)])
  let s2 = make_profile([0.02, 0.3, -0.2, 0.05], 211, [])
  let s3 = make_profile([-0.01, 0.55, -0.4, 0.08], 909, [])
  let multi = @bio.multi_pcf([s1, s2, s3])
  println(
    "\n4. multi_pcf on 3 samples: \{multi.n_intervals} shared segments, gamma scaled by n",
  )
  println("   sample MADs: \{multi.sds[0]}  \{multi.sds[1]}  \{multi.sds[2]}")
  println("   probes       n    sample1   sample2   sample3")
  for seg in multi.segments {
    println(
      "   \{field_int(seg.start, 3)}-\{field_int(seg.start + seg.length - 1, 3)} \{field_int(seg.length, 4)}    \{seg.means[0]}  \{seg.means[1]}  \{seg.means[2]}",
    )
  }

  println("\nDone.")
}
