///|
/// Tests for squidpy-inspired spatial metrics: nhood_enrichment,
/// co_occurrence and ripley K/L/G.

///|
test "scanpy_nhood_enrichment" {
  // 3 A cells mutually connected, 3 B cells mutually connected, no cross edges
  let clusters = ["A", "A", "A", "B", "B", "B"]
  let w : Array[Array[Double]] = Array(capacity=6)
  for _i in 0..<6 {
    w.push(Array::make(6, 0.0))
  }
  for i in 0..<3 {
    for j in 0..<3 {
      if i != j {
        w[i][j] = 1.0
        w[j][i] = 1.0
      }
    }
  }
  for i in 3..<6 {
    for j in 3..<6 {
      if i != j {
        w[i][j] = 1.0
        w[j][i] = 1.0
      }
    }
  }
  let res = @src.scanpy_nhood_enrichment(clusters, w)
  assert_eq(res.clusters.length(), 2)
  // raw counts: each A cell has 2 A neighbors -> 6 ordered A-A pairs
  let aa = match res.clusters[0] == "A" {
    true => res.counts[0][0]
    false => res.counts[1][1]
  }
  assert_eq(aa, 6.0)
  // cross-cluster count is zero and z-score strongly negative
  let ab_counts = res.counts[0][1]
  assert_eq(ab_counts, 0.0)
  let ab_z = res.zscores[0][1]
  assert_true(ab_z < -1.0)
  // same-cluster z-score positive
  let aa_i = if res.clusters[0] == "A" { 0 } else { 1 }
  assert_true(res.zscores[aa_i][aa_i] > 1.0)
}

///|
test "scanpy_co_occurrence" {
  // two blocks at distance 1 within, 100 across
  let clusters = ["A", "A", "A", "B", "B", "B"]
  let d : Array[Array[Double]] = Array(capacity=6)
  for _i in 0..<6 {
    d.push(Array::make(6, 0.0))
  }
  for i in 0..<6 {
    for j in 0..<6 {
      if i != j {
        d[i][j] = if (i < 3) == (j < 3) { 1.0 } else { 100.0 }
      }
    }
  }
  let res = @src.scanpy_co_occurrence(clusters, d, n_bands=2)
  assert_eq(res.intervals.length(), 2)
  // band 1 (0,50]: only same-cluster neighbors -> A-A ratio > 1, A-B = 0
  assert_true(res.ratios[0][0][0] > 1.0)
  assert_eq(res.ratios[0][0][1], 0.0)
  // band 2 (50,100]: only cross-cluster neighbors -> A-B ratio > 1
  assert_true(res.ratios[1][0][1] > 1.0)
  assert_eq(res.ratios[1][0][0], 0.0)
}

///|
test "scanpy_ripley_k" {
  // 3x3 grid of points in a 10x10 area
  let coords : Array[Array[Double]] = Array(capacity=9)
  for gx in 0..<3 {
    for gy in 0..<3 {
      coords.push([gx.to_double() * 5.0, gy.to_double() * 5.0])
    }
  }
  let res = @src.scanpy_ripley_k(coords, radius_max=71.0, n_steps=5)
  assert_eq(res.radii.length(), 5)
  // radii and K increase; G increases to 1 at large radius
  for s in 1..<5 {
    assert_true(res.radii[s] > res.radii[s - 1])
    assert_true(res.k[s] >= res.k[s - 1])
    assert_true(res.g[s] >= res.g[s - 1])
  }
  assert_eq(res.g[4], 1.0)
  // K(0) = 0 -> L(0) = 0
  assert_eq(res.l[0], 0.0)
  // at t=71 every ordered pair (72) is inside: K = area / (n(n-1)) * 72
  assert_eq(res.k[4], 100.0 / (9.0 * 8.0) * 72.0)
}
