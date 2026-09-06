///|
/// Tests for the minet mutual-information network module (port of
/// Bioconductor minet: Meyer, Lafitte & Bontempi, BMC Bioinformatics 2008).
/// Values are hand-computed from deterministic synthetic data; no R
/// runtime is required.

///|
fn minet_close(a : Double, b : Double, tol? : Double = 1.0e-9) -> Bool {
  (a - b).abs() < tol
}

///|
test "minet_discretize_methods" {
  // equalwidth: 0..9 split into 5 bins of width 1.8 -> 1,1,2,2,3,3,4,4,5,5
  let ew : Array[Double] = [0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0]
  let dew = @src.minet_discretize([ew], disc_method="equalwidth", nbin=5)
  assert_eq(dew[0], [1, 1, 2, 2, 3, 3, 4, 4, 5, 5])
  // equalfreq on a shuffled 0..9: each value v has rank v, so the bin
  // follows floor(v*5/10)+1 regardless of position
  let shuffled : Array[Double] = [9.0, 0.0, 8.0, 1.0, 7.0, 2.0, 6.0, 3.0, 5.0, 4.0]
  let def_ = @src.minet_discretize([shuffled], disc_method="equalfreq", nbin=5)
  assert_eq(def_[0], [5, 1, 5, 1, 4, 2, 4, 2, 3, 3])
  // default nbin = floor(sqrt(n_samples)): 25 samples -> 5 bins covering 1..5
  let v25 : Array[Double] = []
  let mut i = 0
  while i < 25 {
    v25.push(i.to_double())
    i = i + 1
  }
  let d25 = @src.minet_discretize([v25])
  let mut lo = 999999
  let mut hi = 0
  for b in d25[0] {
    if b < lo {
      lo = b
    }
    if b > hi {
      hi = b
    }
  }
  assert_eq(lo, 1)
  assert_eq(hi, 5)
  // ties must never straddle a bin boundary
  let ties : Array[Double] = [0.0, 0.0, 1.0, 1.0, 2.0, 2.0, 3.0, 3.0, 4.0, 4.0]
  let dt = @src.minet_discretize([ties], nbin=3)
  let mut k = 0
  while k < 5 {
    assert_eq(dt[0][2 * k], dt[0][2 * k + 1])
    k = k + 1
  }
  // nbin > number of samples raises
  let too_many = try {
    ignore(@src.minet_discretize([[0.0, 1.0]], nbin=5))
    false
  } catch {
    SeqError(_) => true
  }
  assert_true(too_many)
  // unknown discretization method raises
  let bad_disc = try {
    ignore(@src.minet_discretize([[0.0, 1.0]], disc_method="nope", nbin=2))
    false
  } catch {
    SeqError(_) => true
  }
  assert_true(bad_disc)
  // ragged matrix raises
  let ragged = try {
    ignore(@src.minet_discretize([[0.0, 1.0], [0.0, 1.0, 2.0]], nbin=2))
    false
  } catch {
    SeqError(_) => true
  }
  assert_true(ragged)
  // empty matrix raises
  let empty = try {
    ignore(@src.minet_discretize([]))
    false
  } catch {
    SeqError(_) => true
  }
  assert_true(empty)
}

///|
test "minet_entropy_estimators" {
  // perfectly balanced binary variable: H = 1 bit for every estimator
  let bal = [1, 1, 2, 2]
  assert_true(minet_close(@src.minet_entropy(bal, estimator="empirical"), 1.0))
  assert_true(minet_close(@src.minet_entropy(bal, estimator="mi.mm"), 1.125))
  assert_true(minet_close(@src.minet_entropy(bal, estimator="mi.shrink"), 1.0))
  assert_true(minet_close(@src.minet_entropy(bal, estimator="mi.sg"), 1.0))
  // 3:1 split
  let skew = [1, 1, 1, 2]
  assert_true(minet_close(@src.minet_entropy(skew, estimator="empirical"), 0.811278, tol=1.0e-6))
  assert_true(minet_close(@src.minet_entropy(skew, estimator="mi.mm"), 0.936278, tol=1.0e-6))
  // skew + shrinkage collapses to the uniform target (lambda = 1)
  assert_true(minet_close(@src.minet_entropy(skew, estimator="mi.shrink"), 1.0))
  // Schuermann-Grassberger beta prior
  assert_true(minet_close(@src.minet_entropy(skew, estimator="mi.sg"), 0.881291, tol=1.0e-6))
  // four equally likely levels
  let four = [1, 2, 3, 4]
  assert_true(minet_close(@src.minet_entropy(four, estimator="empirical"), 2.0))
  assert_true(minet_close(@src.minet_entropy(four, estimator="mi.mm"), 2.375))
  // unknown estimator raises
  let bad = try {
    ignore(@src.minet_entropy(bal, estimator="bogus"))
    false
  } catch {
    SeqError(_) => true
  }
  assert_true(bad)
  // empty variable raises
  let no_data = try {
    ignore(@src.minet_entropy([]))
    false
  } catch {
    SeqError(_) => true
  }
  assert_true(no_data)
}

///|
test "minet_mutual_information_values" {
  let x = [1, 1, 2, 2]
  let y = [1, 2, 1, 2]
  // independent balanced variables: joint is uniform over 4 cells -> MI 0
  assert_true(minet_close(@src.minet_mutual_information(x, y, estimator="empirical"), 0.0))
  // identical variable: MI(X;X) = H(X) = 1 bit empirical, 1.125 with mm
  assert_true(minet_close(@src.minet_mutual_information(x, x, estimator="empirical"), 1.0))
  assert_true(minet_close(@src.minet_mutual_information(x, x, estimator="mi.mm"), 1.125))
  // relabeled perfect copy is equally informative
  let y2 = [3, 3, 7, 7]
  assert_true(minet_close(@src.minet_mutual_information(x, y2, estimator="empirical"), 1.0))
  // constant variable carries no information
  let cst = [1, 1, 1, 1]
  assert_true(minet_close(@src.minet_mutual_information(x, cst, estimator="empirical"), 0.0))
  // mismatched lengths raise
  let bad = try {
    ignore(@src.minet_mutual_information(x, [1, 2, 3]))
    false
  } catch {
    SeqError(_) => true
  }
  assert_true(bad)
}

///|
test "minet_network_methods_hand_matrix" {
  // hand-built MIM: diagonal = entropy, edges AB = .9 (strong),
  // BC = .8 (strong), AC = .2 (indirect, weakest by far)
  let mim : Array[Array[Double]] = [
    [2.0, 0.9, 0.2],
    [0.9, 2.0, 0.8],
    [0.2, 0.8, 2.0],
  ]
  // ARACNE: DPI prunes the indirect AC edge (0.2 < 0.8*(1-0.1))
  let ar = @src.minet_aracne(mim, eps=0.1)
  assert_true(ar[0][2] == 0.0)
  assert_true(ar[2][0] == 0.0)
  assert_true(minet_close(ar[0][1], 0.9))
  assert_true(minet_close(ar[1][2], 0.8))
  // with a huge tolerance nothing is pruned
  let ar_tol = @src.minet_aracne(mim, eps=3.0)
  assert_true(minet_close(ar_tol[0][2], 0.2))
  // CLR background z-scores (sd over n-1 background edges divides by
  // n-2 as in R): AB = 1.0, BC = 1/sqrt(2), AC = 0
  let clr = @src.minet_clr(mim)
  assert_true(minet_close(clr[0][1], 1.0, tol=1.0e-9))
  assert_true(minet_close(clr[1][2], 0.707107, tol=1.0e-6))
  assert_true(clr[0][2] == 0.0)
  // MRNET: mRMR selects the two direct edges for every target
  let mr = @src.minet_mrnet(mim)
  assert_true(minet_close(mr[0][1], 0.9))
  assert_true(minet_close(mr[1][2], 0.8))
  assert_true(mr[0][2] == 0.0)
  // the MIM itself is never modified
  assert_true(minet_close(mim[0][2], 0.2))
}

///|
/// Synthetic chain of four variables over 20 samples:
///   G1 ~ balanced binary; G2 = G1 with samples 0,1,2 flipped;
///   G3 = G2 with sample 10 flipped; G4 independent (alternating).
/// Direct links are G1-G2 and G2-G3; G1-G3 is indirect and should
/// be pruned by DPI / downweighted by CLR / skipped by mRMR.
fn minet_chain_genes() -> Array[Array[Int]] {
  let g1 : Array[Int] = []
  let g2 : Array[Int] = []
  let g3 : Array[Int] = []
  let g4 : Array[Int] = []
  let mut i = 0
  while i < 20 {
    let v1 = if i < 10 { 1 } else { 2 }
    // flips 0,1,2 in g2 (all v1 == 1 there)
    let v2 = if i == 0 || i == 1 || i == 2 { 2 } else { v1 }
    // flip sample 10 in g3 (v2 == 2 there)
    let v3 = if i == 10 { 1 } else { v2 }
    let v4 = if i % 2 == 0 { 1 } else { 2 }
    g1.push(v1)
    g2.push(v2)
    g3.push(v3)
    g4.push(v4)
    i = i + 1
  }
  [g1, g2, g3, g4]
}

///|
test "minet_chain_recovery" {
  let genes = minet_chain_genes()
  let mim = @src.minet_build_mim(genes, estimator="empirical")
  // MIM is symmetric and the diagonal holds each entropy
  assert_true(minet_close(mim[0][0], 1.0))
  assert_true(minet_close(mim[0][1], mim[1][0]))
  // hand-computed MI values for this chain
  assert_true(minet_close(mim[0][1], 0.4934, tol=1.0e-3))
  assert_true(minet_close(mim[1][2], 0.7166, tol=1.0e-3))
  assert_true(minet_close(mim[0][2], 0.2958, tol=1.0e-3))
  // G4 is independent of G1 and G3 by construction
  assert_true(minet_close(mim[0][3], 0.0))
  assert_true(minet_close(mim[2][3], 0.0))
  // ARACNE removes the indirect G1-G3 edge; direct edges survive
  let ar = @src.minet_aracne(mim, eps=0.1)
  assert_true(ar[0][2] == 0.0)
  assert_true(ar[0][1] > 0.0)
  assert_true(ar[1][2] > 0.0)
  // G4 keeps at most a negligible residual edge
  assert_true(ar[1][3] < 0.05)
  assert_true(ar[0][3] == 0.0)
  assert_true(ar[2][3] == 0.0)
  // CLR ranks the indirect edge below both direct edges
  let clr = @src.minet_clr(mim)
  assert_true(clr[0][2] < clr[0][1])
  assert_true(clr[0][2] < clr[1][2])
  assert_true(clr[0][1] > 0.0)
  assert_true(clr[1][2] > 0.0)
  // MRNET never selects the indirect edge (redundancy exceeds relevance)
  let mr = @src.minet_mrnet(mim)
  assert_true(mr[0][2] == 0.0)
  assert_true(mr[0][1] > 0.0)
  assert_true(mr[1][2] > 0.0)
  // ragged discretized matrix raises
  let bad = try {
    ignore(@src.minet_build_mim([[1, 1, 2, 2], [1, 2]], estimator="empirical"))
    false
  } catch {
    SeqError(_) => true
  }
  assert_true(bad)
}

///|
test "minet_pipeline_end_to_end" {
  // continuous expression of the same chain: two well-separated levels
  // so equal-frequency discretization with 2 bins recovers the pattern
  let genes_int = minet_chain_genes()
  let data : Array[Array[Double]] = []
  let mut g = 0
  while g < genes_int.length() {
    let row : Array[Double] = []
    let mut s = 0
    while s < genes_int[g].length() {
      row.push(if genes_int[g][s] == 1 { 10.0 } else { 20.0 })
      s = s + 1
    }
    data.push(row)
    g = g + 1
  }
  for net_method in ["aracne", "clr", "mrnet"] {
    let net = @src.minet(
      data, net_method=net_method, estimator="empirical", disc="equalwidth", nbin=2
    )
    assert_eq(net.length(), 4)
    assert_eq(net[0].length(), 4)
    assert_true(net[0][1] > 0.0)
    assert_true(net[1][2] > 0.0)
    // indirect edge G1-G3 is absent or strictly weakest
    assert_true(net[0][2] < net[0][1])
    assert_true(net[0][2] < net[1][2])
  }
  // unknown network method raises
  let bad = try {
    ignore(
      @src.minet(data, net_method="santa", estimator="empirical", disc="equalwidth", nbin=2)
    )
    false
  } catch {
    SeqError(_) => true
  }
  assert_true(bad)
}

///|
test "minet_validate_scores" {
  // 4 variables, upper-triangle edges sorted by weight:
  //   AB = .9 (true), AD = .85 (false), BC = .8 (true),
  //   AC = .2 (false), BD = 0, CD = 0
  let net : Array[Array[Double]] = [
    [0.0, 0.9, 0.2, 0.85],
    [0.9, 0.0, 0.8, 0.0],
    [0.2, 0.8, 0.0, 0.0],
    [0.85, 0.0, 0.0, 0.0],
  ]
  // reference: only AB and BC are true edges
  let reference : Array[Array[Int]] = [
    [0, 1, 0, 0],
    [1, 0, 1, 0],
    [0, 1, 0, 0],
    [0, 0, 0, 0],
  ]
  let v = @src.minet_validate(net, reference)
  // initial point + one point per edge
  assert_eq(v.thresholds.length(), 7)
  // the curve reaches full recall and full FPR
  assert_true(minet_close(v.recall[v.recall.length() - 1], 1.0))
  assert_true(minet_close(v.fpr[v.fpr.length() - 1], 1.0))
  // hand-computed areas: AUROC = 0.875, AUPR = 19/24 ~= 0.79167
  assert_true(minet_close(v.auc_roc, 0.875, tol=1.0e-9))
  assert_true(minet_close(v.auc_pr, 0.791667, tol=1.0e-5))
  // best F-score at precision 2/3, recall 1
  assert_true(minet_close(v.max_fscore, 0.8, tol=1.0e-9))
  // dimension mismatch raises
  let bad = try {
    ignore(@src.minet_validate(net, [[0, 1], [1, 0]]))
    false
  } catch {
    SeqError(_) => true
  }
  assert_true(bad)
}
