///|
/// Tests for the Bioconductor apcluster port (affinity propagation).
///
/// Expected exemplars / net similarities come from a numpy mirror script
/// implementing the same message-passing dynamics (median preference,
/// lam = 0.9, convits = 100) over the same deterministic synthetic data.

///|
fn ap_approx(a : Double, b : Double, tol : Double) -> Bool {
  (a - b).abs() <= tol
}

///|
/// Four tight blobs of six points each (deterministic jitter, shared with
/// the numpy reference script).
fn ap_blobs() -> Array[Array[Double]] {
  [
    [-0.1254598811526375, 0.4507143064099162],
    [0.2319939418114051, 0.0986584841970366],
    [-0.3439813595575635, -0.3440054796637974],
    [-0.4419163878318005, 0.3661761457749352],
    [0.1011150117432088, 0.2080725777960455],
    [-0.4794155057041976, 0.4699098521619943],
    [10.3324426408004211, 9.7123391106782755],
    [9.6818249672071008, 9.6834045098534336],
    [9.8042422429595373, 10.0247564316322375],
    [9.9319450186421161, 9.7912291401980411],
    [10.1118528947223787, 9.6394938606520419],
    [9.7921446485352188, 9.8663618432936921],
    [-0.0439300157829641, 10.2851759613930138],
    [-0.3003262178416403, 10.0142344384136113],
    [0.0924145688620425, 9.5464504127199969],
    [0.1075448519014384, 9.6705241236872919],
    [-0.4349484070147205, 10.4488855372533340],
    [0.4656320330745594, 10.3083973481164612],
    [9.8046137691733701, -0.4023278859936161],
    [10.1842330265121568, -0.0598475062603987],
    [9.6220382348447782, -0.0048230898887298],
    [9.5343885211152184, 0.4093204020787821],
    [9.7587799816000178, 0.1625222843539820],
    [9.8117110760894111, 0.0200680211778108],
  ]
}

///|
test "apcluster_neg_dist_mat" {
  // squared Euclidean distance: (3-0)^2 + (4-0)^2 = 25
  let x = [[0.0, 0.0], [3.0, 4.0]]
  let s = @src.ap_neg_dist_mat(x)
  assert_eq(s[0][0], 0.0)
  assert_eq(s[1][1], 0.0)
  assert_true(ap_approx(s[0][1], -25.0, 1.0e-12))
  assert_true(ap_approx(s[1][0], -25.0, 1.0e-12))
  // city block with r = 1
  let s1 = @src.ap_neg_dist_mat(x, r=1.0)
  assert_true(ap_approx(s1[0][1], -7.0, 1.0e-12))
}

///|
test "apcluster_exp_sim_mat" {
  let x = [[0.0, 0.0], [3.0, 4.0]]
  // default w = median of the only pairwise distance^r = 25 -> exp(-1)
  let s = @src.ap_exp_sim_mat(x)
  assert_true(ap_approx(s[0][1], 0.36787944117144233, 1.0e-12))
  assert_eq(s[0][0], 1.0)
  // explicit width
  let s2 = @src.ap_exp_sim_mat(x, w=Some(50.0))
  assert_true(ap_approx(s2[0][1], 0.6065306597126334, 1.0e-12))
}

///|
test "apcluster_lin_sim_mat" {
  let x = [[0.0, 0.0], [3.0, 4.0]]
  // default w = max pairwise distance^r = 25 -> 1 - 25/25 = 0
  let s = @src.ap_lin_sim_mat(x)
  assert_true(ap_approx(s[0][1], 0.0, 1.0e-12))
  assert_eq(s[0][0], 1.0)
  let s2 = @src.ap_lin_sim_mat(x, w=Some(50.0))
  assert_true(ap_approx(s2[0][1], 0.5, 1.0e-12))
}

///|
test "apcluster_cor_sim_mat" {
  // perfectly proportional rows (Pearson == 1), perfectly reversed (== -1)
  let x = [
    [1.0, 2.0, 3.0, 4.0, 5.0],
    [2.0, 4.0, 6.0, 8.0, 10.0],
    [5.0, 4.0, 3.0, 2.0, 1.0],
  ]
  let sp = @src.ap_cor_sim_mat(x)
  assert_true(ap_approx(sp[0][1], 1.0, 1.0e-12))
  assert_true(ap_approx(sp[0][2], -1.0, 1.0e-12))
  assert_eq(sp[0][0], 1.0)
  // Spearman of monotone rows is also 1; reversed-monotone -1
  let ss = @src.ap_cor_sim_mat(x, method="spearman")
  assert_true(ap_approx(ss[0][1], 1.0, 1.0e-12))
  assert_true(ap_approx(ss[1][2], -1.0, 1.0e-12))
  // Spearman with ties: tied ranks average (R ties.method="average")
  let t = [
    [1.0, 2.0, 2.0, 3.0],
    [10.0, 20.0, 20.0, 30.0],
  ]
  let st = @src.ap_cor_sim_mat(t, method="spearman")
  assert_true(ap_approx(st[0][1], 1.0, 1.0e-12))
  // constant rows cannot be correlated (cor() would be NA)
  let c = [[1.0, 1.0, 1.0], [1.0, 2.0, 3.0]]
  let raised = try {
    let _ = @src.ap_cor_sim_mat(c)
    false
  } catch {
    ApError(_) => true
  }
  assert_true(raised)
  // unknown correlation method
  let raised2 = try {
    let _ = @src.ap_cor_sim_mat([[1.0, 2.0], [3.0, 4.0]], method="kendall")
    false
  } catch {
    ApError(_) => true
  }
  assert_true(raised2)
}

///|
test "apcluster_default_four_blobs" {
  let x = ap_blobs()
  let res = @src.apcluster(sim=@src.ap_neg_dist_mat(x))
  assert_true(res.converged())
  assert_eq(res.exemplars(), [0, 9, 12, 23])
  // every blob maps to a single exemplar
  let expected_labels = [
    0, 0, 0, 0, 0, 0,
    9, 9, 9, 9, 9, 9,
    12, 12, 12, 12, 12, 12,
    23, 23, 23, 23, 23, 23,
  ]
  assert_eq(res.labels(), expected_labels)
  assert_eq(res.num_clusters(), 4)
  // cluster sizes, exemplar always its own first member
  for members in res.clusters() {
    assert_eq(members.length(), 6)
  }
  // reference values from the numpy mirror
  assert_true(ap_approx(res.preferences()[0], -97.6990940249, 1.0e-6))
  assert_true(ap_approx(res.dpsim(), -3.8222613499, 1.0e-6))
  assert_true(ap_approx(res.expref(), -390.7963760995, 1.0e-6))
  assert_true(ap_approx(res.netsim(), -394.6186374494, 1.0e-6))
  // objective identity: netsim == dpsim + expref
  assert_true(
    ap_approx(res.netsim(), res.dpsim() + res.expref(), 1.0e-9),
  )
}

///|
test "apcluster_k_exact_and_prune" {
  let x = ap_blobs()
  let s = @src.ap_neg_dist_mat(x)
  // K = 4 at the median preference (first run already exact)
  let r4 = @src.apcluster_k(sim=s, k=4)
  assert_eq(r4.exemplars(), [0, 9, 12, 23])
  assert_true(ap_approx(r4.netsim(), -394.6186374494, 1.0e-6))
  // K = 5 reachable by raising the preference
  let r5 = @src.apcluster_k(sim=s, k=5)
  assert_eq(r5.num_clusters(), 5)
  assert_eq(r5.exemplars(), [0, 9, 12, 15, 23])
  assert_true(ap_approx(r5.netsim(), -6.6892086223, 1.0e-6))
  // K = 3 / K = 2 are unreachable by bisection (solution collapses), so
  // the closest run is pruned by merging closest exemplar pairs
  let r3 = @src.apcluster_k(sim=s, k=3)
  assert_eq(r3.num_clusters(), 3)
  assert_eq(r3.exemplars(), [0, 9, 12])
  let sizes3 = r3.clusters().map(fn(c : Array[Int]) -> Int { c.length() })
  assert_eq(sizes3, [6, 12, 6])
  assert_true(ap_approx(r3.netsim(), -869.8089505034, 1.0e-6))
  let r2 = @src.apcluster_k(sim=s, k=2)
  assert_eq(r2.num_clusters(), 2)
  assert_eq(r2.exemplars(), [0, 9])
  let sizes2 = r2.clusters().map(fn(c : Array[Int]) -> Int { c.length() })
  assert_eq(sizes2, [12, 12])
  assert_true(ap_approx(r2.netsim(), -1324.2021357240, 1.0e-6))
  // K above the number of points is rejected
  let raised = try {
    let _ = @src.apcluster_k(sim=s, k=25)
    false
  } catch {
    ApError(_) => true
  }
  assert_true(raised)
}

///|
test "apcluster_one_dim_two_groups" {
  // two tight groups on a line -> one exemplar each (the middle point)
  let y = [
    [0.0], [1.0], [2.0], [100.0], [101.0], [102.0],
  ]
  let res = @src.apcluster(sim=@src.ap_neg_dist_mat(y))
  assert_true(res.converged())
  assert_eq(res.exemplars(), [1, 4])
  assert_eq(res.labels(), [1, 1, 1, 4, 4, 4])
  assert_true(ap_approx(res.preferences()[0], -4804.0, 1.0e-6))
  assert_true(ap_approx(res.dpsim(), -4.0, 1.0e-6))
  assert_true(ap_approx(res.expref(), -9608.0, 1.0e-6))
  assert_true(ap_approx(res.netsim(), -9612.0, 1.0e-6))
}

///|
test "apcluster_extreme_preferences" {
  let y = [
    [0.0], [1.0], [2.0], [100.0], [101.0], [102.0],
  ]
  let s = @src.ap_neg_dist_mat(y)
  // preference above every similarity: every point becomes its own exemplar
  let hi = @src.apcluster(sim=s, p=Some(10.0))
  assert_eq(hi.num_clusters(), 6)
  for i in 0..<6 {
    assert_eq(hi.labels()[i], i)
  }
  // very low preference: a single cluster
  let lo = @src.apcluster(sim=s, p=Some(-100000.0))
  assert_eq(lo.num_clusters(), 1)
  assert_eq(lo.exemplars(), [3])
  assert_true(ap_approx(lo.expref(), -100000.0, 1.0e-6))
  assert_true(ap_approx(lo.netsim(), -129410.0, 1.0e-6))
  // per-point preference vector overrides the scalar one
  let prefs = Array::make(6, -100000.0)
  prefs[0] = 10.0
  let vec = @src.apcluster(sim=s, preferences=Some(prefs))
  assert_eq(vec.num_clusters(), 1)
}

///|
test "apcluster_input_validation" {
  // non-square similarity matrix
  let bad = [[0.0, -1.0], [-1.0, 0.0, 1.0]]
  let r1 = try {
    let _ = @src.apcluster(sim=bad)
    false
  } catch {
    ApError(_) => true
  }
  assert_true(r1)
  // damping out of range
  let ok = [[0.0, -1.0], [-1.0, 0.0]]
  let r2 = try {
    let _ = @src.apcluster(sim=ok, lam=0.4)
    false
  } catch {
    ApError(_) => true
  }
  assert_true(r2)
  // per-point preferences of wrong length
  let r3 = try {
    let _ = @src.apcluster(sim=ok, preferences=Some([0.0]))
    false
  } catch {
    ApError(_) => true
  }
  assert_true(r3)
  // empty data matrix for similarity construction
  let r4 = try {
    let _ = @src.ap_neg_dist_mat([])
    false
  } catch {
    ApError(_) => true
  }
  assert_true(r4)
}

///|
test "apcluster_noise_runs_reproducible" {
  @src.apcluster_seed(1234)
  let x = ap_blobs()
  let s = @src.ap_neg_dist_mat(x)
  let r1 = @src.apcluster(sim=s, nonoise=5)
  @src.apcluster_seed(1234)
  let r2 = @src.apcluster(sim=s, nonoise=5)
  // same seed -> identical runs (deterministic noise)
  assert_eq(r1.exemplars(), r2.exemplars())
  assert_eq(r1.labels(), r2.labels())
  assert_eq(r1.num_clusters(), r2.num_clusters())
  // four well-separated blobs survive the small perturbations
  assert_eq(r1.num_clusters(), 4)
}
