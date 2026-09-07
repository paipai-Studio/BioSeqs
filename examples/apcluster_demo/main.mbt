///|
/// Bioconductor apcluster Demo - Affinity Propagation Clustering.
///
/// Demonstrates the port of the Bioconductor `apcluster` package
/// (Frey & Dueck, Science 2007): four similarity functions
/// (negDistMat / expSimMat / linSimMat / corSimMat), the default median-
/// preference run, the fixed-cluster-count apclusterK bisection, the
/// preference controls, and the deterministic noise runs.

///|
/// Four tight 2D blobs of six points each (deterministic jitter).
fn demo_blobs() -> Array[Array[Double]] {
  [
    [-0.12545988, 0.45071431], [0.23199394, 0.09865848],
    [-0.34398136, -0.34400548], [-0.44191639, 0.36617615],
    [0.10111501, 0.20807258], [-0.47941551, 0.46990985],
    [10.33244264, 9.71233911], [9.68182497, 9.68340451],
    [9.80424224, 10.02475643], [9.93194502, 9.79122914],
    [10.11185289, 9.63949386], [9.79214465, 9.86636184],
    [-0.04393002, 10.28517596], [-0.30032622, 10.01423444],
    [0.09241457, 9.54645041], [0.10754485, 9.67052412],
    [-0.43494841, 10.44888554], [0.46563203, 10.30839735],
    [9.80461377, -0.40232789], [10.18423303, -0.05984751],
    [9.62203823, -0.00482309], [9.53438852, 0.40932040],
    [9.75877998, 0.16252228], [9.81171108, 0.02006802],
  ]
}

///|
fn cluster_sizes(res : @src.ApResult) -> Array[Int] {
  let out : Array[Int] = []
  for members in res.clusters() {
    out.push(members.length())
  }
  out
}

///|
fn print_result(title : String, res : @src.ApResult) -> Unit {
  println("--- \{title} ---")
  println("exemplars (0-based): \{res.exemplars()}")
  println("cluster sizes:       \{cluster_sizes(res)}")
  println("iters:               \{res.iters()}  converged: \{res.converged()}")
  println("netsim: \{res.netsim()} (dpsim=\{res.dpsim()}, expref=\{res.expref()})")
  println("")
}

///|
fn main raise {
  let x = demo_blobs()

  println("=== 1. negDistMat (negative squared Euclidean) ===")
  let s = @src.ap_neg_dist_mat(x)
  let r = @src.apcluster(sim=s)
  println("median preference: \{r.preferences()[0]}")
  print_result("default apcluster", r)

  println("=== 2. Other distance-based similarities ===")
  let re = @src.apcluster(sim=@src.ap_exp_sim_mat(x))
  print_result("expSimMat (Gaussian kernel)", re)
  let rl = @src.apcluster(sim=@src.ap_lin_sim_mat(x))
  print_result("linSimMat (linear)", rl)

  println("=== 3. apclusterK: target a fixed number of clusters ===")
  let k5 = @src.apcluster_k(sim=s, k=5)
  print_result("apclusterK(k=5)", k5)
  let k3 = @src.apcluster_k(sim=s, k=3)
  println("--- apclusterK(k=3): bisection collapses, prune merge ---")
  println("exemplars: \{k3.exemplars()}  sizes: \{cluster_sizes(k3)}")
  println("")

  println("=== 4. Preference controls cluster granularity ===")
  // high preference: every point wants to be an exemplar
  let hi = @src.apcluster(sim=s, p=Some(10.0))
  println("p=+10   -> clusters: \{hi.num_clusters()}")
  // very low preference: a single consensus exemplar
  let lo = @src.apcluster(sim=s, p=Some(-100000.0))
  println("p=-1e5  -> clusters: \{lo.num_clusters()} exemplar: \{lo.exemplars()}")
  println("")

  println("=== 5. corSimMat: co-expression style grouping ===")
  let z = [
    [1.0, 2.0, 3.0, 4.0, 5.0],
    [1.1, 1.9, 3.1, 3.9, 5.2],
    [5.0, 4.0, 3.0, 2.0, 1.0],
    [5.1, 3.8, 2.9, 2.1, 0.9],
  ]
  let rc = @src.apcluster(sim=@src.ap_cor_sim_mat(z))
  print_result("corSimMat (Pearson)", rc)

  println("=== 6. Noise runs (deterministic with apcluster_seed) ===")
  @src.apcluster_seed(2026)
  let rn = @src.apcluster(sim=s, nonoise=5)
  println("nonoise=5 -> clusters: \{rn.num_clusters()} exemplars: \{rn.exemplars()}")
}
