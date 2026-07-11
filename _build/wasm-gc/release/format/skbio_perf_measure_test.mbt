///|
/// Performance measurement tests using bench clock.
/// Uses @bench.monotonic_clock_start/end for timing.

///|
/// Generate a balanced binary tree Newick string with 2^depth tips.
fn perf_gen_balanced_newick(depth : Int) -> String {
  let buf = StringBuilder::new()
  let counter = [0]
  perf_gen_balanced_rec(buf, depth, counter)
  buf.write_char(';')
  buf.to_string()
}

///|
fn perf_gen_balanced_rec(
  buf : StringBuilder,
  depth : Int,
  counter : Array[Int],
) -> Unit {
  if depth == 0 {
    buf.write_string("tip")
    buf.write_string(counter[0].to_string())
    counter[0] = counter[0] + 1
    buf.write_char(':')
    buf.write_string("0.1")
    return
  }
  buf.write_char('(')
  perf_gen_balanced_rec(buf, depth - 1, counter)
  buf.write_char(',')
  perf_gen_balanced_rec(buf, depth - 1, counter)
  buf.write_char(')')
  buf.write_char(':')
  buf.write_string("0.1")
}

///|
fn perf_gen_caterpillar_newick(n : Int) -> String {
  if n < 2 {
    return "tip0:0.1;"
  }
  let buf = StringBuilder::new()
  for i = 0; i < n - 1; i = i + 1 {
    buf.write_char('(')
  }
  buf.write_string("tip0:0.1,tip1:0.1")
  for i = 2; i < n; i = i + 1 {
    buf.write_char(')')
    buf.write_string("int")
    buf.write_string((i - 2).to_string())
    buf.write_string(":0.1,tip")
    buf.write_string(i.to_string())
    buf.write_string(":1.0")
  }
  buf.write_char(')')
  buf.write_char(';')
  buf.to_string()
}

///|
/// Generate large count vector for testing
fn perf_gen_counts(n : Int, seed : Int) -> Array[Double] {
  let counts = Array::make(n, 0.0)
  let mut s = seed
  for i in 0..<n {
    s = (s * 1103515245 + 12345) & 0x7FFFFFFF
    counts[i] = (s % 100).to_double()
  }
  counts
}

///|
/// Generate a balanced tree with named tips and matching taxa/counts arrays.
/// Returns (tree, taxa, counts) where ~half the taxa are present.
fn perf_gen_faith_data(
  depth : Int,
) -> (TreeNode, Array[String], Array[Double]) raise ParseError {
  let newick = perf_gen_balanced_newick(depth)
  let tree = TreeNode::from_newick(newick)
  let n_tips = tree.count_tips()
  let taxa = Array::new()
  let counts = Array::new()
  let mut s = 7
  for i in 0..<n_tips {
    taxa.push("tip" + i.to_string())
    s = (s * 1103515245 + 12345) & 0x7FFFFFFF
    // ~50% present (odd values -> present), with varied abundances
    if s % 2 == 1 {
      counts.push((s % 50 + 1).to_double())
    } else {
      counts.push(0.0)
    }
  }
  (tree, taxa, counts)
}

///|

///|
/// Report time in microseconds per call. total is in nanoseconds.
fn report_us(name : String, total : Double, iters : Int) -> Unit {
  let per_call = total / iters.to_double() / 1000.0
  println(
    name +
    " x" +
    iters.to_string() +
    ": " +
    (total / 1000000.0).to_string() +
    "ms (" +
    per_call.to_string() +
    " us/call)",
  )
}

///|
/// Report time in milliseconds per call. total is in nanoseconds.
fn report_ms(name : String, total : Double, iters : Int) -> Unit {
  let per_call = total / iters.to_double() / 1000000.0
  println(
    name +
    " x" +
    iters.to_string() +
    ": " +
    (total / 1000000.0).to_string() +
    "ms (" +
    per_call.to_string() +
    " ms/call)",
  )
}

///|
test "perf_measure_diversity_alpha" {
  let counts = perf_gen_counts(1000, 42)
  let iters = 1000

  let t0 = @bench.monotonic_clock_start()
  for i = 0; i < iters; i = i + 1 {
    let _ = shannon(counts)
  }
  report_us("  shannon", @bench.monotonic_clock_end(t0), iters)

  let t0 = @bench.monotonic_clock_start()
  for i = 0; i < iters; i = i + 1 {
    let _ = simpson(counts)
  }
  report_us("  simpson", @bench.monotonic_clock_end(t0), iters)

  let t0 = @bench.monotonic_clock_start()
  for i = 0; i < iters; i = i + 1 {
    let _ = chao1(counts)
  }
  report_us("  chao1", @bench.monotonic_clock_end(t0), iters)

  let t0 = @bench.monotonic_clock_start()
  for i = 0; i < iters; i = i + 1 {
    let _ = ace_default(counts)
  }
  report_us("  ace", @bench.monotonic_clock_end(t0), iters)

  let t0 = @bench.monotonic_clock_start()
  for i = 0; i < iters; i = i + 1 {
    let _ = observed_otus(counts)
  }
  report_us("  observed_otus", @bench.monotonic_clock_end(t0), iters)

  // faith_pd on a 256-tip tree (x100, smaller iter count due to heavier cost)
  let (faith_tree, faith_taxa, faith_counts) = perf_gen_faith_data(8)
  let faith_iters = 100
  let t0 = @bench.monotonic_clock_start()
  for i = 0; i < faith_iters; i = i + 1 {
    let _ = faith_pd(faith_counts, faith_taxa, faith_tree)
  }
  report_us("  faith_pd_256", @bench.monotonic_clock_end(t0), faith_iters)

  assert_true(true)
}

///|
test "perf_measure_diversity_beta" {
  let x = perf_gen_counts(1000, 42)
  let y = perf_gen_counts(1000, 137)
  let iters = 1000

  let t0 = @bench.monotonic_clock_start()
  for i = 0; i < iters; i = i + 1 {
    let _ = braycurtis(x, y)
  }
  report_us("  braycurtis", @bench.monotonic_clock_end(t0), iters)

  let t0 = @bench.monotonic_clock_start()
  for i = 0; i < iters; i = i + 1 {
    let _ = jaccard(x, y)
  }
  report_us("  jaccard", @bench.monotonic_clock_end(t0), iters)

  let t0 = @bench.monotonic_clock_start()
  for i = 0; i < iters; i = i + 1 {
    let _ = euclidean(x, y)
  }
  report_us("  euclidean", @bench.monotonic_clock_end(t0), iters)

  let t0 = @bench.monotonic_clock_start()
  for i = 0; i < iters; i = i + 1 {
    let _ = canberra(x, y)
  }
  report_us("  canberra", @bench.monotonic_clock_end(t0), iters)

  assert_true(true)
}

///|
test "perf_measure_tree_ops" {
  let newick_1024 = perf_gen_balanced_newick(10)
  let newick_4096 = perf_gen_balanced_newick(12)
  let cat_500 = perf_gen_caterpillar_newick(500)

  // Parse 4096-tip tree x10
  let t0 = @bench.monotonic_clock_start()
  for i = 0; i < 10; i = i + 1 {
    let _ = TreeNode::from_newick(newick_4096)
  }
  report_ms("  parse_4096", @bench.monotonic_clock_end(t0), 10)

  // to_newick x10
  let tree_4096 = TreeNode::from_newick(newick_4096)
  let t0 = @bench.monotonic_clock_start()
  for i = 0; i < 10; i = i + 1 {
    let _ = tree_4096.to_newick()
  }
  report_ms("  to_newick_4096", @bench.monotonic_clock_end(t0), 10)

  // Preorder x10
  let t0 = @bench.monotonic_clock_start()
  for i = 0; i < 10; i = i + 1 {
    let _ = tree_4096.preorder()
  }
  report_ms("  preorder_4096", @bench.monotonic_clock_end(t0), 10)

  // Levelorder x10
  let t0 = @bench.monotonic_clock_start()
  for i = 0; i < 10; i = i + 1 {
    let _ = tree_4096.levelorder()
  }
  report_ms("  levelorder_4096", @bench.monotonic_clock_end(t0), 10)

  // LCA on 1024-tip tree x100
  let tree_1024 = TreeNode::from_newick(newick_1024)
  let names = ["tip100", "tip200"]
  let t0 = @bench.monotonic_clock_start()
  for i = 0; i < 100; i = i + 1 {
    let _ = tree_1024.lca(names)
  }
  report_us("  lca_1024", @bench.monotonic_clock_end(t0), 100)

  // Distance on 1024-tip tree x100
  let t0 = @bench.monotonic_clock_start()
  for i = 0; i < 100; i = i + 1 {
    let _ = tree_1024.distance("tip100", "tip200")
  }
  report_us("  distance_1024", @bench.monotonic_clock_end(t0), 100)

  // Find on 4096-tip tree x100
  let t0 = @bench.monotonic_clock_start()
  for i = 0; i < 100; i = i + 1 {
    let _ = tree_4096.find("tip4095")
  }
  report_us("  find_4096", @bench.monotonic_clock_end(t0), 100)

  // Caterpillar LCA x100
  let cat_tree = TreeNode::from_newick(cat_500)
  let cat_names = ["tip0", "tip498"]
  let t0 = @bench.monotonic_clock_start()
  for i = 0; i < 100; i = i + 1 {
    let _ = cat_tree.lca(cat_names)
  }
  report_us("  cat_lca_500", @bench.monotonic_clock_end(t0), 100)

  assert_true(true)
}

///|
/// Generate a test distance matrix of size n.
fn perf_gen_distance_matrix(n : Int, seed : Int) -> DistanceMatrix {
  let data = Array::make(n * n, 0.0)
  let mut s = seed
  let mut i = 0
  while i < n {
    let mut j = i + 1
    while j < n {
      s = (s * 1103515245 + 12345) & 0x7FFFFFFF
      let v = (s % 100).to_double()
      data[i * n + j] = v
      data[j * n + i] = v
      j = j + 1
    }
    i = i + 1
  }
  let ids : Array[String] = Array::make(n, "")
  let mut i = 0
  while i < n {
    ids[i] = "s" + i.to_string()
    i = i + 1
  }
  DistanceMatrix::from_flat(data, n, ids)
}

///|
/// Generate groups array with k groups of n samples.
fn perf_gen_groups(n : Int, k : Int) -> Array[Int] {
  let groups = Array::make(n, 0)
  let mut i = 0
  while i < n {
    groups[i] = i % k
    i = i + 1
  }
  groups
}

///|
test "perf_measure_stats" {
  // P-value adjustment x1000
  let pvals = [
    0.001, 0.005, 0.01, 0.02, 0.03, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45,
    0.5, 0.55, 0.6, 0.65, 0.7, 0.75,
  ]
  let iters = 1000

  let t0 = @bench.monotonic_clock_start()
  for i = 0; i < iters; i = i + 1 {
    let _ = p_adjust_bonferroni(pvals)
  }
  report_us("  p_adjust_bonferroni", @bench.monotonic_clock_end(t0), iters)

  let t0 = @bench.monotonic_clock_start()
  for i = 0; i < iters; i = i + 1 {
    let _ = p_adjust_holm(pvals)
  }
  report_us("  p_adjust_holm", @bench.monotonic_clock_end(t0), iters)

  let t0 = @bench.monotonic_clock_start()
  for i = 0; i < iters; i = i + 1 {
    let _ = p_adjust_fdr_bh(pvals)
  }
  report_us("  p_adjust_fdr_bh", @bench.monotonic_clock_end(t0), iters)

  // Subsample counts x1000
  let counts = perf_gen_counts(500, 42)
  let t0 = @bench.monotonic_clock_start()
  for i = 0; i < iters; i = i + 1 {
    let _ = subsample_counts(counts, 200, true, i)
  }
  report_us("  subsample_replace", @bench.monotonic_clock_end(t0), iters)

  let t0 = @bench.monotonic_clock_start()
  for i = 0; i < iters; i = i + 1 {
    let _ = subsample_counts(counts, 200, false, i)
  }
  report_us("  subsample_no_replace", @bench.monotonic_clock_end(t0), iters)

  // Permutation test x100
  let s1 = perf_gen_counts(50, 7)
  let s2 = perf_gen_counts(50, 13)
  let perm_iters = 100
  let t0 = @bench.monotonic_clock_start()
  for i = 0; i < perm_iters; i = i + 1 {
    let _ = permutation_test(s1, s2, diff_of_means, 500, i)
  }
  report_us(
    "  permutation_test_500",
    @bench.monotonic_clock_end(t0),
    perm_iters,
  )

  // ANOSIM x100 on 20-sample matrix
  let dm20 = perf_gen_distance_matrix(20, 42)
  let groups20 = perf_gen_groups(20, 2)
  let t0 = @bench.monotonic_clock_start()
  for i = 0; i < perm_iters; i = i + 1 {
    let _ = anosim(dm20, groups20, 200, i)
  }
  report_us("  anosim_20x200", @bench.monotonic_clock_end(t0), perm_iters)

  // PERMANOVA x100 on 20-sample matrix
  let t0 = @bench.monotonic_clock_start()
  for i = 0; i < perm_iters; i = i + 1 {
    let _ = permanova(dm20, groups20, 200, i)
  }
  report_us("  permanova_20x200", @bench.monotonic_clock_end(t0), perm_iters)

  // Confidence bound x1000
  let t0 = @bench.monotonic_clock_start()
  for i = 0; i < iters; i = i + 1 {
    let _ = confidence_bound(5.0, 2.0, 100, 0.95)
  }
  report_us("  confidence_bound", @bench.monotonic_clock_end(t0), iters)

  assert_true(true)
}
