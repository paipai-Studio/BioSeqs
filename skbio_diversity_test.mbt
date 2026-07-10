///|
/// Tests for scikit-bio Diversity module.

///|
test "observed_count_basic" {
  let counts = [1.0, 2.0, 0.0, 3.0, 0.0, 4.0]
  assert_eq(observed_otus(counts), 4)
}

///|
test "observed_count_all_zero" {
  let counts = [0.0, 0.0, 0.0]
  assert_eq(observed_otus(counts), 0)
}

///|
test "shannon_basic" {
  let counts = [1.0, 1.0, 1.0, 1.0]
  let h = shannon(counts)
  assert_true((h - 1.3862943611198906).abs() < 0.0001)
}

///|
test "shannon_single_species" {
  let counts = [10.0, 0.0, 0.0]
  let h = shannon(counts)
  assert_true(h < 0.0001)
}

///|
test "shannon_empty" {
  let counts = Array::make(0, 0.0)
  let h = shannon(counts)
  assert_true(h == 0.0)
}

///|
test "shannon_base2" {
  let counts = [1.0, 1.0, 1.0, 1.0]
  let h = shannon_base(counts, 2.0)
  assert_true((h - 2.0).abs() < 0.0001)
}

///|
test "simpson_basic" {
  let counts = [1.0, 1.0, 1.0, 1.0]
  let d = simpson(counts)
  // D = 1 - sum(p_i^2) = 1 - 4*(0.25^2) = 1 - 0.25 = 0.75
  assert_true((d - 0.75).abs() < 0.0001)
}

///|
test "simpson_single_species" {
  let counts = [10.0]
  let d = simpson(counts)
  assert_true(d.abs() < 0.0001)
}

///|
test "inv_simpson_basic" {
  let counts = [1.0, 1.0, 1.0, 1.0]
  let inv_d = inv_simpson(counts)
  // 1/D = 1 / 0.25 = 4.0
  assert_true((inv_d - 4.0).abs() < 0.0001)
}

///|
test "pielou_e_basic" {
  let counts = [1.0, 1.0, 1.0, 1.0]
  let e = pielou_e(counts)
  // Maximum evenness = 1.0
  assert_true((e - 1.0).abs() < 0.0001)
}

///|
test "pielou_e_single" {
  let counts = [10.0]
  let e = pielou_e(counts)
  assert_true(e == 0.0)
}

///|
test "simpson_e_basic" {
  let counts = [1.0, 1.0, 1.0, 1.0]
  let e = simpson_e(counts)
  // 1/D / S = 4 / 4 = 1.0
  assert_true((e - 1.0).abs() < 0.0001)
}

///|
test "heip_e_basic" {
  let counts = [1.0, 1.0, 1.0, 1.0]
  let e = heip_e(counts)
  assert_true((e - 1.0).abs() < 0.0001)
}

///|
test "chao1_simple" {
  let counts = [2.0, 2.0, 2.0, 1.0, 1.0]
  // Bias-corrected Chao1 (scikit-bio default)
  // F1 = 2, F2 = 3, S_obs = 5
  // S_chao1 = 5 + 2*(2-1) / (2*(3+1)) = 5 + 2/8 = 5.25
  let s = chao1(counts)
  assert_true((s - 5.25).abs() < 0.001)
}

///|
test "chao1_no_singletons" {
  let counts = [2.0, 2.0, 2.0]
  let s = chao1(counts)
  // No singletons: S_chao1 = S_obs = 3
  assert_true((s - 3.0).abs() < 0.0001)
}

///|
test "chao1_no_doubletons" {
  let counts = [1.0, 1.0, 1.0]
  let s = chao1(counts)
  // With F1=3, F2=0: bias-corrected Chao1
  // S_chao1 = 3 + 3*2/(2*1) = 3 + 3 = 6
  assert_true(s > 3.0)
}

///|
test "ace_basic" {
  let counts = [1.0, 1.0, 2.0, 2.0, 3.0, 5.0, 10.0]
  let s = ace_default(counts)
  assert_true(s >= 7.0)
}

///|
test "berger_parker_basic" {
  let counts = [1.0, 2.0, 5.0, 2.0]
  // N_max = 5, N = 10, d = 5/10 = 0.5
  let d = berger_parker(counts)
  assert_true((d - 0.5).abs() < 0.0001)
}

///|
test "margalef_basic" {
  let counts = [1.0, 2.0, 3.0, 4.0]
  let d = margalef(counts)
  let expected = 3.0 / @math.ln(10.0)
  assert_true((d - expected).abs() < 0.0001)
}

///|
test "menhinick_basic" {
  let counts = [1.0, 1.0, 1.0, 1.0]
  // S = 4, N = 4, D_Mn = 4 / sqrt(4) = 2.0
  let d = menhinick(counts)
  assert_true((d - 2.0).abs() < 0.0001)
}

///|
test "fisher_alpha_basic" {
  let counts = [10.0, 8.0, 5.0, 3.0, 2.0, 1.0, 1.0]
  let a = fisher_alpha(counts)
  assert_true(a > 0.0)
}

///|
test "braycurtis_identical" {
  let x = [1.0, 2.0, 3.0]
  let y = [1.0, 2.0, 3.0]
  let d = braycurtis(x, y)
  assert_true(d.abs() < 0.0001)
}

///|
test "braycurtis_maximal" {
  let x = [1.0, 0.0]
  let y = [0.0, 1.0]
  let d = braycurtis(x, y)
  assert_true((d - 1.0).abs() < 0.0001)
}

///|
test "braycurtis_basic" {
  let x = [1.0, 2.0, 3.0]
  let y = [2.0, 3.0, 4.0]
  let d = braycurtis(x, y)
  // sum(|x_i - y_i|) = |1-2|+|2-3|+|3-4| = 3
  // sum(x_i + y_i) = 3 + 5 + 7 = 15
  // BC = 3/15 = 0.2
  assert_true((d - 0.2).abs() < 0.0001)
}

///|
test "jaccard_identical" {
  let x = [1.0, 2.0, 3.0]
  let y = [1.0, 2.0, 3.0]
  let d = jaccard(x, y)
  assert_true(d.abs() < 0.0001)
}

///|
test "jaccard_maximal" {
  let x = [1.0, 0.0]
  let y = [0.0, 1.0]
  let d = jaccard(x, y)
  assert_true((d - 1.0).abs() < 0.0001)
}

///|
test "euclidean_basic" {
  let x = [0.0, 0.0]
  let y = [3.0, 4.0]
  let d = euclidean(x, y)
  assert_true((d - 5.0).abs() < 0.0001)
}

///|
test "manhattan_basic" {
  let x = [1.0, 2.0, 3.0]
  let y = [2.0, 4.0, 6.0]
  let d = manhattan(x, y)
  // 1 + 2 + 3 = 6
  assert_true((d - 6.0).abs() < 0.0001)
}

///|
test "chebyshev_basic" {
  let x = [1.0, 2.0, 3.0]
  let y = [2.0, 5.0, 4.0]
  let d = chebyshev(x, y)
  // max(|1-2|, |2-5|, |3-4|) = max(1, 3, 1) = 3
  assert_true((d - 3.0).abs() < 0.0001)
}

///|
test "canberra_basic" {
  let x = [1.0, 2.0, 3.0]
  let y = [2.0, 4.0, 6.0]
  let d = canberra(x, y)
  // 1/3 + 2/6 + 3/9 = 1/3 + 1/3 + 1/3 = 1.0
  assert_true((d - 1.0).abs() < 0.0001)
}

///|
test "cosine_identical" {
  let x = [1.0, 2.0, 3.0]
  let y = [2.0, 4.0, 6.0]
  let d = cosine_distance(x, y)
  assert_true(d.abs() < 0.0001)
}

///|
test "cosine_orthogonal" {
  let x = [1.0, 0.0]
  let y = [0.0, 1.0]
  let d = cosine_distance(x, y)
  assert_true((d - 1.0).abs() < 0.0001)
}

///|
test "faith_pd_simple" {
  let newick = "((A:1.0,B:1.0):2.0,C:3.0):0.0;"
  let tree = TreeNode::from_newick(newick)
  let counts = [1.0, 1.0, 1.0]
  let taxa = ["A", "B", "C"]
  let pd = faith_pd(counts, taxa, tree)
  // Total branch lengths for all 3:
  // A: 1.0 + 2.0 = 3.0 (from tip to root)
  // B: 1.0 + 2.0 = 3.0
  // C: 3.0
  // But shared branches counted once:
  // Root has 0 length
  // Branch to (A,B): 2.0
  // A: 1.0, B: 1.0, C: 3.0
  // PD = 2.0 + 1.0 + 1.0 + 3.0 = 7.0
  assert_true((pd - 7.0).abs() < 0.0001)
}

///|
test "faith_pd_subset" {
  let newick = "((A:1.0,B:1.0):2.0,C:3.0):0.0;"
  let tree = TreeNode::from_newick(newick)
  let counts = [1.0, 0.0, 1.0]
  let taxa = ["A", "B", "C"]
  let pd = faith_pd(counts, taxa, tree)
  // Only A and C present:
  // A: 1.0 + 2.0 = 3.0 (root to A)
  // C: 3.0 (root to C)
  // Shared root (0 length)
  // PD = 2.0 + 1.0 + 3.0 = 6.0
  assert_true((pd - 6.0).abs() < 0.0001)
}

///|
test "faith_pd_single_taxon" {
  let newick = "((A:1.0,B:1.0):2.0,C:3.0):0.0;"
  let tree = TreeNode::from_newick(newick)
  let counts = [1.0, 0.0, 0.0]
  let taxa = ["A", "B", "C"]
  let pd = faith_pd(counts, taxa, tree)
  // Only A: 1.0 + 2.0 = 3.0
  assert_true((pd - 3.0).abs() < 0.0001)
}

///|
test "gini_simpson_basic" {
  let counts = [1.0, 1.0, 1.0, 1.0]
  let gs = gini_simpson(counts)
  // 1 - 0.25 = 0.75
  assert_true((gs - 0.75).abs() < 0.0001)
}

///|
test "dominance_basic" {
  let counts = [1.0, 1.0, 1.0, 1.0]
  let d = dominance(counts)
  // sum(p_i^2) = 0.25
  assert_true((d - 0.25).abs() < 0.0001)
}

///|
test "compare_with_skbio_alpha" {
  // Test data matching Python comparison script: [1, 2, 3, 4, 5, 0, 0, 1, 2, 3]
  let counts = [1.0, 2.0, 3.0, 4.0, 5.0, 0.0, 0.0, 1.0, 2.0, 3.0]

  // scikit-bio values
  assert_eq(observed_otus(counts), 8)
  assert_true((shannon(counts) - 1.9513496063).abs() < 0.0001)
  assert_true((simpson(counts) - 0.8435374150).abs() < 0.0001)
  assert_true((chao1(counts) - 8.3333333333).abs() < 0.0001)
  assert_true((berger_parker(counts) - 0.2380952381).abs() < 0.0001)
  assert_true((dominance(counts) - 0.1564625850).abs() < 0.0001)
}

///|
test "compare_with_skbio_beta" {
  // Test data matching Python comparison script
  let x = [1.0, 2.0, 3.0, 4.0, 5.0, 0.0, 0.0, 1.0, 2.0, 3.0]
  let y = [3.0, 1.0, 4.0, 1.0, 5.0, 9.0, 2.0, 6.0, 5.0, 3.0]

  // scikit-bio values
  assert_true((braycurtis(x, y) - 0.4333333333).abs() < 0.0001)
  assert_true((jaccard(x, y) - 0.2000000000).abs() < 0.0001)
  assert_true((euclidean(x, y) - 11.5758369028).abs() < 0.001)
  assert_true((manhattan(x, y) - 26.0000000000).abs() < 0.0001)
  assert_true((chebyshev(x, y) - 9.0000000000).abs() < 0.0001)
  assert_true((canberra(x, y) - 4.719048).abs() < 0.001)
  assert_true((cosine_distance(x, y) - 0.405915).abs() < 0.001)
}

///|
test "compare_with_skbio_faith_pd" {
  let newick = "((A:1.0,B:1.0):2.0,C:3.0):0.0;"
  let tree = TreeNode::from_newick(newick)
  let taxa = ["A", "B", "C"]

  // All taxa present
  let counts_all = [1.0, 1.0, 1.0]
  let pd_all = faith_pd(counts_all, taxa, tree)
  assert_true((pd_all - 7.0).abs() < 0.0001)

  // Only A and C
  let counts_ac = [1.0, 0.0, 1.0]
  let pd_ac = faith_pd(counts_ac, taxa, tree)
  assert_true((pd_ac - 6.0).abs() < 0.0001)
}
