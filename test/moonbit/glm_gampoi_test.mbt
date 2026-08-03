///|
/// Tests for glmGamPoi differential expression module.

// ---------------------------------------------------------------------------
// Data structure creation and accessors
// ---------------------------------------------------------------------------

test "glm_design_creation" {
  let d = @src.GlmDesign::new(
    ["s1", "s2", "s3"],
    ["Intercept", "Treatment"],
    [[1.0, 0.0], [1.0, 0.0], [1.0, 1.0]],
  )
  assert_eq(d.n_samples(), 3)
  assert_eq(d.n_coefs(), 2)
  assert_eq(d.get(0, 0), 1.0)
  assert_eq(d.get(0, 1), 0.0)
  assert_eq(d.get(2, 0), 1.0)
  assert_eq(d.get(2, 1), 1.0)
}

test "glm_design_intercept_only" {
  let d = @src.GlmDesign::new(
    ["s1", "s2", "s3", "s4"],
    ["Intercept"],
    [[1.0], [1.0], [1.0], [1.0]],
  )
  assert_eq(d.n_samples(), 4)
  assert_eq(d.n_coefs(), 1)
  assert_eq(d.get(0, 0), 1.0)
  assert_eq(d.get(3, 0), 1.0)
}

test "glm_fit_creation" {
  let f = @src.GlmFit::new(
    "GeneA",
    [2.5, 1.2],
    [0.3, 0.4],
    0.15,
    10,
    8.5,
  )
  assert_eq(f.gene, "GeneA")
  assert_eq(f.coefficients.length(), 2)
  assert_eq(f.coefficients[0], 2.5)
  assert_eq(f.coefficients[1], 1.2)
  assert_eq(f.std_errors.length(), 2)
  assert_eq(f.std_errors[0], 0.3)
  assert_eq(f.dispersion, 0.15)
  assert_eq(f.df_residual, 10)
  assert_eq(f.df_dispersion, 8.5)
}

test "glm_fit_single_coef" {
  let f = @src.GlmFit::new(
    "GeneX",
    [3.0],
    [0.5],
    0.2,
    5,
    5.0,
  )
  assert_eq(f.coefficients.length(), 1)
  assert_eq(f.coefficients[0], 3.0)
  assert_eq(f.std_errors[0], 0.5)
}

test "glm_test_result_creation" {
  let r = @src.GlmTestResult::new(
    "GeneA",
    "Treatment",
    1.5,
    0.3,
    5.0,
    0.001,
  )
  assert_eq(r.estimate(), 1.5)
  assert_eq(r.p_value(), 0.001)
  assert_eq(r.adj_p_value(), 0.001)
}

test "glm_test_result_not_significant" {
  let r = @src.GlmTestResult::new(
    "GeneB",
    "Treatment",
    0.1,
    0.5,
    0.2,
    0.8,
  )
  assert_false(r.is_significant)
  assert_eq(r.estimate(), 0.1)
  assert_eq(r.p_value(), 0.8)
}

test "glm_test_result_to_string" {
  let r = @src.GlmTestResult::new(
    "GeneA",
    "Treatment",
    1.5,
    0.3,
    5.0,
    0.001,
  )
  let s = r.to_string()
  assert_true(s.contains("GeneA"))
  assert_true(s.contains("Treatment"))
  assert_true(s.contains("est="))
  assert_true(s.contains("p="))
}

test "glm_test_result_significant_marker" {
  let r = @src.GlmTestResult::new(
    "GeneC",
    "Treatment",
    3.0,
    0.2,
    15.0,
    1.0e-10,
  )
  let s = r.to_string()
  assert_true(s.contains("*"))
}

test "glm_pseudobulk_creation" {
  let pb = @src.GlmPseudobulk::new(
    ["GeneA", "GeneB"],
    ["Control", "Treated"],
    [["s1", "s2"], ["s3", "s4"]],
    [[10.0, 50.0], [20.0, 18.0]],
  )
  assert_eq(pb.n_genes(), 2)
  assert_eq(pb.n_groups(), 2)
}

test "glm_pseudobulk_single_gene" {
  let pb = @src.GlmPseudobulk::new(
    ["GeneX"],
    ["GroupA"],
    [["s1"]],
    [[42.0]],
  )
  assert_eq(pb.n_genes(), 1)
  assert_eq(pb.n_groups(), 1)
}

// ---------------------------------------------------------------------------
// Utility functions
// ---------------------------------------------------------------------------

test "glm_sum_basic" {
  let s = @src.glm_sum([1.0, 2.0, 3.0, 4.0])
  assert_eq(s, 10.0)
}

test "glm_sum_single" {
  assert_eq(@src.glm_sum([5.0]), 5.0)
}

test "glm_sum_empty" {
  assert_eq(@src.glm_sum([]), 0.0)
}

test "glm_sum_negative" {
  assert_eq(@src.glm_sum([-1.0, -2.0, 3.0]), 0.0)
}

test "glm_mean_basic" {
  let m = @src.glm_mean([2.0, 4.0, 6.0])
  assert_eq(m, 4.0)
}

test "glm_mean_single" {
  assert_eq(@src.glm_mean([7.0]), 7.0)
}

test "glm_mean_empty" {
  assert_eq(@src.glm_mean([]), 0.0)
}

test "glm_mean_zeros" {
  assert_eq(@src.glm_mean([0.0, 0.0, 0.0]), 0.0)
}

test "glm_median_odd" {
  let m = @src.glm_median([3.0, 1.0, 2.0])
  assert_eq(m, 2.0)
}

test "glm_median_even" {
  let m = @src.glm_median([1.0, 3.0, 2.0, 4.0])
  assert_eq(m, 2.5)
}

test "glm_median_single" {
  assert_eq(@src.glm_median([5.0]), 5.0)
}

test "glm_median_empty" {
  assert_eq(@src.glm_median([]), 0.0)
}

test "glm_median_already_sorted" {
  let m = @src.glm_median([1.0, 2.0, 3.0, 4.0, 5.0])
  assert_eq(m, 3.0)
}

test "glm_min_positive_basic" {
  let v = @src.glm_min_positive([0.0, 3.0, 2.0, 5.0])
  assert_eq(v, 2.0)
}

test "glm_min_positive_no_positive" {
  let v = @src.glm_min_positive([0.0, -1.0, -5.0])
  assert_eq(v, 0.0)
}

test "glm_min_positive_all_positive" {
  let v = @src.glm_min_positive([10.0, 3.0, 7.0])
  assert_eq(v, 3.0)
}

test "glm_min_positive_single" {
  let v = @src.glm_min_positive([42.0])
  assert_eq(v, 42.0)
}

test "glm_min_positive_empty" {
  let v = @src.glm_min_positive([])
  assert_eq(v, 0.0)
}

// ---------------------------------------------------------------------------
// Linear algebra: solve linear system
// ---------------------------------------------------------------------------

test "glm_solve_linear_system_2x2" {
  // 2x + 3y = 7
  //  x +  y = 3
  // Solution: x = 2, y = 1
  let a = [[2.0, 3.0], [1.0, 1.0]]
  let b = [7.0, 3.0]
  let x = @src.glm_solve_linear_system(a, b, 2)
  assert_true((x[0] - 2.0).abs() < 1.0e-6)
  assert_true((x[1] - 1.0).abs() < 1.0e-6)
}

test "glm_solve_linear_system_3x3" {
  //  x +  y +  z = 6
  // 2x +  y -  z = 1
  //  x -  y +  z = 2
  // Solution: x = 1, y = 2, z = 3
  let a = [[1.0, 1.0, 1.0], [2.0, 1.0, -1.0], [1.0, -1.0, 1.0]]
  let b = [6.0, 1.0, 2.0]
  let x = @src.glm_solve_linear_system(a, b, 3)
  assert_true((x[0] - 1.0).abs() < 1.0e-6)
  assert_true((x[1] - 2.0).abs() < 1.0e-6)
  assert_true((x[2] - 3.0).abs() < 1.0e-6)
}

test "glm_solve_linear_system_identity" {
  // Identity * x = b => x = b
  let a = [[1.0, 0.0], [0.0, 1.0]]
  let b = [5.0, -3.0]
  let x = @src.glm_solve_linear_system(a, b, 2)
  assert_true((x[0] - 5.0).abs() < 1.0e-6)
  assert_true((x[1] - (-3.0)).abs() < 1.0e-6)
}

test "glm_solve_linear_system_1x1" {
  let a = [[3.0]]
  let b = [12.0]
  let x = @src.glm_solve_linear_system(a, b, 1)
  assert_true((x[0] - 4.0).abs() < 1.0e-6)
}

test "glm_solve_linear_system_diagonal" {
  // Diagonal 3x3 system
  // 3x = 9, 2y = 8, 4z = 16
  let a = [[3.0, 0.0, 0.0], [0.0, 2.0, 0.0], [0.0, 0.0, 4.0]]
  let b = [9.0, 8.0, 16.0]
  let x = @src.glm_solve_linear_system(a, b, 3)
  assert_true((x[0] - 3.0).abs() < 1.0e-6)
  assert_true((x[1] - 4.0).abs() < 1.0e-6)
  assert_true((x[2] - 4.0).abs() < 1.0e-6)
}

// ---------------------------------------------------------------------------
// Linear algebra: matrix inverse
// ---------------------------------------------------------------------------

test "glm_invert_matrix_2x2" {
  // A = [[4, 7], [2, 6]]
  // A^{-1} = [[0.6, -0.7], [-0.2, 0.4]]
  let a = [[4.0, 7.0], [2.0, 6.0]]
  let inv = @src.glm_invert_matrix(a, 2)
  assert_true((inv[0][0] - 0.6).abs() < 1.0e-6)
  assert_true((inv[0][1] - (-0.7)).abs() < 1.0e-6)
  assert_true((inv[1][0] - (-0.2)).abs() < 1.0e-6)
  assert_true((inv[1][1] - 0.4).abs() < 1.0e-6)
}

test "glm_invert_matrix_identity_3x3" {
  let a = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]
  let inv = @src.glm_invert_matrix(a, 3)
  for i in 0..<3 {
    for j in 0..<3 {
      let expected = if i == j { 1.0 } else { 0.0 }
      assert_true((inv[i][j] - expected).abs() < 1.0e-6)
    }
  }
}

test "glm_invert_matrix_1x1" {
  let a = [[5.0]]
  let inv = @src.glm_invert_matrix(a, 1)
  assert_true((inv[0][0] - 0.2).abs() < 1.0e-6)
}

test "glm_invert_matrix_verify_product" {
  // A * A^{-1} should be identity
  let a = [[2.0, 1.0], [5.0, 3.0]]
  let inv = @src.glm_invert_matrix(a, 2)
  // Check A * inv ≈ I
  // a[0][0]*inv[0][0] + a[0][1]*inv[1][0] should be 1
  let elem00 = a[0][0] * inv[0][0] + a[0][1] * inv[1][0]
  let elem01 = a[0][0] * inv[0][1] + a[0][1] * inv[1][1]
  let elem10 = a[1][0] * inv[0][0] + a[1][1] * inv[1][0]
  let elem11 = a[1][0] * inv[0][1] + a[1][1] * inv[1][1]
  assert_true((elem00 - 1.0).abs() < 1.0e-4)
  assert_true(elem01.abs() < 1.0e-4)
  assert_true(elem10.abs() < 1.0e-4)
  assert_true((elem11 - 1.0).abs() < 1.0e-4)
}

// ---------------------------------------------------------------------------
// Statistical utilities
// ---------------------------------------------------------------------------

test "glm_phi_standard_normal" {
  // phi(0) ≈ 0.5
  let p0 = @src.glm_phi(0.0)
  assert_true((p0 - 0.5).abs() < 0.01)
  // phi(1.96) ≈ 0.975
  let p1 = @src.glm_phi(1.96)
  assert_true((p1 - 0.975).abs() < 0.01)
  // phi(-1.96) ≈ 0.025
  let p2 = @src.glm_phi(-1.96)
  assert_true((p2 - 0.025).abs() < 0.01)
}

test "glm_phi_symmetry" {
  // phi(x) + phi(-x) = 1
  let x = 1.5
  let p_sum = @src.glm_phi(x) + @src.glm_phi(-x)
  assert_true((p_sum - 1.0).abs() < 1.0e-10)
}

test "glm_phi_non_decreasing" {
  // phi should be non-decreasing
  let p1 = @src.glm_phi(-2.0)
  let p2 = @src.glm_phi(-1.0)
  let p3 = @src.glm_phi(0.0)
  let p4 = @src.glm_phi(1.0)
  let p5 = @src.glm_phi(2.0)
  assert_true(p1 <= p2)
  assert_true(p2 <= p3)
  assert_true(p3 <= p4)
  assert_true(p4 <= p5)
}

test "glm_phi_clamped" {
  // phi should be in [0, 1]
  assert_true(@src.glm_phi(-10.0) >= 0.0)
  assert_true(@src.glm_phi(-10.0) <= 1.0)
  assert_true(@src.glm_phi(10.0) >= 0.0)
  assert_true(@src.glm_phi(10.0) <= 1.0)
}

test "glm_normal_pvalue_two_sided_zero" {
  // z=0 => p=1.0
  let p = @src.glm_normal_pvalue_two_sided(0.0)
  assert_true((p - 1.0).abs() < 0.01)
}

test "glm_normal_pvalue_two_sided_196" {
  // z=1.96 => p≈0.05
  let p = @src.glm_normal_pvalue_two_sided(1.96)
  assert_true((p - 0.05).abs() < 0.01)
}

test "glm_normal_pvalue_two_sided_symmetry" {
  // z and -z should give same p-value
  let p1 = @src.glm_normal_pvalue_two_sided(2.5)
  let p2 = @src.glm_normal_pvalue_two_sided(-2.5)
  assert_true((p1 - p2).abs() < 1.0e-10)
}

test "glm_normal_pvalue_two_sided_in_range" {
  // All p-values should be in [0, 1]
  for z in [-5.0, -2.0, -1.0, 0.0, 1.0, 2.0, 5.0] {
    let p = @src.glm_normal_pvalue_two_sided(z)
    assert_true(p >= 0.0)
    assert_true(p <= 1.0)
  }
}

test "glm_normal_pvalue_two_sided_large_z" {
  // Large |z| should give very small p-value (clamped to 1e-15)
  let p = @src.glm_normal_pvalue_two_sided(10.0)
  assert_true(p <= 1.0e-10)
}

// ---------------------------------------------------------------------------
// Sorting
// ---------------------------------------------------------------------------

test "glm_sort_pairs_basic" {
  let arr = [(3.0, 2), (1.0, 0), (2.0, 1)]
  @src.glm_sort_pairs(arr)
  assert_eq(arr[0].0, 1.0)
  assert_eq(arr[0].1, 0)
  assert_eq(arr[1].0, 2.0)
  assert_eq(arr[1].1, 1)
  assert_eq(arr[2].0, 3.0)
  assert_eq(arr[2].1, 2)
}

test "glm_sort_pairs_already_sorted" {
  let arr = [(1.0, 0), (2.0, 1), (3.0, 2)]
  @src.glm_sort_pairs(arr)
  assert_eq(arr[0].0, 1.0)
  assert_eq(arr[1].0, 2.0)
  assert_eq(arr[2].0, 3.0)
}

test "glm_sort_pairs_single" {
  let arr = [(5.0, 0)]
  @src.glm_sort_pairs(arr)
  assert_eq(arr.length(), 1)
  assert_eq(arr[0].0, 5.0)
}

test "glm_sort_pairs_empty" {
  let arr : Array[(Double, Int)] = []
  @src.glm_sort_pairs(arr)
  assert_eq(arr.length(), 0)
}

// ---------------------------------------------------------------------------
// BH-FDR correction
// ---------------------------------------------------------------------------

test "glm_bh_correct_single" {
  // Single p-value correction returns same value
  let p = @src.glm_bh_correct(0.05)
  assert_eq(p, 0.05)
}

test "glm_bh_correct_array_known_values" {
  // p-values: [0.01, 0.04, 0.03, 0.005]
  // Sorted:   [0.005, 0.01, 0.03, 0.04]
  // ranks:       1      2     3     4
  // BH:  0.005*4/1=0.02, 0.01*4/2=0.02, 0.03*4/3=0.04, 0.04*4/4=0.04
  // cummin from bottom: 0.02, 0.02, 0.04, 0.04
  // Adjusted (in original order): [0.02, 0.04, 0.04, 0.02]
  let pvals = [0.01, 0.04, 0.03, 0.005]
  let adj = @src.glm_bh_correct_array(pvals)
  assert_eq(adj.length(), 4)
  assert_true((adj[0] - 0.02).abs() < 1.0e-6)
  assert_true((adj[1] - 0.04).abs() < 1.0e-6)
  assert_true((adj[2] - 0.04).abs() < 1.0e-6)
  assert_true((adj[3] - 0.02).abs() < 1.0e-6)
}

test "glm_bh_correct_array_all_same" {
  let pvals = [0.05, 0.05, 0.05]
  let adj = @src.glm_bh_correct_array(pvals)
  // All adjusted should be the same
  assert_true((adj[0] - adj[1]).abs() < 1.0e-10)
  assert_true((adj[1] - adj[2]).abs() < 1.0e-10)
}

test "glm_bh_correct_array_single" {
  let adj = @src.glm_bh_correct_array([0.02])
  assert_eq(adj.length(), 1)
  assert_true((adj[0] - 0.02).abs() < 1.0e-6)
}

test "glm_bh_correct_array_empty" {
  let adj = @src.glm_bh_correct_array([])
  assert_eq(adj.length(), 0)
}

test "glm_bh_correct_array_values_in_range" {
  let pvals = [0.001, 0.02, 0.05, 0.1, 0.5]
  let adj = @src.glm_bh_correct_array(pvals)
  for v in adj {
    assert_true(v >= 0.0)
    assert_true(v <= 1.0)
  }
}

test "glm_bh_correct_array_monotonic" {
  let pvals = [0.001, 0.01, 0.03, 0.04]
  let adj = @src.glm_bh_correct_array(pvals)
  let indexed = [
    (pvals[0], 0),
    (pvals[1], 1),
    (pvals[2], 2),
    (pvals[3], 3),
  ]
  @src.glm_sort_pairs(indexed)
  let adj0 = adj[indexed[0].1]
  let adj1 = adj[indexed[1].1]
  let adj2 = adj[indexed[2].1]
  let adj3 = adj[indexed[3].1]
  assert_true(adj0 <= adj1 + 1.0e-10)
  assert_true(adj1 <= adj2 + 1.0e-10)
  assert_true(adj2 <= adj3 + 1.0e-10)
}

// ---------------------------------------------------------------------------
// Size factor calculation
// ---------------------------------------------------------------------------

test "glm_calculate_sf_basic" {
  // 3 genes, 2 samples
  // Sample 0: [10, 20, 30], Sample 1: [20, 40, 60]
  // Sample 1 has twice the counts, so sf[1] ≈ 2 * sf[0]
  let counts = [
    [10.0, 20.0],
    [20.0, 40.0],
    [30.0, 60.0],
  ]
  let sf = @src.glm_calculate_sf(counts)
  assert_eq(sf.length(), 2)
  assert_true(sf[0] > 0.0)
  assert_true(sf[1] > 0.0)
  // sf[1] / sf[0] ≈ 2.0
  assert_true((sf[1] / sf[0] - 2.0).abs() < 0.5)
}

test "glm_calculate_sf_equal_samples" {
  let counts = [
    [10.0, 10.0, 10.0],
    [20.0, 20.0, 20.0],
  ]
  let sf = @src.glm_calculate_sf(counts)
  assert_eq(sf.length(), 3)
  // All samples have same counts, so sf should be equal
  assert_true((sf[0] - sf[1]).abs() < 1.0e-6)
  assert_true((sf[1] - sf[2]).abs() < 1.0e-6)
}

test "glm_calculate_sf_single_sample" {
  let counts = [[5.0], [10.0], [15.0]]
  let sf = @src.glm_calculate_sf(counts)
  assert_eq(sf.length(), 1)
  assert_true(sf[0] > 0.0)
}

test "glm_calculate_sf_empty" {
  let sf = @src.glm_calculate_sf([])
  assert_eq(sf.length(), 0)
}

test "glm_calculate_sf_all_zeros" {
  let counts = [
    [0.0, 0.0],
    [0.0, 0.0],
  ]
  let sf = @src.glm_calculate_sf(counts)
  assert_eq(sf.length(), 2)
  // Should fallback to 1.0 for zero-count samples
  assert_true(sf[0] > 0.0)
  assert_true(sf[1] > 0.0)
}

test "glm_calculate_sf_genes_with_zeros" {
  let counts = [
    [10.0, 0.0],
    [0.0, 20.0],
    [30.0, 30.0],
  ]
  let sf = @src.glm_calculate_sf(counts)
  assert_eq(sf.length(), 2)
  assert_true(sf[0] > 0.0)
  assert_true(sf[1] > 0.0)
}

// ---------------------------------------------------------------------------
// Pseudobulk aggregation
// ---------------------------------------------------------------------------

test "glm_pseudobulk_basic" {
  // 2 genes, 4 samples (2 control, 2 treated)
  let counts = [
    [10.0, 12.0, 50.0, 48.0],
    [20.0, 22.0, 18.0, 20.0],
  ]
  let gene_names = ["GeneA", "GeneB"]
  let group_labels = ["control", "control", "treated", "treated"]
  let pb = @src.glm_pseudobulk(counts, gene_names, group_labels)
  assert_eq(pb.n_genes(), 2)
  assert_eq(pb.n_groups(), 2)
}

test "glm_pseudobulk_single_group" {
  let counts = [[10.0, 20.0, 30.0]]
  let gene_names = ["GeneX"]
  let group_labels = ["A", "A", "A"]
  let pb = @src.glm_pseudobulk(counts, gene_names, group_labels)
  assert_eq(pb.n_genes(), 1)
  assert_eq(pb.n_groups(), 1)
}

test "glm_pseudobulk_single_sample_per_group" {
  let counts = [
    [5.0, 10.0, 15.0],
  ]
  let gene_names = ["GeneX"]
  let group_labels = ["A", "B", "C"]
  let pb = @src.glm_pseudobulk(counts, gene_names, group_labels)
  assert_eq(pb.n_genes(), 1)
  assert_eq(pb.n_groups(), 3)
}

test "glm_pseudobulk_empty" {
  let counts: Array[Array[Double]] = []
  let gene_names: Array[String] = []
  let group_labels: Array[String] = []
  let pb = @src.glm_pseudobulk(counts, gene_names, group_labels)
  assert_eq(pb.n_genes(), 0)
  assert_eq(pb.n_groups(), 0)
}

// ---------------------------------------------------------------------------
// Initial estimates and dispersion
// ---------------------------------------------------------------------------

test "glm_initial_estimates_constant_counts" {
  let counts = [10.0, 10.0, 10.0, 10.0]
  let d = @src.GlmDesign::new(
    ["s1", "s2", "s3", "s4"],
    ["Intercept"],
    [[1.0], [1.0], [1.0], [1.0]],
  )
  let coefs = @src.glm_initial_estimates(counts, d)
  assert_eq(coefs.length(), 1)
  // With constant counts and intercept-only design,
  // the intercept should be approximately ln(10) ≈ 2.3026
  let expected_ln = @math.ln(10.0)
  assert_true((coefs[0] - expected_ln).abs() < 0.5)
}

test "glm_initial_estimates_two_groups" {
  let counts = [10.0, 12.0, 50.0, 48.0]
  let d = @src.GlmDesign::new(
    ["s1", "s2", "s3", "s4"],
    ["Intercept", "Treatment"],
    [[1.0, 0.0], [1.0, 0.0], [1.0, 1.0], [1.0, 1.0]],
  )
  let coefs = @src.glm_initial_estimates(counts, d)
  assert_eq(coefs.length(), 2)
  // Treatment coefficient should be positive (treated > control)
  assert_true(coefs[1] > 0.0)
}

test "glm_initial_estimates_with_zeros" {
  let counts = [0.0, 10.0, 20.0, 30.0]
  let d = @src.GlmDesign::new(
    ["s1", "s2", "s3", "s4"],
    ["Intercept"],
    [[1.0], [1.0], [1.0], [1.0]],
  )
  let coefs = @src.glm_initial_estimates(counts, d)
  assert_eq(coefs.length(), 1)
  assert_true(coefs[0] > 0.0)
}

test "glm_estimate_dispersion_low" {
  let counts = [10.0, 11.0, 10.5, 9.5]
  let d = @src.GlmDesign::new(
    ["s1", "s2", "s3", "s4"],
    ["Intercept"],
    [[1.0], [1.0], [1.0], [1.0]],
  )
  let coefs = @src.glm_initial_estimates(counts, d)
  let disp = @src.glm_estimate_dispersion(counts, d, coefs)
  // Low dispersion for tightly clustered counts
  assert_true(disp >= 0.0)
  assert_true(disp < 5.0)
}

test "glm_estimate_dispersion_high" {
  let counts = [1.0, 50.0, 2.0, 48.0]
  let d = @src.GlmDesign::new(
    ["s1", "s2", "s3", "s4"],
    ["Intercept"],
    [[1.0], [1.0], [1.0], [1.0]],
  )
  let coefs = @src.glm_initial_estimates(counts, d)
  let disp = @src.glm_estimate_dispersion(counts, d, coefs)
  // High dispersion for highly variable counts
  assert_true(disp >= 0.0)
}

test "glm_estimate_dispersion_all_same" {
  let counts = [20.0, 20.0, 20.0, 20.0]
  let d = @src.GlmDesign::new(
    ["s1", "s2", "s3", "s4"],
    ["Intercept"],
    [[1.0], [1.0], [1.0], [1.0]],
  )
  let coefs = @src.glm_initial_estimates(counts, d)
  let disp = @src.glm_estimate_dispersion(counts, d, coefs)
  // Zero variance => dispersion should be very small or 0
  assert_true(disp >= 0.0)
}

// ---------------------------------------------------------------------------
// IWLCS fit
// ---------------------------------------------------------------------------

test "glm_iwlcs_fit_produces_coefs" {
  let counts = [10.0, 12.0, 50.0, 48.0]
  let d = @src.GlmDesign::new(
    ["s1", "s2", "s3", "s4"],
    ["Intercept", "Treatment"],
    [[1.0, 0.0], [1.0, 0.0], [1.0, 1.0], [1.0, 1.0]],
  )
  let initial = @src.glm_initial_estimates(counts, d)
  let coefs = @src.glm_iwlcs_fit(counts, d, 0.1, initial)
  assert_eq(coefs.length(), 2)
  // Treatment effect should remain positive
  assert_true(coefs[1] > 0.0)
}

test "glm_iwlcs_fit_converges" {
  let counts = [10.0, 10.0, 10.0, 10.0]
  let d = @src.GlmDesign::new(
    ["s1", "s2", "s3", "s4"],
    ["Intercept"],
    [[1.0], [1.0], [1.0], [1.0]],
  )
  let initial = @src.glm_initial_estimates(counts, d)
  let coefs = @src.glm_iwlcs_fit(counts, d, 0.1, initial)
  assert_eq(coefs.length(), 1)
  // For constant data, fitted values should match approximately
  let fitted = @math.exp(coefs[0])
  assert_true((fitted - 10.0).abs() < 1.0)
}

// ---------------------------------------------------------------------------
// Standard error computation
// ---------------------------------------------------------------------------

test "glm_compute_se_positive" {
  let d = @src.GlmDesign::new(
    ["s1", "s2", "s3", "s4"],
    ["Intercept", "Treatment"],
    [[1.0, 0.0], [1.0, 0.0], [1.0, 1.0], [1.0, 1.0]],
  )
  let coefs = [2.3, 1.5]
  let sf = [1.0, 1.0, 1.0, 1.0]
  let se = @src.glm_compute_se(d, 0.1, coefs, sf)
  assert_eq(se.length(), 2)
  assert_true(se[0] > 0.0)
  assert_true(se[1] > 0.0)
}

test "glm_compute_se_single_coef" {
  let d = @src.GlmDesign::new(
    ["s1", "s2", "s3"],
    ["Intercept"],
    [[1.0], [1.0], [1.0]],
  )
  let coefs = [2.5]
  let sf = [1.0, 1.0, 1.0]
  let se = @src.glm_compute_se(d, 0.2, coefs, sf)
  assert_eq(se.length(), 1)
  assert_true(se[0] > 0.0)
}

// ---------------------------------------------------------------------------
// Single-gene fitting
// ---------------------------------------------------------------------------

test "glm_fit_one_gene_basic" {
  let counts = [10.0, 12.0, 50.0, 48.0]
  let d = @src.GlmDesign::new(
    ["s1", "s2", "s3", "s4"],
    ["Intercept", "Treatment"],
    [[1.0, 0.0], [1.0, 0.0], [1.0, 1.0], [1.0, 1.0]],
  )
  let sf = [1.0, 1.0, 1.0, 1.0]
  let fit = @src.glm_fit_one_gene(counts, d, sf, "GeneA")
  assert_eq(fit.gene, "GeneA")
  assert_eq(fit.coefficients.length(), 2)
  assert_eq(fit.std_errors.length(), 2)
  assert_true(fit.dispersion >= 0.0)
  assert_eq(fit.df_residual, 2)
  assert_true(fit.df_dispersion > 0.0)
}

test "glm_fit_one_gene_treatment_positive" {
  let counts = [10.0, 12.0, 50.0, 48.0]
  let d = @src.GlmDesign::new(
    ["s1", "s2", "s3", "s4"],
    ["Intercept", "Treatment"],
    [[1.0, 0.0], [1.0, 0.0], [1.0, 1.0], [1.0, 1.0]],
  )
  let sf = [1.0, 1.0, 1.0, 1.0]
  let fit = @src.glm_fit_one_gene(counts, d, sf, "GeneA")
  // Treatment coefficient should be positive (treated > control)
  assert_true(fit.coefficients[1] > 0.0)
}

test "glm_fit_one_gene_no_treatment_effect" {
  let counts = [20.0, 22.0, 18.0, 20.0]
  let d = @src.GlmDesign::new(
    ["s1", "s2", "s3", "s4"],
    ["Intercept", "Treatment"],
    [[1.0, 0.0], [1.0, 0.0], [1.0, 1.0], [1.0, 1.0]],
  )
  let sf = [1.0, 1.0, 1.0, 1.0]
  let fit = @src.glm_fit_one_gene(counts, d, sf, "GeneB")
  // Treatment coefficient should be near zero
  assert_true(fit.coefficients[1].abs() < 1.0)
}

// ---------------------------------------------------------------------------
// Differential expression testing
// ---------------------------------------------------------------------------

test "glm_test_de_basic" {
  // 2 genes, 4 samples
  let counts = [
    [10.0, 12.0, 50.0, 48.0],
    [20.0, 22.0, 18.0, 20.0],
  ]
  let gene_names = ["GeneA", "GeneB"]
  let d = @src.GlmDesign::new(
    ["s1", "s2", "s3", "s4"],
    ["Intercept", "Treatment"],
    [[1.0, 0.0], [1.0, 0.0], [1.0, 1.0], [1.0, 1.0]],
  )
  let sf = [1.0, 1.0, 1.0, 1.0]
  let results = @src.glm_test_de(counts, gene_names, d, sf, 1)
  assert_eq(results.length(), 2)
  assert_eq(results[0].gene, "GeneA")
  assert_eq(results[1].gene, "GeneB")
}

test "glm_test_de_pvalues_in_range" {
  let counts = [
    [10.0, 12.0, 50.0, 48.0],
    [20.0, 22.0, 18.0, 20.0],
  ]
  let gene_names = ["GeneA", "GeneB"]
  let d = @src.GlmDesign::new(
    ["s1", "s2", "s3", "s4"],
    ["Intercept", "Treatment"],
    [[1.0, 0.0], [1.0, 0.0], [1.0, 1.0], [1.0, 1.0]],
  )
  let sf = [1.0, 1.0, 1.0, 1.0]
  let results = @src.glm_test_de(counts, gene_names, d, sf, 1)
  for r in results {
    assert_true(r.p_value() >= 0.0)
    assert_true(r.p_value() <= 1.0)
    assert_true(r.adj_p_value() >= 0.0)
    assert_true(r.adj_p_value() <= 1.0)
  }
}

test "glm_test_de_detects_differential" {
  // GeneA has strong treatment effect, GeneB does not
  let counts = [
    [5.0, 6.0, 80.0, 75.0],
    [20.0, 22.0, 18.0, 20.0],
  ]
  let gene_names = ["GeneA", "GeneB"]
  let d = @src.GlmDesign::new(
    ["s1", "s2", "s3", "s4"],
    ["Intercept", "Treatment"],
    [[1.0, 0.0], [1.0, 0.0], [1.0, 1.0], [1.0, 1.0]],
  )
  let sf = [1.0, 1.0, 1.0, 1.0]
  let results = @src.glm_test_de(counts, gene_names, d, sf, 1)
  // GeneA should have larger absolute z-value
  let z_a = results[0].z_value.abs()
  let z_b = results[1].z_value.abs()
  assert_true(z_a > z_b)
}

test "glm_test_de_estimate_sign" {
  let counts = [
    [10.0, 12.0, 50.0, 48.0],
    [50.0, 48.0, 10.0, 12.0],
  ]
  let gene_names = ["GeneUp", "GeneDown"]
  let d = @src.GlmDesign::new(
    ["s1", "s2", "s3", "s4"],
    ["Intercept", "Treatment"],
    [[1.0, 0.0], [1.0, 0.0], [1.0, 1.0], [1.0, 1.0]],
  )
  let sf = [1.0, 1.0, 1.0, 1.0]
  let results = @src.glm_test_de(counts, gene_names, d, sf, 1)
  // GeneUp is up-regulated in treated (positive estimate)
  assert_true(results[0].estimate() > 0.0)
  // GeneDown is down-regulated in treated (negative estimate)
  assert_true(results[1].estimate() < 0.0)
}

test "glm_test_de_single_gene" {
  let counts = [[10.0, 12.0, 50.0, 48.0]]
  let gene_names = ["GeneA"]
  let d = @src.GlmDesign::new(
    ["s1", "s2", "s3", "s4"],
    ["Intercept", "Treatment"],
    [[1.0, 0.0], [1.0, 0.0], [1.0, 1.0], [1.0, 1.0]],
  )
  let sf = [1.0, 1.0, 1.0, 1.0]
  let results = @src.glm_test_de(counts, gene_names, d, sf, 1)
  assert_eq(results.length(), 1)
  assert_eq(results[0].gene, "GeneA")
}

// ---------------------------------------------------------------------------
// Significant filtering
// ---------------------------------------------------------------------------

test "glm_significant_filters_correctly" {
  let r1 = @src.GlmTestResult::new("GeneA", "Treatment", 3.0, 0.5, 6.0, 1.0e-10)
  let r2 = @src.GlmTestResult::new("GeneB", "Treatment", 0.1, 0.3, 0.3, 0.8)
  let r3 = @src.GlmTestResult::new("GeneC", "Treatment", 2.0, 0.4, 5.0, 0.001)
  let results = [r1, r2, r3]
  let sig = @src.glm_significant(results, 0.05)
  // r1 and r3 have small p-values, r2 does not
  assert_true(sig.length() >= 1)
  // Check that non-significant gene is excluded
  for s in sig {
    if s.gene == "GeneB" {
      assert_false(true)
    }
  }
}

test "glm_significant_empty_input" {
  let sig = @src.glm_significant([], 0.05)
  assert_eq(sig.length(), 0)
}

test "glm_significant_no_passing" {
  let r1 = @src.GlmTestResult::new("GeneA", "Treatment", 0.0, 0.5, 0.0, 0.5)
  let r2 = @src.GlmTestResult::new("GeneB", "Treatment", 0.1, 0.3, 0.3, 0.8)
  let results = [r1, r2]
  let sig = @src.glm_significant(results, 0.0)
  assert_eq(sig.length(), 0)
}

test "glm_significant_all_pass" {
  let r1 = @src.GlmTestResult::new("GeneA", "Treatment", 5.0, 0.5, 10.0, 1.0e-10)
  let r2 = @src.GlmTestResult::new("GeneB", "Treatment", 4.0, 0.3, 13.0, 1.0e-8)
  let results = [r1, r2]
  let sig = @src.glm_significant(results, 1.0)
  assert_eq(sig.length(), 2)
}

// ---------------------------------------------------------------------------
// Edge cases and integration
// ---------------------------------------------------------------------------

test "glm_design_single_sample" {
  let d = @src.GlmDesign::new(
    ["s1"],
    ["Intercept"],
    [[1.0]],
  )
  assert_eq(d.n_samples(), 1)
  assert_eq(d.n_coefs(), 1)
  assert_eq(d.get(0, 0), 1.0)
}

test "glm_design_many_coefs" {
  let d = @src.GlmDesign::new(
    ["s1", "s2"],
    ["A", "B", "C", "D"],
    [[1.0, 0.0, 0.0, 0.0], [0.0, 1.0, 0.0, 0.0]],
  )
  assert_eq(d.n_samples(), 2)
  assert_eq(d.n_coefs(), 4)
  assert_eq(d.get(0, 0), 1.0)
  assert_eq(d.get(1, 1), 1.0)
}

test "glm_fit_constant_counts" {
  let counts = [20.0, 20.0, 20.0, 20.0]
  let d = @src.GlmDesign::new(
    ["s1", "s2", "s3", "s4"],
    ["Intercept"],
    [[1.0], [1.0], [1.0], [1.0]],
  )
  let sf = [1.0, 1.0, 1.0, 1.0]
  let fit = @src.glm_fit_one_gene(counts, d, sf, "GeneConst")
  assert_eq(fit.gene, "GeneConst")
  // For constant counts, dispersion should be near zero
  assert_true(fit.dispersion >= 0.0)
}

test "glm_test_de_with_size_factors" {
  let counts = [
    [10.0, 20.0, 50.0, 100.0],
    [20.0, 40.0, 18.0, 36.0],
  ]
  let gene_names = ["GeneA", "GeneB"]
  let d = @src.GlmDesign::new(
    ["s1", "s2", "s3", "s4"],
    ["Intercept", "Treatment"],
    [[1.0, 0.0], [1.0, 0.0], [1.0, 1.0], [1.0, 1.0]],
  )
  let sf = [1.0, 2.0, 1.0, 2.0]
  let results = @src.glm_test_de(counts, gene_names, d, sf, 1)
  assert_eq(results.length(), 2)
}

test "glm_full_pipeline" {
  // Simulate a full differential expression analysis pipeline
  let counts = [
    [10.0, 12.0, 50.0, 48.0],
    [20.0, 22.0, 18.0, 20.0],
    [5.0, 6.0, 45.0, 52.0],
  ]
  let gene_names = ["GeneA", "GeneB", "GeneC"]

  // Step 1: Calculate size factors
  let sf = @src.glm_calculate_sf(counts)
  assert_eq(sf.length(), 4)
  for s in sf {
    assert_true(s > 0.0)
  }

  // Step 2: Create design matrix
  let d = @src.GlmDesign::new(
    ["s1", "s2", "s3", "s4"],
    ["Intercept", "Treatment"],
    [[1.0, 0.0], [1.0, 0.0], [1.0, 1.0], [1.0, 1.0]],
  )

  // Step 3: Run DE test
  let results = @src.glm_test_de(counts, gene_names, d, sf, 1)
  assert_eq(results.length(), 3)

  // Step 4: Filter significant
  let sig = @src.glm_significant(results, 0.05)
  // Significant genes (if any) should have valid adjusted p-values
  for r in sig {
    assert_true(r.adj_p_value() < 0.05)
  }

  // Step 5: Verify all results have valid statistics
  for r in results {
    assert_true(r.p_value() >= 0.0 && r.p_value() <= 1.0)
    assert_true(r.adj_p_value() >= 0.0 && r.adj_p_value() <= 1.0)
  }
}

test "glm_phi_05" {
  // phi(0.5) should be approximately 0.6915
  let p = @src.glm_phi(0.5)
  assert_true((p - 0.6915).abs() < 0.01)
}

test "glm_phi_2" {
  // phi(2.0) should be approximately 0.9772
  let p = @src.glm_phi(2.0)
  assert_true((p - 0.9772).abs() < 0.01)
}

test "glm_phi_negative" {
  // phi(-1.0) should be approximately 0.1587
  let p = @src.glm_phi(-1.0)
  assert_true((p - 0.1587).abs() < 0.01)
}

test "glm_median_two_elements" {
  let m = @src.glm_median([3.0, 7.0])
  assert_eq(m, 5.0)
}

test "glm_median_with_duplicates" {
  let m = @src.glm_median([1.0, 2.0, 2.0, 3.0, 3.0])
  assert_eq(m, 2.0)
}

test "glm_min_positive_all_zero" {
  let v = @src.glm_min_positive([0.0, 0.0, 0.0])
  assert_eq(v, 0.0)
}

test "glm_calculate_sf_two_genes" {
  let counts = [
    [100.0, 10.0],
    [200.0, 20.0],
  ]
  let sf = @src.glm_calculate_sf(counts)
  assert_eq(sf.length(), 2)
  // Sample 0 has 10x more counts than sample 1
  // sf[0] / sf[1] should be approximately 10
  assert_true((sf[0] / sf[1] - 10.0).abs() < 2.0)
}

test "glm_test_de_all_zeros" {
  let counts = [
    [0.0, 0.0, 0.0, 0.0],
    [0.0, 0.0, 0.0, 0.0],
  ]
  let gene_names = ["GeneA", "GeneB"]
  let d = @src.GlmDesign::new(
    ["s1", "s2", "s3", "s4"],
    ["Intercept", "Treatment"],
    [[1.0, 0.0], [1.0, 0.0], [1.0, 1.0], [1.0, 1.0]],
  )
  let sf = [1.0, 1.0, 1.0, 1.0]
  let results = @src.glm_test_de(counts, gene_names, d, sf, 1)
  assert_eq(results.length(), 2)
  // All p-values should be valid (not NaN)
  for r in results {
    assert_true(r.p_value() >= 0.0)
    assert_true(r.p_value() <= 1.0)
  }
}

test "glm_bh_correct_array_already_corrected" {
  // When all p-values are very small, BH should keep them small
  let pvals = [1.0e-10, 1.0e-8, 1.0e-6, 1.0e-4]
  let adj = @src.glm_bh_correct_array(pvals)
  for v in adj {
    assert_true(v >= 0.0)
    assert_true(v <= 1.0)
  }
  // The smallest p-value should remain small
  assert_true(adj[0] <= 1.0e-4)
}