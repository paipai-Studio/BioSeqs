///|
/// Tests for Bioconductor HTSFilter module - RNA-seq count filtering.

test "cpm_single_value" {
  let cpm = @src.hts_filter_single_cpm(1000.0, 1000000.0)
  assert_eq(cpm, 1000.0)
}

test "cpm_zero_library" {
  let cpm = @src.hts_filter_single_cpm(1000.0, 0.0)
  assert_eq(cpm, 0.0)
}

test "cpm_matrix_basic" {
  let counts = [[100.0, 200.0], [300.0, 400.0]]
  let libs = [1000000.0, 2000000.0]
  let cpm = @src.hts_filter_cpm(counts, libs)
  assert_eq(cpm.length(), 2)
  assert_eq(cpm[0].length(), 2)
  // Gene 0, sample 0: 100/1e6 * 1e6 = 100
  assert_true((cpm[0][0] - 100.0).abs() < 0.01)
}

test "library_sizes" {
  let counts = [[10.0, 20.0, 30.0], [40.0, 50.0, 60.0]]
  let libs = @src.hts_filter_library_sizes(counts)
  assert_eq(libs.length(), 3)
  assert_true((libs[0] - 50.0).abs() < 0.01)
  assert_true((libs[1] - 70.0).abs() < 0.01)
  assert_true((libs[2] - 90.0).abs() < 0.01)
}

test "filter_basic" {
  let (counts, groups, _) = @src.hts_filter_sample_data()
  let result = @src.hts_filter(counts, groups, 100.0, 2)
  assert_eq(result.n_genes_input, 20)
  assert_true(result.n_genes_kept > 0)
  assert_true(result.n_genes_removed >= 0)
  assert_eq(result.keep.length(), 20)
}

test "filter_strict_threshold" {
  let (counts, groups, _) = @src.hts_filter_sample_data()
  // Very high threshold - should remove most genes
  let result = @src.hts_filter(counts, groups, 10000.0, 3)
  assert_true(result.n_genes_kept <= result.n_genes_input)
}

test "filter_keep_mask" {
  let (counts, groups, _) = @src.hts_filter_sample_data()
  let result = @src.hts_filter(counts, groups, 100.0, 2)
  // Keep mask length should match input
  assert_eq(result.keep.length(), result.n_genes_input)
  // Each element should be boolean
  let mut i = 0
  let mut true_count = 0
  while i < result.keep.length() {
    if result.keep[i] { true_count = true_count + 1 }
    i = i + 1
  }
  assert_eq(true_count, result.n_genes_kept)
}

test "filter_apply" {
  let counts = [[10.0, 20.0], [30.0, 40.0], [50.0, 60.0]]
  let keep = [true, false, true]
  let filtered = @src.hts_filter_apply(counts, keep)
  assert_eq(filtered.length(), 2)
  assert_eq(filtered[0][0], 10.0)
  assert_eq(filtered[1][0], 50.0)
}

test "filter_apply_names" {
  let names = ["A", "B", "C"]
  let keep = [true, false, true]
  let filtered = @src.hts_filter_apply_names(names, keep)
  assert_eq(filtered.length(), 2)
  assert_eq(filtered[0], "A")
  assert_eq(filtered[1], "C")
}

test "filter_retention_rate" {
  let (counts, groups, _) = @src.hts_filter_sample_data()
  let result = @src.hts_filter(counts, groups, 100.0, 2)
  let rate = @src.hts_filter_retention_rate(result)
  assert_true(rate >= 0.0 && rate <= 1.0)
}

test "filter_summary" {
  let (counts, groups, _) = @src.hts_filter_sample_data()
  let result = @src.hts_filter(counts, groups, 100.0, 2)
  let summary = @src.hts_filter_summary(result)
  assert_true(summary.length() > 0)
}

test "filter_empty_counts" {
  let result = @src.hts_filter([], [0, 0, 1, 1], 1.0, 1)
  assert_eq(result.n_genes_input, 0)
  assert_eq(result.n_genes_kept, 0)
  assert_eq(result.n_genes_removed, 0)
}

test "filter_single_gene" {
  let counts = [[500.0, 600.0, 550.0, 100.0, 120.0, 110.0]]
  let groups = [0, 0, 0, 1, 1, 1]
  let result = @src.hts_filter(counts, groups, 100.0, 2)
  assert_eq(result.n_genes_input, 1)
  assert_true(result.n_genes_kept <= 1)
}

test "filter_get_keep" {
  let (counts, groups, _) = @src.hts_filter_sample_data()
  let result = @src.hts_filter(counts, groups, 100.0, 2)
  let keep = @src.hts_filter_get_keep(result)
  assert_eq(keep.length(), result.keep.length())
}

test "cpm_zeros" {
  let counts = [[0.0, 0.0], [0.0, 0.0]]
  let libs = [1000000.0, 2000000.0]
  let cpm = @src.hts_filter_cpm(counts, libs)
  assert_eq(cpm[0][0], 0.0)
  assert_eq(cpm[1][1], 0.0)
}
