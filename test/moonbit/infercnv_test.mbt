///|
/// Test file for infercnv module.

// ---------------------------------------------------------------------------
// Data structure basics
// ---------------------------------------------------------------------------

test "infercnv_gene_position_creation" {
  let gp = @src.gene_position("TP53", "chr17", 7668402, 7687490)
  assert_eq(gp.gene_id, "TP53")
  assert_eq(gp.chromosome, "chr17")
  assert_eq(gp.start, 7668402)
  assert_eq(gp.end, 7687490)
}

test "infercnv_ordered_genes_natural_chromosome_sort" {
  let genes = [
    @src.gene_position("G_chr10_a", "chr10", 100, 200),
    @src.gene_position("G_chr1_a", "chr1", 300, 400),
    @src.gene_position("G_chr2_a", "chr2", 100, 200),
    @src.gene_position("G_chr1_b", "chr1", 100, 200),
    @src.gene_position("G_chrX_a", "chrX", 100, 200),
  ]
  let og = @src.ordered_genes(genes)
  assert_eq(og.n_genes, 5)
  // chr1 (b within chr1 first, then a) -> chr2 -> chr10 -> chrX
  assert_eq(og.gene_order[0].gene_id, "G_chr1_b")
  assert_eq(og.gene_order[1].gene_id, "G_chr1_a")
  assert_eq(og.gene_order[2].gene_id, "G_chr2_a")
  assert_eq(og.gene_order[3].gene_id, "G_chr10_a")
  assert_eq(og.gene_order[4].gene_id, "G_chrX_a")
}

test "infercnv_exclude_chromosomes_filters_sex_and_mito" {
  let genes = [
    @src.gene_position("A", "chr1", 1, 100),
    @src.gene_position("B", "chrX", 1, 100),
    @src.gene_position("C", "chr2", 1, 100),
    @src.gene_position("D", "chrY", 1, 100),
    @src.gene_position("E", "chrM", 1, 100),
  ]
  let og = @src.ordered_genes(genes)
  let filtered = @src.exclude_chromosomes(og, ["chrX", "chrY", "chrM"])
  assert_eq(filtered.n_genes, 2)
  let ids = filtered.gene_order.map(fn (g) -> String { g.gene_id })
  assert_eq(ids, ["A", "C"])
}

// ---------------------------------------------------------------------------
// Log normalization
// ---------------------------------------------------------------------------

test "infercnv_log_normalize_counts_basic" {
  let raw = [
    [1.0, 2.0, 7.0],    // total = 10 -> per gene * 10000
  ]
  let norm = @src.log_normalize_counts(raw, target_sum=10000.0)
  // first row: 1000, 2000, 7000; log2(x+1)
  // just check it's positive and reasonable
  assert_true(norm[0][0] > 9.0 && norm[0][0] < 10.0)
  assert_true(norm[0][1] > 10.0 && norm[0][1] < 11.0)
  assert_true(norm[0][2] > 12.0 && norm[0][2] < 13.0)
}

test "infercnv_log_normalize_counts_zero_cell_preserved" {
  let raw = [[0.0, 0.0, 0.0]]
  let norm = @src.log_normalize_counts(raw)
  assert_eq(norm[0], [0.0, 0.0, 0.0])
}

// ---------------------------------------------------------------------------
// Reference method helpers
// ---------------------------------------------------------------------------

test "infercnv_ref_method_constructors" {
  match @src.ref_method_global_mean() {
    @src.ReferenceMethod::GlobalMean => assert_true(true)
    _ => assert_true(false)
  }
  match @src.ref_method_reference_categories() {
    @src.ReferenceMethod::ReferenceCategories => assert_true(true)
    _ => assert_true(false)
  }
  let cust = @src.ref_method_custom([1.0, 2.0, 3.0])
  match cust {
    @src.ReferenceMethod::Custom(arr) => assert_eq(arr.length(), 3)
    _ => assert_true(false)
  }
}

test "infercnv_default_params_sensible" {
  let p = @src.default_cnv_params()
  assert_eq(p.window_size, 100)
  assert_eq(p.lfc_cap, 1.5)
  assert_eq(p.noise_threshold, 0.2)
}

// ---------------------------------------------------------------------------
// Synthetic data generator
// ---------------------------------------------------------------------------

test "infercnv_sample_data_shape" {
  let (input, tumour_cats, normal_cats) = @src.infercnv_sample_data(
    n_tumour=20, n_normal=10, n_chr=3, n_genes_per_chr=20,
  )
  assert_eq(input.n_cells, 30)
  assert_eq(input.n_genes, 60)
  assert_eq(input.annotations.length(), 30)
  assert_eq(tumour_cats, ["Tumour"])
  assert_eq(normal_cats, ["Normal"])
  // first 20 cells are tumour
  assert_eq(input.annotations[0].category, "Tumour")
  assert_eq(input.annotations[29].category, "Normal")
}

// ---------------------------------------------------------------------------
// End-to-end pipeline with sample data and all reference strategies
// ---------------------------------------------------------------------------

test "infercnv_run_pipeline_global_mean" {
  let (input, _tumour, _normal) = @src.infercnv_sample_data(
    n_tumour=8, n_normal=6, n_chr=3, n_genes_per_chr=30, seed=1,
  )
  let params = @src.make_cnv_params(
    11, 1.5, 0.05, @src.ref_method_global_mean(), [],
  )
  let res = @src.run_infercnv(input, params=params, excluded_chromosomes=[])
  assert_eq(res.n_cells, 14)
  assert_eq(res.n_genes, 90)
  assert_eq(res.cnv_matrix.length(), 14)
  assert_eq(res.cnv_matrix[0].length(), 90)
  // per_cluster should have Tumour / Normal
  assert_eq(res.per_cluster_cnv_score.length(), 2)
  // With GlobalMean, tumour CNV signal should be detectable (non-zero)
  let t_score = res.cluster_score("Tumour")
  let n_score = res.cluster_score("Normal")
  assert_true(t_score > 0.0)
  assert_true(n_score >= 0.0)
}

test "infercnv_run_pipeline_reference_categories" {
  let (input, _tumour, normal) = @src.infercnv_sample_data(
    n_tumour=10, n_normal=8, n_chr=2, n_genes_per_chr=20, seed=7,
  )
  let params = @src.make_cnv_params(
    9, 1.5, 0.05, @src.ref_method_reference_categories(), normal,
  )
  let res = @src.run_infercnv(input, params=params, excluded_chromosomes=[])
  assert_eq(res.n_cells, 18)
  assert_eq(res.n_genes, 40)
  let t_score = res.cluster_score("Tumour")
  let n_score = res.cluster_score("Normal")
  // tumour score must exceed normal (simulated deletions + amplifications)
  assert_true(t_score > n_score * 1.2)
}

test "infercnv_run_pipeline_custom_reference" {
  let (input, _t, _n) = @src.infercnv_sample_data(
    n_tumour=5, n_normal=4, n_chr=2, n_genes_per_chr=16, seed=99,
  )
  // flat reference of 3.0 for all genes (for testing the path only)
  let ref : Array[Double] = []
  for _i in 0..<input.n_genes { ref.push(3.0) }
  let params = @src.make_cnv_params(
    100, 1.5, 0.0, @src.ref_method_custom(ref), [],
  )
  let res = @src.run_infercnv(input, params=params, excluded_chromosomes=[])
  assert_eq(res.n_genes, input.n_genes)
  assert_eq(res.per_cell_cnv_score.length(), input.n_cells)
}

// ---------------------------------------------------------------------------
// CNVResult methods
// ---------------------------------------------------------------------------

test "infercnv_result_summary_nonempty" {
  let (input, _, _) = @src.infercnv_sample_data(
    n_tumour=6, n_normal=4, n_chr=2, n_genes_per_chr=20, seed=3,
  )
  let res = @src.run_infercnv(input, excluded_chromosomes=[])
  let s = res.summary()
  assert_true(s.contains("cells="))
  assert_true(s.contains("genes="))
  assert_true(s.contains("mean_score="))
  assert_true(s.contains("top_cluster="))
}

test "infercnv_result_at_out_of_range_is_zero" {
  let (input, _, _) = @src.infercnv_sample_data(n_tumour=3, n_normal=2, n_chr=2, n_genes_per_chr=10)
  let res = @src.run_infercnv(input, excluded_chromosomes=[])
  assert_eq(res.at(-1, 0), 0.0)
  assert_eq(res.at(99, 0), 0.0)
  assert_eq(res.at(0, 99999), 0.0)
}

test "infercnv_result_chromosome_list_matches_input" {
  let (input, _, _) = @src.infercnv_sample_data(
    n_tumour=4, n_normal=3, n_chr=3, n_genes_per_chr=12,
  )
  let res = @src.run_infercnv(input, excluded_chromosomes=[])
  // 3 autosomes in synthetic data -> 3 boundaries
  assert_eq(res.chromosome_list().length(), 3)
  let (chr1, s1, e1) = res.chromosome_list()[0]
  assert_eq(chr1, "chr1")
  assert_eq(s1, 0)
  assert_eq(e1, 11)   // 12 genes per chr, indexed 0..11
}

test "infercnv_predict_tumour_cells_flags_tumour_like" {
  let (input, _t, _n) = @src.infercnv_sample_data(
    n_tumour=15, n_normal=10, n_chr=3, n_genes_per_chr=30, seed=42,
  )
  let params = @src.make_cnv_params(
    11, 1.5, 0.05, @src.ref_method_reference_categories(), ["Normal"],
  )
  let res = @src.run_infercnv(input, params=params, excluded_chromosomes=[])
  let predicted = res.predict_tumour_cells("Normal", threshold_factor=1.3)
  // Majority of the tumour cells (indices 0..14) should be flagged
  let mut tumour_flagged = 0
  for idx in predicted {
    if idx < 15 { tumour_flagged = tumour_flagged + 1 }
  }
  // At least 10/15 tumour cells called
  assert_true(tumour_flagged >= 10)
  // Very few normal cells (indices 15..24) should be flagged
  let normal_false = predicted.length() - tumour_flagged
  assert_true(normal_false <= 3)
}
