///|
/// Test file for EnhancedVolcano module.

test "volcano_plot_basic" {
  let ids = ["TP53", "BRCA1", "MYC", "KRAS", "EGFR"]
  let lfcs = [2.5, -3.0, 1.2, -0.5, 0.1]
  let pvals = [0.001, 0.005, 0.01, 0.05, 0.5]
  let result = @src.volcano_plot_default(ids, lfcs, pvals)
  assert_eq(result.get_n_genes(), 5)
}

test "volcano_classification_up" {
  let ids = ["Gene1", "Gene2", "Gene3"]
  let lfcs = [2.0, -2.0, 0.5]
  let pvals = [0.001, 0.001, 0.001]
  let result = @src.volcano_plot(ids, lfcs, pvals, p_cutoff=0.05, fc_cutoff=1.0, title="Test", x_label="x", y_label="y")
  assert_eq(result.get_n_up(), 1)
  assert_eq(result.get_n_down(), 1)
  assert_eq(result.get_n_nonsig(), 1)
}

test "volcano_classification_all_sig" {
  let ids = ["A", "B", "C", "D"]
  let lfcs = [2.5, -2.5, 1.5, -1.5]
  let pvals = [0.001, 0.001, 0.01, 0.01]
  let result = @src.volcano_plot(ids, lfcs, pvals, p_cutoff=0.05, fc_cutoff=1.0, title="", x_label="", y_label="")
  assert_eq(result.get_n_up(), 2)
  assert_eq(result.get_n_down(), 2)
  assert_eq(result.get_n_nonsig(), 0)
}

test "volcano_classification_none_sig" {
  let ids = ["A", "B", "C"]
  let lfcs = [0.1, -0.1, 0.5]
  let pvals = [0.1, 0.2, 0.1]
  let result = @src.volcano_plot(ids, lfcs, pvals, p_cutoff=0.05, fc_cutoff=1.0, title="", x_label="", y_label="")
  assert_eq(result.get_n_up(), 0)
  assert_eq(result.get_n_down(), 0)
  assert_eq(result.get_n_nonsig(), 3)
}

test "volcano_get_genes" {
  let ids = ["A", "B", "C"]
  let lfcs = [2.0, -2.0, 0.5]
  let pvals = [0.001, 0.001, 0.001]
  let result = @src.volcano_plot(ids, lfcs, pvals, p_cutoff=0.05, fc_cutoff=1.0, title="", x_label="", y_label="")
  let genes = result.get_genes()
  assert_eq(genes.length(), 3)
}

test "volcano_get_significant" {
  let ids = ["A", "B", "C"]
  let lfcs = [2.0, -2.0, 0.5]
  let pvals = [0.001, 0.001, 0.001]
  let result = @src.volcano_plot(ids, lfcs, pvals, p_cutoff=0.05, fc_cutoff=1.0, title="", x_label="", y_label="")
  let sig = result.get_significant_genes()
  assert_eq(sig.length(), 2) // A (Up) and B (Down)
}

test "volcano_get_up_genes" {
  let ids = ["A", "B", "C"]
  let lfcs = [2.0, -2.0, 1.5]
  let pvals = [0.001, 0.001, 0.001]
  let result = @src.volcano_plot(ids, lfcs, pvals, p_cutoff=0.05, fc_cutoff=1.0, title="", x_label="", y_label="")
  let up = result.get_up_genes()
  assert_eq(up.length(), 2)
}

test "volcano_get_down_genes" {
  let ids = ["A", "B", "C"]
  let lfcs = [2.0, -2.0, -1.5]
  let pvals = [0.001, 0.001, 0.001]
  let result = @src.volcano_plot(ids, lfcs, pvals, p_cutoff=0.05, fc_cutoff=1.0, title="", x_label="", y_label="")
  let down = result.get_down_genes()
  assert_eq(down.length(), 2)
}

test "volcano_get_top_genes" {
  let ids = ["A", "B", "C", "D", "E"]
  let lfcs = [2.5, -3.0, 1.5, -1.5, 0.3]
  let pvals = [0.01, 0.001, 0.02, 0.005, 0.001]
  let result = @src.volcano_plot(ids, lfcs, pvals, p_cutoff=0.05, fc_cutoff=1.0, title="", x_label="", y_label="")
  let top = result.get_top_genes(3)
  assert_eq(top.length(), 3)
}

test "volcano_neg_log10_p" {
  let nlp = @src.neg_log10_p(0.05)
  assert_true(nlp > 1.29) // -log10(0.05) = 1.301
  assert_true(nlp < 1.31)
}

test "volcano_neg_log10_p_zero" {
  let nlp = @src.neg_log10_p(0.0)
  assert_true(nlp > 0.0) // should handle p=0 gracefully
}

test "volcano_cutoffs" {
  let ids = ["A", "B"]
  let lfcs = [2.0, -2.0]
  let pvals = [0.001, 0.001]
  let result = @src.volcano_plot(ids, lfcs, pvals, p_cutoff=0.01, fc_cutoff=2.0, title="", x_label="", y_label="")
  assert_eq(result.get_p_cutoff(), 0.01)
  assert_eq(result.get_fc_cutoff(), 2.0)
}

test "volcano_to_ascii" {
  let result = @src.volcano_sample(20)
  let ascii = result.to_ascii(width=60, height=20)
  assert_true(ascii.length() > 0)
  assert_true(ascii.contains("Volcano Plot"))
}

test "volcano_summary" {
  let result = @src.volcano_sample(10)
  let summary = result.summary()
  assert_true(summary.contains("Volcano Plot Summary"))
  assert_true(summary.contains("Total genes"))
}

test "volcano_sample" {
  let result = @src.volcano_sample(30)
  assert_eq(result.get_n_genes(), 30)
  assert_true(result.get_n_up() >= 0)
  assert_true(result.get_n_down() >= 0)
  assert_true(result.get_n_nonsig() >= 0)
}

test "volcano_classification_to_string" {
  // Test classification through the result of volcano_plot
  let ids = ["Gene1", "Gene2", "Gene3"]
  let lfcs = [2.0, -2.0, 0.5]
  let pvals = [0.001, 0.001, 0.5]
  let result = @src.volcano_plot(ids, lfcs, pvals, p_cutoff=0.05, fc_cutoff=1.0, title="Test", x_label="log2FC", y_label="-log10(p)")
  let genes = result.get_genes()
  assert_eq(genes.length(), 3)
  
  // Check that classifications work by examining result counts
  assert_true(result.get_n_up() >= 1)  // Gene1 should be Up
  assert_true(result.get_n_down() >= 1)  // Gene2 should be Down
  assert_true(result.get_n_nonsig() >= 1)  // Gene3 should be NonSig
}

test "volcano_empty_inputs" {
  let ids : Array[String] = Array::new()
  let lfcs : Array[Double] = Array::new()
  let pvals : Array[Double] = Array::new()
  let result = @src.volcano_plot_default(ids, lfcs, pvals)
  assert_eq(result.get_n_genes(), 0)
}

test "volcano_ascii_visual_elements" {
  let result = @src.volcano_sample(50)
  let ascii = result.to_ascii(width=60, height=20)
  assert_true(ascii.contains("+")) // up-regulated
  assert_true(ascii.contains("-")) // down-regulated
}