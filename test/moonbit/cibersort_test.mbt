///|
/// Tests for Bioconductor CIBERSORT module - immune cell deconvolution
/// via non-negative least squares (NNLS).

// ---------------------------------------------------------------------------
// Data structure creation and accessors
// ---------------------------------------------------------------------------

test "cib_signature_matrix_creation" {
  let sig = @src.CibSignatureMatrix::new(
    ["G1", "G2", "G3"],
    ["T cells", "B cells"],
    [[1.0, 0.0], [0.5, 0.5], [0.0, 1.0]],
  )
  assert_eq(sig.gene_names().length(), 3)
  assert_eq(sig.cell_types().length(), 2)
  assert_eq(sig.cell_types()[0], "T cells")
  assert_eq(sig.matrix().length(), 3)
  assert_eq(sig.matrix()[0][0], 1.0)
  assert_eq(sig.matrix()[2][1], 1.0)
}

test "cib_mixture_matrix_creation" {
  let mix = @src.CibMixtureMatrix::new(
    ["G1", "G2"],
    ["S1", "S2"],
    [[1.0, 2.0], [3.0, 4.0]],
  )
  assert_eq(mix.gene_names().length(), 2)
  assert_eq(mix.sample_names().length(), 2)
  assert_eq(mix.sample_names()[0], "S1")
  assert_eq(mix.matrix()[1][0], 3.0)
}

test "cib_result_accessors" {
  // Use cib_run to get a real result and test accessors
  let sig = @src.cib_default_signature()
  let mix = @src.cib_sample_mixture()
  let result = @src.cib_run(sig, mix)
  let r0 = result.results()[0]
  assert_eq(r0.sample_id(), "Sample1")
  assert_eq(r0.fractions().length(), 10)
  assert_eq(r0.cell_types().length(), 10)
  assert_true(r0.pearson_r() > 0.0 && r0.pearson_r() <= 1.0)
  assert_true(r0.rmse() >= 0.0)
}

// ---------------------------------------------------------------------------
// Built-in signature and marker genes
// ---------------------------------------------------------------------------

test "cib_cell_type_names_count" {
  let names = @src.cib_cell_type_names()
  assert_eq(names.length(), 10)
  assert_eq(names[0], "B cells")
  assert_eq(names[9], "Mast cells")
}

test "cib_marker_genes_count" {
  let genes = @src.cib_marker_genes()
  // 10 cell types × 4 markers each = 40 genes
  assert_eq(genes.length(), 40)
  assert_eq(genes[0], "CD19") // B cell marker
  assert_eq(genes[4], "CD8A") // T cell CD8 marker
}

test "cib_default_signature_dimensions" {
  let sig = @src.cib_default_signature()
  assert_eq(sig.gene_names().length(), 40)
  assert_eq(sig.cell_types().length(), 10)
  assert_eq(sig.matrix().length(), 40)
  assert_eq(sig.matrix()[0].length(), 10)
  // B cell markers should be highly expressed in B cells (col 0)
  assert_true(sig.matrix()[0][0] > sig.matrix()[0][1])
}

// ---------------------------------------------------------------------------
// NNLS deconvolution correctness
// ---------------------------------------------------------------------------

test "cib_run_single_cell_type_recovery" {
  // Build a simple mixture: pure B cells (column 0 of signature)
  let sig = @src.cib_default_signature()
  let gene_names = sig.gene_names().copy()
  // Mixture = exactly the B-cell column + small noise
  let matrix : Array[Array[Double]] = []
  for gi in 0..<gene_names.length() {
    let val = sig.matrix()[gi][0] * 1.0 + 0.1 // B cells only
    matrix.push([val])
  }
  let mix = @src.CibMixtureMatrix::new(gene_names, ["Pure_B"], matrix)
  let result = @src.cib_run(sig, mix)
  assert_eq(result.sample_names().length(), 1)
  assert_eq(result.cell_types().length(), 10)
  let r0 = result.results()[0]
  // B-cell fraction should be the dominant component
  let b_frac = r0.get_fraction("B cells")
  assert_true(b_frac > 0.5)
  // Pearson correlation should be high (good fit)
  assert_true(r0.pearson_r() > 0.9)
}

test "cib_run_mixed_composition_recovery" {
  // Use the built-in sample mixture generator
  let sig = @src.cib_default_signature()
  let mix = @src.cib_sample_mixture()
  let result = @src.cib_run(sig, mix)
  assert_eq(result.sample_names().length(), 3)
  // For each sample, fractions should sum to ~1.0
  for r in result.results() {
    let mut sum = 0.0
    for f in r.fractions() {
      sum = sum + f
    }
    assert_true((sum - 1.0).abs() < 0.01)
    // All fractions should be non-negative
    for f in r.fractions() {
      assert_true(f >= -1.0e-9)
    }
  }
}

test "cib_run_sample1_dominant_t_cell_cd8" {
  // Sample1 was constructed as 60% T cells CD8 + 30% B cells + 10% NK cells
  let sig = @src.cib_default_signature()
  let mix = @src.cib_sample_mixture()
  let result = @src.cib_run(sig, mix)
  let s1 = result.results()[0]
  let t_cd8 = s1.get_fraction("T cells CD8")
  let b = s1.get_fraction("B cells")
  let nk = s1.get_fraction("NK cells")
  // T cells CD8 should be the largest among the three
  assert_true(t_cd8 > b)
  assert_true(t_cd8 > nk)
}

test "cib_run_sample2_dominant_monocytes" {
  // Sample2: 50% Monocytes + 30% Macrophages M1 + 20% Macrophages M2
  let sig = @src.cib_default_signature()
  let mix = @src.cib_sample_mixture()
  let result = @src.cib_run(sig, mix)
  let s2 = result.results()[1]
  // Find max fraction cell type
  let mut max_idx = 0
  let mut max_val = -1.0
  for i in 0..<s2.fractions().length() {
    if s2.fractions()[i] > max_val {
      max_val = s2.fractions()[i]
      max_idx = i
    }
  }
  // The top cell type should be Monocytes, Macrophages M1, or Macrophages M2
  let top_type = s2.cell_types()[max_idx]
  assert_true(
    top_type == "Monocytes" ||
    top_type == "Macrophages M1" ||
    top_type == "Macrophages M2",
  )
}

test "cib_run_pearson_high_for_clean_data" {
  // For synthetic mixtures constructed from the signature, fit should be good
  let sig = @src.cib_default_signature()
  let mix = @src.cib_sample_mixture()
  let result = @src.cib_run(sig, mix)
  for r in result.results() {
    assert_true(r.pearson_r() > 0.85)
  }
}

test "cib_run_handles_partial_gene_overlap" {
  // Mixture has only a subset of signature genes
  let sig = @src.cib_default_signature()
  // Take only first 10 genes from the mixture
  let sub_genes : Array[String] = []
  let sub_matrix : Array[Array[Double]] = []
  for gi in 0..<10 {
    sub_genes.push(sig.gene_names()[gi])
    sub_matrix.push([sig.matrix()[gi][0]]) // pure B cells
  }
  let mix = @src.CibMixtureMatrix::new(sub_genes, ["Sub_B"], sub_matrix)
  let result = @src.cib_run(sig, mix)
  assert_eq(result.sample_names().length(), 1)
  let r0 = result.results()[0]
  // Should still identify B cells as dominant
  assert_true(r0.get_fraction("B cells") > 0.3)
}

test "cib_run_empty_mixture_genes" {
  // Mixture has no genes in common with the signature -> all zeros
  let sig = @src.cib_default_signature()
  let mix = @src.CibMixtureMatrix::new(
    ["NONEXIST1", "NONEXIST2"],
    ["Empty"],
    [[1.0], [2.0]],
  )
  let result = @src.cib_run(sig, mix)
  assert_eq(result.sample_names().length(), 1)
  let r0 = result.results()[0]
  for f in r0.fractions() {
    assert_true(f >= -1.0e-9)
  }
}

// ---------------------------------------------------------------------------
// to_string formatting
// ---------------------------------------------------------------------------

test "cib_deconvolution_to_string" {
  let sig = @src.cib_default_signature()
  let mix = @src.cib_sample_mixture()
  let result = @src.cib_run(sig, mix)
  let s = result.to_string()
  assert_true(s.contains("CIBERSORT"))
  assert_true(s.contains("Sample1"))
  assert_true(s.contains("Pearson_r"))
}

// ---------------------------------------------------------------------------
// NNLS edge cases
// ---------------------------------------------------------------------------

test "cib_get_fraction_unknown_type" {
  // Use cib_run to get a real result and test get_fraction
  let sig = @src.cib_default_signature()
  let mix = @src.cib_sample_mixture()
  let result = @src.cib_run(sig, mix)
  let r0 = result.results()[0]
  // Unknown cell type returns 0.0
  assert_true(r0.get_fraction("Unknown").abs() < 1.0e-12)
  // Known cell type returns a non-negative fraction
  let b_frac = r0.get_fraction("B cells")
  assert_true(b_frac >= 0.0)
}

test "cib_run_with_custom_tolerance" {
  let sig = @src.cib_default_signature()
  let mix = @src.cib_sample_mixture()
  let result = @src.cib_run(sig, mix, max_iters=500, tol=0.0001)
  assert_eq(result.sample_names().length(), 3)
  // With fewer iterations, results should still be reasonable
  for r in result.results() {
    let mut sum = 0.0
    for f in r.fractions() {
      sum = sum + f
    }
    assert_true((sum - 1.0).abs() < 0.05)
  }
}
