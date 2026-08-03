///|
/// Tests for Bioconductor stageR - two-stage hypothesis testing.

///|
test "stage_r_config_default" {
  let config = @src.StageRConfig::new()
  assert_eq(config.alpha, 0.05)
  assert_eq(config.pScreenAdjusted, false)
}

///|
test "stage_r_bh_adjust_basic" {
  let pvals = [0.01, 0.04, 0.03, 0.20, 0.05]
  let adj = @src.stage_r_bh_adjust(pvals)
  assert_eq(adj.length(), 5)
  // All adjusted values should be >= original and <= 1.
  let mut i = 0
  while i < pvals.length() {
    assert_true(adj[i] >= pvals[i] - 1.0e-10)
    assert_true(adj[i] <= 1.0 + 1.0e-10)
    i = i + 1
  }
}

///|
test "stage_r_bh_adjust_sorted" {
  let pvals = [0.001, 0.01, 0.02, 0.03, 0.04, 0.05, 0.1, 0.2, 0.5, 0.8]
  let adj = @src.stage_r_bh_adjust(pvals)
  // Smallest p-value should get the smallest adjusted value.
  assert_true(adj[0] <= adj[1])
  assert_true(adj[1] <= adj[2])
}

///|
test "stage_r_bh_adjust_empty" {
  let adj = @src.stage_r_bh_adjust([])
  assert_eq(adj.length(), 0)
}

///|
test "stage_r_bh_adjust_single" {
  let adj = @src.stage_r_bh_adjust([0.05])
  assert_eq(adj.length(), 1)
  assert_true((adj[0] - 0.05).abs() < 1.0e-10)
}

///|
test "stage_r_bh_adjust_all_significant" {
  let pvals = [0.001, 0.002, 0.003]
  let adj = @src.stage_r_bh_adjust(pvals)
  // All should be significant at alpha=0.05.
  let mut i = 0
  while i < adj.length() {
    assert_true(adj[i] <= 0.05 + 1.0e-10)
    i = i + 1
  }
}

///|
test "stage_r_bh_adjust_none_significant" {
  let pvals = [0.5, 0.6, 0.7, 0.8, 0.9]
  let adj = @src.stage_r_bh_adjust(pvals)
  // None should be significant at alpha=0.05.
  let mut i = 0
  while i < adj.length() {
    assert_true(adj[i] > 0.05)
    i = i + 1
  }
}

///|
test "stage_r_simes_basic" {
  // Simes: min(3*p_(i)/i) for sorted p-values
  let pvals = [0.01, 0.02, 0.03]
  let simes = @src.simes_aggregate(pvals)
  // Sorted: 0.01, 0.02, 0.03
  // Values: 3*0.01/1=0.03, 3*0.02/2=0.03, 3*0.03/3=0.03
  // Min = 0.03
  assert_true((simes - 0.03).abs() < 1.0e-10)
}

///|
test "stage_r_simes_unsorted" {
  let pvals = [0.03, 0.01, 0.02]
  let simes = @src.simes_aggregate(pvals)
  // Should be the same as sorted.
  assert_true((simes - 0.03).abs() < 1.0e-10)
}

///|
test "stage_r_simes_empty" {
  let simes = @src.simes_aggregate([])
  assert_true((simes - 1.0).abs() < 1.0e-10)
}

///|
test "stage_r_simes_single" {
  let simes = @src.simes_aggregate([0.05])
  assert_true((simes - 0.05).abs() < 1.0e-10)
}

///|
test "stage_r_simes_caps_at_one" {
  let simes = @src.simes_aggregate([0.5, 0.6, 0.7])
  // Values: 3*0.5/1=1.5, 3*0.6/2=0.9, 3*0.7/3=0.7
  // Min = 0.7, but capped at 1.0
  assert_true((simes - 0.7).abs() < 1.0e-10)
}

///|
test "stage_r_adjustment_holm" {
  let (pScreen, pConfirmation) = @src.stage_r_sample()
  let config = @src.StageRConfig::new(alpha=0.05, pScreenAdjusted=false, stageMethod=@src.stage_r_method_holm())
  let result = @src.stage_wise_adjustment(pScreen, pConfirmation, config)
  assert_eq(result.nGenes, 10)
  assert_eq(result.nHypotheses, 3)
  assert_eq(result.pAdjScreen.length(), 10)
  assert_eq(result.pAdjConfirmation.length(), 10)
}

///|
test "stage_r_adjustment_none" {
  let (pScreen, pConfirmation) = @src.stage_r_sample()
  let config = @src.StageRConfig::new(alpha=0.05, pScreenAdjusted=false, stageMethod=@src.stage_r_method_none())
  let result = @src.stage_wise_adjustment(pScreen, pConfirmation, config)
  assert_eq(result.nGenes, 10)
}

///|
test "stage_r_adjustment_dte" {
  let (pScreen, pConfirmation) = @src.stage_r_sample()
  let config = @src.StageRConfig::new(alpha=0.05, pScreenAdjusted=false, stageMethod=@src.stage_r_method_dte())
  let result = @src.stage_wise_adjustment(pScreen, pConfirmation, config)
  assert_eq(result.nGenes, 10)
  assert_eq(result.nHypotheses, 3)
}

///|
test "stage_r_adjustment_dtu" {
  let (pScreen, pConfirmation) = @src.stage_r_sample()
  let config = @src.StageRConfig::new(alpha=0.05, pScreenAdjusted=false, stageMethod=@src.stage_r_method_dtu())
  let result = @src.stage_wise_adjustment(pScreen, pConfirmation, config)
  assert_eq(result.nGenes, 10)
  assert_eq(result.nHypotheses, 3)
}

///|
test "stage_r_alpha_adjusted" {
  let (pScreen, pConfirmation) = @src.stage_r_sample()
  let config = @src.StageRConfig::new(alpha=0.05, pScreenAdjusted=false, stageMethod=@src.stage_r_method_holm())
  let result = @src.stage_wise_adjustment(pScreen, pConfirmation, config)
  // alphaAdjusted = (R/G) * alpha
  let expected = result.nSignificantGenes.to_double() / result.nGenes.to_double() * config.alpha
  assert_true((result.alphaAdjusted - expected).abs() < 1.0e-10)
}

///|
test "stage_r_significant_genes" {
  let (pScreen, pConfirmation) = @src.stage_r_sample()
  let config = @src.StageRConfig::new(alpha=0.05, pScreenAdjusted=false, stageMethod=@src.stage_r_method_holm())
  let result = @src.stage_wise_adjustment(pScreen, pConfirmation, config)
  let sigGenes = @src.get_significant_genes(result)
  assert_eq(sigGenes.length(), result.nSignificantGenes)
  // All indices should be valid.
  let mut i = 0
  while i < sigGenes.length() {
    assert_true(sigGenes[i] >= 0)
    assert_true(sigGenes[i] < result.nGenes)
    i = i + 1
  }
}

///|
test "stage_r_significant_hypotheses" {
  let (pScreen, pConfirmation) = @src.stage_r_sample()
  let config = @src.StageRConfig::new(alpha=0.05, pScreenAdjusted=false, stageMethod=@src.stage_r_method_holm())
  let result = @src.stage_wise_adjustment(pScreen, pConfirmation, config)
  let sigHyps = @src.get_significant_hypotheses(result)
  // All pairs should have valid indices.
  let mut i = 0
  while i < sigHyps.length() {
    let (g, h) = sigHyps[i]
    assert_true(g >= 0 && g < result.nGenes)
    assert_true(h >= 0 && h < result.nHypotheses)
    i = i + 1
  }
}

///|
test "stage_r_get_results" {
  let (pScreen, pConfirmation) = @src.stage_r_sample()
  let config = @src.StageRConfig::new(alpha=0.05, pScreenAdjusted=false, stageMethod=@src.stage_r_method_holm())
  let result = @src.stage_wise_adjustment(pScreen, pConfirmation, config)
  let mat = @src.get_results(result)
  assert_eq(mat.length(), result.nGenes)
  // Each row should have nHypotheses + 1 columns.
  let mut i = 0
  while i < mat.length() {
    assert_eq(mat[i].length(), result.nHypotheses + 1)
    i = i + 1
  }
}

///|
test "stage_r_non_significant_genes_na" {
  let (pScreen, pConfirmation) = @src.stage_r_sample()
  let config = @src.StageRConfig::new(alpha=0.001, pScreenAdjusted=false, stageMethod=@src.stage_r_method_holm())
  let result = @src.stage_wise_adjustment(pScreen, pConfirmation, config)
  // Non-significant genes should have -1.0 in confirmation p-values.
  let mut i = 0
  while i < result.nGenes {
    if result.pAdjScreen[i] > config.alpha {
      let mut j = 0
      while j < result.nHypotheses {
        assert_eq(result.pAdjConfirmation[i][j], -1.0)
        j = j + 1
      }
    }
    i = i + 1
  }
}

///|
test "stage_r_pscreen_adjusted" {
  let (pScreen, pConfirmation) = @src.stage_r_sample()
  // Pre-adjust the screening p-values.
  let adjScreen = @src.stage_r_bh_adjust(pScreen)
  let config = @src.StageRConfig::new(alpha=0.05, pScreenAdjusted=true, stageMethod=@src.stage_r_method_holm())
  let result = @src.stage_wise_adjustment(adjScreen, pConfirmation, config)
  // pAdjScreen should be the same as input.
  let mut i = 0
  while i < pScreen.length() {
    assert_true((result.pAdjScreen[i] - adjScreen[i]).abs() < 1.0e-10)
    i = i + 1
  }
}

///|
test "stage_r_confirmation_rescaled" {
  let (pScreen, pConfirmation) = @src.stage_r_sample()
  let config = @src.StageRConfig::new(alpha=0.05, pScreenAdjusted=false, stageMethod=@src.stage_r_method_holm())
  let result = @src.stage_wise_adjustment(pScreen, pConfirmation, config)
  // For significant genes, confirmation p-values should be rescaled.
  let mut i = 0
  while i < result.nGenes {
    if result.pAdjScreen[i] <= config.alpha {
      let mut j = 0
      while j < result.nHypotheses {
        // Rescaled values should be >= 0.
        assert_true(result.pAdjConfirmation[i][j] >= 0.0)
        // And <= 1.
        assert_true(result.pAdjConfirmation[i][j] <= 1.0)
        j = j + 1
      }
    }
    i = i + 1
  }
}

///|
test "stage_r_summary" {
  let (pScreen, pConfirmation) = @src.stage_r_sample()
  let config = @src.StageRConfig::new()
  let result = @src.stage_wise_adjustment(pScreen, pConfirmation, config)
  let summary = @src.stage_r_summary(result)
  assert_true(summary.length() > 0)
}

///|
test "stage_r_sample_data" {
  let (pScreen, pConfirmation) = @src.stage_r_sample()
  assert_eq(pScreen.length(), 10)
  assert_eq(pConfirmation.length(), 10)
  assert_eq(pConfirmation[0].length(), 3)
}
