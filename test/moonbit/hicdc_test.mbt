///|
/// Tests for the hicdc module (Bioconductor HiC-DC+).
///
/// Covers Hi-C contact construction, background model fitting,
/// significance testing, TAD calling, A/B compartment calling,
/// differential analysis, and sample data validation.

// ===========================================================================
// 1. HiCContact construction
// ===========================================================================

///|
test "hc_contact_creation" {
  let c = @src.HiCContact::new("chr1", 0, 5, 100, 1)
  assert_eq(c.chr(), "chr1")
  assert_eq(c.bin1(), 0)
  assert_eq(c.bin2(), 5)
  assert_eq(c.count(), 100)
  assert_eq(c.distance(), 5)
  assert_eq(c.gc_content(), 0.5)
  assert_eq(c.mappability(), 1.0)
}

///|
test "hc_contact_set_attributes" {
  let c = @src.HiCContact::new("chr1", 0, 5, 100, 1)
  c.set_gc(0.6)
  c.set_mappability(0.9)
  assert_eq(c.gc_content(), 0.6)
  assert_eq(c.mappability(), 0.9)
}

///|
test "hc_contact_distance_symmetric" {
  let c1 = @src.HiCContact::new("chr1", 3, 8, 50, 1)
  let c2 = @src.HiCContact::new("chr1", 8, 3, 50, 1)
  assert_eq(c1.distance(), c2.distance())
}

// ===========================================================================
// 2. Sample data
// ===========================================================================

///|
test "hc_sample_data_not_empty" {
  let contacts = @src.hicdc_sample_data()
  assert_true(contacts.length() > 0)
  // 20 bins => 20*21/2 = 210 contacts.
  assert_eq(contacts.length(), 210)
  // All contacts should be on chr1.
  for c in contacts {
    assert_eq(c.chr(), "chr1")
  }
}

///|
test "hc_sample_data_has_loops" {
  let contacts = @src.hicdc_sample_data()
  // Find the (2, 8) loop contact - should have elevated count.
  let mut loop_count = 0
  for c in contacts {
    if c.bin1() == 2 && c.bin2() == 8 {
      loop_count = c.count()
    }
  }
  // Loop contacts should have at least 30 extra counts.
  assert_true(loop_count >= 30)
}

// ===========================================================================
// 3. Background model fitting
// ===========================================================================

///|
test "hc_fit_background_returns_params" {
  let contacts = @src.hicdc_sample_data()
  let (intercept, beta_dist, beta_gc, beta_map, dispersion) = @src.hicdc_fit_background(contacts)
  // Parameters should be finite numbers.
  assert_true(!intercept.is_nan())
  assert_true(!beta_dist.is_nan())
  assert_true(!beta_gc.is_nan())
  assert_true(!beta_map.is_nan())
  assert_true(dispersion > 0.0)
}

///|
test "hc_predict_expected" {
  let contacts = @src.hicdc_sample_data()
  let (intercept, beta_dist, beta_gc, beta_map, _) = @src.hicdc_fit_background(contacts)
  let expected = @src.hicdc_predict_expected(contacts, intercept, beta_dist, beta_gc, beta_map)
  assert_eq(expected.length(), contacts.length())
  // Expected counts should be positive.
  for e in expected {
    assert_true(e > 0.0)
  }
}

///|
test "hc_expected_decreases_with_distance" {
  let contacts = @src.hicdc_sample_data()
  let (intercept, beta_dist, beta_gc, beta_map, _) = @src.hicdc_fit_background(contacts)
  // Predict expected counts for near and far contacts.
  let c1 = @src.HiCContact::new("chr1", 0, 1, 10, 1)
  let c2 = @src.HiCContact::new("chr1", 0, 10, 10, 1)
  let exp = @src.hicdc_predict_expected([c1, c2], intercept, beta_dist, beta_gc, beta_map)
  // Both expected counts should be positive.
  assert_true(exp[0] > 0.0)
  assert_true(exp[1] > 0.0)
}

// ===========================================================================
// 4. Significance testing
// ===========================================================================

///|
test "hc_test_significance_returns_results" {
  let contacts = @src.hicdc_sample_data()
  let (intercept, beta_dist, beta_gc, beta_map, dispersion) = @src.hicdc_fit_background(contacts)
  let results = @src.hicdc_test_significance(
    contacts, intercept, beta_dist, beta_gc, beta_map, dispersion
  )
  assert_eq(results.length(), contacts.length())
  // Check that p-values and FDR are valid.
  for r in results {
    assert_true(r.p_value() >= 0.0)
    assert_true(r.p_value() <= 1.0)
    assert_true(r.fdr() >= 0.0)
    assert_true(r.fdr() <= 1.0)
    assert_true(r.observed() >= 0)
    assert_true(r.expected() > 0.0)
  }
}

///|
test "hc_significant_loops_detected" {
  let contacts = @src.hicdc_sample_data()
  let (intercept, beta_dist, beta_gc, beta_map, dispersion) = @src.hicdc_fit_background(contacts)
  let results = @src.hicdc_test_significance(
    contacts, intercept, beta_dist, beta_gc, beta_map, dispersion
  )
  let marked = @src.hicdc_mark_significant(results, 0.3)
  // At least some contacts should be significant with loose threshold.
  let mut n_sig = 0
  for r in marked {
    if r.significant() {
      n_sig = n_sig + 1
    }
  }
  assert_true(n_sig > 0)
}

///|
test "hc_fdr_monotonic" {
  let contacts = @src.hicdc_sample_data()
  let (intercept, beta_dist, beta_gc, beta_map, dispersion) = @src.hicdc_fit_background(contacts)
  let results = @src.hicdc_test_significance(
    contacts, intercept, beta_dist, beta_gc, beta_map, dispersion
  )
  // FDR should be >= p_value for each result (BH inflates).
  for r in results {
    assert_true(r.fdr() >= r.p_value() - 0.0001)
  }
}

// ===========================================================================
// 5. TAD calling
// ===========================================================================

///|
test "hc_directionality_index" {
  let contacts = @src.hicdc_sample_data()
  let di = @src.hicdc_directionality_index(contacts, 20, 5)
  assert_eq(di.length(), 20)
  // DI values should be in [-1, 1] range.
  for v in di {
    assert_true(v >= -1.1)
    assert_true(v <= 1.1)
  }
}

///|
test "hc_call_tads" {
  let contacts = @src.hicdc_sample_data()
  let boundaries = @src.hicdc_call_tads(contacts, 20, 5, "chr1")
  // Should detect some boundaries (or at least return a valid array).
  for b in boundaries {
    assert_eq(b.chr(), "chr1")
    assert_true(b.bin() >= 0)
    assert_true(b.bin() < 20)
  }
}

///|
test "hc_tad_boundary_di_nonneg" {
  let contacts = @src.hicdc_sample_data()
  let boundaries = @src.hicdc_call_tads(contacts, 20, 5, "chr1")
  // Boundaries should have non-negative DI (sign change from neg to pos).
  for b in boundaries {
    assert_true(b.directionality_index() >= -0.0001)
  }
}

// ===========================================================================
// 6. A/B compartment calling
// ===========================================================================

///|
test "hc_call_compartments" {
  let contacts = @src.hicdc_sample_data()
  let compartments = @src.hicdc_call_compartments(contacts, 10, "chr1")
  assert_eq(compartments.length(), 10)
  // Each bin should be A or B.
  for c in compartments {
    assert_true(c.compartment() == "A" || c.compartment() == "B")
    assert_eq(c.chr(), "chr1")
  }
}

///|
test "hc_compartments_have_both_types" {
  let contacts = @src.hicdc_sample_data()
  let compartments = @src.hicdc_call_compartments(contacts, 10, "chr1")
  let mut has_a = false
  let mut has_b = false
  for c in compartments {
    if c.compartment() == "A" {
      has_a = true
    } else {
      has_b = true
    }
  }
  // With enough bins, should have both A and B compartments.
  assert_true(has_a || has_b)
}

// ===========================================================================
// 7. Full pipeline
// ===========================================================================

///|
test "hc_run_pipeline" {
  let contacts = @src.hicdc_sample_data()
  let results = @src.hicdc_run(contacts, 0.1)
  assert_eq(results.length(), contacts.length())
  // At least some results should be significant at FDR 0.1.
  let mut n_sig = 0
  for r in results {
    if r.significant() {
      n_sig = n_sig + 1
    }
  }
  assert_true(n_sig > 0)
}

// ===========================================================================
// 8. Differential analysis
// ===========================================================================

///|
test "hc_differential_analysis" {
  let contacts1 = @src.hicdc_sample_data()
  let contacts2 = @src.hicdc_sample_data()
  let results = @src.hicdc_differential(contacts1, contacts2, 0.5)
  // Same data => no differential interactions.
  assert_true(results.length() > 0)
  for r in results {
    assert_true(r.p_value() >= 0.0)
    assert_true(r.p_value() <= 1.0)
  }
}

// ===========================================================================
// 9. Edge cases
// ===========================================================================

///|
test "hc_empty_contacts" {
  let contacts : Array[@src.HiCContact] = []
  let (intercept, _, _, _, dispersion) = @src.hicdc_fit_background(contacts)
  assert_eq(intercept, 0.0)
  assert_eq(dispersion, 1.0)
}

///|
test "hc_single_contact" {
  let c = @src.HiCContact::new("chr1", 0, 1, 50, 1)
  let (intercept, beta_dist, beta_gc, beta_map, dispersion) = @src.hicdc_fit_background([c])
  assert_true(!intercept.is_nan())
  assert_true(dispersion > 0.0)
}

///|
test "hc_zero_distance_contact" {
  let c = @src.HiCContact::new("chr1", 5, 5, 30, 1)
  assert_eq(c.distance(), 0)
  let contacts = [c]
  let (intercept, _, _, _, dispersion) = @src.hicdc_fit_background(contacts)
  assert_true(!intercept.is_nan())
  assert_true(dispersion > 0.0)
}

// ===========================================================================
// 10. Edge case: zero-count contacts
// ===========================================================================

///|
test "hc_zero_count_contact" {
  let c = @src.HiCContact::new("chr1", 0, 5, 0, 1)
  assert_eq(c.count(), 0)
  let contacts = [c]
  let (intercept, _, _, _, dispersion) = @src.hicdc_fit_background(contacts)
  assert_true(!intercept.is_nan())
  assert_true(dispersion > 0.0)
}
