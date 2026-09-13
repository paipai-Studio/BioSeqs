///|
/// Tests for the DOSE module (src/dose.mbt).
///
/// All expected values are frozen from Python 3.12 + scipy 1.18.1.
/// See /tmp/gen_dose_expected.py for the generator script.

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

///|
/// Near-equality with absolute + relative tolerance.
fn dose_near(a : Double, b : Double) -> Bool {
  if a.is_nan() && b.is_nan() {
    return true
  }
  let diff = if a - b < 0.0 { -(a - b) } else { a - b }
  let denom = if a.abs() > b.abs() { a.abs() } else { b.abs() }
  diff <= 1.0e-10 || diff <= denom * 1.0e-12
}

// ---------------------------------------------------------------------------
// Fisher exact tests — hand-verified 2x2 tables (frozen from scipy 1.18.1)
// ---------------------------------------------------------------------------

///|
test "dose_fisher_basic_two_sided_and_greater" {
  // a=4, b=6, c=16, d=974
  let r2 = @src.fisher_exact_two_sided(4, 6, 16, 974)
  let rg = @src.fisher_exact_greater(4, 6, 16, 974)
  // scipy: two=2.273022846873337e-05, greater=2.273022846873337e-05 (same, 1-sided extreme)
  assert_true(dose_near(r2.p_value, 2.273022846873337e-05))
  assert_true(dose_near(rg.p_value, 2.273022846873337e-05))
  assert_true(dose_near(r2.odds_ratio, 40.583333333333336))
  assert_true(dose_near(rg.odds_ratio, 40.583333333333336))
}

///|
test "dose_fisher_basic_less" {
  let rl = @src.fisher_exact_less(4, 6, 16, 974)
  // scipy: less=9.999995554761869e-01
  assert_true(dose_near(rl.p_value, 9.999995554761869e-01))
}

///|
test "dose_fisher_medium_table" {
  // a=45, b=255, c=755, d=3945
  let r2 = @src.fisher_exact_two_sided(45, 255, 755, 3945)
  let rg = @src.fisher_exact_greater(45, 255, 755, 3945)
  let rl = @src.fisher_exact_less(45, 255, 755, 3945)
  // scipy: two=6.848326864305199e-01, greater=7.113565655338840e-01,
  //        less=3.474054301567979e-01, OR=0.922088040514219
  assert_true(dose_near(r2.p_value, 6.848326864305199e-01))
  assert_true(dose_near(rg.p_value, 7.113565655338840e-01))
  assert_true(dose_near(rl.p_value, 3.474054301567979e-01))
  assert_true(dose_near(r2.odds_ratio, 0.922088040514219))
}

///|
test "dose_fisher_balanced_table" {
  // a=5, b=5, c=5, d=5 → OR=1.0, two-sided p=1.0
  let r2 = @src.fisher_exact_two_sided(5, 5, 5, 5)
  let rg = @src.fisher_exact_greater(5, 5, 5, 5)
  let rl = @src.fisher_exact_less(5, 5, 5, 5)
  assert_true(dose_near(r2.p_value, 1.0))
  assert_true(dose_near(rg.p_value, 6.718591006516703e-01))
  assert_true(dose_near(rl.p_value, 6.718591006516703e-01))
  assert_true(dose_near(r2.odds_ratio, 1.0))
}

///|
test "dose_fisher_zero_enrichment" {
  // a=0, b=10, c=10, d=20 → OR=0 (then Haldane correction)
  let r2 = @src.fisher_exact_two_sided(0, 10, 10, 20)
  let rg = @src.fisher_exact_greater(0, 10, 10, 20)
  let rl = @src.fisher_exact_less(0, 10, 10, 20)
  // scipy: two=4.283217137132048e-02, greater=1.0, less=3.544463143858929e-02
  assert_true(dose_near(r2.p_value, 4.283217137132048e-02))
  assert_true(dose_near(rg.p_value, 1.0))
  assert_true(dose_near(rl.p_value, 3.544463143858929e-02))
  // OR with Haldane correction: (0.5 * 20.5) / (10.5 * 10.5) = 10.25 / 110.25 ≈ 0.09297
  assert_true(r2.odds_ratio > 0.0 && r2.odds_ratio < 0.2)
}

///|
test "dose_fisher_one_sided_unilateral" {
  // Extreme imbalance: a=3, b=100, c=0, d=1000
  // Here two-sided should equal greater because only one tail has mass.
  let r2 = @src.fisher_exact_two_sided(3, 100, 0, 1000)
  let rg = @src.fisher_exact_greater(3, 100, 0, 1000)
  assert_true(dose_near(r2.p_value, rg.p_value))
  assert_true(dose_near(rg.p_value, 7.928925629120149e-04))
}

// ---------------------------------------------------------------------------
// Hypergeometric enrichment test
// ---------------------------------------------------------------------------

///|
test "dose_hypergeometric_small_enrichment" {
  // k=12 overlap, query=30, term=40, bg=200
  // Table: a=12, b=18, c=28, d=142
  let p_h = @src.dose_hypergeometric_test(12, 30, 40, 200)
  // scipy hypergeom: sf(11, 200, 40, 30) = 5.017999300851878e-03
  assert_true(dose_near(p_h, 5.017999300851878e-03))
}

///|
test "dose_hypergeometric_medium" {
  // k=10, query=50, term=30, bg=500 — strong enrichment
  let p_h = @src.dose_hypergeometric_test(10, 50, 30, 500)
  assert_true(dose_near(p_h, 2.639650131469328e-04))
}

///|
test "dose_hypergeometric_depletion" {
  // k=2, query=100, term=20, bg=200 — depletion (FE=0.2)
  let p_h = @src.dose_hypergeometric_test(2, 100, 20, 200)
  // scipy hypergeom: sf(1, 200, 20, 100) = 9.999914661348041e-01
  assert_true(dose_near(p_h, 9.999914661348041e-01))
}

///|
test "dose_hypergeometric_zero_overlap" {
  let p_h = @src.dose_hypergeometric_test(0, 20, 10, 100)
  assert_true(dose_near(p_h, 1.0))
}

///|
test "dose_hyper_equals_fisher_greater" {
  let p_h = @src.dose_hypergeometric_test(12, 30, 40, 200)
  let p_f = @src.fisher_exact_greater(12, 18, 28, 142).p_value
  assert_true(dose_near(p_h, p_f))
}

// ---------------------------------------------------------------------------
// BH-FDR correction (original order preserved)
// ---------------------------------------------------------------------------

///|
test "dose_bh_fdr_standard" {
  // Input: [0.04, 0.03, 0.02, 0.01, 0.005]
  // Expected (preserving input order): [0.04, 0.0375, 0.033333333333333, 0.025, 0.025]
  let bh = @src.dose_bh_fdr([0.04, 0.03, 0.02, 0.01, 0.005])
  assert_true(bh.length() == 5)
  assert_true(dose_near(bh[0], 0.040000000000000))
  assert_true(dose_near(bh[1], 0.037500000000000))
  assert_true(dose_near(bh[2], 0.033333333333333))
  assert_true(dose_near(bh[3], 0.025000000000000))
  assert_true(dose_near(bh[4], 0.025000000000000))
  // Monotonicity check
  assert_true(bh[0] >= bh[1] - 1.0e-14)
  assert_true(bh[1] >= bh[2] - 1.0e-14)
  assert_true(bh[2] >= bh[3] - 1.0e-14)
  assert_true(bh[3] >= bh[4] - 1.0e-14)
}

///|
test "dose_bh_fdr_all_high_pvalues" {
  // Input: [0.9, 0.95, 0.99]
  // Sorted: [0.9, 0.95, 0.99], m=3
  // rank 3: 0.99*3/3 = 0.99, min_so_far=0.99, out[2]=0.99
  // rank 2: 0.95*3/2 = 1.425 → clamp to 1.0, but min_so_far=min(0.99, 1.0)=0.99, out[1]=0.99
  // rank 1: 0.9*3/1 = 2.7 → clamp to 1.0, but min_so_far=min(0.99, 1.0)=0.99, out[0]=0.99
  let bh = @src.dose_bh_fdr([0.9, 0.95, 0.99])
  assert_true(bh.length() == 3)
  assert_true(dose_near(bh[0], 0.99))
  assert_true(dose_near(bh[1], 0.99))
  assert_true(dose_near(bh[2], 0.99))
}

///|
test "dose_bh_fdr_single_small" {
  let bh = @src.dose_bh_fdr([1.0e-5, 2.0e-5, 3.0e-5, 0.01])
  assert_true(bh.length() == 4)
  assert_true(dose_near(bh[0], 0.000040000000000))
  assert_true(dose_near(bh[1], 0.000040000000000))
  assert_true(dose_near(bh[2], 0.000040000000000))
  assert_true(dose_near(bh[3], 0.010000000000000))
}

///|
test "dose_bh_fdr_single_pvalue" {
  let bh = @src.dose_bh_fdr([0.5])
  assert_true(bh.length() == 1)
  assert_true(dose_near(bh[0], 0.5))
}

///|
test "dose_bh_fdr_empty" {
  let bh = @src.dose_bh_fdr([])
  assert_true(bh.length() == 0)
}

// ---------------------------------------------------------------------------
// Bonferroni correction
// ---------------------------------------------------------------------------

///|
test "dose_bonferroni_standard" {
  let bonf = @src.dose_bonferroni([0.04, 0.03, 0.02, 0.01, 0.005])
  assert_true(bonf.length() == 5)
  assert_true(dose_near(bonf[0], 0.2))
  assert_true(dose_near(bonf[1], 0.15))
  assert_true(dose_near(bonf[2], 0.1))
  assert_true(dose_near(bonf[3], 0.05))
  assert_true(dose_near(bonf[4], 0.025))
}

///|
test "dose_bonferroni_clamped" {
  let bonf = @src.dose_bonferroni([0.9, 0.95, 0.99])
  assert_true(bonf.length() == 3)
  assert_true(dose_near(bonf[0], 1.0))
  assert_true(dose_near(bonf[1], 1.0))
  assert_true(dose_near(bonf[2], 1.0))
}

///|
test "dose_bonferroni_single" {
  let bonf = @src.dose_bonferroni([0.5])
  assert_true(bonf.length() == 1)
  assert_true(dose_near(bonf[0], 0.5))
}

///|
test "dose_bonferroni_empty" {
  let bonf = @src.dose_bonferroni([])
  assert_true(bonf.length() == 0)
}

// ---------------------------------------------------------------------------
// Fold enrichment / odds ratio
// ---------------------------------------------------------------------------

///|
test "dose_fold_enrichment_standard" {
  // k=12, m=30, M=40, N=200 → 12*200/(30*40) = 2400/1200 = 2.0
  assert_true(dose_near(@src.dose_fold_enrichment(12, 30, 40, 200), 2.0))
  // k=10, m=50, M=30, N=500 → 10*500/(50*30) = 5000/1500 = 3.3333...
  assert_true(dose_near(@src.dose_fold_enrichment(10, 50, 30, 500), 3.33333333333333))
}

///|
test "dose_fold_enrichment_zero" {
  assert_true(dose_near(@src.dose_fold_enrichment(0, 20, 10, 100), 0.0))
}

// ---------------------------------------------------------------------------
// Batch enrichment (small tables to keep log_factorial fast)
// ---------------------------------------------------------------------------

///|
/// Build a tiny 3-term database and verify enrich_disease_fisher.
/// Background size = 200 so log_factorial stays fast.
fn mk_tiny_test_case() -> (Array[String], Map[String, Array[String]], Array[String]) {
  let query : Array[String] = ["G01", "G02", "G03", "G04", "G05", "G06", "G07", "G08",
                               "G09", "G10", "G11", "G12", "G13", "G14", "G15"]
  let db : Map[String, Array[String]] = Map([])
  // TERM_X: k=12 overlap (G01..G12), term size 20
  let tx : Array[String] = ["G01","G02","G03","G04","G05","G06","G07","G08",
                             "G09","G10","G11","G12","T01","T02","T03","T04",
                             "T05","T06","T07","T08"]
  db["TERM_X"] = tx
  // TERM_Y: k=5 overlap (G01..G05), term size 25
  let ty : Array[String] = ["G01","G02","G03","G04","G05","T11","T12","T13","T14","T15",
                             "T16","T17","T18","T19","T20","T21","T22","T23","T24","T25",
                             "T26","T27","T28","T29","T30"]
  db["TERM_Y"] = ty
  // TERM_Z: k=0 overlap, term size 10
  let tz : Array[String] = ["S01","S02","S03","S04","S05","S06","S07","S08","S09","S10"]
  db["TERM_Z"] = tz
  // Background: 200 genes including all term genes + query genes
  let bg : Array[String] = Array::new()
  let all_term_keys = db.keys().collect()
  let mut i = 0
  while i < all_term_keys.length() {
    let t = db[all_term_keys[i]]
    let mut j = 0
    while j < t.length() {
      bg.push(t[j])
      j = j + 1
    }
    i = i + 1
  }
  // Add query genes explicitly
  let mut j = 0
  while j < query.length() {
    bg.push(query[j])
    j = j + 1
  }
  // Pad to 200
  i = 1
  while bg.length() < 200 {
    bg.push("B" + i.to_string())
    i = i + 1
  }
  (query, db, bg)
}

///|
test "dose_enrichment_tiny_smoke" {
  let (query, db, bg) = mk_tiny_test_case()
  let result = @src.enrich_disease_fisher(query, db, bg)

  assert_true(result.query_count == 15)
  assert_true(result.background_count == bg.length())
  assert_true(result.n_terms == 3)
  assert_true(result.terms.length() == 3)

  // Verify all 3 terms exist
  let mut found_x = false
  let mut found_y = false
  let mut found_z = false
  let mut i = 0
  while i < result.terms.length() {
    let t = result.terms[i]
    if t.term_id == "TERM_X" { found_x = true }
    if t.term_id == "TERM_Y" { found_y = true }
    if t.term_id == "TERM_Z" { found_z = true }
    i = i + 1
  }
  assert_true(found_x && found_y && found_z)

  // Verify TERM_X has correct features
  i = 0
  while i < result.terms.length() {
    let t = result.terms[i]
    if t.term_id == "TERM_X" {
      assert_true(t.overlap_count == 12)
      assert_true(t.gene_count == 20)
      let expected_fe = 12.0 * 200.0 / (15.0 * 20.0)  // = 8.0
      assert_true(dose_near(t.fold_enrichment, expected_fe))
    }
    i = i + 1
  }
}

///|
test "dose_enrichment_sorted_by_p_then_fe" {
  let (query, db, bg) = mk_tiny_test_case()
  let result = @src.enrich_disease_fisher(query, db, bg)

  // TERM_X (k=12) should be top — smallest p-value
  let top = result.terms[0]
  assert_true(top.term_id == "TERM_X")

  // TERM_Z (k=0) should be bottom — largest p-value
  let bottom = result.terms[result.terms.length() - 1]
  assert_true(bottom.term_id == "TERM_Z")

  // Verify adjusted p-values are monotonic
  let mut ok = true
  let mut i = 1
  while i < result.terms.length() {
    if result.terms[i].adjusted_p_value < result.terms[i - 1].adjusted_p_value - 1.0e-14 {
      ok = false
    }
    i = i + 1
  }
  assert_true(ok)
}

// ---------------------------------------------------------------------------
// DOSEResult accessors
// ---------------------------------------------------------------------------

///|
test "dose_filter_by_adjp_below_cutoff" {
  let (query, db, bg) = mk_tiny_test_case()
  let result = @src.enrich_disease_fisher(query, db, bg)
  let filtered = @src.dose_filter_by_adjp(result, 2.0)
  // cutoff=2.0 picks all 7 terms because all adj_p ≤ 1.0
  assert_true(filtered.terms.length() == 3)
}

///|
test "dose_filter_by_adjp_strict_less" {
  // R convention uses strict < — verify a cutoff of 0 picks nothing
  let (query, db, bg) = mk_tiny_test_case()
  let result = @src.enrich_disease_fisher(query, db, bg)
  let filtered = @src.dose_filter_by_adjp(result, 0.0)
  assert_true(filtered.terms.length() == 0)
}

///|
test "dose_top_terms_basic" {
  let (query, db, bg) = mk_tiny_test_case()
  let result = @src.enrich_disease_fisher(query, db, bg)
  let top = @src.dose_top_terms(result, 2)
  assert_true(top.terms.length() == 2)
}

///|
test "dose_term_ids_accessor" {
  let (query, db, bg) = mk_tiny_test_case()
  let result = @src.enrich_disease_fisher(query, db, bg)
  let ids = @src.dose_term_ids(result)
  assert_true(ids.length() == 3)
}

///|
test "dose_adjusted_pvalues_accessor" {
  let (query, db, bg) = mk_tiny_test_case()
  let result = @src.enrich_disease_fisher(query, db, bg)
  let adj = @src.dose_adjusted_pvalues(result)
  assert_true(adj.length() == 3)
}

///|
test "dose_summarize_non_empty" {
  let (query, db, bg) = mk_tiny_test_case()
  let result = @src.enrich_disease_fisher(query, db, bg)
  let s = @src.dose_summarize(result)
  assert_true(s.length() > 0)
}
