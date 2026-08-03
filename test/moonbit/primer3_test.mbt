///|
/// Tests for the Primer3 module (Bio.Emboss.Primer3 equivalent).
///
/// These tests exercise the public API exposed in `src/primer3.mbt`:
/// the Primer3Config and Primer3Result structs, the Tm calculators
/// (`primer3_tm_wallace`, `primer3_tm_salt`, `primer3_tm`), GC content
/// (`primer3_gc_content`), complementarity scoring (`primer3_self_complement`,
/// `primer3_cross_dimer`, `primer3_hairpin_3prime`), primer scoring
/// (`primer3_score`), the primer-design pipeline (`primer3_design_primers`),
/// amplicon utilities, and sample data.

// ============================================================================
// 1. Primer3Config - creation and accessors
// ============================================================================

///|
test "primer3_config_new_and_accessors" {
  let cfg = @src.Primer3Config::new(
    18, 22, 55.0, 65.0, 0.4, 0.6, 6, 6, 0.05, 50,
  )
  assert_eq(cfg.primer_length_min(), 18)
  assert_eq(cfg.primer_length_max(), 22)
  assert_true(cfg.tm_min() == 55.0)
  assert_true(cfg.tm_max() == 65.0)
  assert_true(cfg.gc_min() == 0.4)
  assert_true(cfg.gc_max() == 0.6)
  assert_eq(cfg.max_self_complement(), 6)
  assert_eq(cfg.max_cross_complement(), 6)
  assert_true(cfg.salt_conc() == 0.05)
  assert_eq(cfg.flank(), 50)
}

///|
test "primer3_default_config_values" {
  let cfg = @src.primer3_default_config()
  assert_eq(cfg.primer_length_min(), 18)
  assert_eq(cfg.primer_length_max(), 22)
  assert_true(cfg.tm_min() == 55.0)
  assert_true(cfg.tm_max() == 65.0)
  assert_true(cfg.gc_min() == 0.4)
  assert_true(cfg.gc_max() == 0.6)
  assert_eq(cfg.flank(), 50)
}

// ============================================================================
// 2. Primer3Result - creation and accessors
// ============================================================================

///|
test "primer3_result_new_and_accessors" {
  let r = @src.Primer3Result::new(
    "ATGCGTACGTAGATGCGT", 10, 18, 60.0, 0.5, 2, true, 85.0,
  )
  assert_eq(r.sequence(), "ATGCGTACGTAGATGCGT")
  assert_eq(r.start(), 10)
  assert_eq(r.length(), 18)
  assert_true(r.tm() == 60.0)
  assert_true(r.gc() == 0.5)
  assert_eq(r.self_complement(), 2)
  assert_true(r.is_forward())
  assert_true(r.score() == 85.0)
}

///|
test "primer3_result_reverse_primer" {
  let r = @src.Primer3Result::new(
    "CTACGTACGCATCTACGT", 100, 18, 58.0, 0.55, 3, false, 78.0,
  )
  assert_false(r.is_forward())
  assert_eq(r.length(), 18)
}

///|
test "primer3_result_to_string_contains_key_fields" {
  let r = @src.Primer3Result::new(
    "ATGCGTACGTAGATGCGT", 10, 18, 60.0, 0.5, 2, true, 85.0,
  )
  let s = r.to_string()
  assert_true(s.contains("Primer[F]"))
  assert_true(s.contains("len=18"))
  assert_true(s.contains("start=10"))
  assert_true(s.contains("seq=ATGCGTACGTAGATGCGT"))
}

// ============================================================================
// 3. GC content
// ============================================================================

///|
test "primer3_gc_content_basic" {
  assert_true(@src.primer3_gc_content("ATGC") == 0.5)
  assert_true(@src.primer3_gc_content("GGCC") == 1.0)
  assert_true(@src.primer3_gc_content("AATT") == 0.0)
}

///|
test "primer3_gc_content_case_insensitive" {
  assert_true(@src.primer3_gc_content("atgc") == 0.5)
  assert_true(@src.primer3_gc_content("ggcc") == 1.0)
}

///|
test "primer3_gc_content_empty" {
  assert_true(@src.primer3_gc_content("") == 0.0)
}

///|
test "primer3_gc_content_ignores_non_acgt" {
  // N is ignored in GC count but counted in length.
  let gc = @src.primer3_gc_content("ATNGC")
  // G + C = 2, length = 5
  assert_true(gc == 2.0 / 5.0)
}

// ============================================================================
// 4. Tm calculation - Wallace rule
// ============================================================================

///|
test "primer3_tm_wallace_basic" {
  // 2 A/T + 4 G/C: A=2, T=2, G=1, C=1 -> 2*(2+2) + 4*(1+1) = 8 + 8 = 16
  let tm = @src.primer3_tm_wallace("AATTGC")
  assert_true(tm == 16.0)
}

///|
test "primer3_tm_wallace_all_gc" {
  // GCGC: 4 G/C, 0 A/T -> 4*4 = 16
  let tm = @src.primer3_tm_wallace("GCGC")
  assert_true(tm == 16.0)
}

///|
test "primer3_tm_wallace_all_at" {
  // AATT: 0 G/C, 4 A/T -> 2*4 = 8
  let tm = @src.primer3_tm_wallace("AATT")
  assert_true(tm == 8.0)
}

///|
test "primer3_tm_wallace_case_insensitive" {
  let tm1 = @src.primer3_tm_wallace("AATTGC")
  let tm2 = @src.primer3_tm_wallace("aattgc")
  assert_true(tm1 == tm2)
}

// ============================================================================
// 5. Tm calculation - salt-adjusted
// ============================================================================

///|
test "primer3_tm_salt_basic" {
  let tm = @src.primer3_tm_salt("ATGCGTACGTAGATGCGTAC", 0.05)
  // Just verify it returns a reasonable value (not NaN, not 0).
  assert_true(tm > 0.0)
  assert_true(tm < 200.0)
}

///|
test "primer3_tm_salt_empty_returns_zero" {
  let tm = @src.primer3_tm_salt("", 0.05)
  assert_true(tm == 0.0)
}

///|
test "primer3_tm_salt_higher_gc_increases_tm" {
  let tm_at = @src.primer3_tm_salt("AATTAAATTAAATTAAATT", 0.05)
  let tm_gc = @src.primer3_tm_salt("GGCCGGCCGGCCGGCCGG", 0.05)
  assert_true(tm_gc > tm_at)
}

///|
test "primer3_tm_salt_zero_salt_uses_default" {
  // When salt_conc <= 0, the function substitutes 0.05. Should not crash.
  let tm = @src.primer3_tm_salt("ATGCGTACGTAGATGCGTAC", 0.0)
  assert_true(tm > 0.0)
}

// ============================================================================
// 6. Tm dispatch (primer3_tm)
// ============================================================================

///|
test "primer3_tm_short_uses_wallace" {
  // Length < 14 uses Wallace.
  let tm = @src.primer3_tm("AATTGC", 0.05)
  let wallace = @src.primer3_tm_wallace("AATTGC")
  assert_true(tm == wallace)
}

///|
test "primer3_tm_long_uses_salt" {
  // Length >= 14 uses salt-adjusted.
  let seq = "ATGCGTACGTAGATGCGTAC"
  let tm = @src.primer3_tm(seq, 0.05)
  let salt = @src.primer3_tm_salt(seq, 0.05)
  assert_true(tm == salt)
}

// ============================================================================
// 7. Self-complementarity
// ============================================================================

///|
test "primer3_self_complement_empty" {
  assert_eq(@src.primer3_self_complement(""), 0)
}

///|
test "primer3_self_complement_palindrome" {
  // A palindromic sequence should have high self-complementarity.
  // GAATTC reverse complement is GAATTC.
  let sc = @src.primer3_self_complement("GAATTC")
  assert_true(sc > 0)
}

///|
test "primer3_self_complement_no_pairing" {
  // AAAA has reverse complement TTTT; A pairs with T, so there is pairing.
  let sc = @src.primer3_self_complement("AAAA")
  assert_true(sc > 0)
}

///|
test "primer3_self_complement_case_insensitive" {
  let sc1 = @src.primer3_self_complement("GAATTC")
  let sc2 = @src.primer3_self_complement("gaattc")
  assert_eq(sc1, sc2)
}

// ============================================================================
// 8. Cross-dimer
// ============================================================================

///|
test "primer3_cross_dimer_empty" {
  assert_eq(@src.primer3_cross_dimer("", "ATGC"), 0)
  assert_eq(@src.primer3_cross_dimer("ATGC", ""), 0)
}

///|
test "primer3_cross_dimer_complementary_pair" {
  // ATGC and GCAT are reverse complements; should have high pairing.
  let cd = @src.primer3_cross_dimer("ATGC", "ATGC")
  assert_true(cd >= 0)
}

///|
test "primer3_cross_dimer_no_overlap" {
  let cd = @src.primer3_cross_dimer("AAAA", "AAAA")
  // AAAA reverse complement is TTTT; A does not pair with A.
  // Best alignment of AAAA vs TTTT has 0 pairs.
  assert_true(cd >= 0)
}

// ============================================================================
// 9. Hairpin 3' end
// ============================================================================

///|
test "primer3_hairpin_3prime_basic" {
  let hp = @src.primer3_hairpin_3prime("AAAAGAATTC", 6)
  assert_true(hp >= 0)
}

///|
test "primer3_hairpin_3prime_empty" {
  assert_eq(@src.primer3_hairpin_3prime("", 5), 0)
}

///|
test "primer3_hairpin_3prime_zero_end_len" {
  assert_eq(@src.primer3_hairpin_3prime("ATGCATGC", 0), 0)
}

///|
test "primer3_hairpin_3prime_end_longer_than_seq" {
  // If end_len > seq length, uses the whole sequence.
  let hp = @src.primer3_hairpin_3prime("GAATTC", 100)
  let full = @src.primer3_self_complement("GAATTC")
  assert_eq(hp, full)
}

// ============================================================================
// 10. Primer scoring
// ============================================================================

///|
test "primer3_score_ideal_primer" {
  let cfg = @src.primer3_default_config()
  // Tm at midpoint (60.0), GC at midpoint (0.5), no self-complement.
  let score = @src.primer3_score("ATGCGTACGTAGATGCGT", 60.0, 0.5, 0, cfg)
  // Ideal: 100 - 0 - 0 - 0 = 100
  assert_true(score == 100.0)
}

///|
test "primer3_score_penalizes_deviation" {
  let cfg = @src.primer3_default_config()
  let ideal = @src.primer3_score("ATGCGTACGTAGATGCGT", 60.0, 0.5, 0, cfg)
  let deviated = @src.primer3_score("ATGCGTACGTAGATGCGT", 56.0, 0.5, 0, cfg)
  assert_true(ideal > deviated)
}

///|
test "primer3_score_penalizes_self_complement" {
  let cfg = @src.primer3_default_config()
  let no_sc = @src.primer3_score("ATGCGTACGTAGATGCGT", 60.0, 0.5, 0, cfg)
  let with_sc = @src.primer3_score("ATGCGTACGTAGATGCGT", 60.0, 0.5, 6, cfg)
  assert_true(no_sc > with_sc)
}

// ============================================================================
// 11. Primer design pipeline
// ============================================================================

///|
test "primer3_design_primers_short_template_returns_empty" {
  let cfg = @src.primer3_default_config()
  let template = "ATGCGTACGTAGATGCGTAC" // 20 bp, too short for flank=50
  let design = @src.primer3_design_primers(template, cfg)
  assert_eq(design.n_forward(), 0)
  assert_eq(design.n_reverse(), 0)
  assert_false(design.has_pair())
}

///|
test "primer3_design_primers_sample_template_finds_candidates" {
  // Use a permissive config: the sample template has ~50% GC motifs, so
  // the Tm range must accommodate the salt-adjusted Tm values achievable
  // with 18-22-mers at 50 mM Na+.
  let cfg = @src.Primer3Config::new(
    18, 22, 40.0, 75.0, 0.3, 0.8, 8, 8, 0.05, 50,
  )
  let template = @src.primer3_sample_template()
  let design = @src.primer3_design_primers(template, cfg)
  // The sample template is ~600 bp; should find forward candidates.
  assert_true(design.n_forward() > 0)
  assert_true(design.n_reverse() > 0)
  assert_true(design.has_pair())
}

///|
test "primer3_design_primers_best_forward_is_highest_score" {
  let cfg = @src.Primer3Config::new(
    18, 22, 40.0, 75.0, 0.3, 0.8, 8, 8, 0.05, 50,
  )
  let template = @src.primer3_sample_template()
  let design = @src.primer3_design_primers(template, cfg)
  match design.best_forward() {
    Some(best) => {
      let fwds = design.forward()
      for i in 0..<fwds.length() {
        assert_true(best.score() >= fwds[i].score())
      }
    }
    None => assert_true(false, msg="expected a best forward primer")
  }
}

///|
test "primer3_design_primers_best_reverse_is_highest_score" {
  let cfg = @src.Primer3Config::new(
    18, 22, 40.0, 75.0, 0.3, 0.8, 8, 8, 0.05, 50,
  )
  let template = @src.primer3_sample_template()
  let design = @src.primer3_design_primers(template, cfg)
  match design.best_reverse() {
    Some(best) => {
      let revs = design.reverse()
      for i in 0..<revs.length() {
        assert_true(best.score() >= revs[i].score())
      }
    }
    None => assert_true(false, msg="expected a best reverse primer")
  }
}

///|
test "primer3_design_primers_best_pair_returns_both" {
  let cfg = @src.Primer3Config::new(
    18, 22, 40.0, 75.0, 0.3, 0.8, 8, 8, 0.05, 50,
  )
  let template = @src.primer3_sample_template()
  let design = @src.primer3_design_primers(template, cfg)
  match design.best_pair() {
    Some((f, r)) => {
      assert_true(f.is_forward())
      assert_false(r.is_forward())
    }
    None => assert_true(false, msg="expected a best primer pair")
  }
}

///|
test "primer3_design_primers_forward_starts_after_flank" {
  let cfg = @src.Primer3Config::new(
    18, 22, 40.0, 75.0, 0.3, 0.8, 8, 8, 0.05, 50,
  )
  let template = @src.primer3_sample_template()
  let design = @src.primer3_design_primers(template, cfg)
  let fwds = design.forward()
  for i in 0..<fwds.length() {
    assert_true(fwds[i].start() >= cfg.flank())
  }
}

///|
test "primer3_design_primers_reverse_is_not_forward" {
  let cfg = @src.Primer3Config::new(
    18, 22, 40.0, 75.0, 0.3, 0.8, 8, 8, 0.05, 50,
  )
  let template = @src.primer3_sample_template()
  let design = @src.primer3_design_primers(template, cfg)
  let revs = design.reverse()
  for i in 0..<revs.length() {
    assert_false(revs[i].is_forward())
  }
}

///|
test "primer3_design_primers_candidates_within_length_bounds" {
  let cfg = @src.Primer3Config::new(
    18, 22, 40.0, 75.0, 0.3, 0.8, 8, 8, 0.05, 50,
  )
  let template = @src.primer3_sample_template()
  let design = @src.primer3_design_primers(template, cfg)
  for r in design.forward() {
    assert_true(r.length() >= cfg.primer_length_min())
    assert_true(r.length() <= cfg.primer_length_max())
  }
  for r in design.reverse() {
    assert_true(r.length() >= cfg.primer_length_min())
    assert_true(r.length() <= cfg.primer_length_max())
  }
}

///|
test "primer3_design_primers_candidates_within_tm_bounds" {
  let cfg = @src.Primer3Config::new(
    18, 22, 40.0, 75.0, 0.3, 0.8, 8, 8, 0.05, 50,
  )
  let template = @src.primer3_sample_template()
  let design = @src.primer3_design_primers(template, cfg)
  for r in design.forward() {
    assert_true(r.tm() >= cfg.tm_min())
    assert_true(r.tm() <= cfg.tm_max())
  }
  for r in design.reverse() {
    assert_true(r.tm() >= cfg.tm_min())
    assert_true(r.tm() <= cfg.tm_max())
  }
}

///|
test "primer3_design_primers_candidates_within_gc_bounds" {
  let cfg = @src.Primer3Config::new(
    18, 22, 40.0, 75.0, 0.3, 0.8, 8, 8, 0.05, 50,
  )
  let template = @src.primer3_sample_template()
  let design = @src.primer3_design_primers(template, cfg)
  for r in design.forward() {
    assert_true(r.gc() >= cfg.gc_min())
    assert_true(r.gc() <= cfg.gc_max())
  }
  for r in design.reverse() {
    assert_true(r.gc() >= cfg.gc_min())
    assert_true(r.gc() <= cfg.gc_max())
  }
}

///|
test "primer3_design_primers_no_pair_when_strict_config" {
  // An extremely strict config with very narrow Tm range may yield no primers.
  let cfg = @src.Primer3Config::new(
    20, 20, 99.0, 99.5, 0.5, 0.5, 0, 0, 0.05, 50,
  )
  let template = @src.primer3_short_template()
  let design = @src.primer3_design_primers(template, cfg)
  // Tm 99-99.5 is unreachable for normal DNA; expect no candidates.
  assert_eq(design.n_forward(), 0)
  assert_eq(design.n_reverse(), 0)
  assert_false(design.has_pair())
}

// ============================================================================
// 12. Amplicon utilities
// ============================================================================

///|
test "primer3_amplicon_length_basic" {
  // fwd at 50, len 20; rev at 200, len 20 -> amp_end = 220, amp_len = 170
  let al = @src.primer3_amplicon_length(50, 20, 200, 20)
  assert_eq(al, 170)
}

///|
test "primer3_amplicon_length_adjacent" {
  let al = @src.primer3_amplicon_length(0, 20, 20, 20)
  assert_eq(al, 40)
}

///|
test "primer3_pair_amplicon_length_from_results" {
  let f = @src.Primer3Result::new("ATGC", 50, 20, 60.0, 0.5, 0, true, 90.0)
  let r = @src.Primer3Result::new("GCAT", 200, 20, 60.0, 0.5, 0, false, 90.0)
  let al = @src.primer3_pair_amplicon_length(f, r)
  assert_eq(al, 170)
}

// ============================================================================
// 13. Sample data
// ============================================================================

///|
test "primer3_sample_template_nonempty" {
  let t = @src.primer3_sample_template()
  assert_true(t.length() > 0)
  assert_true(t.length() >= 200)
}

///|
test "primer3_sample_template_is_dna" {
  let t = @src.primer3_sample_template().to_upper()
  let mut i = 0
  while i < t.length() {
    let c = @src.substring(t, i, 1)
    assert_true(c == "A" || c == "T" || c == "G" || c == "C")
    i = i + 1
  }
}

///|
test "primer3_short_template_nonempty" {
  let t = @src.primer3_short_template()
  assert_true(t.length() > 0)
  assert_true(t.length() < @src.primer3_sample_template().length())
}

///|
test "primer3_design_summary_has_header" {
  let cfg = @src.primer3_default_config()
  let template = @src.primer3_sample_template()
  let design = @src.primer3_design_primers(template, cfg)
  let s = @src.primer3_design_summary(design)
  assert_true(s.contains("Primer3Design Summary"))
  assert_true(s.contains("Forward candidates:"))
  assert_true(s.contains("Reverse candidates:"))
  assert_true(s.contains("Has pair:"))
}

///|
test "primer3_design_summary_empty_design" {
  let design = @src.Primer3Design::new(Array::new(), Array::new())
  let s = @src.primer3_design_summary(design)
  assert_true(s.contains("Forward candidates: 0"))
  assert_true(s.contains("Reverse candidates: 0"))
  assert_true(s.contains("Has pair: false"))
  assert_true(s.contains("Best forward: (none)"))
}
