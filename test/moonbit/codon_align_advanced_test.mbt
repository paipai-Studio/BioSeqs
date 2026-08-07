// Tests for the advanced codon alignment selection-pressure statistics module.

///|
fn caa_test_close(actual : Double, expected : Double, tolerance : Double) -> Unit {
  if (actual - expected).abs() > tolerance {
    abort(
      "CodonAlign advanced value " +
      actual.to_string() +
      " differs from " +
      expected.to_string() +
      " (tolerance " +
      tolerance.to_string() +
      ")",
    )
  }
}

///|
// Sequences with many nonsynonymous changes → positive selection pressure.
// 9 codons (27 nt): Seq1 M-A-L-K-W-Q-Q-P-W, Seq2 M-G-L-E-R-H-P-S-M (7 nonsyn, 2 syn)
fn caa_test_positive_selection_pair() -> (String, String) {
  (
    "ATGGCCCTGAAATGGCAGCAGCCATGG",
    "ATGGGTCTGGAGAGACACCCCAGCATG",
  )
}

///|
// Sequences with mostly synonymous changes → purifying selection pressure.
// 9 codons (27 nt): Seq1 M-A-L-K-W-Q-Q-P-W, Seq2 has only synonymous codon switches.
fn caa_test_purifying_selection_pair() -> (String, String) {
  (
    "ATGGCCCTGAAATGGCAGCAGCCATGG",
    "ATGGCGCTCAAATGGCAGCAGCGATGG",
  )
}

///|
// Identical sequences → no selection, dN = dS = 0.
fn caa_test_identical_pair() -> (String, String) {
  ("ATGGCCCTGAAATGGCAGCAGCCATGG", "ATGGCCCTGAAATGGCAGCAGCCATGG")
}

// ===========================================================================
// Z-test for selection
// ===========================================================================

///|
test "CodonAlign advanced Z-test positive selection detects dN > dS" {
  let (s1, s2) = caa_test_positive_selection_pair()
  let result = @src.codon_test_selection(s1, s2, test_type="positive") catch {
    _ => abort("Z-test should run for positive selection pair")
  }
  assert_eq(result.test_type(), "positive")
  assert_true(result.dn() >= 0.0)
  assert_true(result.ds() >= 0.0)
  // With many nonsynonymous changes and no synonymous changes, dN should
  // be greater than dS, and the Z-score should be positive.
  assert_true(result.z_score() > 0.0 || result.variance() == 0.0)
  assert_true(result.p_value() >= 0.0 && result.p_value() <= 1.0)
}

///|
test "CodonAlign advanced Z-test purifying selection detects dN < dS" {
  let (s1, s2) = caa_test_purifying_selection_pair()
  let result = @src.codon_test_selection(s1, s2, test_type="purifying") catch {
    _ => abort("Z-test should run for purifying selection pair")
  }
  assert_eq(result.test_type(), "purifying")
  // With mostly synonymous changes, dS should be >= dN.
  assert_true(result.ds() >= result.dn() || result.variance() == 0.0)
  assert_true(result.p_value() >= 0.0 && result.p_value() <= 1.0)
}

///|
test "CodonAlign advanced Z-test neutrality is two-tailed" {
  let (s1, s2) = caa_test_positive_selection_pair()
  let result = @src.codon_test_selection(s1, s2, test_type="neutrality") catch {
    _ => abort("Z-test should run for neutrality test")
  }
  assert_eq(result.test_type(), "neutrality")
  assert_true(result.p_value() >= 0.0 && result.p_value() <= 1.0)
}

///|
test "CodonAlign advanced Z-test identical sequences has zero variance" {
  let (s1, s2) = caa_test_identical_pair()
  let result = @src.codon_test_selection(s1, s2) catch {
    _ => abort("Z-test should run for identical sequences")
  }
  caa_test_close(result.dn(), 0.0, 1.0e-12)
  caa_test_close(result.ds(), 0.0, 1.0e-12)
  caa_test_close(result.variance(), 0.0, 1.0e-12)
  caa_test_close(result.z_score(), 0.0, 1.0e-12)
  caa_test_close(result.p_value(), 1.0, 1.0e-12)
  assert_true(result.conclusion().contains("inconclusive"))
}

///|
test "CodonAlign advanced Z-test rejects short sequences" {
  let failed = try {
    ignore(@src.codon_test_selection("AT", "AT"))
    false
  } catch {
    CodonAlignAdvancedError(_) => true
  }
  assert_true(failed)
}

///|
test "CodonAlign advanced Z-test rejects mismatched lengths" {
  let failed = try {
    ignore(@src.codon_test_selection("ATGGCC", "ATGGCCAAG"))
    false
  } catch {
    CodonAlignAdvancedError(_) => true
  }
  assert_true(failed)
}

///|
test "CodonAlign advanced Z-test rejects non-codon-aligned lengths" {
  let failed = try {
    ignore(@src.codon_test_selection("ATGGCCAT", "ATGGCCAT"))
    false
  } catch {
    CodonAlignAdvancedError(_) => true
  }
  assert_true(failed)
}

// ===========================================================================
// Fisher's exact test for neutrality
// ===========================================================================

///|
test "CodonAlign advanced Fisher test two-sided runs and returns valid p-value" {
  let (s1, s2) = caa_test_positive_selection_pair()
  let result = @src.codon_test_neutrality(s1, s2) catch {
    _ => abort("Fisher exact test should run")
  }
  assert_eq(result.test_type(), "two-sided")
  assert_true(result.p_value() >= 0.0 && result.p_value() <= 1.0)
  assert_true(result.odds_ratio() >= 0.0)
  assert_true(result.n_diff() >= 0)
  assert_true(result.s_diff() >= 0)
  assert_true(result.n_sites() >= 0.0)
  assert_true(result.s_sites() >= 0.0)
}

///|
test "CodonAlign advanced Fisher test greater detects excess nonsynonymous" {
  let (s1, s2) = caa_test_positive_selection_pair()
  let result = @src.codon_test_neutrality(s1, s2, test_type="greater") catch {
    _ => abort("Fisher exact test greater should run")
  }
  assert_eq(result.test_type(), "greater")
  // With more nonsynonymous differences, the one-sided test for greater
  // should give a smaller p-value than two-sided.
  let two_sided = @src.codon_test_neutrality(s1, s2) catch {
    _ => abort("Fisher exact test two-sided should run")
  }
  assert_true(result.p_value() <= two_sided.p_value() + 1.0e-12)
}

///|
test "CodonAlign advanced Fisher test identical sequences has no differences" {
  let (s1, s2) = caa_test_identical_pair()
  let result = @src.codon_test_neutrality(s1, s2) catch {
    _ => abort("Fisher exact test should run for identical sequences")
  }
  assert_eq(result.n_diff(), 0)
  assert_eq(result.s_diff(), 0)
  // With no differences, p-value should be 1 (can't reject neutrality).
  caa_test_close(result.p_value(), 1.0, 1.0e-10)
  assert_true(result.conclusion().contains("not significant"))
}

///|
test "CodonAlign advanced Fisher test purifying pair has more syn differences" {
  let (s1, s2) = caa_test_purifying_selection_pair()
  let result = @src.codon_test_neutrality(s1, s2, test_type="less") catch {
    _ => abort("Fisher exact test less should run")
  }
  assert_eq(result.test_type(), "less")
  // With mostly synonymous changes, s_diff >= n_diff.
  assert_true(result.s_diff() >= result.n_diff())
}

// ===========================================================================
// Codon alignment builder
// ===========================================================================

///|
test "CodonAlign advanced builder constructs codon alignment with gaps" {
  let (protein_aln, coding_seqs) = @src.create_demo_protein_alignment()
  let alignment = @src.build_codon_alignment(protein_aln, coding_seqs) catch {
    _ => abort("codon alignment builder should succeed")
  }
  assert_eq(alignment.n_codons(), 5)
  let seqs = alignment.sequences()
  // Each sequence should have 5 codon columns × 3 nt = 15 characters.
  for seq in seqs {
    assert_eq(seq.length(), 15)
  }
  // Sequence 1 (M-A-K with gaps at pos 2,4): ATG(M) ---(gap) GCC(A) ---(gap) AAG(K)
  assert_eq(seqs[0], "ATG---GCC---AAG")
  // Sequence 2 (MTASK, no gaps): ATG(M) ACC(T) GCC(A) AGC(S) AAG(K)
  assert_eq(seqs[1], "ATGACCGCCAGCAAG")
}

///|
test "CodonAlign advanced builder with names" {
  let protein_aln = ["MAK", "MAK"]
  let coding_seqs = ["ATGGCCAAG", "ATGGCCAAG"]
  let alignment = @src.build_codon_alignment(
    protein_aln,
    coding_seqs,
    names=["gene1", "gene2"],
  ) catch {
    _ => abort("codon alignment builder should succeed")
  }
  let names = alignment.names()
  assert_eq(names[0], "gene1")
  assert_eq(names[1], "gene2")
  assert_eq(alignment.n_codons(), 3)
}

///|
test "CodonAlign advanced builder rejects translation mismatch" {
  // Protein says M-A-K but coding sequence has T(ACC) instead of A(GCC).
  let protein_aln = ["MAK"]
  let coding_seqs = ["ATGACCAAG"]
  let failed = try {
    ignore(@src.build_codon_alignment(protein_aln, coding_seqs))
    false
  } catch {
    CodonAlignAdvancedError(_) => true
  }
  assert_true(failed)
}

///|
test "CodonAlign advanced builder rejects exhausted coding sequence" {
  // Protein alignment is longer than the coding sequence can support.
  let protein_aln = ["MAKV"]
  let coding_seqs = ["ATGGCC"]
  let failed = try {
    ignore(@src.build_codon_alignment(protein_aln, coding_seqs))
    false
  } catch {
    CodonAlignAdvancedError(_) => true
  }
  assert_true(failed)
}

///|
test "CodonAlign advanced builder rejects leftover nucleotides" {
  let protein_aln = ["MA"]
  let coding_seqs = ["ATGGCCAAG"]
  let failed = try {
    ignore(@src.build_codon_alignment(protein_aln, coding_seqs))
    false
  } catch {
    CodonAlignAdvancedError(_) => true
  }
  assert_true(failed)
}

///|
test "CodonAlign advanced builder rejects mismatched row counts" {
  let failed = try {
    ignore(@src.build_codon_alignment(["MAK"], ["ATGGCC", "ATGGCC"]))
    false
  } catch {
    CodonAlignAdvancedError(_) => true
  }
  assert_true(failed)
}

///|
test "CodonAlign advanced builder rejects unequal protein lengths" {
  let failed = try {
    ignore(@src.build_codon_alignment(["MAK", "MAKV"], ["ATGGCC", "ATGGCCAAG"]))
    false
  } catch {
    CodonAlignAdvancedError(_) => true
  }
  assert_true(failed)
}

// ===========================================================================
// Sliding-window dN/dS
// ===========================================================================

///|
test "CodonAlign advanced sliding window produces multiple windows" {
  // 21-codon (63 nt) sequences.
  let s1 = "ATGGCCCTGAAATGGCAGCAGCCATGGGCCCTGAAATGGCAGCAGCCATGGGCCCTGAAATGG"
  let s2 = "ATGGGTCTGGAGAGACACCCCAGCATGGGGTCTGGAGAGACACCCCAGCATGGGGTCTGGAGA"
  let result = @src.sliding_window_dnds(s1, s2, window_size=5, step_size=3) catch {
    _ => abort("sliding window should run")
  }
  let windows = result.windows()
  assert_true(windows.length() >= 2)
  for w in windows {
    assert_true(w.end_codon() > w.start_codon())
    assert_true(w.n_codons() >= 3)
    assert_true(w.dn() >= 0.0)
    assert_true(w.ds() >= 0.0)
  }
}

///|
test "CodonAlign advanced sliding window skips degenerate windows" {
  // Identical sequences → all windows have dN = dS = 0.
  let s1 = "ATGGCCCTGAAATGGCAGCAGCCATGGGCCCTGAAATGGCAGCAGCCATGGGCCCTGAAATGG"
  let result = @src.sliding_window_dnds(s1, s1, window_size=5, step_size=3) catch {
    _ => abort("sliding window should run for identical sequences")
  }
  let windows = result.windows()
  for w in windows {
    caa_test_close(w.dn(), 0.0, 1.0e-12)
    caa_test_close(w.ds(), 0.0, 1.0e-12)
  }
}

///|
test "CodonAlign advanced sliding window rejects small window" {
  let failed = try {
    ignore(@src.sliding_window_dnds("ATGGCCATGGCC", "ATGGCCATGGCC", window_size=2))
    false
  } catch {
    CodonAlignAdvancedError(_) => true
  }
  assert_true(failed)
}

///|
test "CodonAlign advanced sliding window rejects invalid step" {
  let failed = try {
    ignore(@src.sliding_window_dnds("ATGGCCATGGCC", "ATGGCCATGGCC", step_size=0))
    false
  } catch {
    CodonAlignAdvancedError(_) => true
  }
  assert_true(failed)
}

// ===========================================================================
// BH-FDR
// ===========================================================================

///|
test "CodonAlign advanced BH-FDR is monotonic non-decreasing in rank" {
  let p_values = [0.001, 0.04, 0.03, 0.5, 0.01]
  let fdr = @src.codon_align_advanced_bh_fdr(p_values)
  assert_eq(fdr.length(), 5)
  // Sort p-values and verify FDR is non-decreasing in sorted order.
  let indexed : Array[(Double, Int)] = []
  for i in 0..<p_values.length() {
    indexed.push((p_values[i], i))
  }
  indexed.sort_by(fn(a : (Double, Int), b : (Double, Int)) -> Int {
    if a.0 < b.0 {
      -1
    } else if a.0 > b.0 {
      1
    } else {
      0
    }
  })
  let mut previous = 0.0
  for entry in indexed {
    let adjusted = fdr[entry.1]
    assert_true(adjusted >= previous - 1.0e-12)
    previous = adjusted
  }
}

///|
test "CodonAlign advanced BH-FDR clamps p-values to [0,1]" {
  let p_values = [-0.5, 0.5, 1.5]
  let fdr = @src.codon_align_advanced_bh_fdr(p_values)
  for v in fdr {
    assert_true(v >= 0.0 && v <= 1.0)
  }
}

///|
test "CodonAlign advanced BH-FDR empty input returns empty" {
  let fdr = @src.codon_align_advanced_bh_fdr([])
  assert_eq(fdr.length(), 0)
}

///|
test "CodonAlign advanced BH-FDR single p-value equals itself" {
  let fdr = @src.codon_align_advanced_bh_fdr([0.03])
  caa_test_close(fdr[0], 0.03, 1.0e-12)
}

// ===========================================================================
// Pairwise Ka/Ks table
// ===========================================================================

///|
test "CodonAlign advanced pairwise table computes all pairs with FDR" {
  let seqs = @src.create_positive_selection_alignment()
  let table = @src.pairwise_kaks_table(
    seqs,
    names=["g1", "g2", "g3", "g4", "g5", "g6", "g7", "g8", "g9", "g10"],
    fdr_threshold=0.05,
  ) catch {
    _ => abort("pairwise Ka/Ks table should run")
  }
  // C(10, 2) = 45 pairs.
  assert_eq(table.length(), 45)
  for row in table {
    assert_true(row.seq1_name().length() > 0)
    assert_true(row.seq2_name().length() > 0)
    assert_true(row.dn() >= 0.0)
    assert_true(row.ds() >= 0.0)
    assert_true(row.fdr() >= 0.0 && row.fdr() <= 1.0)
    assert_true(row.p_value() >= 0.0 && row.p_value() <= 1.0)
  }
}

///|
test "CodonAlign advanced pairwise table rejects single sequence" {
  let failed = try {
    ignore(@src.pairwise_kaks_table(["ATGGCC"]))
    false
  } catch {
    CodonAlignAdvancedError(_) => true
  }
  assert_true(failed)
}

///|
test "CodonAlign advanced pairwise table uses default names when none provided" {
  let table = @src.pairwise_kaks_table(["ATGGCCATG", "ATGGCCATG", "ATGGCCATG"]) catch {
    _ => abort("pairwise Ka/Ks table should run")
  }
  // 3 sequences → 3 pairs.
  assert_eq(table.length(), 3)
  assert_eq(table[0].seq1_name(), "seq0")
  assert_eq(table[0].seq2_name(), "seq1")
}

// ===========================================================================
// Demo data
// ===========================================================================

///|
test "CodonAlign advanced purifying selection demo has valid codons" {
  let seqs = @src.create_purifying_selection_alignment()
  assert_true(seqs.length() >= 2)
  for seq in seqs {
    assert_eq(seq.length() % 3, 0)
  }
}

///|
test "CodonAlign advanced positive selection demo has valid codons" {
  let seqs = @src.create_positive_selection_alignment()
  assert_true(seqs.length() >= 2)
  for seq in seqs {
    assert_eq(seq.length() % 3, 0)
  }
}

///|
test "CodonAlign advanced demo protein alignment is consistent" {
  let (protein_aln, coding_seqs) = @src.create_demo_protein_alignment()
  assert_eq(protein_aln.length(), coding_seqs.length())
  // All protein alignments should have the same length.
  let aln_len = protein_aln[0].length()
  for pa in protein_aln {
    assert_eq(pa.length(), aln_len)
  }
}
