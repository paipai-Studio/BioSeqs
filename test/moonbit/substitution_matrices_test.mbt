///|
/// Test file for substitution_matrices module.
/// Implements Bio.Align.substitution_matrices equivalent functionality.

// ---------------------------------------------------------------------------
// Alphabets
// ---------------------------------------------------------------------------

test "submat_protein_alphabet_has_20_letters" {
  let alpha = @src.subs_protein_alphabet()
  assert_eq(alpha.length(), 20)
  assert_eq(alpha[0], "A")
  assert_eq(alpha[19], "Y")
}

///|
test "submat_nucleotide_alphabet_has_4_letters" {
  let alpha = @src.subs_nucleotide_alphabet()
  assert_eq(alpha.length(), 4)
  assert_eq(alpha, ["A", "C", "G", "T"])
}

// ---------------------------------------------------------------------------
// ArrayData structure
// ---------------------------------------------------------------------------

///|
test "submat_array_data_from_dict_basic" {
  let pairs = [
    ("A", "A", 5.0),
    ("A", "C", -1.0),
    ("C", "A", -1.0),
    ("C", "C", 9.0),
  ]
  let ad = @src.array_data_from_dict(pairs, ["A", "C"], ["A", "C"])
  assert_eq(ad.rows, ["A", "C"])
  assert_eq(ad.cols, ["A", "C"])
  assert_true((ad.get(0, 0) - 5.0).abs() < 1.0e-10)
  assert_true((ad.get(0, 1) - -1.0).abs() < 1.0e-10)
  assert_true((ad.get(1, 1) - 9.0).abs() < 1.0e-10)
}

///|
test "submat_array_data_get_by_label" {
  let pairs = [("A", "B", 3.0), ("B", "A", 7.0)]
  let ad = @src.array_data_from_dict(pairs, ["A", "B"], ["A", "B"])
  assert_true((ad.get_by_label("A", "B") - 3.0).abs() < 1.0e-10)
  assert_true((ad.get_by_label("B", "A") - 7.0).abs() < 1.0e-10)
  // Missing label returns 0.0
  assert_true(ad.get_by_label("X", "Y").abs() < 1.0e-10)
}

///|
test "submat_array_data_set_updates_value" {
  let pairs : Array[(String, String, Double)] = []
  let ad = @src.array_data_from_dict(pairs, ["A", "B"], ["A", "B"])
  ad.set(0, 1, 4.5)
  assert_true((ad.get(0, 1) - 4.5).abs() < 1.0e-10)
  // Out-of-bounds set is silently ignored
  ad.set(99, 99, 1.0)
  assert_true(true)
}

///|
test "submat_array_data_get_out_of_bounds_returns_zero" {
  let pairs : Array[(String, String, Double)] = []
  let ad = @src.array_data_from_dict(pairs, ["A"], ["A"])
  assert_true(ad.get(-1, 0).abs() < 1.0e-10)
  assert_true(ad.get(0, 99).abs() < 1.0e-10)
}

// ---------------------------------------------------------------------------
// Built-in BLOSUM matrices
// ---------------------------------------------------------------------------

///|
test "submat_blosum62_known_scores" {
  let m = @src.subs_blosum62_matrix()
  assert_eq(m.name(), "BLOSUM62")
  assert_eq(m.n_letters(), 20)
  assert_true((m.get_score("A", "A") - 4.0).abs() < 1.0e-10)
  assert_true((m.get_score("R", "R") - 5.0).abs() < 1.0e-10)
  assert_true((m.get_score("A", "R") - -1.0).abs() < 1.0e-10)
  assert_true((m.get_score("W", "W") - 11.0).abs() < 1.0e-10)
  assert_true((m.get_score("C", "C") - 9.0).abs() < 1.0e-10)
}

///|
test "submat_blosum62_is_symmetric" {
  let m = @src.subs_blosum62_matrix()
  let alpha = m.alphabet()
  for i in 0..<alpha.length() {
    for j in 0..<alpha.length() {
      let s_ij = m.get_score(alpha[i], alpha[j])
      let s_ji = m.get_score(alpha[j], alpha[i])
      assert_true((s_ij - s_ji).abs() < 1.0e-10)
    }
  }
}

///|
test "submat_blosum45_known_scores" {
  let m = @src.subs_blosum45_matrix()
  assert_eq(m.name(), "BLOSUM45")
  assert_true((m.get_score("A", "A") - 5.0).abs() < 1.0e-10)
  assert_true((m.get_score("R", "R") - 6.0).abs() < 1.0e-10)
  assert_true((m.get_score("W", "W") - 14.0).abs() < 1.0e-10)
  assert_true((m.get_score("C", "C") - 12.0).abs() < 1.0e-10)
}

///|
test "submat_blosum80_known_scores" {
  let m = @src.subs_blosum80_matrix()
  assert_eq(m.name(), "BLOSUM80")
  assert_true((m.get_score("A", "A") - 7.0).abs() < 1.0e-10)
  assert_true((m.get_score("R", "R") - 8.0).abs() < 1.0e-10)
  assert_true((m.get_score("W", "W") - 15.0).abs() < 1.0e-10)
}

///|
test "submat_blosum90_known_scores" {
  let m = @src.subs_blosum90_matrix()
  assert_eq(m.name(), "BLOSUM90")
  assert_true((m.get_score("A", "A") - 7.0).abs() < 1.0e-10)
  assert_true((m.get_score("R", "R") - 9.0).abs() < 1.0e-10)
  assert_true((m.get_score("C", "C") - 14.0).abs() < 1.0e-10)
}

// ---------------------------------------------------------------------------
// Built-in PAM matrices
// ---------------------------------------------------------------------------

///|
test "submat_pam30_known_scores" {
  let m = @src.subs_pam30_matrix()
  assert_eq(m.name(), "PAM30")
  assert_true((m.get_score("A", "A") - 6.0).abs() < 1.0e-10)
  assert_true((m.get_score("C", "C") - 12.0).abs() < 1.0e-10)
  assert_true((m.get_score("W", "W") - 14.0).abs() < 1.0e-10)
}

///|
test "submat_pam70_known_scores" {
  let m = @src.subs_pam70_matrix()
  assert_eq(m.name(), "PAM70")
  assert_true((m.get_score("A", "A") - 5.0).abs() < 1.0e-10)
  assert_true((m.get_score("C", "C") - 10.0).abs() < 1.0e-10)
  assert_true((m.get_score("W", "W") - 13.0).abs() < 1.0e-10)
}

///|
test "submat_pam250_known_scores" {
  let m = @src.subs_pam250_matrix()
  assert_eq(m.name(), "PAM250")
  assert_true((m.get_score("A", "A") - 2.0).abs() < 1.0e-10)
  assert_true((m.get_score("C", "C") - 12.0).abs() < 1.0e-10)
  assert_true((m.get_score("W", "W") - 17.0).abs() < 1.0e-10)
}

// ---------------------------------------------------------------------------
// Built-in nucleotide matrix
// ---------------------------------------------------------------------------

///|
test "submat_nuc44_match_mismatch" {
  let m = @src.subs_nuc44_matrix()
  assert_eq(m.name(), "NUC4.4")
  assert_eq(m.n_letters(), 4)
  assert_true((m.get_score("A", "A") - 5.0).abs() < 1.0e-10)
  assert_true((m.get_score("A", "C") - -4.0).abs() < 1.0e-10)
  assert_true((m.get_score("G", "G") - 5.0).abs() < 1.0e-10)
  assert_true((m.get_score("T", "T") - 5.0).abs() < 1.0e-10)
  assert_true((m.get_score("A", "T") - -4.0).abs() < 1.0e-10)
}

// ---------------------------------------------------------------------------
// SubsMatrix construction and methods
// ---------------------------------------------------------------------------

///|
test "submat_subs_matrix_from_pairs" {
  let alpha = ["A", "C"]
  let pairs = [
    ("A", "A", 5.0),
    ("A", "C", -2.0),
    ("C", "A", -2.0),
    ("C", "C", 9.0),
  ]
  let m = @src.subs_matrix_from_pairs("TEST", alpha, pairs)
  assert_eq(m.name(), "TEST")
  assert_eq(m.n_letters(), 2)
  assert_true((m.get_score("A", "A") - 5.0).abs() < 1.0e-10)
  assert_true((m.get_score("C", "C") - 9.0).abs() < 1.0e-10)
  assert_true((m.get_score("A", "C") - -2.0).abs() < 1.0e-10)
}

///|
test "submat_get_score_idx" {
  let m = @src.subs_blosum62_matrix()
  // A=0, R=14 in alphabet ACDEFGHIKLMNPQRSTVWY
  assert_true((m.get_score_idx(0, 0) - 4.0).abs() < 1.0e-10)
  // A-R = -1
  assert_true((m.get_score_idx(0, 14) - -1.0).abs() < 1.0e-10)
}

///|
test "submat_select_submatrix" {
  let m = @src.subs_blosum62_matrix()
  let sub = m.select(["A", "C", "D"])
  assert_eq(sub.n_letters(), 3)
  // A-A=4, C-C=9, D-D=6, A-C=0
  assert_true((sub.get_score("A", "A") - 4.0).abs() < 1.0e-10)
  assert_true((sub.get_score("C", "C") - 9.0).abs() < 1.0e-10)
  assert_true((sub.get_score("D", "D") - 6.0).abs() < 1.0e-10)
  assert_true((sub.get_score("A", "C") - 0.0).abs() < 1.0e-10)
}

///|
test "submat_to_table_string_has_name" {
  let m = @src.subs_blosum62_matrix()
  let s = m.to_table_string()
  assert_true(s.contains("BLOSUM62"))
  // Should contain all 20 letters as header
  assert_true(s.contains("A"))
  assert_true(s.contains("Y"))
}

// ---------------------------------------------------------------------------
// Matrix registry
// ---------------------------------------------------------------------------

///|
test "submat_registry_initialize_and_list" {
  @src.subs_initialize_registry()
  let names = @src.list_matrices()
  // Should include all 8 built-in matrices
  assert_true(names.length() >= 8)
  let mut has_blosum62 = false
  let mut has_pam250 = false
  let mut has_nuc44 = false
  for n in names {
    if n == "BLOSUM62" {
      has_blosum62 = true
    }
    if n == "PAM250" {
      has_pam250 = true
    }
    if n == "NUC4.4" {
      has_nuc44 = true
    }
  }
  assert_true(has_blosum62)
  assert_true(has_pam250)
  assert_true(has_nuc44)
}

///|
test "submat_registry_load_existing_matrix" {
  @src.subs_initialize_registry()
  let opt = @src.load_matrix("BLOSUM62")
  assert_true(opt.is_some())
  let m = opt.unwrap()
  assert_eq(m.name(), "BLOSUM62")
  assert_true((m.get_score("A", "A") - 4.0).abs() < 1.0e-10)
}

///|
test "submat_registry_load_missing_returns_none" {
  @src.subs_initialize_registry()
  let opt = @src.load_matrix("NONEXISTENT_MATRIX_XYZ")
  assert_true(opt.is_none())
}

///|
test "submat_registry_register_custom_matrix" {
  let alpha = ["A", "C"]
  let pairs = [
    ("A", "A", 1.0),
    ("A", "C", -1.0),
    ("C", "A", -1.0),
    ("C", "C", 1.0),
  ]
  let custom = @src.subs_matrix_from_pairs("CUSTOM_TEST_MATRIX", alpha, pairs)
  @src.register_matrix("CUSTOM_TEST_MATRIX", custom)
  let opt = @src.load_matrix("CUSTOM_TEST_MATRIX")
  assert_true(opt.is_some())
  let loaded = opt.unwrap()
  assert_true((loaded.get_score("A", "A") - 1.0).abs() < 1.0e-10)
}

///|
test "submat_registry_register_overwrites_existing" {
  let alpha = ["A"]
  let pairs1 = [("A", "A", 5.0)]
  let pairs2 = [("A", "A", 7.0)]
  let m1 = @src.subs_matrix_from_pairs("OVERWRITE_TEST", alpha, pairs1)
  let m2 = @src.subs_matrix_from_pairs("OVERWRITE_TEST", alpha, pairs2)
  @src.register_matrix("OVERWRITE_TEST", m1)
  @src.register_matrix("OVERWRITE_TEST", m2)
  let loaded = @src.load_matrix("OVERWRITE_TEST").unwrap()
  assert_true((loaded.get_score("A", "A") - 7.0).abs() < 1.0e-10)
}

// ---------------------------------------------------------------------------
// Frequency matrix calculation
// ---------------------------------------------------------------------------

///|
test "submat_calculate_frequency_matrix_identical_sequences" {
  let alignment = ["ACGT", "ACGT"]
  let alpha = @src.subs_nucleotide_alphabet()
  let freq = @src.calculate_frequency_matrix(alignment, alpha)
  // 2 sequences, 4 positions, all matches
  // Each match contributes 2 to the diagonal (symmetric)
  // A-A=2, C-C=2, G-G=2, T-T=2
  assert_true((freq.get(0, 0) - 2.0).abs() < 1.0e-10) // A
  assert_true((freq.get(1, 1) - 2.0).abs() < 1.0e-10) // C
  assert_true((freq.get(2, 2) - 2.0).abs() < 1.0e-10) // G
  assert_true((freq.get(3, 3) - 2.0).abs() < 1.0e-10) // T
  // Off-diagonal should be 0
  assert_true(freq.get(0, 1).abs() < 1.0e-10)
}

///|
test "submat_calculate_frequency_matrix_with_substitutions" {
  let alignment = ["ACGT", "AGGT"] // position 1: C->G
  let alpha = @src.subs_nucleotide_alphabet()
  let freq = @src.calculate_frequency_matrix(alignment, alpha)
  // Position 0: A-A -> freq[0][0] += 2
  // Position 1: C-G -> freq[1][2] += 1, freq[2][1] += 1
  // Position 2: G-G -> freq[2][2] += 2
  // Position 3: T-T -> freq[3][3] += 2
  assert_true((freq.get(0, 0) - 2.0).abs() < 1.0e-10)
  assert_true((freq.get(1, 2) - 1.0).abs() < 1.0e-10)
  assert_true((freq.get(2, 1) - 1.0).abs() < 1.0e-10)
  assert_true((freq.get(2, 2) - 2.0).abs() < 1.0e-10)
  assert_true((freq.get(3, 3) - 2.0).abs() < 1.0e-10)
}

///|
test "submat_calculate_frequency_matrix_skips_gaps" {
  let alignment = ["A-A", "ACA"] // position 1: gap in seq1
  let alpha = @src.subs_nucleotide_alphabet()
  let freq = @src.calculate_frequency_matrix(alignment, alpha)
  // Position 0: A-A -> freq[0][0] += 2
  // Position 1: gap-C -> skipped
  // Position 2: A-A -> freq[0][0] += 2
  assert_true((freq.get(0, 0) - 4.0).abs() < 1.0e-10)
  // No C-C pairs
  assert_true(freq.get(1, 1).abs() < 1.0e-10)
}

///|
test "submat_calculate_frequency_matrix_single_sequence_returns_zeros" {
  let alignment = ["ACGT"]
  let alpha = @src.subs_nucleotide_alphabet()
  let freq = @src.calculate_frequency_matrix(alignment, alpha)
  for i in 0..<4 {
    for j in 0..<4 {
      assert_true(freq.get(i, j).abs() < 1.0e-10)
    }
  }
}

// ---------------------------------------------------------------------------
// Substitution matrix calculation (log-odds)
// ---------------------------------------------------------------------------

///|
test "submat_calculate_substitution_matrix_identical_seqs" {
  let alignment = ["AC", "AC"]
  let alpha = ["A", "C"]
  let m = @src.calculate_substitution_matrix(alignment, alpha, scale=2.0)
  // A-A: q=2/4=0.5, p(A)=2/4=0.5, exp=0.25
  // score = 2 * log2(0.5/0.25) = 2 * log2(2) = 2 * 1 = 2
  assert_true((m.get_score("A", "A") - 2.0).abs() < 1.0e-10)
  assert_true((m.get_score("C", "C") - 2.0).abs() < 1.0e-10)
  // Off-diagonal: observed=0 -> -999
  assert_true((m.get_score("A", "C") - -999.0).abs() < 1.0e-10)
}

///|
test "submat_calculate_substitution_matrix_default_scale" {
  let alignment = ["AA", "AA"]
  let alpha = ["A"]
  let m = @src.calculate_substitution_matrix(alignment, alpha)
  // Default scale = 2.0
  // q(A,A) = 4/4 = 1.0, p(A) = 4/4 = 1.0
  // exp = 1.0, log2(1.0/1.0) = 0, score = 0
  assert_true((m.get_score("A", "A") - 0.0).abs() < 1.0e-10)
}

// ---------------------------------------------------------------------------
// Shannon entropy
// ---------------------------------------------------------------------------

///|
test "submat_shannon_entropy_uniform_distribution" {
  // 2x2 matrix with all equal values: H = log2(4) = 2
  let pairs = [
    ("A", "A", 1.0),
    ("A", "C", 1.0),
    ("C", "A", 1.0),
    ("C", "C", 1.0),
  ]
  let freq = @src.array_data_from_dict(pairs, ["A", "C"], ["A", "C"])
  let h = @src.subs_shannon_entropy(freq)
  assert_true((h - 2.0).abs() < 1.0e-10)
}

///|
test "submat_shannon_entropy_certain_distribution" {
  // All mass on one cell: H = 0
  let pairs = [("A", "A", 4.0)]
  let freq = @src.array_data_from_dict(pairs, ["A", "C"], ["A", "C"])
  let h = @src.subs_shannon_entropy(freq)
  assert_true(h.abs() < 1.0e-10)
}

///|
test "submat_shannon_entropy_empty_matrix" {
  let pairs : Array[(String, String, Double)] = []
  let freq = @src.array_data_from_dict(pairs, ["A"], ["A"])
  let h = @src.subs_shannon_entropy(freq)
  assert_true(h.abs() < 1.0e-10)
}

///|
test "submat_shannon_entropy_two_equal_cells" {
  // Two non-zero cells with equal weight: H = 1
  let pairs = [("A", "A", 1.0), ("C", "C", 1.0)]
  let freq = @src.array_data_from_dict(pairs, ["A", "C"], ["A", "C"])
  let h = @src.subs_shannon_entropy(freq)
  assert_true((h - 1.0).abs() < 1.0e-10)
}

// ---------------------------------------------------------------------------
// Relative entropy (KL divergence)
// ---------------------------------------------------------------------------

///|
test "submat_relative_entropy_uniform_is_zero" {
  // When all letters equally frequent, observed == expected, KL = 0
  let pairs = [
    ("A", "A", 1.0),
    ("A", "C", 1.0),
    ("C", "A", 1.0),
    ("C", "C", 1.0),
  ]
  let freq = @src.array_data_from_dict(pairs, ["A", "C"], ["A", "C"])
  let kl = @src.subs_relative_entropy(freq)
  assert_true(kl.abs() < 1.0e-10)
}

///|
test "submat_relative_entropy_identical_pairs_zero" {
  // Only A-A pairs: q(A,A)=1, p(A)=1, exp=1, KL = 1 * log2(1/1) = 0
  let pairs = [("A", "A", 4.0)]
  let freq = @src.array_data_from_dict(pairs, ["A", "C"], ["A", "C"])
  let kl = @src.subs_relative_entropy(freq)
  assert_true(kl.abs() < 1.0e-10)
}

///|
test "submat_relative_entropy_empty_matrix" {
  let pairs : Array[(String, String, Double)] = []
  let freq = @src.array_data_from_dict(pairs, ["A"], ["A"])
  let kl = @src.subs_relative_entropy(freq)
  assert_true(kl.abs() < 1.0e-10)
}

// ---------------------------------------------------------------------------
// NCBI matrix parsing
// ---------------------------------------------------------------------------

///|
test "submat_parse_ncbi_matrix_simple" {
  let content =
    #|# Test matrix
    #|   A  C
    #|A  1 -1
    #|C -1  1
  let opt = @src.parse_ncbi_matrix(content)
  assert_true(opt.is_some())
  let m = opt.unwrap()
  assert_true((m.get_score("A", "A") - 1.0).abs() < 1.0e-10)
  assert_true((m.get_score("A", "C") - -1.0).abs() < 1.0e-10)
  assert_true((m.get_score("C", "C") - 1.0).abs() < 1.0e-10)
}

///|
test "submat_parse_ncbi_matrix_with_comments" {
  let content =
    #|# Generated by NCBI
    #|# First line: header
    #|   A  C  G  T
    #|A  5 -4 -4 -4
    #|C -4  5 -4 -4
    #|G -4 -4  5 -4
    #|T -4 -4 -4  5
  let opt = @src.parse_ncbi_matrix(content)
  assert_true(opt.is_some())
  let m = opt.unwrap()
  assert_true((m.get_score("A", "A") - 5.0).abs() < 1.0e-10)
  assert_true((m.get_score("C", "C") - 5.0).abs() < 1.0e-10)
  assert_true((m.get_score("A", "C") - -4.0).abs() < 1.0e-10)
  assert_true((m.get_score("G", "T") - -4.0).abs() < 1.0e-10)
}

///|
test "submat_parse_ncbi_matrix_empty_returns_none" {
  let content = "# only comments\n# nothing else\n"
  let opt = @src.parse_ncbi_matrix(content)
  assert_true(opt.is_none())
}

///|
test "submat_parse_ncbi_matrix_handles_extra_whitespace" {
  let content = "A    C\nA   2   -3\nC  -3    2\n"
  let opt = @src.parse_ncbi_matrix(content)
  assert_true(opt.is_some())
  let m = opt.unwrap()
  assert_true((m.get_score("A", "A") - 2.0).abs() < 1.0e-10)
  assert_true((m.get_score("C", "C") - 2.0).abs() < 1.0e-10)
  assert_true((m.get_score("A", "C") - -3.0).abs() < 1.0e-10)
}

// ---------------------------------------------------------------------------
// Matrix correlation
// ---------------------------------------------------------------------------

///|
test "submat_matrix_correlation_identical_matrices" {
  let m1 = @src.subs_blosum62_matrix()
  let m2 = @src.subs_blosum62_matrix()
  let corr = @src.matrix_correlation(m1, m2)
  // Identical matrices should have correlation = 1.0
  assert_true((corr - 1.0).abs() < 1.0e-9)
}

///|
test "submat_matrix_correlation_blosum45_vs_blosum62" {
  let m45 = @src.subs_blosum45_matrix()
  let m62 = @src.subs_blosum62_matrix()
  let corr = @src.matrix_correlation(m45, m62)
  // BLOSUM45 and BLOSUM62 should be highly correlated
  assert_true(corr > 0.85)
  assert_true(corr <= 1.0)
}

///|
test "submat_matrix_correlation_disjoint_alphabets" {
  let alpha1 = ["A"]
  let alpha2 = ["C"]
  let m1 = @src.subs_matrix_from_pairs("M1", alpha1, [("A", "A", 5.0)])
  let m2 = @src.subs_matrix_from_pairs("M2", alpha2, [("C", "C", 9.0)])
  let corr = @src.matrix_correlation(m1, m2)
  // No common letters -> correlation 0
  assert_true(corr.abs() < 1.0e-10)
}

// ---------------------------------------------------------------------------
// Biopython 1.88 canonical matrix registry (30 matrices)
// ---------------------------------------------------------------------------

///|
test "submat_registry_contains_all_30_biopython_matrices" {
  @src.subs_initialize_registry()
  let names = @src.list_matrices()
  // 30 canonical Biopython names plus the legacy "NUC4.4" alias
  assert_true(names.length() >= 31)
  let expected = [
    "BENNER22", "BENNER6", "BENNER74", "BLASTN", "BLASTP", "BLOSUM45", "BLOSUM50",
    "BLOSUM62", "BLOSUM80", "BLOSUM90", "DAYHOFF", "FENG", "GENETIC", "GONNET1992",
    "HOXD70", "JOHNSON", "JONES", "LEVIN", "MCLACHLAN", "MDM78", "MEGABLAST", "NUC.4.4",
    "PAM250", "PAM30", "PAM70", "RAO", "RISLER", "SCHNEIDER", "STR", "TRANS",
  ]
  for name in expected {
    let m = @src.load_matrix(name)
    assert_true(m.is_some(), msg=name)
  }
}

///|
test "submat_canonical_alphabet_sizes" {
  @src.subs_initialize_registry()
  let size24 = [
    "BLOSUM45", "BLOSUM50", "BLOSUM62", "BLOSUM80", "BLOSUM90", "PAM30", "PAM70",
    "PAM250",
  ]
  for name in size24 {
    let m = @src.load_matrix(name).unwrap()
    assert_eq(m.n_letters(), 24, msg=name)
    // canonical 24-letter alphabet includes B, Z, X and '*'
    assert_true(m.has_letter("B"))
    assert_true(m.has_letter("Z"))
    assert_true(m.has_letter("X"))
    assert_true(m.has_letter("*"))
  }
  let size15 = ["BLASTN", "MEGABLAST", "NUC.4.4"]
  for name in size15 {
    let m = @src.load_matrix(name).unwrap()
    assert_eq(m.n_letters(), 15, msg=name)
  }
  assert_eq(@src.load_matrix("BLASTP").unwrap().n_letters(), 27)
  assert_eq(@src.load_matrix("SCHNEIDER").unwrap().n_letters(), 64)
  let size4 = ["HOXD70", "TRANS"]
  for name in size4 {
    assert_eq(@src.load_matrix(name).unwrap().n_letters(), 4, msg=name)
  }
  let size20 = [
    "BENNER22", "BENNER6", "BENNER74", "DAYHOFF", "FENG", "GENETIC", "GONNET1992",
    "JOHNSON", "JONES", "LEVIN", "MCLACHLAN", "MDM78", "RAO", "RISLER", "STR",
  ]
  for name in size20 {
    assert_eq(@src.load_matrix(name).unwrap().n_letters(), 20, msg=name)
  }
  // Legacy simplified matrix is still 4 letters
  assert_eq(@src.subs_nuc44_matrix().n_letters(), 4)
}

///|
test "submat_canonical_blosum62_scores" {
  @src.subs_initialize_registry()
  // The registered BLOSUM62 is now the authoritative 24-letter matrix,
  // with the exact Biopython alphabet order.
  let m = @src.load_matrix("BLOSUM62").unwrap()
  assert_eq(m.alphabet().length(), 24)
  assert_eq(m.alphabet(), [
    "A", "R", "N", "D", "C", "Q", "E", "G", "H", "I", "L", "K", "M", "F", "P", "S",
    "T", "W", "Y", "V", "B", "Z", "X", "*",
  ])
  assert_true((m.get_score("A", "A") - 4.0).abs() < 1.0e-10)
  assert_true((m.get_score("W", "W") - 11.0).abs() < 1.0e-10)
  assert_true((m.get_score("X", "X") - -1.0).abs() < 1.0e-10)
  assert_true((m.get_score("*", "*") - 1.0).abs() < 1.0e-10)
  assert_true((m.get_score("A", "*") - -4.0).abs() < 1.0e-10)
  // ambiguous residues: B pairs with D/N, Z pairs with E/Q
  assert_true((m.get_score("B", "B") - 4.0).abs() < 1.0e-10)
  assert_true((m.get_score("B", "D") - 4.0).abs() < 1.0e-10)
  assert_true((m.get_score("Z", "Z") - 4.0).abs() < 1.0e-10)
  assert_true((m.get_score("Z", "E") - 4.0).abs() < 1.0e-10)
  assert_true(!m.has_letter("J"))
}

///|
test "submat_canonical_nuc44_iupac_scores" {
  @src.subs_initialize_registry()
  let m = @src.load_matrix("NUC.4.4").unwrap()
  assert_eq(m.alphabet(), [
    "A", "T", "G", "C", "S", "W", "R", "Y", "K", "M", "B", "V", "H", "D", "N",
  ])
  assert_true((m.get_score("A", "A") - 5.0).abs() < 1.0e-10)
  assert_true((m.get_score("A", "C") - -4.0).abs() < 1.0e-10)
  assert_true((m.get_score("A", "R") - 1.0).abs() < 1.0e-10)
  assert_true((m.get_score("S", "S") - -1.0).abs() < 1.0e-10)
  assert_true((m.get_score("S", "C") - 1.0).abs() < 1.0e-10)
  assert_true((m.get_score("N", "A") - -2.0).abs() < 1.0e-10)
}

///|
test "submat_canonical_dna_and_codon_matrix_scores" {
  @src.subs_initialize_registry()
  let hox = @src.load_matrix("HOXD70").unwrap()
  assert_eq(hox.alphabet(), ["A", "C", "G", "T"])
  assert_true((hox.get_score("A", "A") - 91.0).abs() < 1.0e-10)
  assert_true((hox.get_score("A", "C") - -114.0).abs() < 1.0e-10)
  let trans = @src.load_matrix("TRANS").unwrap()
  assert_eq(trans.alphabet(), ["A", "T", "C", "G"])
  // SCHNEIDER is a 64-codon scoring matrix
  let sch = @src.load_matrix("SCHNEIDER").unwrap()
  assert_eq(sch.alphabet()[0], "AAA")
  assert_true((sch.get_score("AAA", "AAA") - 11.6).abs() < 1.0e-10)
  assert_true((sch.get_score("AAA", "AAC") - -2.7).abs() < 1.0e-10)
  // BLASTP uses the full A-Z + '*' alphabet
  let blastp = @src.load_matrix("BLASTP").unwrap()
  assert_true((blastp.get_score("A", "A") - 4.0).abs() < 1.0e-10)
  assert_true(blastp.has_letter("J"))
}

///|
test "submat_canonical_matrices_are_symmetric" {
  @src.subs_initialize_registry()
  for name in ["DAYHOFF", "GONNET1992", "PAM250", "BLASTN", "BLASTP", "HOXD70"] {
    let m = @src.load_matrix(name).unwrap()
    let alpha = m.alphabet()
    for a in alpha {
      for b in alpha {
        assert_true(
          (m.get_score(a, b) - m.get_score(b, a)).abs() < 1.0e-9,
          msg=name,
        )
      }
    }
  }
}

// ---------------------------------------------------------------------------
// End-to-end: derive matrix from alignment, then score
// ---------------------------------------------------------------------------

///|
test "submat_end_to_end_derive_matrix_and_score" {
  // Build a small alignment where A and C never co-occur at the same position
  let alignment = ["AAA", "AAA", "CCC", "CCC"]
  let alpha = ["A", "C"]
  let m = @src.calculate_substitution_matrix(alignment, alpha, scale=2.0)
  // Pairs:
  // (0,1): AAA vs AAA -> 3x A-A -> freq[A][A] += 6 (3*2 symmetric)
  // (0,2): AAA vs CCC -> 3x A-C -> freq[A][C] += 3, freq[C][A] += 3
  // (0,3): AAA vs CCC -> 3x A-C -> freq[A][C] += 3, freq[C][A] += 3
  // (1,2): AAA vs CCC -> 3x A-C -> freq[A][C] += 3, freq[C][A] += 3
  // (1,3): AAA vs CCC -> 3x A-C -> freq[A][C] += 3, freq[C][A] += 3
  // (2,3): CCC vs CCC -> 3x C-C -> freq[C][C] += 6
  // Total: A-A=6, C-C=6, A-C=12, C-A=12 -> total=36
  // letter counts: A=18, C=18, p(A)=p(C)=0.5
  // q(A,A)=6/36=1/6, expected=0.25, ratio=2/3 < 1 -> negative diagonal
  // q(A,C)=12/36=1/3, expected=0.25, ratio=4/3 > 1 -> positive off-diagonal
  let aa = m.get_score("A", "A")
  let cc = m.get_score("C", "C")
  let ac = m.get_score("A", "C")
  // Diagonal is negative (matches less frequent than expected due to cross pairs)
  assert_true(aa < 0.0)
  assert_true(cc < 0.0)
  // Off-diagonal is positive (substitutions more frequent than expected)
  assert_true(ac > 0.0)
  // Symmetric
  assert_true((aa - cc).abs() < 1.0e-9)
  assert_true((ac - m.get_score("C", "A")).abs() < 1.0e-9)
}
