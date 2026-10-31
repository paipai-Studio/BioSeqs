///|
test "tree_construction_distance_calculator" {
  let calc = @src.DistanceCalculator::new()

  let seq1 = "ACGTACGTACGT"
  let seq2 = "ACGTACGTACGT"
  let dist_identical = calc.pairwise(seq1, seq2)
  assert_true(dist_identical >= 0.0 && dist_identical < 0.001)

  let seq3 = "ACGTACGTACGT"
  let seq4 = "ACGTACGTATGT"
  let dist_diff = calc.pairwise(seq3, seq4)
  assert_true(dist_diff > 0.0)
}

///|
test "tree_construction_upgma" {
  let tree_cons = @src.TreeConstructor::new(tree_method="upgma")

  let dist_matrix = [
    [0.0, 5.0, 9.0, 9.0, 8.0],
    [5.0, 0.0, 10.0, 10.0, 9.0],
    [9.0, 10.0, 0.0, 7.0, 3.0],
    [9.0, 10.0, 7.0, 0.0, 4.0],
    [8.0, 9.0, 3.0, 4.0, 0.0],
  ]
  let names = ["A", "B", "C", "D", "E"]

  let tree = tree_cons.build_tree(dist_matrix, names)
  assert_true(tree.clades.length() == 2)
}

///|
test "tree_construction_nj" {
  let tree_cons = @src.TreeConstructor::new(tree_method="nj")

  let dist_matrix = [
    [0.0, 5.0, 9.0, 9.0, 8.0],
    [5.0, 0.0, 10.0, 10.0, 9.0],
    [9.0, 10.0, 0.0, 7.0, 3.0],
    [9.0, 10.0, 7.0, 0.0, 4.0],
    [8.0, 9.0, 3.0, 4.0, 0.0],
  ]
  let names = ["A", "B", "C", "D", "E"]

  let tree = tree_cons.build_tree(dist_matrix, names)
  assert_true(tree.clades.length() >= 2)
}

///|
test "tree_construction_wpgma" {
  let tree_cons = @src.TreeConstructor::new(tree_method="wpgma")

  let dist_matrix = [
    [0.0, 2.0, 6.0, 6.0],
    [2.0, 0.0, 6.0, 6.0],
    [6.0, 6.0, 0.0, 2.0],
    [6.0, 6.0, 2.0, 0.0],
  ]
  let names = ["W", "X", "Y", "Z"]

  let tree = tree_cons.build_tree(dist_matrix, names)
  assert_true(tree.clades.length() == 2)
}

///|
test "tree_construction_jukes_cantor" {
  let calc = @src.DistanceCalculator::new(model="jukes_cantor")

  let seq1 = "ACGTACGTACGTACGTACGT"
  let seq2 = "ACGTACGTACGTACGTACGT"
  let dist = calc.pairwise(seq1, seq2)
  assert_true(dist >= 0.0 && dist < 0.01)
}

///|
test "tree_construction_kimura" {
  let calc = @src.DistanceCalculator::new(model="kimura")

  let seq1 = "AAAAAAAAAA"
  let seq2 = "AGAAAAGAAA"
  let dist = calc.pairwise(seq1, seq2)
  assert_true(dist >= 0.0)
}

///|
test "tree_construction_distance_matrix" {
  let calc = @src.DistanceCalculator::new()
  let seqs = [("seq1", "AAAA"), ("seq2", "AAAT"), ("seq3", "AATT")]
  let matrix = calc.get_distance(seqs)
  assert_eq(matrix.length(), 3)
  assert_eq(matrix[0].length(), 3)
  assert_true(matrix[0][0] == 0.0)
  assert_true(matrix[0][1] >= 0.0)
}

// ===========================================================================
// Biopython 1.88 DistanceCalculator compatibility tests.
// All expected values were computed with the local Biopython 1.88
// Bio.Phylo.TreeConstruction.DistanceCalculator (TreeConstruction/msa.phy
// doc-test alignment and additional cases).
// ===========================================================================

///|
fn tree_dist_close(a : Double, b : Double) -> Bool {
  (a - b).abs() < 1.0e-9
}

///|
fn tree_dist_close6(a : Double, b : Double) -> Bool {
  (a - b).abs() < 1.0e-6
}

///|
fn tree_dist_raises(model : String, s1 : String, s2 : String) -> String {
  try {
    let calc = @src.DistanceCalculator::new(model~)
    ignore(calc.pairwise(s1, s2))
    "NO_ERROR"
  } catch {
    SeqError(msg) => msg
  }
}

///|
test "tree_dist_identity_msa_phylip_matches_biopython" {
  // Biopython doctest: TreeConstruction/msa.phy
  let seqs = [
    ("Alpha", "AACGTGGCCACAT"),
    ("Beta", "AAGGTCGCCACAC"),
    ("Gamma", "CAGTTCGCCACAA"),
    ("Delta", "GAGATTTCCGCCT"),
    ("Epsilon", "GAGATCTCCGCCC"),
  ]
  let expected = [
    [0.0, 0.230769, 0.384615, 0.538462, 0.615385],
    [0.230769, 0.0, 0.230769, 0.538462, 0.384615],
    [0.384615, 0.230769, 0.0, 0.538462, 0.461538],
    [0.538462, 0.538462, 0.538462, 0.0, 0.153846],
    [0.615385, 0.384615, 0.461538, 0.153846, 0.0],
  ]
  let calc = @src.DistanceCalculator::new(model="identity")
  let matrix = calc.get_distance(seqs)
  for i in 0..<5 {
    for j in 0..<5 {
      assert_true(tree_dist_close6(matrix[i][j], expected[i][j]))
    }
  }
  // Identity skips no letters: gaps are compared like any character
  assert_true(tree_dist_close(calc.pairwise("A--A", "ACCA"), 0.5))
}

///|
test "tree_dist_blosum62_msa_phylip_matches_biopython" {
  let seqs = [
    ("Alpha", "AACGTGGCCACAT"),
    ("Beta", "AAGGTCGCCACAC"),
    ("Gamma", "CAGTTCGCCACAA"),
    ("Delta", "GAGATTTCCGCCT"),
    ("Epsilon", "GAGATCTCCGCCC"),
  ]
  let expected = [
    [0.0, 0.369048, 0.493976, 0.585366, 0.7],
    [0.369048, 0.0, 0.25, 0.547619, 0.355556],
    [0.493976, 0.25, 0.0, 0.566265, 0.488889],
    [0.585366, 0.547619, 0.566265, 0.0, 0.222222],
    [0.7, 0.355556, 0.488889, 0.222222, 0.0],
  ]
  let calc = @src.DistanceCalculator::new(model="blosum62")
  let matrix = calc.get_distance(seqs)
  for i in 0..<5 {
    for j in 0..<5 {
      assert_true(tree_dist_close6(matrix[i][j], expected[i][j]))
    }
  }
}

///|
test "tree_dist_identity_uses_len_seq1_as_max_score" {
  let calc = @src.DistanceCalculator::new(model="identity")
  // max_score is len(seq1) even when seq2 is shorter (zip semantics)
  assert_true(tree_dist_close(calc.pairwise("AAAA", "AA"), 0.5))
  // empty seq1 -> distance 1
  assert_true(tree_dist_close(calc.pairwise("", "AA"), 1.0))
}

///|
test "tree_dist_dna_matrix_models_match_biopython" {
  let blastn = @src.DistanceCalculator::new(model="blastn")
  assert_true(tree_dist_close(blastn.pairwise("ACGT", "ACGT"), 0.0))
  assert_true(tree_dist_close(blastn.pairwise("ACGT", "ACGA"), 0.45))
  assert_true(tree_dist_close(blastn.pairwise("ACGT", "TCGA"), 0.9))
  // "blastn" is the canonical NUC.4.4 IUPAC matrix
  let nuc44 = @src.DistanceCalculator::new(model="NUC.4.4")
  assert_true(
    tree_dist_close(
      nuc44.pairwise("ACGT", "TCGA"),
      blastn.pairwise("ACGT", "TCGA"),
    ),
  )
  let hoxd70 = @src.DistanceCalculator::new(model="hoxd70")
  assert_true(
    tree_dist_close(hoxd70.pairwise("ACGT", "ACGA"), 0.5602094240837696),
  )
  let trans = @src.DistanceCalculator::new(model="trans")
  assert_true(tree_dist_close(trans.pairwise("ACGT", "TCGA"), 0.5))
  let megablast = @src.DistanceCalculator::new(model="megablast")
  assert_true(tree_dist_close(megablast.pairwise("ACGT", "ACGA"), 0.75))
}

///|
test "tree_dist_protein_matrix_models_match_biopython" {
  let blosum62 = @src.DistanceCalculator::new(model="blosum62")
  assert_true(
    tree_dist_close(blosum62.pairwise("ARNDCQ", "ARNDCE"), 0.08571428571428574),
  )
  let pam250 = @src.DistanceCalculator::new(model="pam250")
  assert_true(
    tree_dist_close(
      pam250.pairwise("ACDEFGHIKLMNPQRSTVWY", "ACDEFGHIKLMNPQRSTVWF"),
      0.025423728813559365,
    ),
  )
  let dayhoff = @src.DistanceCalculator::new(model="dayhoff")
  assert_true(
    tree_dist_close(dayhoff.pairwise("ACDEFGH", "ACDEFGW"), 0.37735849056603776),
  )
  let gonnet = @src.DistanceCalculator::new(model="gonnet1992")
  assert_true(tree_dist_close(gonnet.pairwise("ACW", "ACW"), 0.0))
  assert_true(
    tree_dist_close(gonnet.pairwise("ACW", "ATW"), 0.4270462633451958),
  )
  let blastp = @src.DistanceCalculator::new(model="blastp")
  assert_true(
    tree_dist_close(
      blastp.pairwise("ACDEFGHIKLMNPQRSTVWY", "ACDEFGHIKLMNPQRSTVWB"),
      0.08620689655172409,
    ),
  )
  let benner22 = @src.DistanceCalculator::new(model="benner22")
  assert_true(
    tree_dist_close(benner22.pairwise("ACW", "ATW"), 0.44480519480519487),
  )
  let jones = @src.DistanceCalculator::new(model="jones")
  assert_true(tree_dist_close(jones.pairwise("ACW", "ATW"), 0.4285714285714285))
}

///|
test "tree_dist_matrix_models_skip_gap_and_star" {
  let calc = @src.DistanceCalculator::new(model="blosum62")
  // Columns containing '-' or '*' are ignored entirely
  assert_true(tree_dist_close(calc.pairwise("AC-DE", "ACWDE"), 0.0))
  assert_true(tree_dist_close(calc.pairwise("AC*DE", "ACWDE"), 0.0))
  // Every column skipped -> max_score 0 -> distance 1
  assert_true(tree_dist_close(calc.pairwise("----", "----"), 1.0))
  assert_true(tree_dist_close(calc.pairwise("****", "WWWW"), 1.0))
  // Empty sequence -> distance 1
  assert_true(tree_dist_close(calc.pairwise("", ""), 1.0))
}

///|
test "tree_dist_custom_skip_letters" {
  // Explicit skip_letters overrides defaults (max_score stays len(seq1))
  let calc = @src.DistanceCalculator::new(
    model="identity",
    skip_letters=Some(["#"]),
  )
  assert_true(tree_dist_close(calc.pairwise("A#A", "A#A"), 1.0 / 3.0))
}

///|
test "tree_dist_iupac_matrix_can_exceed_one" {
  // Biopython does not clamp distances; low/negative IUPAC self-scores
  // can give distances above 1.
  let blastn = @src.DistanceCalculator::new(model="blastn")
  assert_true(tree_dist_close(blastn.pairwise("ATGC", "SSWW"), 1.8))
}

///|
test "tree_dist_bad_letter_raises_seq_error" {
  let msg1 = tree_dist_raises("blosum62", "ACWJ", "ACWY")
  assert_true(msg1.contains("Bad letter 'J'"))
  assert_true(msg1.contains("position '3'"))
  let msg2 = tree_dist_raises("blastn", "ACGX", "ACGA")
  assert_true(msg2.contains("Bad letter 'X'"))
  // A bad letter in a column skipped via the other sequence does not raise
  let calc = @src.DistanceCalculator::new(model="blosum62")
  assert_true(tree_dist_close(calc.pairwise("AC-J", "AC--"), 0.0))
}

///|
test "tree_dist_unknown_model_raises_seq_error" {
  let msg = tree_dist_raises("not_a_model", "AC", "AC")
  assert_true(msg.contains("Model not supported. Available models:"))
  assert_true(msg.contains("blosum62"))
}

///|
test "tree_dist_model_name_lists_match_biopython" {
  let dna = @src.distance_dna_models()
  let protein = @src.distance_protein_models()
  let models = @src.distance_models()
  let expected_dna = [
    "benner22", "benner6", "benner74", "blastn", "dayhoff", "feng", "genetic", "gonnet1992",
    "hoxd70", "johnson", "jones", "levin", "mclachlan", "mdm78", "megablast", "blastn",
    "rao", "risler", "schneider", "str", "trans",
  ]
  let expected_protein = [
    "blastp", "blosum45", "blosum50", "blosum62", "blosum80", "blosum90", "pam250",
    "pam30", "pam70",
  ]
  assert_eq(dna, expected_dna)
  assert_eq(protein, expected_protein)
  assert_eq(models.length(), 31)
  assert_eq(models[0], "identity")
  assert_eq(models[models.length() - 1], "pam70")
}

///|
test "tree_dist_model_name_and_jc_aliases" {
  let calc = @src.DistanceCalculator::new(model="blosum62")
  assert_eq(calc.model_name(), "blosum62")
  // MoonBit extension models (not in Biopython 1.88) keep working under
  // all supported aliases.
  let jc1 = @src.DistanceCalculator::new(model="jc69")
  let jc2 = @src.DistanceCalculator::new(model="jukes-cantor")
  let jc3 = @src.DistanceCalculator::new(model="jukes_cantor")
  let s1 = "ACGTACGT"
  let s2 = "ACGTACGA"
  assert_true(tree_dist_close(jc1.pairwise(s1, s2), jc2.pairwise(s1, s2)))
  assert_true(tree_dist_close(jc2.pairwise(s1, s2), jc3.pairwise(s1, s2)))
  let k2p = @src.DistanceCalculator::new(model="k2p")
  let kimura = @src.DistanceCalculator::new(model="kimura")
  assert_true(tree_dist_close(k2p.pairwise(s1, s2), kimura.pairwise(s1, s2)))
}

///|
test "tree_dist_get_distance_is_symmetric" {
  let calc = @src.DistanceCalculator::new(model="blosum62")
  let seqs = [("a", "ARNDCQ"), ("b", "ARNDCE"), ("c", "ARNDCC")]
  let matrix = calc.get_distance(seqs)
  for i in 0..<3 {
    for j in 0..<3 {
      assert_true(tree_dist_close(matrix[i][j], matrix[j][i]))
    }
    assert_true(tree_dist_close(matrix[i][i], 0.0))
  }
}
