///|
/// Tests for the Biopython Bio.Align MSA statistics port:
/// MultipleSeqAlignment.substitutions (symmetric weighted substitution
/// count matrix), Alignment.frequencies (per-column letter counts with
/// a '-' gap row) and the relative pair frequencies q_ij / e_ij.
/// Reference values are the exact numbers printed in the Biopython
/// documentation (plus hand-computed weighted cases).

///|
fn afreq_close(a : Double, b : Double, tol? : Double = 1.0e-9) -> Bool {
  (a - b).abs() < tol
}

///|
fn afreq_val(m : Map[String, Double], key : String) -> Double {
  match m.get(key) {
    Some(v) => v
    None => abort("missing map key " + key)
  }
}

///|
fn afreq_mat_close(
  a : Array[Array[Double]],
  b : Array[Array[Double]],
  tol? : Double = 1.0e-9,
) -> Bool {
  if a.length() != b.length() {
    return false
  }
  let mut i = 0
  while i < a.length() {
    if a[i].length() != b[i].length() {
      return false
    }
    let mut j = 0
    while j < a[i].length() {
      if not(afreq_close(a[i][j], b[i][j], tol~)) {
        return false
      }
      j = j + 1
    }
    i = i + 1
  }
  true
}

///|
/// Biopython doc example for MultipleSeqAlignment.substitutions:
///   ACGT / A--A / ACGT / TTTC  ->
///          A    C    G    T
///   A    3.0  0.5  0.0  2.5
///   C    0.5  1.0  0.0  2.0
///   G    0.0  0.0  1.0  1.0
///   T    2.5  2.0  1.0  1.0
test "alignment_substitutions_biopython_doc" {
  let rows : Array[String] = ["ACGT", "A--A", "ACGT", "TTTC"]
  let m = @src.alignment_substitution_matrix(rows) catch {
    _ => abort("doc example substitution matrix should build")
  }
  assert_eq(m.letters, ["A", "C", "G", "T"])
  let expected : Array[Array[Double]] = [
    [3.0, 0.5, 0.0, 2.5],
    [0.5, 1.0, 0.0, 2.0],
    [0.0, 0.0, 1.0, 1.0],
    [2.5, 2.0, 1.0, 1.0],
  ]
  assert_true(afreq_mat_close(m.matrix, expected))
  // symmetry
  assert_true(afreq_close(m.count_of("A", "T"), m.count_of("T", "A")))
  assert_true(afreq_close(m.count_of("C", "T"), m.count_of("T", "C")))
  // total mass: unordered pair contributions over non-gap columns
  assert_true(afreq_close(m.total(), 18.0))
}

///|
test "alignment_substitutions_simple_pair" {
  // Two sequences ACGT vs ACGA: A-A, C-C, G-G matches, plus one
  // A/T pair split symmetrically as 0.5/0.5.
  let m = @src.alignment_substitution_matrix(["ACGT", "ACGA"]) catch {
    _ => abort("simple substitution matrix should build")
  }
  let expected : Array[Array[Double]] = [
    [1.0, 0.0, 0.0, 0.5],
    [0.0, 1.0, 0.0, 0.0],
    [0.0, 0.0, 1.0, 0.0],
    [0.5, 0.0, 0.0, 0.0],
  ]
  assert_true(afreq_mat_close(m.matrix, expected))
  assert_true(afreq_close(m.total(), 4.0))
}

///|
/// Weighted doc example: GTATC (weight 0.5), AT--C (0.8), CTGTC (1.0).
/// First column pairs: (A,G) = 0.5*0.8 = 0.4, (C,G) = 0.5*1.0 = 0.5,
/// (A,C) = 0.8*1.0 = 0.8; after symmetrization each unordered pair is
/// split over the two symmetric cells.
test "alignment_substitutions_weighted" {
  let rows : Array[String] = ["GTATC", "AT--C", "CTGTC"]
  let weights : Array[Double] = [0.5, 0.8, 1.0]
  let m = @src.alignment_substitution_matrix(rows, weights~) catch {
    _ => abort("weighted substitution matrix should build")
  }
  let expected : Array[Array[Double]] = [
    [0.0, 0.4, 0.45, 0.0],
    [0.4, 1.7, 0.25, 0.0],
    [0.45, 0.25, 0.0, 0.0],
    [0.0, 0.0, 0.0, 2.2],
  ]
  assert_true(afreq_mat_close(m.matrix, expected))
  // first column unordered contributions
  assert_true(afreq_close(m.count_of("A", "G"), 0.45))
  assert_true(afreq_close(m.count_of("G", "C"), 0.25))
  assert_true(afreq_close(m.count_of("A", "C"), 0.4))
  // total weighted pair mass: 1.7 + 1.7 + 0.5 + 0.5 + 1.7 = 6.1
  assert_true(afreq_close(m.total(), 6.1))
  // diagonal holds total self-pair weight per residue
  assert_true(afreq_close(m.count_of("T", "T"), 2.2))
  assert_true(afreq_close(m.count_of("C", "C"), 1.7))
}

///|
/// Per-column counts for the four-sequence doc MSA. Gaps appear in the
/// '-' row; raw counts are integers (unit weights).
test "alignment_letter_counts_doc" {
  let rows : Array[String] = ["ACGT", "A--A", "ACGT", "TTTC"]
  let counts = @src.alignment_letter_counts(rows) catch {
    _ => abort("letter counts should build")
  }
  assert_eq(counts.letters, ["-", "A", "C", "G", "T"])
  let expected : Array[Array[Double]] = [
    [0.0, 1.0, 1.0, 0.0], // '-'
    [3.0, 0.0, 0.0, 1.0], // A
    [0.0, 2.0, 0.0, 1.0], // C
    [0.0, 0.0, 2.0, 0.0], // G
    [1.0, 1.0, 1.0, 2.0], // T
  ]
  assert_true(afreq_mat_close(counts.counts, expected))
  assert_true(afreq_close(counts.count_of("A", 0), 3.0))
  assert_true(afreq_close(counts.count_of("-", 1), 1.0))
  assert_true(afreq_close(counts.count_of("X", 0), 0.0))
  // column-normalized frequencies: every column sums to 1
  let freq = counts.column_frequencies()
  let mut j = 0
  while j < counts.ncols {
    let mut s = 0.0
    let mut k = 0
    while k < counts.letters.length() {
      s += freq[k][j]
      k = k + 1
    }
    assert_true(afreq_close(s, 1.0))
    j = j + 1
  }
  // spot checks matching the column composition
  assert_true(afreq_close(freq[1][0], 0.75)) // A at col 0
  assert_true(afreq_close(freq[4][0], 0.25)) // T at col 0
  assert_true(afreq_close(freq[0][1], 0.25)) // gap at col 1
  assert_true(afreq_close(freq[2][1], 0.5)) // C at col 1
}

///|
/// Biopython doc example for Alignment.frequencies on a global pairwise
/// alignment: -GACCTG vs CGATC-G (7 aligned columns, gaps as '-').
test "alignment_letter_counts_global_pair" {
  let rows : Array[String] = ["-GACCTG", "CGATC-G"]
  let counts = @src.alignment_letter_counts(rows) catch {
    _ => abort("pairwise letter counts should build")
  }
  assert_eq(counts.letters, ["-", "A", "C", "G", "T"])
  let expected : Array[Array[Double]] = [
    [1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0], // '-'
    [0.0, 0.0, 2.0, 0.0, 0.0, 0.0, 0.0], // A
    [1.0, 0.0, 0.0, 1.0, 2.0, 0.0, 0.0], // C
    [0.0, 2.0, 0.0, 0.0, 0.0, 0.0, 2.0], // G
    [0.0, 0.0, 0.0, 1.0, 0.0, 1.0, 0.0], // T
  ]
  assert_true(afreq_mat_close(counts.counts, expected))
  // '.' gaps are reported under '-' as well
  let dot = @src.alignment_letter_counts(["A.T", "A-G"]) catch {
    _ => abort("dot-gap letter counts should build")
  }
  assert_true(afreq_close(dot.count_of("-", 1), 2.0))
}

///|
/// A local (ungapped) alignment has no '-' row; every column sums to the
/// number of sequences.
test "alignment_letter_counts_local_no_gap" {
  let rows : Array[String] = ["GACCTG", "GATCTG"]
  let counts = @src.alignment_letter_counts(rows) catch {
    _ => abort("local letter counts should build")
  }
  assert_eq(counts.letters, ["A", "C", "G", "T"])
  let expected : Array[Array[Double]] = [
    [0.0, 2.0, 0.0, 0.0, 0.0, 0.0],
    [0.0, 0.0, 1.0, 2.0, 0.0, 0.0],
    [2.0, 0.0, 0.0, 0.0, 0.0, 2.0],
    [0.0, 0.0, 1.0, 0.0, 2.0, 0.0],
  ]
  assert_true(afreq_mat_close(counts.counts, expected))
  let freq = counts.column_frequencies()
  assert_true(afreq_close(freq[2][0], 1.0)) // G at col 0
  assert_true(afreq_close(freq[1][2], 0.5)) // C at col 2
  assert_true(afreq_close(freq[3][2], 0.5)) // T at col 2
  // weighted column counts
  let wcounts = @src.alignment_letter_counts(rows, weights=[0.5, 1.5]) catch {
    _ => abort("weighted local letter counts should build")
  }
  assert_true(afreq_close(wcounts.count_of("G", 0), 2.0))
  assert_true(afreq_close(wcounts.count_of("C", 2), 0.5))
  assert_true(afreq_close(wcounts.count_of("T", 2), 1.5))
}

///|
/// Background letter frequencies: column-count based (gaps excluded) and
/// pair-count based (BLOSUM-style row sums of q_ij).
test "alignment_background_frequencies" {
  let rows : Array[String] = ["ACGT", "A--A", "ACGT", "TTTC"]
  let counts = @src.alignment_letter_counts(rows) catch {
    _ => abort("letter counts should build")
  }
  let col_bg = counts.background_frequencies()
  // 14 non-gap residues: A=4, C=3, G=2, T=5
  assert_true(afreq_close(afreq_val(col_bg, "A"), 4.0 / 14.0))
  assert_true(afreq_close(afreq_val(col_bg, "C"), 3.0 / 14.0))
  assert_true(afreq_close(afreq_val(col_bg, "G"), 2.0 / 14.0))
  assert_true(afreq_close(afreq_val(col_bg, "T"), 5.0 / 14.0))
  assert_true(col_bg.get("-").is_none())

  let m = @src.alignment_substitution_matrix(rows) catch {
    _ => abort("substitution matrix should build")
  }
  let pair_bg = m.background_frequencies()
  // q row sums: A = 6/18, C = 3.5/18, G = 2/18, T = 6.5/18
  assert_true(afreq_close(afreq_val(pair_bg, "A"), 6.0 / 18.0))
  assert_true(afreq_close(afreq_val(pair_bg, "C"), 3.5 / 18.0))
  assert_true(afreq_close(afreq_val(pair_bg, "G"), 2.0 / 18.0))
  assert_true(afreq_close(afreq_val(pair_bg, "T"), 6.5 / 18.0))
  let mut psum = 0.0
  for kv in pair_bg {
    psum += kv.1
  }
  assert_true(afreq_close(psum, 1.0))
}

///|
/// Relative pair frequencies q_ij / e_ij with the MLE background and
/// with a user-supplied uniform background.
test "alignment_relative_frequencies" {
  let rows : Array[String] = ["ACGT", "A--A", "ACGT", "TTTC"]
  let m = @src.alignment_substitution_matrix(rows) catch {
    _ => abort("substitution matrix should build")
  }
  let rel = m.relative_frequencies() catch {
    _ => abort("relative frequencies should build")
  }
  // letters A,C,G,T; hand values:
  //   AA: (3/18)/(1/3)^2      = 1.5
  //   GG: (1/18)/(1/9)^2      = 4.5
  //   AG: 0
  //   TT: (1/18)/(13/36)^2    = 72/169
  assert_true(afreq_close(rel[0][0], 1.5))
  assert_true(afreq_close(rel[2][2], 4.5))
  assert_true(afreq_close(rel[0][2], 0.0))
  assert_true(afreq_close(rel[3][3], 72.0 / 169.0, tol=1.0e-9))
  // symmetric
  assert_true(afreq_close(rel[0][3], rel[3][0]))

  let uniform : Map[String, Double] = {
    "A": 0.25,
    "C": 0.25,
    "G": 0.25,
    "T": 0.25,
  }
  let relu = m.relative_frequencies_expected(uniform) catch {
    _ => abort("relative frequencies with uniform background should build")
  }
  // e_ij = 1/16: AA = (3/18)*16 = 8/3, GG = (1/18)*16 = 8/9
  assert_true(afreq_close(relu[0][0], 8.0 / 3.0))
  assert_true(afreq_close(relu[2][2], 8.0 / 9.0))
  assert_true(afreq_close(relu[0][2], 0.0))
}

///|
test "alignment_frequency_errors" {
  // empty alignment
  let empty = try {
    ignore(@src.alignment_substitution_matrix([]))
    false
  } catch {
    AlignmentFrequenciesError(_) => true
  }
  assert_true(empty)
  // ragged rows
  let ragged = try {
    ignore(@src.alignment_letter_counts(["ACG", "ACGT"]))
    false
  } catch {
    AlignmentFrequenciesError(_) => true
  }
  assert_true(ragged)
  // empty row
  let empty_row = try {
    ignore(@src.alignment_letter_counts(["", "ACGT"]))
    false
  } catch {
    AlignmentFrequenciesError(_) => true
  }
  assert_true(empty_row)
  // weight length mismatch
  let bad_weights = try {
    ignore(@src.alignment_substitution_matrix(["ACGT", "ACGA"], weights=[0.5]))
    false
  } catch {
    AlignmentFrequenciesError(_) => true
  }
  assert_true(bad_weights)
  // unknown residue lookup
  let bad_residue = try {
    let m = @src.alignment_substitution_matrix(["ACGT", "ACGA"])
    ignore(m.count_of("A", "Z"))
    false
  } catch {
    AlignmentFrequenciesError(_) => true
  }
  assert_true(bad_residue)
  // missing expected background key
  let missing_bg = try {
    let m = @src.alignment_substitution_matrix(["ACGT", "ACGA"])
    let bg : Map[String, Double] = { "A": 0.25, "C": 0.25, "G": 0.25 }
    ignore(m.relative_frequencies_expected(bg))
    false
  } catch {
    AlignmentFrequenciesError(_) => true
  }
  assert_true(missing_bg)
  // non-positive expected frequency
  let bad_bg = try {
    let m = @src.alignment_substitution_matrix(["ACGT", "ACGA"])
    let bg : Map[String, Double] = Map([
      ("A", 0.25),
      ("C", 0.25),
      ("G", 0.25),
      ("T", 0.0),
    ])
    ignore(m.relative_frequencies_expected(bg))
    false
  } catch {
    AlignmentFrequenciesError(_) => true
  }
  assert_true(bad_bg)
}
