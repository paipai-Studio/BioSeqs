///|
test "AbstractAlignment::new construction" {
  let sequences = ["ATCGATCG", "ATCGATCG"]
  let identifiers = ["seq1", "seq2"]
  let aln = @src.AbstractAlignment::new(
    sequences=sequences,
    identifiers=identifiers,
    alignment_type=@src.align_abstract_type_from_str("nucleotide"),
  )
  assert_eq(aln.abstract_n_seqs(), 2)
  assert_eq(aln.abstract_alignment_length(), 8)
}

///|
test "AbstractAlignment::validate success" {
  let sequences = ["ATCGATCG", "ATCGATCG"]
  let identifiers = ["seq1", "seq2"]
  let aln = @src.AbstractAlignment::new(
    sequences=sequences,
    identifiers=identifiers,
    alignment_type=@src.align_abstract_type_from_str("nucleotide"),
  )
  let (valid, msg) = aln.validate()
  assert_true(valid)
  assert_true(msg.length() > 0)
}

///|
test "AbstractAlignment::validate different lengths" {
  let sequences = ["ATCG", "ATCGATCG"]
  let identifiers = ["seq1", "seq2"]
  let aln = @src.AbstractAlignment::new(
    sequences=sequences,
    identifiers=identifiers,
    alignment_type=@src.align_abstract_type_from_str("nucleotide"),
  )
  let (valid, _msg) = aln.validate()
  assert_false(valid)
}

///|
test "AbstractAlignment::validate empty alignment" {
  let sequences : Array[String] = []
  let identifiers : Array[String] = []
  let aln = @src.AbstractAlignment::new(
    sequences=sequences,
    identifiers=identifiers,
    alignment_type=@src.align_abstract_type_from_str("nucleotide"),
  )
  let (valid, msg) = aln.validate()
  assert_false(valid)
  assert_true(msg.contains("No sequences"))
}

///|
test "AbstractAlignment::validate invalid nucleotide character" {
  let sequences = ["ATCGNTCG", "ATCGXTCG"]
  let identifiers = ["seq1", "seq2"]
  let aln = @src.AbstractAlignment::new(
    sequences=sequences,
    identifiers=identifiers,
    alignment_type=@src.align_abstract_type_from_str("nucleotide"),
  )
  let (valid, _msg) = aln.validate()
  assert_false(valid)
}

///|
test "AbstractAlignment::validate protein" {
  let sequences = ["ALIGMENT", "ALIGMENT"]
  let identifiers = ["seq1", "seq2"]
  let aln = @src.AbstractAlignment::new(
    sequences=sequences,
    identifiers=identifiers,
    alignment_type=@src.align_abstract_type_from_str("protein"),
  )
  let (valid, _msg) = aln.validate()
  assert_true(valid)
}

///|
test "AbstractAlignment::validate protein invalid char" {
  let sequences = ["ALIGMENT", "ALIGM1NT"]
  let identifiers = ["seq1", "seq2"]
  let aln = @src.AbstractAlignment::new(
    sequences=sequences,
    identifiers=identifiers,
    alignment_type=@src.align_abstract_type_from_str("protein"),
  )
  let (valid, _msg) = aln.validate()
  assert_false(valid)
}

///|
test "abstract_column_stats identical columns" {
  let sequences = ["AAAA", "AAAA"]
  let identifiers = ["seq1", "seq2"]
  let aln = @src.AbstractAlignment::new(
    sequences=sequences,
    identifiers=identifiers,
    alignment_type=@src.align_abstract_type_from_str("nucleotide"),
  )
  let stats = @src.abstract_column_stats(aln, 0)
  let diff_cons = (stats.conservation - 1.0).abs()
  assert_true(diff_cons < 1.0e-10)
  let gap_val = (stats.gap_fraction).abs()
  assert_true(gap_val < 1.0e-10)
  assert_eq(stats.n_unique_chars, 1)
}

///|
test "abstract_column_stats with gaps" {
  let sequences = ["A-CG", "AC-G"]
  let identifiers = ["seq1", "seq2"]
  let aln = @src.AbstractAlignment::new(
    sequences=sequences,
    identifiers=identifiers,
    alignment_type=@src.align_abstract_type_from_str("nucleotide"),
  )
  let stats = @src.abstract_column_stats(aln, 1)
  assert_true(stats.gap_fraction > 0.0)
}

///|
test "abstract_column_stats diversity" {
  let sequences = ["ACGT", "TGCA"]
  let identifiers = ["seq1", "seq2"]
  let aln = @src.AbstractAlignment::new(
    sequences=sequences,
    identifiers=identifiers,
    alignment_type=@src.align_abstract_type_from_str("nucleotide"),
  )
  let stats = @src.abstract_column_stats(aln, 0)
  assert_true(stats.diversity > 0.0)
}

///|
test "abstract_column_stats out of bounds" {
  let sequences = ["ACGT", "TGCA"]
  let identifiers = ["seq1", "seq2"]
  let aln = @src.AbstractAlignment::new(
    sequences=sequences,
    identifiers=identifiers,
    alignment_type=@src.align_abstract_type_from_str("nucleotide"),
  )
  let stats = @src.abstract_column_stats(aln, 100)
  let cons_val = (stats.conservation).abs()
  assert_true(cons_val < 1.0e-10)
  assert_eq(stats.n_unique_chars, 0)
}

///|
test "abstract_identity_matrix identical sequences" {
  let sequences = ["ATCGATCG", "ATCGATCG"]
  let identifiers = ["seq1", "seq2"]
  let aln = @src.AbstractAlignment::new(
    sequences=sequences,
    identifiers=identifiers,
    alignment_type=@src.align_abstract_type_from_str("nucleotide"),
  )
  let matrix = @src.abstract_identity_matrix(aln)
  assert_eq(matrix.length(), 2)
  assert_eq(matrix[0].length(), 2)
  let d00 = matrix[0][0]
  let d11 = matrix[1][1]
  let d01 = matrix[0][1]
  assert_true((d00 - 1.0).abs() < 1.0e-10)
  assert_true((d11 - 1.0).abs() < 1.0e-10)
  assert_true((d01 - 1.0).abs() < 1.0e-10)
}

///|
test "abstract_identity_matrix different sequences" {
  let sequences = ["AAAA", "TTTT"]
  let identifiers = ["seq1", "seq2"]
  let aln = @src.AbstractAlignment::new(
    sequences=sequences,
    identifiers=identifiers,
    alignment_type=@src.align_abstract_type_from_str("nucleotide"),
  )
  let matrix = @src.abstract_identity_matrix(aln)
  assert_true((matrix[0][1]).abs() < 1.0e-10)
}

///|
test "abstract_identity_matrix partial match" {
  let sequences = ["ACGT", "ACGA"]
  let identifiers = ["seq1", "seq2"]
  let aln = @src.AbstractAlignment::new(
    sequences=sequences,
    identifiers=identifiers,
    alignment_type=@src.align_abstract_type_from_str("nucleotide"),
  )
  let matrix = @src.abstract_identity_matrix(aln)
  assert_true((matrix[0][1] - 0.75).abs() < 1.0e-10)
}

///|
test "abstract_overall_identity" {
  let sequences = ["ACGT", "ACGA", "ACGT"]
  let identifiers = ["seq1", "seq2", "seq3"]
  let aln = @src.AbstractAlignment::new(
    sequences=sequences,
    identifiers=identifiers,
    alignment_type=@src.align_abstract_type_from_str("nucleotide"),
  )
  let identity = @src.abstract_overall_identity(aln)
  assert_true(identity > 0.0)
  assert_true(identity <= 1.0)
}

///|
test "abstract_overall_identity single sequence" {
  let sequences = ["ACGT"]
  let identifiers = ["seq1"]
  let aln = @src.AbstractAlignment::new(
    sequences=sequences,
    identifiers=identifiers,
    alignment_type=@src.align_abstract_type_from_str("nucleotide"),
  )
  let identity = @src.abstract_overall_identity(aln)
  assert_true((identity - 1.0).abs() < 1.0e-10)
}

///|
test "abstract_coverage" {
  let sequences = ["ACGT", "A-GT"]
  let identifiers = ["seq1", "seq2"]
  let aln = @src.AbstractAlignment::new(
    sequences=sequences,
    identifiers=identifiers,
    alignment_type=@src.align_abstract_type_from_str("nucleotide"),
  )
  let coverage = @src.abstract_coverage(aln)
  assert_true(coverage > 0.0)
  assert_true(coverage < 1.0)
}

///|
test "abstract_coverage no gaps" {
  let sequences = ["ACGT", "ACGT"]
  let identifiers = ["seq1", "seq2"]
  let aln = @src.AbstractAlignment::new(
    sequences=sequences,
    identifiers=identifiers,
    alignment_type=@src.align_abstract_type_from_str("nucleotide"),
  )
  let coverage = @src.abstract_coverage(aln)
  assert_true((coverage - 1.0).abs() < 1.0e-10)
}

///|
test "AbstractAlignment::abstract_get_column" {
  let sequences = ["ACGT", "TGCA"]
  let identifiers = ["seq1", "seq2"]
  let aln = @src.AbstractAlignment::new(
    sequences=sequences,
    identifiers=identifiers,
    alignment_type=@src.align_abstract_type_from_str("nucleotide"),
  )
  let col = aln.abstract_get_column(0)
  assert_eq(col.length(), 2)
}

///|
test "AbstractAlignment::abstract_get_column out of bounds" {
  let sequences = ["ACGT"]
  let identifiers = ["seq1"]
  let aln = @src.AbstractAlignment::new(
    sequences=sequences,
    identifiers=identifiers,
    alignment_type=@src.align_abstract_type_from_str("nucleotide"),
  )
  let col = aln.abstract_get_column(10)
  assert_eq(col.length(), 0)
}

///|
test "AbstractAlignment::abstract_get_seq" {
  let sequences = ["ATCGATCG"]
  let identifiers = ["seq1"]
  let aln = @src.AbstractAlignment::new(
    sequences=sequences,
    identifiers=identifiers,
    alignment_type=@src.align_abstract_type_from_str("nucleotide"),
  )
  assert_eq(aln.abstract_get_seq(0), "ATCGATCG")
}

///|
test "AbstractAlignment::abstract_get_seq out of bounds" {
  let sequences = ["ATCG"]
  let identifiers = ["seq1"]
  let aln = @src.AbstractAlignment::new(
    sequences=sequences,
    identifiers=identifiers,
    alignment_type=@src.align_abstract_type_from_str("nucleotide"),
  )
  assert_eq(aln.abstract_get_seq(5), "")
}

///|
test "AbstractAlignment::abstract_get_id" {
  let sequences = ["ATCG"]
  let identifiers = ["gene1"]
  let aln = @src.AbstractAlignment::new(
    sequences=sequences,
    identifiers=identifiers,
    alignment_type=@src.align_abstract_type_from_str("nucleotide"),
  )
  assert_eq(aln.abstract_get_id(0), "gene1")
}

///|
test "AbstractAlignment::abstract_type" {
  let sequences = ["ATCG"]
  let identifiers = ["seq1"]
  let aln = @src.AbstractAlignment::new(
    sequences=sequences,
    identifiers=identifiers,
    alignment_type=@src.align_abstract_type_from_str("nucleotide"),
  )
  let t = aln.abstract_type()
  assert_true(t is @src.AlignAbstractType::Nucleotide)
}

///|
test "abstract_variable_sites" {
  let sequences = ["ACGT", "ACGA"]
  let identifiers = ["seq1", "seq2"]
  let aln = @src.AbstractAlignment::new(
    sequences=sequences,
    identifiers=identifiers,
    alignment_type=@src.align_abstract_type_from_str("nucleotide"),
  )
  let num_sites = @src.abstract_variable_sites(aln)
  assert_eq(num_sites, 1)
}

///|
test "abstract_variable_sites none" {
  let sequences = ["ACGT", "ACGT"]
  let identifiers = ["seq1", "seq2"]
  let aln = @src.AbstractAlignment::new(
    sequences=sequences,
    identifiers=identifiers,
    alignment_type=@src.align_abstract_type_from_str("nucleotide"),
  )
  let num_sites = @src.abstract_variable_sites(aln)
  assert_eq(num_sites, 0)
}

///|
test "abstract_parsimony_sites" {
  let sequences = ["ACGT", "ACGA", "TCGT"]
  let identifiers = ["seq1", "seq2", "seq3"]
  let aln = @src.AbstractAlignment::new(
    sequences=sequences,
    identifiers=identifiers,
    alignment_type=@src.align_abstract_type_from_str("nucleotide"),
  )
  let parsimony = @src.abstract_parsimony_sites(aln)
  assert_true(parsimony >= 0)
}

///|
test "abstract_consensus_sequence" {
  let sequences = ["ACGT", "ACGA", "ACGT"]
  let identifiers = ["seq1", "seq2", "seq3"]
  let aln = @src.AbstractAlignment::new(
    sequences=sequences,
    identifiers=identifiers,
    alignment_type=@src.align_abstract_type_from_str("nucleotide"),
  )
  let consensus = @src.abstract_consensus_sequence(aln)
  assert_eq(consensus.length(), 4)
}

///|
test "abstract_filter_gaps" {
  let sequences = ["A-CGT", "A-CGA"]
  let identifiers = ["seq1", "seq2"]
  let aln = @src.AbstractAlignment::new(
    sequences=sequences,
    identifiers=identifiers,
    alignment_type=@src.align_abstract_type_from_str("nucleotide"),
  )
  let filtered = @src.abstract_filter_gaps(aln, 0.5)
  assert_true(filtered.abstract_alignment_length() <= aln.abstract_alignment_length())
}

///|
test "abstract_trim" {
  let sequences = ["ACGTACGT", "ACGTACGT"]
  let identifiers = ["seq1", "seq2"]
  let aln = @src.AbstractAlignment::new(
    sequences=sequences,
    identifiers=identifiers,
    alignment_type=@src.align_abstract_type_from_str("nucleotide"),
  )
  let trimmed = @src.abstract_trim(aln, 2, 5)
  assert_eq(trimmed.abstract_alignment_length(), 4)
  assert_eq(trimmed.abstract_get_seq(0), "GTAC")
}

///|
test "abstract_distance_matrix" {
  let sequences = ["ACGT", "ACGA"]
  let identifiers = ["seq1", "seq2"]
  let aln = @src.AbstractAlignment::new(
    sequences=sequences,
    identifiers=identifiers,
    alignment_type=@src.align_abstract_type_from_str("nucleotide"),
  )
  let matrix = @src.abstract_distance_matrix(aln)
  assert_eq(matrix.length(), 2)
  assert_true((matrix[0][0]).abs() < 1.0e-10)
  assert_true((matrix[1][1]).abs() < 1.0e-10)
}

///|
test "align_abstract_type_from_str nucleotide" {
  let t = @src.align_abstract_type_from_str("nucleotide")
  assert_true(t is @src.AlignAbstractType::Nucleotide)
}

///|
test "align_abstract_type_from_str protein" {
  let t = @src.align_abstract_type_from_str("protein")
  assert_true(t is @src.AlignAbstractType::Protein)
}

///|
test "align_abstract_type_from_str generic" {
  let t = @src.align_abstract_type_from_str("unknown")
  assert_true(t is @src.AlignAbstractType::Generic)
}

///|
test "AbstractAlignment::multiple sequences" {
  let sequences = ["AAAA", "CCCC", "GGGG", "TTTT"]
  let identifiers = ["s1", "s2", "s3", "s4"]
  let aln = @src.AbstractAlignment::new(
    sequences=sequences,
    identifiers=identifiers,
    alignment_type=@src.align_abstract_type_from_str("nucleotide"),
  )
  assert_eq(aln.abstract_n_seqs(), 4)
  assert_eq(aln.abstract_alignment_length(), 4)
}

///|
test "abstract_singleton_sites" {
  let sequences = ["ACGT", "ACGT", "ACGA"]
  let identifiers = ["seq1", "seq2", "seq3"]
  let aln = @src.AbstractAlignment::new(
    sequences=sequences,
    identifiers=identifiers,
    alignment_type=@src.align_abstract_type_from_str("nucleotide"),
  )
  let singleton = @src.abstract_singleton_sites(aln)
  assert_true(singleton >= 0)
}

///|
test "abstract_summary" {
  let sequences = ["ACGT", "ACGT"]
  let identifiers = ["seq1", "seq2"]
  let aln = @src.AbstractAlignment::new(
    sequences=sequences,
    identifiers=identifiers,
    alignment_type=@src.align_abstract_type_from_str("nucleotide"),
  )
  let summary = @src.abstract_summary(aln)
  assert_true(summary.contains("Sequences"))
  assert_true(summary.contains("Length"))
}

///|
test "AbstractAlignment::single sequence identity" {
  let sequences = ["ACGT"]
  let identifiers = ["seq1"]
  let aln = @src.AbstractAlignment::new(
    sequences=sequences,
    identifiers=identifiers,
    alignment_type=@src.align_abstract_type_from_str("nucleotide"),
  )
  let matrix = @src.abstract_identity_matrix(aln)
  assert_eq(matrix.length(), 1)
  assert_true((matrix[0][0] - 1.0).abs() < 1.0e-10)
}
