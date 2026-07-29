///|
test "letter_frequency_dna_basic" {
  let seq = "AGCT"
  let freq = @src.letter_frequency_dna(seq)
  assert_eq(freq.a, 1)
  assert_eq(freq.c, 1)
  assert_eq(freq.g, 1)
  assert_eq(freq.t, 1)
  assert_eq(freq.other, 0)
}

///|
test "letter_frequency_dna_with_other" {
  let seq = "AGCTN"
  let freq = @src.letter_frequency_dna(seq)
  assert_eq(freq.a, 1)
  assert_eq(freq.c, 1)
  assert_eq(freq.g, 1)
  assert_eq(freq.t, 1)
  assert_eq(freq.n, 1)
}

///|
test "letter_frequency_aa_basic" {
  let seq = "ACDEFGHIKLMNPQRSTVWY"
  let freq = @src.letter_frequency_aa(seq)
  assert_eq(freq.a, 1)
  assert_eq(freq.c, 1)
  assert_eq(freq.e, 1)
  assert_eq(freq.other, 0)
}

///|
test "oligonucleotide_frequency_mono" {
  let seq = "AATCG"
  let freq = @src.oligonucleotide_frequency(seq, 1, false)
  assert_eq(freq["A"], 2.0)
  assert_eq(freq["T"], 1.0)
  assert_eq(freq["C"], 1.0)
  assert_eq(freq["G"], 1.0)
}

///|
test "oligonucleotide_frequency_di" {
  let seq = "AATCG"
  let freq = @src.oligonucleotide_frequency(seq, 2, false)
  assert_eq(freq["AA"], 1.0)
  assert_eq(freq["AT"], 1.0)
  assert_eq(freq["TC"], 1.0)
  assert_eq(freq["CG"], 1.0)
}

///|
test "oligonucleotide_frequency_tri" {
  let seq = "AATCG"
  let freq = @src.oligonucleotide_frequency(seq, 3, false)
  assert_eq(freq["AAT"], 1.0)
  assert_eq(freq["ATC"], 1.0)
  assert_eq(freq["TCG"], 1.0)
}

///|
test "rscu_basic" {
  let seq = "ATGTGCTGA"
  let rscu = @src.rscu(seq)
  assert_eq(rscu.length() > 0, true)
}

///|
test "gc_content_by_position_basic" {
  let seq = "AGCT"
  let gc = @src.gc_content_by_position(seq, 1)
  assert_eq(gc.length(), 4)
  assert_eq(gc[0], 0.0)
  assert_eq(gc[1], 1.0)
  assert_eq(gc[2], 1.0)
  assert_eq(gc[3], 0.0)
}

///|
test "shannon_entropy_basic" {
  let entropy = @src.shannon_entropy("AAAA")
  assert_eq(entropy, 0.0)
}

///|
test "shannon_entropy_max" {
  let entropy = @src.shannon_entropy("AGCT")
  assert_eq(entropy >= 1.9 && entropy <= 2.1, true)
}

///|
test "dust_complexity_basic" {
  let seq = "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"
  let score = @src.dust_complexity(seq)
  assert_eq(score > 0, true)
}

///|
test "iupac_match_exact" {
  let matches = @src.iupac_match("A", "A")
  assert_eq(matches.length() > 0, true)
}

///|
test "iupac_match_degenerate" {
  let matches = @src.iupac_match("N", "A")
  assert_eq(matches.length() > 0, true)
}

///|
test "tm_nn_basic" {
  let tm = @src.tm_nn("AGCT")
  assert_eq(tm > 0, true)
}

///|
test "nedit_at_basic" {
  let edits = @src.nedit_at("AGCT", "AGCT", true)
  assert_eq(edits, 0)
}

///|
test "nedit_at_with_mismatch" {
  let edits = @src.nedit_at("AGCT", "AGTT", true)
  assert_eq(edits, 1)
}

///|
test "percent_identity_basic" {
  let pid = @src.percent_identity("AGCT", "AGCT")
  assert_eq(pid, 100.0)
}

///|
test "percent_identity_partial" {
  let pid = @src.percent_identity("AGCT", "AGTT")
  assert_eq(pid, 75.0)
}

///|
test "iupac_reverse_complement_basic" {
  let rc = @src.iupac_reverse_complement("AGCT")
  assert_eq(rc, "AGCT")
}

///|
test "iupac_reverse_complement_degenerate" {
  let rc = @src.iupac_reverse_complement("N")
  assert_eq(rc, "N")
}

///|
test "iupac_reverse_complement_mixed" {
  let rc = @src.iupac_reverse_complement("ACGTMRWSYKVHDBN")
  assert_eq(rc, "NVHDBMRSWYKACGT")
}

///|
test "biostrings_translate_basic" {
  let protein = @src.biostrings_translate("ATGGCTTAA")
  assert_eq(protein, "MA*")
}

///|
test "biostrings_translate_frame" {
  let protein = @src.biostrings_translate("ATGGCTTAA", frame=1)
  assert_eq(protein, "WL")
}

///|
test "biostrings_translate_lowercase" {
  let protein = @src.biostrings_translate("atggcttaa")
  assert_eq(protein, "MA*")
}

///|
test "biostrings_translate_empty" {
  let protein = @src.biostrings_translate("")
  assert_eq(protein, "")
}

///|
test "biostrings_reverse_complement_basic" {
  let rc = @src.biostrings_reverse_complement("ATGC")
  assert_eq(rc, "GCAT")
}

///|
test "biostrings_reverse_complement_iupac" {
  let rc = @src.biostrings_reverse_complement("AR")
  assert_eq(rc, "YT")
}

///|
test "biostrings_reverse_complement_empty" {
  let rc = @src.biostrings_reverse_complement("")
  assert_eq(rc, "")
}

///|
test "letter_frequency_matrix_basic" {
  let matrix = @src.letter_frequency_matrix("ACGT")
  let a_counts = matrix.get_or_default("A", Array::make(4, 0))
  assert_eq(a_counts[0], 1)
}

///|
test "biostrings_consensus_sequence_basic" {
  let alignment = ["ACGT", "ACGT", "TCGA"]
  let consensus = @src.biostrings_consensus_sequence(alignment)
  assert_eq(consensus, "ACGT")
}

///|
test "biostrings_consensus_sequence_empty" {
  let alignment : Array[String] = []
  let consensus = @src.biostrings_consensus_sequence(alignment)
  assert_eq(consensus, "")
}

///|
test "expected_matches_basic" {
  let n = @src.expected_matches("ACGT", 12, 0.5)
  assert_true(n > 0.0)
}

///|
test "sequence_complexity_basic" {
  let entropy = @src.sequence_complexity("ACGTACGTACGT", 2)
  assert_true(entropy > 0.0)
}
