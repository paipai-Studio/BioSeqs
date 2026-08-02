///|
/// Tests for the AlignACE module (Bio.Motifs.AlignAce equivalent).
///
/// These tests exercise the public API exposed in `src/alignace.mbt`:
/// the AlignAceMotif / AlignAceSite / AlignAceRecord structs, the parser
/// (`alignace_parse`), the serializer (`alignace_write`), the PWM and
/// information content calculations, and utility functions.

// ============================================================================
// 1. AlignAceSite - construction
// ============================================================================

///|
test "alignace_site_new" {
  let s = @src.AlignAceSite::new(
    sequence_id="seq1",
    position=10,
    strand="+",
    sequence="TGACTCGATA",
  )
  assert_eq(s.sequence_id, "seq1")
  assert_eq(s.position, 10)
  assert_eq(s.strand, "+")
  assert_eq(s.sequence, "TGACTCGATA")
}

///|
test "alignace_site_new_negative_strand" {
  let s = @src.AlignAceSite::new(
    sequence_id="seq2",
    position=25,
    strand="-",
    sequence="ATATCGAGTC",
  )
  assert_eq(s.sequence_id, "seq2")
  assert_eq(s.position, 25)
  assert_eq(s.strand, "-")
}

// ============================================================================
// 2. AlignAceMotif - construction and consensus
// ============================================================================

///|
test "alignace_motif_new_basic" {
  let m = @src.AlignAceMotif::new(
    [
      [8, 0, 0, 0],
      [0, 8, 0, 0],
      [0, 0, 8, 0],
      [0, 0, 0, 8],
    ],
  )
  assert_eq(m.width, 4)
  // Total sites = sum of column 0 = 8 + 0 + 0 + 0 = 8
  assert_eq(m.num_sites, 8)
  // Consensus should be ACGT (each column has max count for that base)
  assert_eq(m.consensus, "ACGT")
}

///|
test "alignace_motif_new_consensus_tgac" {
  // Motif where T is most common at each position
  let m = @src.AlignAceMotif::new(
    [
      [1, 0, 0, 0],
      [0, 1, 0, 0],
      [0, 0, 1, 0],
      [9, 9, 9, 9],
    ],
  )
  assert_eq(m.width, 4)
  assert_eq(m.num_sites, 1 + 0 + 0 + 9)
  assert_eq(m.consensus, "TTTT")
}

///|
test "alignace_motif_consensus_ties_default_to_a" {
  // When counts are tied, A (first row) wins by default
  let m = @src.AlignAceMotif::new(
    [
      [5, 0, 5],
      [5, 0, 0],
      [0, 5, 0],
      [0, 5, 0],
    ],
  )
  // Position 0: A=5, C=5 -> A wins (first checked)
  // Position 1: G=5, T=5 -> G wins (checked before T)
  // Position 2: A=5 -> A
  assert_eq(m.consensus, "AGA")
}

///|
test "alignace_motif_motif_consensus" {
  let m = @src.AlignAceMotif::new(
    [
      [8, 0, 0, 0],
      [0, 8, 0, 0],
      [0, 0, 8, 0],
      [0, 0, 0, 8],
    ],
  )
  assert_eq(@src.alignace_motif_consensus(m), "ACGT")
}

// ============================================================================
// 3. alignace_to_pwm - position weight matrix
// ============================================================================

///|
test "alignace_to_pwm_uniform_distribution" {
  let m = @src.AlignAceMotif::new(
    [
      [2, 2],
      [2, 2],
      [2, 2],
      [2, 2],
    ],
  )
  let pwm = @src.alignace_to_pwm(m)
  assert_eq(pwm.length(), 4)
  assert_eq(pwm[0].length(), 2)
  // Each frequency should be 0.25 (uniform: 2 out of 8)
  for row in 0..<4 {
    for col in 0..<2 {
      assert_true(pwm[row][col] > 0.24 && pwm[row][col] < 0.26)
    }
  }
}

///|
test "alignace_to_pwm_pure_motif" {
  // Motif where only one base appears at each position
  let m = @src.AlignAceMotif::new(
    [
      [10, 0],
      [0, 0],
      [0, 10],
      [0, 0],
    ],
  )
  let pwm = @src.alignace_to_pwm(m)
  // Column 0: A=1.0, C=0.0, G=0.0, T=0.0
  assert_true(pwm[0][0] > 0.99)
  assert_true(pwm[1][0] < 0.01)
  assert_true(pwm[2][0] < 0.01)
  assert_true(pwm[3][0] < 0.01)
  // Column 1: A=0.0, C=0.0, G=1.0, T=0.0
  assert_true(pwm[0][1] < 0.01)
  assert_true(pwm[2][1] > 0.99)
}

///|
test "alignace_to_pwm_zero_column_uses_default" {
  // Column with all zeros should default to 0.25
  let m = @src.AlignAceMotif::new(
    [
      [0, 5],
      [0, 0],
      [0, 0],
      [0, 0],
    ],
  )
  let pwm = @src.alignace_to_pwm(m)
  // Column 0: all zeros, should be 0.25 each
  assert_true(pwm[0][0] > 0.24 && pwm[0][0] < 0.26)
  assert_true(pwm[1][0] > 0.24 && pwm[1][0] < 0.26)
  assert_true(pwm[2][0] > 0.24 && pwm[2][0] < 0.26)
  assert_true(pwm[3][0] > 0.24 && pwm[3][0] < 0.26)
}

// ============================================================================
// 4. alignace_information_content
// ============================================================================

///|
test "alignace_information_content_uniform_is_zero" {
  // Uniform distribution has 0 information content
  let m = @src.AlignAceMotif::new(
    [
      [2, 2],
      [2, 2],
      [2, 2],
      [2, 2],
    ],
  )
  let ic = @src.alignace_information_content(m)
  assert_eq(ic.length(), 2)
  // IC should be approximately 0 for uniform distribution
  assert_true(ic[0] < 0.01)
  assert_true(ic[0] > -0.01)
}

///|
test "alignace_information_content_pure_is_two_bits" {
  // Pure motif (only one base) has IC = 2 bits
  let m = @src.AlignAceMotif::new(
    [
      [10, 0],
      [0, 0],
      [0, 10],
      [0, 0],
    ],
  )
  let ic = @src.alignace_information_content(m)
  // Column 0: only A, IC should be ~2 bits
  assert_true(ic[0] > 1.9)
  // Column 1: only G, IC should be ~2 bits
  assert_true(ic[1] > 1.9)
}

///|
test "alignace_total_ic_uniform_is_zero" {
  let m = @src.AlignAceMotif::new(
    [
      [2, 2, 2],
      [2, 2, 2],
      [2, 2, 2],
      [2, 2, 2],
    ],
  )
  let total = @src.alignace_total_ic(m)
  assert_true(total < 0.01)
  assert_true(total > -0.01)
}

///|
test "alignace_total_ic_pure_is_width_times_two" {
  let m = @src.AlignAceMotif::new(
    [
      [10, 10, 10],
      [0, 0, 0],
      [0, 0, 0],
      [0, 0, 0],
    ],
  )
  let total = @src.alignace_total_ic(m)
  // 3 columns × 2 bits = 6 bits
  assert_true(total > 5.9)
}

// ============================================================================
// 5. alignace_parse - basic parsing
// ============================================================================

///|
test "alignace_parse_empty_returns_empty_record" {
  let record = @src.alignace_parse("")
  assert_eq(@src.alignace_num_motifs(record), 0)
}

///|
test "alignace_parse_extracts_version" {
  let text = "# AlignACE version 4.0\n# Program: AlignACE\n#\n"
  let record = @src.alignace_parse(text)
  assert_eq(record.version, "4.0")
}

///|
test "alignace_parse_extracts_command" {
  let text = "# AlignACE version 4.0\n# Program: AlignACE\n# Command: alignace -i=input.fa\n#\n"
  let record = @src.alignace_parse(text)
  assert_eq(record.command, "alignace -i=input.fa")
}

///|
test "alignace_parse_extracts_parameters" {
  let text = "# AlignACE version 4.0\n# Program: AlignACE\n# Command: alignace\n#\n# Parameters used:\n# numcols=10\n# expect=10\n#\n"
  let record = @src.alignace_parse(text)
  assert_eq(record.parameters["numcols"], "10")
  assert_eq(record.parameters["expect"], "10")
}

// ============================================================================
// 6. alignace_parse - motif matrix parsing
// ============================================================================

///|
test "alignace_parse_extracts_motif_matrix" {
  let text = "# AlignACE version 4.0\n# Program: AlignACE\n# Command: alignace\n#\n# Parameters used:\n# numcols=4\n#\ni\tA\tC\tG\tT\n\t0\t1\t2\t3\nA\t8\t0\t0\t0\nC\t0\t8\t0\t0\nG\t0\t0\t8\t0\nT\t0\t0\t0\t8\n#\n"
  let record = @src.alignace_parse(text)
  assert_eq(@src.alignace_num_motifs(record), 1)
  let m = @src.alignace_get_motif(record, 0)
  assert_true(m is Some(_))
  match m {
    Some(motif) => {
      assert_eq(motif.width, 4)
      assert_eq(motif.consensus, "ACGT")
      assert_eq(motif.num_sites, 8)
    }
    None => assert_true(false)
  }
}

///|
test "alignace_parse_extracts_multiple_motifs" {
  let text = "# AlignACE version 4.0\n#\ni\tA\tC\n\t0\t1\nA\t8\t0\nC\t0\t8\nG\t0\t0\nT\t0\t0\n#\ni\tG\tT\n\t0\t1\nA\t0\t0\nC\t0\t0\nG\t8\t0\nT\t0\t8\n#\n"
  let record = @src.alignace_parse(text)
  assert_eq(@src.alignace_num_motifs(record), 2)
}

// ============================================================================
// 7. alignace_parse - site parsing
// ============================================================================

///|
test "alignace_parse_extracts_sites" {
  let text = "# AlignACE version 4.0\n#\ni\tA\tC\n\t0\t1\nA\t8\t0\nC\t0\t8\nG\t0\t0\nT\t0\t0\n#\n# Sites:\n# seq1\t10\t+\tAC\n# seq2\t25\t+\tAC\n#\n"
  let record = @src.alignace_parse(text)
  assert_eq(@src.alignace_num_motifs(record), 1)
  let m = @src.alignace_get_motif(record, 0)
  match m {
    Some(motif) => assert_eq(motif.sites.length(), 2)
    None => assert_true(false)
  }
}

// ============================================================================
// 8. alignace_write - serialization
// ============================================================================

///|
test "alignace_write_includes_header" {
  let record = @src.alignace_sample()
  let text = @src.alignace_write(record)
  assert_true(text.contains("AlignACE version"))
  assert_true(text.contains("Program: AlignACE"))
  assert_true(text.contains("Command:"))
}

///|
test "alignace_write_includes_parameters" {
  let record = @src.alignace_sample()
  let text = @src.alignace_write(record)
  assert_true(text.contains("Parameters used:"))
  assert_true(text.contains("numcols="))
  assert_true(text.contains("expect="))
}

///|
test "alignace_write_includes_motifs" {
  let record = @src.alignace_sample()
  let text = @src.alignace_write(record)
  // The "i" line starts each motif matrix
  assert_true(text.contains("i\t"))
  // Base rows
  assert_true(text.contains("A\t"))
  assert_true(text.contains("C\t"))
  assert_true(text.contains("G\t"))
  assert_true(text.contains("T\t"))
}

// ============================================================================
// 9. alignace_write / alignace_parse round trip
// ============================================================================

///|
test "alignace_write_parse_round_trip_version" {
  let record = @src.alignace_sample()
  let text = @src.alignace_write(record)
  let parsed = @src.alignace_parse(text)
  assert_eq(parsed.version, record.version)
}

///|
test "alignace_write_parse_round_trip_command" {
  let record = @src.alignace_sample()
  let text = @src.alignace_write(record)
  let parsed = @src.alignace_parse(text)
  assert_eq(parsed.command, record.command)
}

///|
test "alignace_write_parse_round_trip_parameters" {
  let record = @src.alignace_sample()
  let text = @src.alignace_write(record)
  let parsed = @src.alignace_parse(text)
  assert_eq(parsed.parameters["numcols"], record.parameters["numcols"])
  assert_eq(parsed.parameters["expect"], record.parameters["expect"])
}

///|
test "alignace_write_parse_round_trip_motif_count" {
  let record = @src.alignace_sample()
  let text = @src.alignace_write(record)
  let parsed = @src.alignace_parse(text)
  assert_eq(@src.alignace_num_motifs(parsed), @src.alignace_num_motifs(record))
}

///|
test "alignace_write_parse_round_trip_motif_width" {
  let record = @src.alignace_sample()
  let text = @src.alignace_write(record)
  let parsed = @src.alignace_parse(text)
  let m1 = @src.alignace_get_motif(record, 0)
  let m2 = @src.alignace_get_motif(parsed, 0)
  match (m1, m2) {
    (Some(a), Some(b)) => assert_eq(b.width, a.width)
    _ => assert_true(false)
  }
}

// ============================================================================
// 10. alignace_sample - sample data
// ============================================================================

///|
test "alignace_sample_has_motifs" {
  let record = @src.alignace_sample()
  assert_true(@src.alignace_num_motifs(record) >= 1)
}

///|
test "alignace_sample_version" {
  let record = @src.alignace_sample()
  assert_eq(record.version, "4.0")
}

///|
test "alignace_sample_has_parameters" {
  let record = @src.alignace_sample()
  assert_true(record.parameters.size() >= 1)
}

///|
test "alignace_sample_first_motif_has_consensus" {
  let record = @src.alignace_sample()
  let m = @src.alignace_get_motif(record, 0)
  match m {
    Some(motif) => assert_true(motif.consensus.length() > 0)
    None => assert_true(false)
  }
}

///|
test "alignace_sample_first_motif_has_sites" {
  let record = @src.alignace_sample()
  let m = @src.alignace_get_motif(record, 0)
  match m {
    Some(motif) => assert_true(motif.sites.length() > 0)
    None => assert_true(false)
  }
}

// ============================================================================
// 11. alignace_num_motifs / alignace_get_motif
// ============================================================================

///|
test "alignace_num_motifs_empty" {
  let record = @src.AlignAceRecord::new()
  assert_eq(@src.alignace_num_motifs(record), 0)
}

///|
test "alignace_get_motif_out_of_range" {
  let record = @src.alignace_sample()
  assert_true(@src.alignace_get_motif(record, 999) is None)
}

///|
test "alignace_get_motif_negative_index" {
  let record = @src.alignace_sample()
  assert_true(@src.alignace_get_motif(record, -1) is None)
}

// ============================================================================
// 12. alignace_summary
// ============================================================================

///|
test "alignace_summary_contains_key_info" {
  let record = @src.alignace_sample()
  let summary = @src.alignace_summary(record)
  assert_true(summary.contains("AlignACE Record Summary:"))
  assert_true(summary.contains("Version:"))
  assert_true(summary.contains("Motifs:"))
  assert_true(summary.contains("Consensus:"))
}

///|
test "alignace_summary_empty_record" {
  let record = @src.AlignAceRecord::new()
  let summary = @src.alignace_summary(record)
  assert_true(summary.contains("AlignACE Record Summary:"))
  assert_true(summary.contains("Motifs: 0"))
}

// ============================================================================
// 13. Edge cases and integration
// ============================================================================

///|
test "alignace_motif_width_matches_matrix" {
  let m = @src.AlignAceMotif::new(
    [
      [1, 2, 3, 4, 5],
      [0, 0, 0, 0, 0],
      [0, 0, 0, 0, 0],
      [0, 0, 0, 0, 0],
    ],
  )
  assert_eq(m.width, 5)
  assert_eq(m.consensus, "AAAAA")
}

///|
test "alignace_motif_empty_matrix" {
  let m = @src.AlignAceMotif::new([])
  assert_eq(m.width, 0)
  assert_eq(m.consensus, "")
  assert_eq(m.num_sites, 0)
}

///|
test "alignace_information_content_length_matches_width" {
  let m = @src.AlignAceMotif::new(
    [
      [5, 5, 5, 5, 5],
      [0, 0, 0, 0, 0],
      [0, 0, 0, 0, 0],
      [5, 5, 5, 5, 5],
    ],
  )
  let ic = @src.alignace_information_content(m)
  assert_eq(ic.length(), 5)
}

///|
test "alignace_to_pwm_dimensions" {
  let m = @src.AlignAceMotif::new(
    [
      [1, 2, 3],
      [0, 0, 0],
      [0, 0, 0],
      [0, 0, 0],
    ],
  )
  let pwm = @src.alignace_to_pwm(m)
  // PWM should have 4 rows (one per base) and width columns
  assert_eq(pwm.length(), 4)
  for row in 0..<4 {
    assert_eq(pwm[row].length(), 3)
  }
}

///|
test "alignace_parse_motif_with_consensus_line" {
  // The "i" line should contain the consensus letters
  let text = "# AlignACE version 4.0\n#\ni\tT\tG\tA\tC\n\t0\t1\t2\t3\nA\t1\t0\t8\t0\nC\t0\t0\t0\t1\nG\t0\t8\t0\t0\nT\t7\t0\t0\t8\n#\n"
  let record = @src.alignace_parse(text)
  assert_eq(@src.alignace_num_motifs(record), 1)
  let m = @src.alignace_get_motif(record, 0)
  match m {
    Some(motif) => {
      assert_eq(motif.width, 4)
      // Column 0: T=7 is max -> T
      // Column 1: G=8 is max -> G
      // Column 2: A=8 is max -> A
      // Column 3: T=8 is max -> C
      assert_eq(motif.consensus, "TGAT")
    }
    None => assert_true(false)
  }
}
