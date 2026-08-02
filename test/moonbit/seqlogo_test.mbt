///|
/// seqLogo module tests.
/// Covers PWM construction, information content, logo computation,
/// rendering, consensus extraction, and sample data validation.

// ============================================================================
// 1. SeqLogoPwm creation (new, width, get, to_string, alphabet)
// ============================================================================

///|
test "seqlogo_pwm_new_basic" {
  let matrix = [
    [0.25, 0.25],
    [0.25, 0.25],
    [0.25, 0.25],
    [0.25, 0.25],
  ]
  let pwm = @src.SeqLogoPwm::new(matrix)
  assert_eq(pwm.width(), 2)
  assert_eq(pwm.alphabet_size(), 4)
}

///|
test "seqlogo_pwm_new_get_values" {
  let matrix = [
    [0.9, 0.1],
    [0.03, 0.1],
    [0.04, 0.7],
    [0.03, 0.1],
  ]
  let pwm = @src.SeqLogoPwm::new(matrix)
  assert_true((pwm.get(0, 0) - 0.9).abs() < 0.001)
  assert_true((pwm.get(0, 1) - 0.1).abs() < 0.001)
  assert_true((pwm.get(2, 1) - 0.7).abs() < 0.001)
  assert_true((pwm.get(3, 0) - 0.03).abs() < 0.001)
}

///|
test "seqlogo_pwm_new_get_out_of_bounds" {
  let matrix = [[0.25], [0.25], [0.25], [0.25]]
  let pwm = @src.SeqLogoPwm::new(matrix)
  assert_true((pwm.get(-1, 0)).abs() < 0.001)
  assert_true((pwm.get(4, 0)).abs() < 0.001)
  assert_true((pwm.get(0, -1)).abs() < 0.001)
  assert_true((pwm.get(0, 100)).abs() < 0.001)
}

///|
test "seqlogo_pwm_to_string" {
  let matrix = [[0.25], [0.25], [0.25], [0.25]]
  let pwm = @src.SeqLogoPwm::new(matrix)
  let s = pwm.to_string()
  assert_true(s.contains("SeqLogoPwm"))
  assert_true(s.contains("width=1"))
  assert_true(s.contains("A:"))
  assert_true(s.contains("C:"))
  assert_true(s.contains("G:"))
  assert_true(s.contains("T:"))
}

///|
test "seqlogo_pwm_alphabet" {
  let matrix = [[0.25], [0.25], [0.25], [0.25]]
  let pwm = @src.SeqLogoPwm::new(matrix)
  let alpha = pwm.alphabet()
  assert_eq(alpha.length(), 4)
  assert_eq(alpha[0], "A")
  assert_eq(alpha[1], "C")
  assert_eq(alpha[2], "G")
  assert_eq(alpha[3], "T")
}

// ============================================================================
// 2. SeqLogoPwm::from_sequences
// ============================================================================

///|
test "seqlogo_from_sequences_basic" {
  let seqs = ["ACGT", "ACGT"]
  let pwm = @src.SeqLogoPwm::from_sequences(seqs)
  assert_eq(pwm.width(), 4)
  assert_eq(pwm.alphabet_size(), 4)
  // Position 0: all A -> freq 1.0
  assert_true((pwm.get(0, 0) - 1.0).abs() < 0.001)
  // Position 1: all C -> freq 1.0
  assert_true((pwm.get(1, 1) - 1.0).abs() < 0.001)
  // Position 2: all G -> freq 1.0
  assert_true((pwm.get(2, 2) - 1.0).abs() < 0.001)
  // Position 3: all T -> freq 1.0
  assert_true((pwm.get(3, 3) - 1.0).abs() < 0.001)
  // Off-diagonal frequencies are 0
  assert_true((pwm.get(0, 1)).abs() < 0.001)
  assert_true((pwm.get(1, 0)).abs() < 0.001)
}

///|
test "seqlogo_from_sequences_mixed" {
  let seqs = ["AC", "AG", "AT", "AA"]
  let pwm = @src.SeqLogoPwm::from_sequences(seqs)
  assert_eq(pwm.width(), 2)
  // Position 0: all A
  assert_true((pwm.get(0, 0) - 1.0).abs() < 0.001)
  // Position 1: A=0.25, C=0.25, G=0.25, T=0.25
  assert_true((pwm.get(0, 1) - 0.25).abs() < 0.001)
  assert_true((pwm.get(1, 1) - 0.25).abs() < 0.001)
  assert_true((pwm.get(2, 1) - 0.25).abs() < 0.001)
  assert_true((pwm.get(3, 1) - 0.25).abs() < 0.001)
}

///|
test "seqlogo_from_sequences_lowercase" {
  let seqs = ["acgt"]
  let pwm = @src.SeqLogoPwm::from_sequences(seqs)
  assert_eq(pwm.width(), 4)
  assert_true((pwm.get(0, 0) - 1.0).abs() < 0.001)
  assert_true((pwm.get(1, 1) - 1.0).abs() < 0.001)
  assert_true((pwm.get(2, 2) - 1.0).abs() < 0.001)
  assert_true((pwm.get(3, 3) - 1.0).abs() < 0.001)
}

// ============================================================================
// 3. SeqLogoLetter creation and accessors
// ============================================================================

///|
test "seqlogo_letter_new_and_accessors" {
  let letter = @src.SeqLogoLetter::new('A', 1.5, "#00CC00")
  assert_eq(letter.letter(), 'A')
  assert_true((letter.height() - 1.5).abs() < 0.001)
  assert_eq(letter.color(), "#00CC00")
}

///|
test "seqlogo_letter_zero_height" {
  let letter = @src.SeqLogoLetter::new('T', 0.0, "#CC0000")
  assert_eq(letter.letter(), 'T')
  assert_true((letter.height()).abs() < 0.001)
  assert_eq(letter.color(), "#CC0000")
}

// ============================================================================
// 4. SeqLogoColumn creation, add_letter, accessors
// ============================================================================

///|
test "seqlogo_column_new_and_accessors" {
  let col = @src.SeqLogoColumn::new(3, 1.5)
  assert_eq(col.position(), 3)
  assert_true((col.info_content() - 1.5).abs() < 0.001)
  assert_eq(col.letters().length(), 0)
}

///|
test "seqlogo_column_add_letter" {
  let col = @src.SeqLogoColumn::new(0, 2.0)
  col.add_letter(@src.SeqLogoLetter::new('A', 2.0, "#00CC00"))
  col.add_letter(@src.SeqLogoLetter::new('C', 0.0, "#0000CC"))
  assert_eq(col.letters().length(), 2)
  assert_eq(col.letters()[0].letter(), 'A')
  assert_eq(col.letters()[1].letter(), 'C')
  assert_true((col.letters()[0].height() - 2.0).abs() < 0.001)
}

// ============================================================================
// 5. seqlogo_information_content (uniform background)
// ============================================================================

///|
test "seqlogo_ic_uniform_pwm_zero" {
  // All positions have uniform 0.25 frequency -> IC = 0
  let matrix = [
    [0.25, 0.25],
    [0.25, 0.25],
    [0.25, 0.25],
    [0.25, 0.25],
  ]
  let pwm = @src.SeqLogoPwm::new(matrix)
  let ic = @src.seqlogo_information_content(pwm)
  assert_eq(ic.length(), 2)
  assert_true((ic[0]).abs() < 0.001)
  assert_true((ic[1]).abs() < 0.001)
}

///|
test "seqlogo_ic_conserved_position" {
  // Fully conserved A -> IC = log2(4) = 2.0
  let matrix = [
    [1.0],
    [0.0],
    [0.0],
    [0.0],
  ]
  let pwm = @src.SeqLogoPwm::new(matrix)
  let ic = @src.seqlogo_information_content(pwm)
  assert_eq(ic.length(), 1)
  assert_true((ic[0] - 2.0).abs() < 0.001)
}

///|
test "seqlogo_ic_mixed_position" {
  // A=0.5, C=0.25, G=0.125, T=0.125
  // IC = 0.5*log2(2) + 0.25*log2(1) + 0.125*log2(0.5) + 0.125*log2(0.5)
  //    = 0.5 - 0.125 - 0.125 = 0.25
  let matrix = [
    [0.5],
    [0.25],
    [0.125],
    [0.125],
  ]
  let pwm = @src.SeqLogoPwm::new(matrix)
  let ic = @src.seqlogo_information_content(pwm)
  assert_true((ic[0] - 0.25).abs() < 0.001)
}

// ============================================================================
// 6. seqlogo_information_content_bg (non-uniform background)
// ============================================================================

///|
test "seqlogo_ic_bg_non_uniform" {
  // A=0.5, C=0.5, G=0.0, T=0.0
  // bg = [0.1, 0.4, 0.4, 0.1]
  // IC = 0.5*log2(0.5/0.1) + 0.5*log2(0.5/0.4)
  //    = 0.5*log2(5) + 0.5*log2(1.25)
  //    = 0.5*2.321928 + 0.5*0.321928 = 1.321928
  let matrix = [
    [0.5],
    [0.5],
    [0.0],
    [0.0],
  ]
  let pwm = @src.SeqLogoPwm::new(matrix)
  let ic = @src.seqlogo_information_content_bg(pwm, [0.1, 0.4, 0.4, 0.1])
  assert_true((ic[0] - 1.321928).abs() < 0.001)
}

///|
test "seqlogo_ic_bg_uniform_equals_default" {
  let matrix = [
    [0.9],
    [0.03],
    [0.04],
    [0.03],
  ]
  let pwm = @src.SeqLogoPwm::new(matrix)
  let ic_default = @src.seqlogo_information_content(pwm)
  let ic_bg = @src.seqlogo_information_content_bg(pwm, [
    0.25, 0.25, 0.25, 0.25,
  ])
  assert_eq(ic_default.length(), ic_bg.length())
  assert_true((ic_default[0] - ic_bg[0]).abs() < 0.001)
}

// ============================================================================
// 7. seqlogo_total_information_content
// ============================================================================

///|
test "seqlogo_total_ic_conserved" {
  // Two fully conserved positions -> total IC = 4.0
  let matrix = [
    [1.0, 1.0],
    [0.0, 0.0],
    [0.0, 0.0],
    [0.0, 0.0],
  ]
  let pwm = @src.SeqLogoPwm::new(matrix)
  let total = @src.seqlogo_total_information_content(pwm)
  assert_true((total - 4.0).abs() < 0.001)
}

///|
test "seqlogo_total_ic_uniform_zero" {
  let matrix = [
    [0.25, 0.25, 0.25],
    [0.25, 0.25, 0.25],
    [0.25, 0.25, 0.25],
    [0.25, 0.25, 0.25],
  ]
  let pwm = @src.SeqLogoPwm::new(matrix)
  let total = @src.seqlogo_total_information_content(pwm)
  assert_true(total.abs() < 0.001)
}

// ============================================================================
// 8. seqlogo_max_information_content
// ============================================================================

///|
test "seqlogo_max_ic_dna" {
  // DNA alphabet size 4 -> log2(4) = 2.0
  let max_ic = @src.seqlogo_max_information_content(4)
  assert_true((max_ic - 2.0).abs() < 0.001)
}

///|
test "seqlogo_max_ic_protein" {
  // Protein alphabet size 20 -> log2(20) ~= 4.321928
  let max_ic = @src.seqlogo_max_information_content(20)
  assert_true((max_ic - 4.321928).abs() < 0.001)
}

///|
test "seqlogo_max_ic_edge_cases" {
  // alphabet_size <= 1 -> 0.0
  assert_true((@src.seqlogo_max_information_content(0)).abs() < 0.001)
  assert_true((@src.seqlogo_max_information_content(1)).abs() < 0.001)
  // Binary alphabet -> log2(2) = 1.0
  assert_true((@src.seqlogo_max_information_content(2) - 1.0).abs() < 0.001)
}

// ============================================================================
// 9. seqlogo_compute_logo
// ============================================================================

///|
test "seqlogo_compute_logo_columns_count" {
  let matrix = [
    [0.25, 1.0],
    [0.25, 0.0],
    [0.25, 0.0],
    [0.25, 0.0],
  ]
  let pwm = @src.SeqLogoPwm::new(matrix)
  let logo = @src.seqlogo_compute_logo(pwm)
  assert_eq(logo.width(), 2)
  assert_eq(logo.columns().length(), 2)
}

///|
test "seqlogo_compute_logo_conserved_heights" {
  // Fully conserved A: IC=2.0, A height=2.0, others=0.0
  let matrix = [
    [1.0],
    [0.0],
    [0.0],
    [0.0],
  ]
  let pwm = @src.SeqLogoPwm::new(matrix)
  let logo = @src.seqlogo_compute_logo(pwm)
  let col = logo.columns()[0]
  assert_true((col.info_content() - 2.0).abs() < 0.001)
  assert_eq(col.letters().length(), 4)
  // Letters sorted by frequency descending; A is first
  assert_eq(col.letters()[0].letter(), 'A')
  assert_true((col.letters()[0].height() - 2.0).abs() < 0.001)
  // Other letters have height 0
  assert_true((col.letters()[1].height()).abs() < 0.001)
  assert_true((col.letters()[2].height()).abs() < 0.001)
  assert_true((col.letters()[3].height()).abs() < 0.001)
}

///|
test "seqlogo_compute_logo_uniform_zero_ic" {
  let matrix = [
    [0.25],
    [0.25],
    [0.25],
    [0.25],
  ]
  let pwm = @src.SeqLogoPwm::new(matrix)
  let logo = @src.seqlogo_compute_logo(pwm)
  let col = logo.columns()[0]
  assert_true((col.info_content()).abs() < 0.001)
  // All letter heights are 0 since IC=0
  for letter in col.letters() {
    assert_true((letter.height()).abs() < 0.001)
  }
}

///|
test "seqlogo_compute_logo_total_ic" {
  let matrix = [
    [1.0, 0.25],
    [0.0, 0.25],
    [0.0, 0.25],
    [0.0, 0.25],
  ]
  let pwm = @src.SeqLogoPwm::new(matrix)
  let logo = @src.seqlogo_compute_logo(pwm)
  // Position 0: IC=2.0, Position 1: IC=0.0 -> total = 2.0
  assert_true((logo.total_ic() - 2.0).abs() < 0.001)
}

///|
test "seqlogo_compute_logo_letter_colors" {
  let matrix = [
    [1.0],
    [0.0],
    [0.0],
    [0.0],
  ]
  let pwm = @src.SeqLogoPwm::new(matrix)
  let logo = @src.seqlogo_compute_logo(pwm)
  let col = logo.columns()[0]
  assert_eq(col.letters()[0].color(), "#00CC00") // A is green
  assert_eq(col.letters()[1].color(), "#0000CC") // C is blue
  assert_eq(col.letters()[2].color(), "#FFB300") // G is orange
  assert_eq(col.letters()[3].color(), "#CC0000") // T is red
}

// ============================================================================
// 10. seqlogo_base_color
// ============================================================================

///|
test "seqlogo_base_color_dna" {
  assert_eq(@src.seqlogo_base_color('A'), "#00CC00")
  assert_eq(@src.seqlogo_base_color('C'), "#0000CC")
  assert_eq(@src.seqlogo_base_color('G'), "#FFB300")
  assert_eq(@src.seqlogo_base_color('T'), "#CC0000")
}

///|
test "seqlogo_base_color_unknown" {
  // Unknown bases return black
  assert_eq(@src.seqlogo_base_color('N'), "#000000")
  assert_eq(@src.seqlogo_base_color('X'), "#000000")
  assert_eq(@src.seqlogo_base_color('-'), "#000000")
}

// ============================================================================
// 11. seqlogo_to_ascii
// ============================================================================

///|
test "seqlogo_to_ascii_non_empty" {
  let pwm = @src.seqlogo_sample_pwm()
  let logo = @src.seqlogo_compute_logo(pwm)
  let ascii = @src.seqlogo_to_ascii(logo, 10)
  assert_true(ascii.length() > 0)
}

///|
test "seqlogo_to_ascii_contains_letters" {
  let matrix = [
    [1.0, 0.0, 0.0, 0.0],
    [0.0, 1.0, 0.0, 0.0],
    [0.0, 0.0, 1.0, 0.0],
    [0.0, 0.0, 0.0, 1.0],
  ]
  let pwm = @src.SeqLogoPwm::new(matrix)
  let logo = @src.seqlogo_compute_logo(pwm)
  let ascii = @src.seqlogo_to_ascii(logo, 10)
  // The conserved positions should produce uppercase letters
  assert_true(ascii.contains("A"))
  assert_true(ascii.contains("C"))
  assert_true(ascii.contains("G"))
  assert_true(ascii.contains("T"))
}

///|
test "seqlogo_to_ascii_zero_height" {
  let ascii = @src.seqlogo_to_ascii(@src.seqlogo_compute_logo(@src.SeqLogoPwm::new([
    [0.25],
    [0.25],
    [0.25],
    [0.25],
  ])), 5)
  // IC=0 so no letters rendered, only spaces and newlines
  assert_true(ascii.length() > 0)
  assert_false(ascii.contains("A"))
}

// ============================================================================
// 12. seqlogo_to_text_table
// ============================================================================

///|
test "seqlogo_to_text_table_has_header" {
  let pwm = @src.seqlogo_sample_pwm()
  let logo = @src.seqlogo_compute_logo(pwm)
  let table = @src.seqlogo_to_text_table(logo)
  assert_true(table.length() > 0)
  assert_true(table.contains("position"))
  assert_true(table.contains("base"))
  assert_true(table.contains("height"))
  assert_true(table.contains("IC"))
}

///|
test "seqlogo_to_text_table_tab_separated" {
  let matrix = [
    [1.0],
    [0.0],
    [0.0],
    [0.0],
  ]
  let pwm = @src.SeqLogoPwm::new(matrix)
  let logo = @src.seqlogo_compute_logo(pwm)
  let table = @src.seqlogo_to_text_table(logo)
  // Header line is tab-separated
  assert_true(table.contains("position\tbase\theight\tIC"))
  // Data rows contain tab characters
  let lines = table.split("\n").to_array()
  assert_true(lines.length() >= 2)
}

// ============================================================================
// 13. seqlogo_consensus_sequence
// ============================================================================

///|
test "seqlogo_consensus_basic" {
  let matrix = [
    [1.0, 0.0, 0.0, 0.0],
    [0.0, 1.0, 0.0, 0.0],
    [0.0, 0.0, 1.0, 0.0],
    [0.0, 0.0, 0.0, 1.0],
  ]
  let pwm = @src.SeqLogoPwm::new(matrix)
  let consensus = @src.seqlogo_consensus_sequence(pwm)
  assert_eq(consensus, "ACGT")
}

///|
test "seqlogo_consensus_from_sequences" {
  let seqs = ["AAAA", "AACA", "AAGA", "AATA"]
  let pwm = @src.SeqLogoPwm::from_sequences(seqs)
  let consensus = @src.seqlogo_consensus_sequence(pwm)
  assert_eq(consensus.length(), 4)
  // Position 0: all A
  assert_eq(consensus[0:1].to_string(), "A")
  // Position 2: all A (most frequent)
  assert_eq(consensus[2:3].to_string(), "A")
}

///|
test "seqlogo_consensus_tie_first_base" {
  // When frequencies are tied, the lowest-index base wins
  let matrix = [
    [0.25],
    [0.25],
    [0.25],
    [0.25],
  ]
  let pwm = @src.SeqLogoPwm::new(matrix)
  let consensus = @src.seqlogo_consensus_sequence(pwm)
  assert_eq(consensus, "A")
}

// ============================================================================
// 14. seqlogo_logo_summary
// ============================================================================

///|
test "seqlogo_logo_summary_non_empty" {
  let pwm = @src.seqlogo_sample_pwm()
  let logo = @src.seqlogo_compute_logo(pwm)
  let summary = @src.seqlogo_logo_summary(logo)
  assert_true(summary.length() > 0)
  assert_true(summary.contains("Sequence Logo Summary"))
  assert_true(summary.contains("Total IC"))
  assert_true(summary.contains("Mean IC"))
}

///|
test "seqlogo_logo_summary_width" {
  let matrix = [
    [1.0, 0.25],
    [0.0, 0.25],
    [0.0, 0.25],
    [0.0, 0.25],
  ]
  let pwm = @src.SeqLogoPwm::new(matrix)
  let logo = @src.seqlogo_compute_logo(pwm)
  let summary = @src.seqlogo_logo_summary(logo)
  assert_true(summary.contains("Width: 2"))
}

// ============================================================================
// 15. Sample data validation
// ============================================================================

///|
test "seqlogo_sample_pwm_width" {
  let pwm = @src.seqlogo_sample_pwm()
  assert_eq(pwm.width(), 10)
  assert_eq(pwm.alphabet_size(), 4)
}

///|
test "seqlogo_sample_sequences_count" {
  let seqs = @src.seqlogo_sample_sequences()
  assert_eq(seqs.length(), 20)
}

///|
test "seqlogo_sample_sequences_length" {
  let seqs = @src.seqlogo_sample_sequences()
  for seq in seqs {
    assert_eq(seq.length(), 10)
  }
}

///|
test "seqlogo_sample_pwm_from_sample_sequences" {
  let seqs = @src.seqlogo_sample_sequences()
  let pwm = @src.SeqLogoPwm::from_sequences(seqs)
  assert_eq(pwm.width(), 10)
  // Position 0 of sample sequences is mostly A
  let consensus = @src.seqlogo_consensus_sequence(pwm)
  assert_eq(consensus.length(), 10)
  assert_eq(consensus[0:1].to_string(), "A")
}

// ============================================================================
// 16. Edge cases
// ============================================================================

///|
test "seqlogo_edge_single_position" {
  let matrix = [
    [1.0],
    [0.0],
    [0.0],
    [0.0],
  ]
  let pwm = @src.SeqLogoPwm::new(matrix)
  assert_eq(pwm.width(), 1)
  let ic = @src.seqlogo_information_content(pwm)
  assert_eq(ic.length(), 1)
  assert_true((ic[0] - 2.0).abs() < 0.001)
  let consensus = @src.seqlogo_consensus_sequence(pwm)
  assert_eq(consensus, "A")
}

///|
test "seqlogo_edge_all_same_base" {
  let seqs = ["AAAA", "AAAA"]
  let pwm = @src.SeqLogoPwm::from_sequences(seqs)
  assert_eq(pwm.width(), 4)
  let ic = @src.seqlogo_information_content(pwm)
  // Each position is fully conserved -> IC = 2.0
  for v in ic {
    assert_true((v - 2.0).abs() < 0.001)
  }
  let total = @src.seqlogo_total_information_content(pwm)
  assert_true((total - 8.0).abs() < 0.001)
  assert_eq(@src.seqlogo_consensus_sequence(pwm), "AAAA")
}

///|
test "seqlogo_edge_all_equal_frequencies" {
  let matrix = [
    [0.25, 0.25, 0.25],
    [0.25, 0.25, 0.25],
    [0.25, 0.25, 0.25],
    [0.25, 0.25, 0.25],
  ]
  let pwm = @src.SeqLogoPwm::new(matrix)
  let ic = @src.seqlogo_information_content(pwm)
  for v in ic {
    assert_true(v.abs() < 0.001)
  }
  let logo = @src.seqlogo_compute_logo(pwm)
  assert_true((logo.total_ic()).abs() < 0.001)
}

///|
test "seqlogo_edge_empty_pwm" {
  let pwm = @src.SeqLogoPwm::from_sequences([])
  assert_eq(pwm.width(), 0)
  let ic = @src.seqlogo_information_content(pwm)
  assert_eq(ic.length(), 0)
  let total = @src.seqlogo_total_information_content(pwm)
  assert_true(total.abs() < 0.001)
  let logo = @src.seqlogo_compute_logo(pwm)
  assert_eq(logo.width(), 0)
  assert_eq(logo.columns().length(), 0)
  let consensus = @src.seqlogo_consensus_sequence(pwm)
  assert_eq(consensus, "")
}

///|
test "seqlogo_edge_empty_logo_rendering" {
  let pwm = @src.SeqLogoPwm::from_sequences([])
  let logo = @src.seqlogo_compute_logo(pwm)
  // Rendering an empty logo should not crash
  let ascii = @src.seqlogo_to_ascii(logo, 5)
  assert_true(ascii.length() > 0) // contains newlines
  let table = @src.seqlogo_to_text_table(logo)
  assert_true(table.contains("position")) // header is still present
  let summary = @src.seqlogo_logo_summary(logo)
  assert_true(summary.contains("Sequence Logo Summary"))
}

///|
test "seqlogo_edge_to_ascii_zero_max_height" {
  let pwm = @src.seqlogo_sample_pwm()
  let logo = @src.seqlogo_compute_logo(pwm)
  let ascii = @src.seqlogo_to_ascii(logo, 0)
  assert_eq(ascii, "")
}

///|
test "seqlogo_compute_logo_bg_non_uniform" {
  let matrix = [
    [0.5],
    [0.5],
    [0.0],
    [0.0],
  ]
  let pwm = @src.SeqLogoPwm::new(matrix)
  let logo = @src.seqlogo_compute_logo_bg(pwm, [0.1, 0.4, 0.4, 0.1])
  let col = logo.columns()[0]
  // IC with non-uniform background ~= 1.321928
  assert_true((col.info_content() - 1.321928).abs() < 0.001)
  // Total IC matches column IC
  assert_true((logo.total_ic() - col.info_content()).abs() < 0.001)
}
