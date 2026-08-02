///|
/// Tests for the FreqTable module (Bio.SubsMat.FreqTable equivalent).
///
/// These tests exercise the public API exposed in `src/freq_table.mbt`:
/// the FreqTableType enum, the FreqTable struct and its accessors,
/// the parsers (`freq_table_read_count`, `freq_table_read_freq`),
/// convenience constructors (`freq_table_from_counts`,
/// `freq_table_from_freqs`, `freq_table_from_sequence`),
/// statistical helpers (`freq_table_entropy`,
/// `freq_table_relative_entropy`, `freq_table_normalize`,
/// `freq_table_js_distance`), and the sample data helpers.

// ============================================================================
// 1. FreqTableType enum
// ============================================================================

///|
test "freq_table_count_type_tag" {
  let t = @src.freq_table_count_type()
  ignore(t)
}

///|
test "freq_table_freq_type_tag" {
  let t = @src.freq_table_freq_type()
  ignore(t)
}

// ============================================================================
// 2. FreqTable construction - COUNT
// ============================================================================

///|
test "freq_table_new_count_basic" {
  let dict : Map[String, Double] = Map([], capacity=4)
  dict["A"] = 50.0
  dict["C"] = 37.0
  dict["G"] = 23.0
  dict["T"] = 58.0
  let ft = @src.FreqTable::new(dict, @src.freq_table_count_type(), None)
  assert_eq(ft.size(), 4)
  assert_true(ft.has_counts())
  assert_eq(ft.total_count(), 168)
}

///|
test "freq_table_new_count_frequencies_sum_to_one" {
  let dict : Map[String, Double] = Map([], capacity=4)
  dict["A"] = 50.0
  dict["C"] = 37.0
  dict["G"] = 23.0
  dict["T"] = 58.0
  let ft = @src.FreqTable::new(dict, @src.freq_table_count_type(), None)
  let mut sum = 0.0
  for k in ft.data().keys() {
    sum = sum + ft.data()[k]
  }
  assert_true(sum > 0.999)
  assert_true(sum < 1.001)
}

///|
test "freq_table_new_count_specific_frequency" {
  let dict : Map[String, Double] = Map([], capacity=4)
  dict["A"] = 50.0
  dict["C"] = 37.0
  dict["G"] = 23.0
  dict["T"] = 58.0
  let ft = @src.FreqTable::new(dict, @src.freq_table_count_type(), None)
  // A: 50/168 = 0.2976...
  let freq_a = ft.frequency("A")
  assert_true(freq_a > 0.29)
  assert_true(freq_a < 0.30)
}

///|
test "freq_table_new_count_with_explicit_alphabet" {
  let dict : Map[String, Double] = Map([], capacity=4)
  dict["A"] = 50.0
  dict["C"] = 37.0
  dict["G"] = 23.0
  dict["T"] = 58.0
  let alph = ["A", "C", "G", "T"]
  let ft = @src.FreqTable::new(dict, @src.freq_table_count_type(), Some(alph))
  assert_eq(ft.size(), 4)
  let letters = ft.alphabet()
  assert_eq(letters[0], "A")
  assert_eq(letters[1], "C")
  assert_eq(letters[2], "G")
  assert_eq(letters[3], "T")
}

///|
test "freq_table_new_count_pads_missing_letters" {
  let dict : Map[String, Double] = Map([], capacity=2)
  dict["A"] = 10.0
  dict["T"] = 10.0
  let alph = ["A", "C", "G", "T"]
  let ft = @src.FreqTable::new(dict, @src.freq_table_count_type(), Some(alph))
  assert_eq(ft.size(), 4)
  // C and G should have 0 count / 0 frequency.
  assert_eq(ft.count_of("C"), 0)
  assert_true(ft.frequency("G") == 0.0)
}

// ============================================================================
// 3. FreqTable construction - FREQ
// ============================================================================

///|
test "freq_table_new_freq_basic" {
  let dict : Map[String, Double] = Map([], capacity=4)
  dict["A"] = 0.2976
  dict["C"] = 0.2202
  dict["G"] = 0.1369
  dict["T"] = 0.3452
  let ft = @src.FreqTable::new(dict, @src.freq_table_freq_type(), None)
  assert_eq(ft.size(), 4)
  assert_false(ft.has_counts())
  assert_eq(ft.total_count(), 0)
}

///|
test "freq_table_new_freq_preserves_values" {
  let dict : Map[String, Double] = Map([], capacity=4)
  dict["A"] = 0.2976
  dict["C"] = 0.2202
  dict["G"] = 0.1369
  dict["T"] = 0.3452
  let ft = @src.FreqTable::new(dict, @src.freq_table_freq_type(), None)
  let freq_a = ft.frequency("A")
  assert_true(freq_a == 0.2976)
}

// ============================================================================
// 4. Accessors
// ============================================================================

///|
test "freq_table_frequency_missing_letter_returns_zero" {
  let dict : Map[String, Double] = Map([], capacity=2)
  dict["A"] = 10.0
  dict["T"] = 10.0
  let ft = @src.FreqTable::new(dict, @src.freq_table_count_type(), None)
  assert_true(ft.frequency("N") == 0.0)
}

///|
test "freq_table_count_of_missing_letter_returns_zero" {
  let dict : Map[String, Double] = Map([], capacity=2)
  dict["A"] = 10.0
  dict["T"] = 10.0
  let ft = @src.FreqTable::new(dict, @src.freq_table_count_type(), None)
  assert_eq(ft.count_of("N"), 0)
}

///|
test "freq_table_contains" {
  let dict : Map[String, Double] = Map([], capacity=2)
  dict["A"] = 10.0
  dict["T"] = 10.0
  let ft = @src.FreqTable::new(dict, @src.freq_table_count_type(), None)
  assert_true(ft.contains("A"))
  assert_false(ft.contains("N"))
}

///|
test "freq_table_letters_sorted" {
  let dict : Map[String, Double] = Map([], capacity=4)
  dict["T"] = 10.0
  dict["A"] = 10.0
  dict["G"] = 10.0
  dict["C"] = 10.0
  let ft = @src.FreqTable::new(dict, @src.freq_table_count_type(), None)
  let letters = ft.letters()
  assert_eq(letters[0], "A")
  assert_eq(letters[1], "C")
  assert_eq(letters[2], "G")
  assert_eq(letters[3], "T")
}

// ============================================================================
// 5. to_string
// ============================================================================

///|
test "freq_table_to_string_count_includes_letters" {
  let dict : Map[String, Double] = Map([], capacity=2)
  dict["A"] = 10.0
  dict["T"] = 10.0
  let ft = @src.FreqTable::new(dict, @src.freq_table_count_type(), None)
  let s = ft.to_string()
  assert_true(s.contains("A"))
  assert_true(s.contains("T"))
}

///|
test "freq_table_to_string_freq_no_counts" {
  let dict : Map[String, Double] = Map([], capacity=2)
  dict["A"] = 0.5
  dict["T"] = 0.5
  let ft = @src.FreqTable::new(dict, @src.freq_table_freq_type(), None)
  let s = ft.to_string()
  assert_true(s.contains("A"))
  assert_true(s.contains("T"))
}

// ============================================================================
// 6. Parsers - read_count / read_freq
// ============================================================================

///|
test "freq_table_read_count_basic" {
  let text = "A 50\nC 37\nG 23\nT 58"
  let ft = @src.freq_table_read_count(text)
  assert_eq(ft.size(), 4)
  assert_true(ft.has_counts())
  assert_eq(ft.count_of("A"), 50)
  assert_eq(ft.count_of("C"), 37)
  assert_eq(ft.count_of("G"), 23)
  assert_eq(ft.count_of("T"), 58)
}

///|
test "freq_table_read_count_computes_frequencies" {
  let text = "A 50\nC 37\nG 23\nT 58"
  let ft = @src.freq_table_read_count(text)
  // 50 / 168 = 0.2976...
  let f = ft.frequency("A")
  assert_true(f > 0.29)
  assert_true(f < 0.30)
}

///|
test "freq_table_read_count_handles_spaces_and_tabs" {
  let text = "A\t50 C 37\nG\t23 T 58"
  let ft = @src.freq_table_read_count(text)
  assert_eq(ft.size(), 4)
  assert_eq(ft.count_of("A"), 50)
}

///|
test "freq_table_read_count_empty_text" {
  let ft = @src.freq_table_read_count("")
  assert_eq(ft.size(), 0)
}

///|
test "freq_table_read_freq_basic" {
  let text = "A 0.2976\nC 0.2202\nG 0.1369\nT 0.3452"
  let ft = @src.freq_table_read_freq(text)
  assert_eq(ft.size(), 4)
  assert_false(ft.has_counts())
  // Use approximate comparison for floating-point values.
  let fa = ft.frequency("A")
  assert_true(fa > 0.29 && fa < 0.30)
  let ft2 = ft.frequency("T")
  assert_true(ft2 > 0.34 && ft2 < 0.35)
}

///|
test "freq_table_read_freq_empty_text" {
  let ft = @src.freq_table_read_freq("")
  assert_eq(ft.size(), 0)
}

// ============================================================================
// 7. Convenience constructors
// ============================================================================

///|
test "freq_table_from_counts_basic" {
  let counts : Map[String, Int] = Map([], capacity=4)
  counts["A"] = 50
  counts["C"] = 37
  counts["G"] = 23
  counts["T"] = 58
  let ft = @src.freq_table_from_counts(counts)
  assert_eq(ft.size(), 4)
  assert_true(ft.has_counts())
  assert_eq(ft.count_of("A"), 50)
}

///|
test "freq_table_from_freqs_basic" {
  let freqs : Map[String, Double] = Map([], capacity=4)
  freqs["A"] = 0.25
  freqs["C"] = 0.25
  freqs["G"] = 0.25
  freqs["T"] = 0.25
  let ft = @src.freq_table_from_freqs(freqs)
  assert_eq(ft.size(), 4)
  assert_false(ft.has_counts())
  assert_true(ft.frequency("A") == 0.25)
}

///|
test "freq_table_from_sequence_basic" {
  let ft = @src.freq_table_from_sequence("AATTGGCC", None)
  assert_eq(ft.count_of("A"), 2)
  assert_eq(ft.count_of("T"), 2)
  assert_eq(ft.count_of("G"), 2)
  assert_eq(ft.count_of("C"), 2)
}

///|
test "freq_table_from_sequence_with_alphabet_restricts" {
  let ft = @src.freq_table_from_sequence("AATTGGCCNN", Some(["A", "T", "G", "C"]))
  assert_eq(ft.size(), 4)
  assert_eq(ft.count_of("A"), 2)
  // N is not in alphabet, so it should not be counted.
  assert_false(ft.contains("N"))
}

///|
test "freq_table_from_sequence_empty" {
  let ft = @src.freq_table_from_sequence("", None)
  assert_eq(ft.size(), 0)
}

// ============================================================================
// 8. Statistical helpers - entropy
// ============================================================================

///|
test "freq_table_entropy_uniform_distribution" {
  // Uniform distribution over 4 letters: H = log2(4) = 2.0 bits.
  let dict : Map[String, Double] = Map([], capacity=4)
  dict["A"] = 0.25
  dict["C"] = 0.25
  dict["G"] = 0.25
  dict["T"] = 0.25
  let ft = @src.FreqTable::new(dict, @src.freq_table_freq_type(), None)
  let h = @src.freq_table_entropy(ft)
  assert_true(h > 1.99)
  assert_true(h < 2.01)
}

///|
test "freq_table_entropy_skewed_distribution_lower" {
  // Skewed distribution has lower entropy.
  let dict : Map[String, Double] = Map([], capacity=4)
  dict["A"] = 0.9
  dict["C"] = 0.05
  dict["G"] = 0.03
  dict["T"] = 0.02
  let ft = @src.FreqTable::new(dict, @src.freq_table_freq_type(), None)
  let h = @src.freq_table_entropy(ft)
  assert_true(h < 2.0)
}

///|
test "freq_table_entropy_single_letter_zero" {
  let dict : Map[String, Double] = Map([], capacity=1)
  dict["A"] = 1.0
  let ft = @src.FreqTable::new(dict, @src.freq_table_freq_type(), None)
  let h = @src.freq_table_entropy(ft)
  assert_true(h == 0.0)
}

// ============================================================================
// 9. Statistical helpers - relative entropy (KL divergence)
// ============================================================================

///|
test "freq_table_relative_entropy_identical_tables_zero" {
  let dict : Map[String, Double] = Map([], capacity=4)
  dict["A"] = 0.25
  dict["C"] = 0.25
  dict["G"] = 0.25
  dict["T"] = 0.25
  let ft1 = @src.FreqTable::new(dict, @src.freq_table_freq_type(), None)
  let ft2 = @src.FreqTable::new(dict, @src.freq_table_freq_type(), None)
  let d = @src.freq_table_relative_entropy(ft1, ft2)
  assert_true(d < 0.001)
}

///|
test "freq_table_relative_entropy_different_tables_positive" {
  let dict1 : Map[String, Double] = Map([], capacity=2)
  dict1["A"] = 0.9
  dict1["T"] = 0.1
  let dict2 : Map[String, Double] = Map([], capacity=2)
  dict2["A"] = 0.1
  dict2["T"] = 0.9
  let ft1 = @src.FreqTable::new(dict1, @src.freq_table_freq_type(), None)
  let ft2 = @src.FreqTable::new(dict2, @src.freq_table_freq_type(), None)
  let d = @src.freq_table_relative_entropy(ft1, ft2)
  assert_true(d > 0.0)
}

// ============================================================================
// 10. Normalize
// ============================================================================

///|
test "freq_table_normalize_sums_to_one" {
  let dict : Map[String, Double] = Map([], capacity=3)
  dict["A"] = 2.0
  dict["C"] = 3.0
  dict["G"] = 5.0
  let ft = @src.FreqTable::new(dict, @src.freq_table_freq_type(), None)
  let norm = @src.freq_table_normalize(ft)
  let mut sum = 0.0
  for k in norm.data().keys() {
    sum = sum + norm.data()[k]
  }
  assert_true(sum > 0.999)
  assert_true(sum < 1.001)
}

///|
test "freq_table_normalize_preserves_proportions" {
  let dict : Map[String, Double] = Map([], capacity=3)
  dict["A"] = 2.0
  dict["C"] = 3.0
  dict["G"] = 5.0
  let ft = @src.FreqTable::new(dict, @src.freq_table_freq_type(), None)
  let norm = @src.freq_table_normalize(ft)
  // A should be 0.2, C should be 0.3, G should be 0.5
  assert_true(norm.frequency("A") > 0.19)
  assert_true(norm.frequency("A") < 0.21)
  assert_true(norm.frequency("G") > 0.49)
  assert_true(norm.frequency("G") < 0.51)
}

// ============================================================================
// 11. Jensen-Shannon distance
// ============================================================================

///|
test "freq_table_js_distance_identical_tables_zero" {
  let dict : Map[String, Double] = Map([], capacity=4)
  dict["A"] = 0.25
  dict["C"] = 0.25
  dict["G"] = 0.25
  dict["T"] = 0.25
  let ft1 = @src.FreqTable::new(dict, @src.freq_table_freq_type(), None)
  let ft2 = @src.FreqTable::new(dict, @src.freq_table_freq_type(), None)
  let js = @src.freq_table_js_distance(ft1, ft2)
  assert_true(js < 0.001)
}

///|
test "freq_table_js_distance_different_tables_positive" {
  let dict1 : Map[String, Double] = Map([], capacity=2)
  dict1["A"] = 0.9
  dict1["T"] = 0.1
  let dict2 : Map[String, Double] = Map([], capacity=2)
  dict2["A"] = 0.1
  dict2["T"] = 0.9
  let ft1 = @src.FreqTable::new(dict1, @src.freq_table_freq_type(), None)
  let ft2 = @src.FreqTable::new(dict2, @src.freq_table_freq_type(), None)
  let js = @src.freq_table_js_distance(ft1, ft2)
  assert_true(js > 0.0)
}

///|
test "freq_table_js_distance_symmetric" {
  let dict1 : Map[String, Double] = Map([], capacity=2)
  dict1["A"] = 0.9
  dict1["T"] = 0.1
  let dict2 : Map[String, Double] = Map([], capacity=2)
  dict2["A"] = 0.3
  dict2["T"] = 0.7
  let ft1 = @src.FreqTable::new(dict1, @src.freq_table_freq_type(), None)
  let ft2 = @src.FreqTable::new(dict2, @src.freq_table_freq_type(), None)
  let js12 = @src.freq_table_js_distance(ft1, ft2)
  let js21 = @src.freq_table_js_distance(ft2, ft1)
  assert_true((js12 - js21).abs() < 0.001)
}

// ============================================================================
// 12. Sample data
// ============================================================================

///|
test "freq_table_sample_count_text_has_four_letters" {
  let text = @src.freq_table_sample_count_text()
  let ft = @src.freq_table_read_count(text)
  assert_eq(ft.size(), 4)
  assert_eq(ft.total_count(), 168)
}

///|
test "freq_table_sample_freq_text_has_four_letters" {
  let text = @src.freq_table_sample_freq_text()
  let ft = @src.freq_table_read_freq(text)
  assert_eq(ft.size(), 4)
}

///|
test "freq_table_sample_dna_basic" {
  let ft = @src.freq_table_sample_dna()
  assert_eq(ft.size(), 4)
  assert_true(ft.has_counts())
  assert_eq(ft.count_of("A"), 50)
}

///|
test "freq_table_sample_dna_freq_basic" {
  let ft = @src.freq_table_sample_dna_freq()
  assert_eq(ft.size(), 4)
  assert_false(ft.has_counts())
}

///|
test "freq_table_sample_protein_basic" {
  let ft = @src.freq_table_sample_protein()
  assert_eq(ft.size(), 20)
  assert_true(ft.has_counts())
}

///|
test "freq_table_sample_protein_count_text_has_20_letters" {
  let text = @src.freq_table_sample_protein_count_text()
  let ft = @src.freq_table_read_count(text)
  assert_eq(ft.size(), 20)
  // Total = sum of all counts.
  assert_true(ft.total_count() > 0)
}

///|
test "freq_table_sample_dna_frequencies_sum_to_one" {
  let ft = @src.freq_table_sample_dna()
  let mut sum = 0.0
  for k in ft.data().keys() {
    sum = sum + ft.data()[k]
  }
  assert_true(sum > 0.999)
  assert_true(sum < 1.001)
}
