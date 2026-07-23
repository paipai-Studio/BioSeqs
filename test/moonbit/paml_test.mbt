///|
/// Tests for PAML module.

test "PAMLAlignment creation" {
  let seqs = [("Human", "ATGCCG"), ("Mouse", "ATGCCG")]
  let alignment = @src.PAMLAlignment::new(seqs)
  assert_eq(alignment.sequences.length(), 2)
}

test "PAMLResult creation" {
  let result = @src.PAMLResult::new(-100.0)
  assert_eq(result.ln_likelihood, -100.0)
}

test "DNDSResult creation" {
  let dnds = @src.DNDSResult::new(0.1, 0.5, 0.2)
  assert_eq(dnds.dN, 0.1)
  assert_eq(dnds.dS, 0.5)
  assert_eq(dnds.omega, 0.2)
}

test "paml_calculate_dnds identical sequences" {
  let dnds = @src.paml_calculate_dnds("ATGCCG", "ATGCCG", "Nei-Gojobori")
  assert_true(dnds.omega >= 0.0)
}

test "paml_calculate_dnds different sequences" {
  let dnds = @src.paml_calculate_dnds("ATGCCG", "ATTTTT", "Nei-Gojobori")
  assert_true(dnds.dN >= 0.0)
  assert_true(dnds.dS >= 0.0)
}

test "paml_estimate_parameters" {
  let seqs = [("Seq1", "ATGC"), ("Seq2", "ATGC")]
  let alignment = @src.PAMLAlignment::new(seqs)
  let result = @src.paml_estimate_parameters(alignment)
  assert_true(result.parameters.contains("pi_A"))
}

test "paml_calculate_substitution_matrix" {
  let pi = Map([("pi_A", 0.25), ("pi_T", 0.25), ("pi_C", 0.25), ("pi_G", 0.25)], capacity=4)
  let matrix = @src.paml_calculate_substitution_matrix(2.0, pi)
  assert_eq(matrix.length(), 4)
}

test "paml_create_example_alignment" {
  let alignment = @src.paml_create_example_alignment()
  assert_eq(alignment.sequences.length(), 4)
}

test "paml_run_likelihood" {
  let alignment = @src.paml_create_example_alignment()
  let result = @src.paml_run_likelihood(alignment)
  assert_true(result.dnds_ratios.length() > 0)
}

test "get_standard_codon_table" {
  let table = @src.get_standard_codon_table()
  assert_true(table.contains("ATG"))
}