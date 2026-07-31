///|
/// Tests for MutationalPatterns module.

///|
test "mutation_category_to_string" {
  assert_eq(@src.mutation_c_a().to_string(), "C>A")
  assert_eq(@src.mutation_c_g().to_string(), "C>G")
  assert_eq(@src.mutation_c_t().to_string(), "C>T")
  assert_eq(@src.mutation_t_a().to_string(), "T>A")
  assert_eq(@src.mutation_t_c().to_string(), "T>C")
  assert_eq(@src.mutation_t_g().to_string(), "T>G")
}

///|
test "all_mutation_types_count" {
  let types = @src.all_mutation_types()
  assert_eq(types.length(), 6)
}

///|
test "trinucleotide_contexts_count" {
  let contexts = @src.trinucleotide_contexts()
  assert_eq(contexts.length(), 16)
}

///|
test "channel_labels_count" {
  let labels = @src.get_channel_labels()
  assert_eq(labels.length(), 96)
}

///|
test "get_mutation_type_ct" {
  let mtype = @src.get_mutation_type("C", "T")
  assert_eq(mtype.to_string(), "C>T")
}

///|
test "get_mutation_type_ca" {
  let mtype = @src.get_mutation_type("C", "A")
  assert_eq(mtype.to_string(), "C>A")
}

///|
test "get_mutation_type_purine_ref" {
  // G>A should be converted to C>T (pyrimidine-based)
  let mtype = @src.get_mutation_type("G", "A")
  assert_eq(mtype.to_string(), "C>T")
}

///|
test "normalize_context_pyrimidine" {
  // ACT (middle is C) should stay as ACT
  let norm = @src.normalize_context("ACT")
  assert_eq(norm, "ACT")
}

///|
test "normalize_context_purine" {
  // AGT (middle is G, a purine) should be reverse complemented to ACT
  let norm = @src.normalize_context("AGT")
  assert_eq(norm, "ACT")
}

///|
test "get_channel_index_range" {
  let idx1 = @src.get_channel_index("C", "T", "ACT")
  let idx2 = @src.get_channel_index("C", "A", "CCA")
  assert_true(idx1 >= 0 && idx1 < 96)
  assert_true(idx2 >= 0 && idx2 < 96)
}

///|
test "build_mutation_matrix" {
  let mutations = @src.mutational_patterns_sample_data()
  let sample_ids = ["S1", "S2"]
  let matrix = @src.build_mutation_matrix(mutations, sample_ids)
  assert_eq(matrix.sample_ids.length(), 2)
  assert_eq(matrix.channel_labels.length(), 96)
  // S1 should have 5 mutations
  assert_eq(matrix.total_mutations[0], 5)
  // S2 should have 4 mutations
  assert_eq(matrix.total_mutations[1], 4)
}

///|
test "mutation_matrix_get_count" {
  let mutations = @src.mutational_patterns_sample_data()
  let sample_ids = ["S1", "S2"]
  let matrix = @src.build_mutation_matrix(mutations, sample_ids)
  // Get total counts for sample 0
  let mut total = 0
  let mut i = 0
  while i < 96 {
    total = total + matrix.get_count(i, 0)
    i = i + 1
  }
  assert_eq(total, 5)
}

///|
test "cosine_similarity_identical" {
  let vec = [1.0, 2.0, 3.0]
  let sim = @src.cosine_similarity(vec, vec)
  assert_true(sim > 0.99)
}

///|
test "cosine_similarity_orthogonal" {
  let v1 = [1.0, 0.0]
  let v2 = [0.0, 1.0]
  let sim = @src.cosine_similarity(v1, v2)
  assert_true(sim < 0.01)
}

///|
test "relative_frequencies_sum" {
  let mutations = @src.mutational_patterns_sample_data()
  let sample_ids = ["S1", "S2"]
  let matrix = @src.build_mutation_matrix(mutations, sample_ids)
  let freqs = matrix.to_relative_frequencies()
  // Sum of frequencies for each sample should be ~1.0
  let mut sum0 = 0.0
  let mut sum1 = 0.0
  let mut i = 0
  while i < 96 {
    sum0 = sum0 + freqs[i][0]
    sum1 = sum1 + freqs[i][1]
    i = i + 1
  }
  assert_true(sum0 > 0.99)
  assert_true(sum0 < 1.01)
  assert_true(sum1 > 0.99)
  assert_true(sum1 < 1.01)
}

///|
test "sample_cosine_similarity" {
  let mutations = @src.mutational_patterns_sample_data()
  let sample_ids = ["S1", "S2"]
  let matrix = @src.build_mutation_matrix(mutations, sample_ids)
  let sim = matrix.sample_cosine_similarity(0, 0)
  assert_true(sim > 0.99) // Identical samples
}

///|
test "known_cosmic_signatures" {
  let sigs = @src.known_cosmic_signatures()
  assert_eq(sigs.length(), 2)
  assert_eq(sigs[0].id, "SBS1")
  assert_eq(sigs[1].id, "SBS5")
}

///|
test "fit_signatures" {
  let sigs = @src.known_cosmic_signatures()
  let profile : Array[Double] = Array::new()
  let mut i = 0
  while i < 96 {
    profile.push(1.0 / 96.0)
    i = i + 1
  }
  let contributions = @src.fit_signatures(profile, sigs)
  assert_eq(contributions.length(), 2)
  // Sum should be ~1.0
  let mut total = 0.0
  let mut j = 0
  while j < contributions.length() {
    total = total + contributions[j]
    j = j + 1
  }
  assert_true(total > 0.99)
  assert_true(total < 1.01)
}

///|
test "spectrum_summary" {
  let mutations = @src.mutational_patterns_sample_data()
  let sample_ids = ["S1", "S2"]
  let matrix = @src.build_mutation_matrix(mutations, sample_ids)
  let summary = matrix.spectrum_summary()
  assert_true(summary.length() > 0)
}
