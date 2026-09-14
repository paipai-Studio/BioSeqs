///|
/// Tests for dada2 module.

///|
test "dada2 nt_to_int" {
  assert_eq(@src.nt_to_int('A'), 0)
  assert_eq(@src.nt_to_int('C'), 1)
  assert_eq(@src.nt_to_int('G'), 2)
  assert_eq(@src.nt_to_int('T'), 3)
  assert_eq(@src.nt_to_int('a'), 0)
  assert_eq(@src.nt_to_int('N'), -1)
}

///|
test "dada2 int_to_nt" {
  assert_eq(@src.int_to_nt(0), "A")
  assert_eq(@src.int_to_nt(1), "C")
  assert_eq(@src.int_to_nt(2), "G")
  assert_eq(@src.int_to_nt(3), "T")
  assert_eq(@src.int_to_nt(4), "N")
}

///|
test "dada2 hamming distance" {
  assert_eq(@src.dada_hamming_distance("ACGT", "ACGT"), 0)
  assert_eq(@src.dada_hamming_distance("ACGT", "ACGA"), 1)
  assert_eq(@src.dada_hamming_distance("ACGT", "TGCA"), 4)
  assert_eq(@src.dada_hamming_distance("ACGT", "ACGTA"), -1)
}

///|
test "dada2 dereplicate basic" {
  let seqs = ["ACGT", "ACGT", "ACGT", "TTGC", "TTGC", "GGCC"]
  let derep = @src.dereplicate(seqs, [])
  assert_eq(derep.total_reads, 6)
  assert_eq(derep.uniques.length(), 3)
  // Most abundant first
  assert_eq(derep.uniques[0].sequence, "ACGT")
  assert_eq(derep.uniques[0].abundance, 3)
  assert_eq(derep.uniques[1].sequence, "TTGC")
  assert_eq(derep.uniques[1].abundance, 2)
  assert_eq(derep.uniques[2].sequence, "GGCC")
  assert_eq(derep.uniques[2].abundance, 1)
}

///|
test "dada2 dereplicate with quality" {
  let seqs = ["ACGT", "ACGT", "ACGT"]
  let quals = [[10.0, 20.0, 30.0, 40.0], [20.0, 20.0, 30.0, 40.0], [30.0, 20.0, 30.0, 40.0]]
  let derep = @src.dereplicate(seqs, quals)
  assert_eq(derep.uniques[0].abundance, 3)
  // Average quality: pos0 = (10+20+30)/3 = 20
  assert_true((derep.uniques[0].qual_scores[0] - 20.0).abs() < 0.001)
}

///|
test "dada2 quality profile" {
  let quals = [
    [10.0, 20.0, 30.0],
    [20.0, 30.0, 40.0],
    [30.0, 40.0, 50.0],
  ]
  let profile = @src.quality_profile(quals)
  assert_eq(profile.positions, 3)
  // Mean: pos0 = (10+20+30)/3 = 20
  assert_true((profile.mean[0] - 20.0).abs() < 0.001)
  assert_eq(profile.min[0], 10.0)
  assert_eq(profile.max[0], 30.0)
}

///|
test "dada2 init error model" {
  let model = @src.init_error_model(40)
  assert_eq(model.max_qual, 40)
  // Diagonal (match) rates should be close to 1
  assert_true(model.err_mat[0][20] > 0.99) // AA at qual 20
  // Off-diagonal (mismatch) rates should be small
  assert_true(model.err_mat[1][20] < 0.01) // AC at qual 20
}

///|
test "dada2 transition index" {
  assert_eq(@src.transition_index(0, 0), 0) // AA
  assert_eq(@src.transition_index(0, 1), 1) // AC
  assert_eq(@src.transition_index(3, 3), 15) // TT
  assert_eq(@src.transition_index(1, 2), 6) // CG
}

///|
test "dada2 abundance p-value" {
  // When observed >> expected, p should be very small (real variant)
  let p1 = @src.abundance_p_value(100, 1.0)
  assert_true(p1 < 0.001)
  // When observed << expected, p should be large (consistent with error)
  let p2 = @src.abundance_p_value(1, 100.0)
  assert_true(p2 > 0.99)
  // When observed ≈ expected, p should be moderate
  let p3 = @src.abundance_p_value(10, 10.0)
  assert_true(p3 > 0.3 && p3 < 0.7)
}

///|
test "dada2 reverse complement" {
  assert_eq(@src.dada_reverse_complement("ACGT"), "ACGT")
  assert_eq(@src.dada_reverse_complement("AAAA"), "TTTT")
  assert_eq(@src.dada_reverse_complement("ACGTACGT"), "ACGTACGT")
  assert_eq(@src.dada_reverse_complement("ATCG"), "CGAT")
}

///|
test "dada2 merge pair perfect overlap" {
  // Forward: ACGTACGTACGT (12bp)
  // Reverse: ACGTACGTACGT (12bp, rc = ACGTACGTACGT)
  // Full overlap = 12
  let fwd = "ACGTACGTACGT"
  let rev = "ACGTACGTACGT"
  let merged = @src.merge_pair(fwd, rev, min_overlap=8)
  assert_true(merged.accepted)
  assert_eq(merged.sequence, "ACGTACGTACGT")
}

///|
test "dada2 merge pair partial overlap" {
  // Forward: ACGTACGTACGTAAAA (16bp)
  // Reverse: TTTTACGTACGTACGT (16bp, rc = ACGTACGTACGTAAAA)
  // Overlap = 12
  let fwd = "ACGTACGTACGTAAAA"
  let rev = "TTTTACGTACGTACGT"
  let merged = @src.merge_pair(fwd, rev, min_overlap=8)
  assert_true(merged.accepted)
  assert_true(merged.overlap >= 8)
}

///|
test "dada2 merge pair no overlap" {
  let fwd = "ACGTACGTACGT"
  let rev = "TTTTGGGGCCCC"
  let merged = @src.merge_pair(fwd, rev, min_overlap=12)
  assert_true(!merged.accepted)
}

///|
test "dada2 is bimera true" {
  // Parent A: AAAAAAATTTTTTT
  // Parent B: GGGGGGGCCCCCCC
  // Chimera:  AAAAAAACCCCCCC (first half from A, second half from B)
  let parent_a = "AAAAAAATTTTTTT"
  let parent_b = "GGGGGGGCCCCCCC"
  let chimera = "AAAAAAACCCCCCC"
  assert_true(@src.is_bimera(chimera, parent_a, parent_b))
}

///|
test "dada2 is bimera false" {
  let parent_a = "AAAAAAAAAAAAA"
  let parent_b = "GGGGGGGGGGGGG"
  let not_chimera = "ACGTACGTACGTG"
  assert_true(!@src.is_bimera(not_chimera, parent_a, parent_b))
}

///|
test "dada2 remove chimeras" {
  let parent_a = "AAAAAAATTTTTTT"
  let parent_b = "GGGGGGGCCCCCCC"
  let chimera = "AAAAAAACCCCCCC"
  let asvs = [parent_a, parent_b, chimera]
  let abundances = [100, 100, 10]
  let (clean, clean_abund, flags) = @src.remove_chimeras(asvs, abundances)
  assert_eq(clean.length(), 2)
  assert_true(flags[2]) // chimera flagged
  assert_true(!flags[0]) // parent_a not flagged
  assert_true(!flags[1]) // parent_b not flagged
}

///|
test "dada2 extract kmers" {
  let kmers = @src.extract_kmers("ACGTAC", 3)
  assert_eq(kmers.length(), 4)
  assert_eq(kmers[0], "ACG")
  assert_eq(kmers[1], "CGT")
  assert_eq(kmers[2], "GTA")
  assert_eq(kmers[3], "TAC")
}

///|
test "dada2 extract kmers too short" {
  let kmers = @src.extract_kmers("AC", 4)
  assert_eq(kmers.length(), 0)
}

///|
test "dada2 learn errors" {
  let derep = @src.create_example_derep()
  let model = @src.learn_errors(derep, 40)
  assert_eq(model.max_qual, 40)
  assert_eq(model.n_reads, 12)
}

///|
test "dada2 dada inference" {
  let derep = @src.create_example_derep()
  let model = @src.learn_errors(derep, 40)
  let result = @src.dada_inference(derep, model)
  assert_true(result.n_clusters >= 1)
  assert_eq(result.n_input_reads, 12)
}

///|
test "dada2 make sequence table" {
  let derep = @src.create_example_derep()
  let model = @src.learn_errors(derep, 40)
  let result1 = @src.dada_inference(derep, model)
  let result2 = @src.dada_inference(derep, model)
  let table = @src.make_sequence_table([("sample1", result1), ("sample2", result2)])
  assert_true(table.n_asvs >= 1)
  assert_eq(table.n_samples, 2)
}

///|
test "dada2 seq table sample totals" {
  let derep = @src.create_example_derep()
  let model = @src.learn_errors(derep, 40)
  let result = @src.dada_inference(derep, model)
  let table = @src.make_sequence_table([("s1", result)])
  let totals = @src.seq_table_sample_totals(table)
  assert_eq(totals.length(), 1)
  assert_eq(totals[0], 12)
}

///|
test "dada2 seq table asv totals" {
  let derep = @src.create_example_derep()
  let model = @src.learn_errors(derep, 40)
  let result = @src.dada_inference(derep, model)
  let table = @src.make_sequence_table([("s1", result)])
  let totals = @src.seq_table_asv_totals(table)
  let mut sum = 0
  let mut i = 0
  while i < totals.length() {
    sum = sum + totals[i]
    i = i + 1
  }
  assert_eq(sum, 12)
}

///|
test "dada2 assign taxonomy" {
  let ref_db = @src.create_example_ref_db()
  let asvs = ["ACGTACGTACGTACGTACGTACGTACGTACGT"]
  let results = @src.assign_taxonomy(asvs, ref_db)
  assert_eq(results.length(), 1)
  // Should match E. coli reference
  assert_eq(results[0].taxonomy[0], "Bacteria")
  assert_eq(results[0].taxonomy[6], "coli")
  assert_true(results[0].confidence > 0.5)
}

///|
test "dada2 seq table to relative" {
  let derep = @src.create_example_derep()
  let model = @src.learn_errors(derep, 40)
  let result1 = @src.dada_inference(derep, model)
  let result2 = @src.dada_inference(derep, model)
  let table = @src.make_sequence_table([("s1", result1), ("s2", result2)])
  let rel = @src.seq_table_to_relative(table)
  // Each sample column should sum to ~1.0
  let mut col_sum = 0.0
  let mut ai = 0
  while ai < table.n_asvs {
    col_sum = col_sum + rel[ai][0]
    ai = ai + 1
  }
  assert_true((col_sum - 1.0).abs() < 0.001)
}

///|
test "dada2 filter table" {
  let derep = @src.create_example_derep()
  let model = @src.learn_errors(derep, 40)
  let result = @src.dada_inference(derep, model)
  let table = @src.make_sequence_table([("s1", result)])
  let flags = Array::make(table.n_asvs, false)
  let filtered = @src.filter_table(table, flags)
  assert_eq(filtered.n_asvs, table.n_asvs)
}

///|
test "dada2 merge pairs batch" {
  let fwd = ["ACGTACGTACGTACGT", "TTTTGGGGCCCCAAAA"]
  let rev = ["ACGTACGTACGTACGT", "TTTTGGGGCCCCAAAA"]
  let results = @src.merge_pairs(fwd, rev, min_overlap=8)
  assert_eq(results.length(), 2)
  assert_true(results[0].accepted)
  assert_true(results[1].accepted)
}

///|
test "dada2 expected abundance" {
  let model = @src.init_error_model(40)
  // Same sequence: expected abundance = true_abund (all match rates ~1)
  let exp = @src.expected_abundance("ACGT", "ACGT", 100, model)
  assert_true(exp > 90.0)
  // Different sequence: expected abundance should be very small
  let exp2 = @src.expected_abundance("ACGT", "TGCA", 100, model)
  assert_true(exp2 < 1.0)
}

///|
test "dada2 example derep" {
  let derep = @src.create_example_derep()
  assert_eq(derep.total_reads, 12)
  // Most abundant: ACGTACGTACGTACGT with 5 reads
  assert_eq(derep.uniques[0].sequence, "ACGTACGTACGTACGT")
  assert_eq(derep.uniques[0].abundance, 5)
}