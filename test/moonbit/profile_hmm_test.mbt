///|
/// Tests for profile HMM construction and sequence search (HMMER core).

// ============================================================
// Test fixture: a small conserved protein MSA
// ============================================================

///|
/// A toy MSA of 4 "protein-like" sequences over the small alphabet
/// {A, C, G, T} treated as 4 residues. Using a 4-letter alphabet keeps the
/// fixture tiny while still exercising match/insert/delete logic.
fn prof_test_msa() -> Array[String] {
  // A clean, valid alignment: 4 same-length sequences with two point
  // mutations relative to the consensus "ACGTACGTACGT".
  [
    "ACGTACGTACGT",
    "ACGTACGTACGT",
    "TCGTACGTACGT", // T at position 1 instead of A
    "ACGTACGTACGA", // A at position 12 instead of T
  ]
}

fn prof_dna() -> Array[String] {
  @src.profile_dna_alphabet()
}

// ============================================================
// Construction tests
// ============================================================

test "build_profile_hmm_basic" {
  let msa = prof_test_msa()
  let profile = @src.build_profile_hmm(msa, prof_dna())
  assert_eq(profile.n_seqs(), 4)
  // 12 columns -> should have ~12 match columns (all columns conserved enough)
  assert_true(profile.length() > 0)
}

test "build_profile_hmm_match_emission_sums_to_one" {
  let msa = prof_test_msa()
  let profile = @src.build_profile_hmm(msa, prof_dna())
  let alphabet = prof_dna()
  // Every match state's emissions must sum to ~1.0 (within rounding tolerance).
  let mut k = 1
  while k <= profile.length() {
    let mut total = 0.0
    for res in alphabet {
      total = total + profile.match_emission(k, res)
    }
    assert_true(total > 0.999 && total < 1.001)
    k = k + 1
  }
}

test "build_profile_hmm_consensus_scores_highest" {
  // The consensus sequence (matching most columns) should score higher than
  // unrelated or heavily-mutated sequences.
  let msa = prof_test_msa()
  let profile = @src.build_profile_hmm(msa, prof_dna())

  let consensus = "ACGTACGTACGT" // exactly matches 2 of 4, close to all
  let mutated = "TTTTTTTTTTTT" // very different
  let random = "ATCGGCTAGGCA" // unrelated

  let s_consensus = @src.profile_forward_score(profile, consensus)
  let s_mutated = @src.profile_forward_score(profile, mutated)
  let s_random = @src.profile_forward_score(profile, random)

  assert_true(s_consensus > s_mutated)
  assert_true(s_consensus > s_random)
}

test "profile_forward_score_identical_non_negative" {
  let msa = prof_test_msa()
  let profile = @src.build_profile_hmm(msa, prof_dna())
  // A sequence in the MSA should score positively (better than null model).
  let s = @src.profile_forward_score(profile, "ACGTACGTACGT")
  assert_true(s >= 0.0)
}

// ============================================================
// Viterbi tests
// ============================================================

test "profile_viterbi_returns_alignment" {
  let msa = prof_test_msa()
  let profile = @src.build_profile_hmm(msa, prof_dna())
  let result = @src.profile_viterbi(profile, "ACGTACGTACGT")
  // The aligned target should be non-empty and contain only ACGT and '-'
  let t = result.target_seq
  assert_true(t.length() > 0)
  let mut ok = true
  for i in 0..<t.length() {
    let c = t[i]
    if c != 'A' && c != 'C' && c != 'G' && c != 'T' && c != '-' {
      ok = false
    }
  }
  assert_true(ok)
  // Path should contain state symbols M / I / D
  let p = result.path
  assert_true(p.contains("M") || p.contains("I") || p.contains("D"))
}

test "profile_viterbi_score_consistency" {
  // Viterbi score should be <= forward score (Viterbi is the single best path,
  // forward sums over all paths).
  let msa = prof_test_msa()
  let profile = @src.build_profile_hmm(msa, prof_dna())
  let v = @src.profile_viterbi(profile, "ACGTACGTACGT").score
  let f = @src.profile_forward_score(profile, "ACGTACGTACGT")
  assert_true(v <= f + 0.001)
}

// ============================================================
// Search / ranking tests
// ============================================================

test "hmmer_search_ranks_consensus_first" {
  let msa = prof_test_msa()
  let profile = @src.build_profile_hmm(msa, prof_dna())
  let names = ["consensus", "mutated", "unrelated"]
  let seqs = ["ACGTACGTACGT", "TTTTTTTTTTTT", "ATCGGCTAGGCA"]
  let hits = @src.hmmer_search(profile, names, seqs)
  assert_eq(hits.length(), 3)
  // Sorted descending by score -> consensus should be first
  assert_eq(hits[0].name, "consensus")
  // Check scores are monotonic non-increasing
  assert_true(hits[0].score >= hits[1].score)
  assert_true(hits[1].score >= hits[2].score)
}

test "hmmer_search_evalue_small_for_best" {
  let msa = prof_test_msa()
  let profile = @src.build_profile_hmm(msa, prof_dna())
  let names = ["consensus", "unrelated"]
  let seqs = ["ACGTACGTACGT", "ATCGGCTAGGCA"]
  let hits = @src.hmmer_search(profile, names, seqs)
  // Best hit should have a small E-value; worst a larger one.
  assert_true(hits[0].evalue < hits[1].evalue)
}

test "hmmer_search_empty_input" {
  let msa = prof_test_msa()
  let profile = @src.build_profile_hmm(msa, prof_dna())
  let hits = @src.hmmer_search(profile, [], [])
  assert_eq(hits.length(), 0)
}

// ============================================================
// Background model test
// ============================================================

test "profile_aa_alphabet_length" {
  let aa = @src.profile_aa_alphabet()
  assert_eq(aa.length(), 20)
  let dna = @src.profile_dna_alphabet()
  assert_eq(dna.length(), 4)
}

// ============================================================
// Insert / delete handling: single-column MSA edge cases
// ============================================================

test "build_profile_hmm_single_column" {
  // A single conserved column.
  let msa = ["A", "A", "A"]
  let profile = @src.build_profile_hmm(msa, prof_dna())
  assert_eq(profile.length(), 1)
  assert_eq(profile.n_seqs(), 3)
  // "A" should score higher than "T"
  let sA = @src.profile_forward_score(profile, "A")
  let sT = @src.profile_forward_score(profile, "T")
  assert_true(sA > sT)
}

test "build_profile_hmm_insert_column_detection" {
  // A column that is mostly gaps should be classified as an insert column.
  // 4 sequences, one column has only 1 residue (gap fraction 3/4 = 0.75 > 0.5).
  let msa = [
    "AC-A",
    "AC-A",
    "ACTA",
    "AC-A",
  ]
  let profile = @src.build_profile_hmm(msa, prof_dna())
  // Only 3 columns are match columns (the 3rd is insert)
  assert_eq(profile.length(), 3)
}