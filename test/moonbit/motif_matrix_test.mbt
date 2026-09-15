///|
/// Tests for motif_matrix module (Bio.motifs.matrix).

test "motif_pfm_from_sites basic counts" {
  let alphabet = ["A", "C", "G", "T"]
  let sites = ["ACGT", "ACGT", "ACGA"]
  let pfm = @src.motif_pfm_from_sites(alphabet, sites)
  assert_eq(pfm.alphabet.length(), 4)
  assert_eq(pfm.counts.length(), 4)
  assert_eq(pfm.n_seqs, 3.0)
  // Position 0: A=3, C=0, G=0, T=0
  assert_eq(pfm.counts[0][0], 3.0)
  assert_eq(pfm.counts[0][1], 0.0)
  // Position 2: G=3
  assert_eq(pfm.counts[2][2], 3.0)
}

test "motif_pfm_to_ppm normalization" {
  let alphabet = ["A", "C", "G", "T"]
  let sites = ["ACGT", "ACGT", "ACGA"]
  let pfm = @src.motif_pfm_from_sites(alphabet, sites)
  let ppm = @src.motif_pfm_to_ppm(pfm, pseudocount=0.0)
  // Each position should sum to 1.0
  let mut p = 0
  while p < 4 {
    let sum = ppm.values[p][0] + ppm.values[p][1] + ppm.values[p][2] + ppm.values[p][3]
    assert_true((sum - 1.0).abs() < 1.0e-9)
    p = p + 1
  }
  // Position 0: A=1.0
  assert_true((ppm.values[0][0] - 1.0).abs() < 1.0e-9)
}

test "motif_pfm_to_ppm with pseudocount" {
  let alphabet = ["A", "C", "G", "T"]
  let sites = ["ACGT"]
  let pfm = @src.motif_pfm_from_sites(alphabet, sites)
  let ppm = @src.motif_pfm_to_ppm(pfm, pseudocount=1.0)
  // With pseudocount=1, total = 1 + 4 = 5
  // Position 0: A=(1+1)/5=0.4, others=(0+1)/5=0.2
  assert_true((ppm.values[0][0] - 0.4).abs() < 1.0e-9)
  assert_true((ppm.values[0][1] - 0.2).abs() < 1.0e-9)
}

test "motif_ppm_to_pwm log_odds" {
  let alphabet = ["A", "C", "G", "T"]
  let sites = ["AAAA", "AAAA", "AAAA", "AAAA"]
  let pfm = @src.motif_pfm_from_sites(alphabet, sites)
  let ppm = @src.motif_pfm_to_ppm(pfm, pseudocount=1.0)
  let pwm = @src.motif_ppm_to_pwm(ppm)
  // With uniform background and A dominant, PWM[A] > 0, PWM[others] < 0
  assert_true(pwm.values[0][0] > 0.0)
  assert_true(pwm.values[0][1] < 0.0)
}

test "motif_pwm_to_pssm rounding" {
  let alphabet = ["A", "C", "G", "T"]
  let sites = ["ACGT", "ACGT", "ACGA"]
  let pfm = @src.motif_pfm_from_sites(alphabet, sites)
  let pwm = @src.motif_pfm_to_pwm(pfm, pseudocount=1.0)
  let pssm = @src.motif_pwm_to_pssm(pwm, scale=100.0)
  assert_eq(pssm.values.length(), 4)
  assert_eq(pssm.scale, 100.0)
}

test "motif_build_pssm pipeline" {
  let alphabet = ["A", "C", "G", "T"]
  let sites = ["ACGTT", "ACGTT", "ACGTA"]
  let pssm = @src.motif_build_pssm(alphabet, sites, pseudocount=1.0, scale=100.0)
  assert_eq(pssm.values.length(), 5)
  assert_eq(pssm.alphabet.length(), 4)
}

test "motif_ppm_ic_per_position" {
  let alphabet = ["A", "C", "G", "T"]
  let sites = ["AAAA", "AAAA", "AAAA", "AAAA"]
  let pfm = @src.motif_pfm_from_sites(alphabet, sites)
  let ppm = @src.motif_pfm_to_ppm(pfm, pseudocount=0.0)
  let ic = @src.motif_ppm_ic_per_position(ppm)
  // All-A position with no pseudocount: IC = 2 bits (fully conserved)
  assert_true((ic[0] - 2.0).abs() < 1.0e-9)
}

test "motif_ppm_total_ic" {
  let alphabet = ["A", "C", "G", "T"]
  let sites = ["AAAA", "AAAA", "AAAA", "AAAA"]
  let pfm = @src.motif_pfm_from_sites(alphabet, sites)
  let ppm = @src.motif_pfm_to_ppm(pfm, pseudocount=0.0)
  let total = @src.motif_ppm_total_ic(ppm)
  // 4 positions × 2 bits = 8 bits
  assert_true((total - 8.0).abs() < 1.0e-9)
}

test "motif_pwm_score sequence" {
  let alphabet = ["A", "C", "G", "T"]
  let sites = ["ACGT", "ACGT", "ACGT"]
  let pfm = @src.motif_pfm_from_sites(alphabet, sites)
  let pwm = @src.motif_pfm_to_pwm(pfm, pseudocount=1.0)
  let score = @src.motif_pwm_score(pwm, "ACGT")
  // Matching sequence should have positive score
  assert_true(score > 0.0)
}

test "motif_pssm_score integer" {
  let alphabet = ["A", "C", "G", "T"]
  let sites = ["ACGT", "ACGT", "ACGT"]
  let pssm = @src.motif_build_pssm(alphabet, sites, pseudocount=1.0, scale=100.0)
  let score = @src.motif_pssm_score(pssm, "ACGT")
  assert_true(score > 0)
}

test "motif_pwm_scan threshold" {
  let alphabet = ["A", "C", "G", "T"]
  let sites = ["ACGT", "ACGT", "ACGT"]
  let pfm = @src.motif_pfm_from_sites(alphabet, sites)
  let pwm = @src.motif_pfm_to_pwm(pfm, pseudocount=1.0)
  let hits = @src.motif_pwm_scan(pwm, "TTACGTTT", 0.0)
  assert_true(hits.length() >= 1)
}

test "motif_ppm_consensus" {
  let alphabet = ["A", "C", "G", "T"]
  let sites = ["ACGT", "ACGT", "ACGA"]
  let pfm = @src.motif_pfm_from_sites(alphabet, sites)
  let ppm = @src.motif_pfm_to_ppm(pfm, pseudocount=0.0)
  let cons = @src.motif_ppm_consensus(ppm)
  assert_eq(cons, "ACGT")
}

test "motif_pwm_reverse_complement" {
  let alphabet = ["A", "C", "G", "T"]
  let sites = ["ACGT", "ACGT", "ACGT"]
  let pfm = @src.motif_pfm_from_sites(alphabet, sites)
  let pwm = @src.motif_pfm_to_pwm(pfm, pseudocount=1.0)
  let rc = @src.motif_pwm_reverse_complement(pwm)
  assert_eq(rc.values.length(), 4)
  assert_eq(rc.alphabet.length(), 4)
}

test "motif_matrix_example" {
  let pwm = @src.motif_matrix_example()
  assert_eq(pwm.values.length(), 5)
  assert_eq(pwm.alphabet.length(), 4)
}

test "motif_pfm_to_string format" {
  let alphabet = ["A", "C", "G", "T"]
  let sites = ["ACGT", "ACGT"]
  let pfm = @src.motif_pfm_from_sites(alphabet, sites)
  let s = @src.motif_pfm_to_string(pfm)
  assert_true(s.contains("PFM"))
  assert_true(s.contains("A:"))
}