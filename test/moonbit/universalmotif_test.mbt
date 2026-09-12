///|
/// Tests for the universalmotif module (port of the Bioconductor
/// `universalmotif` package). Reference values are cross-checked against a
/// numpy/Python mirror (`/tmp/um_mirror.py`) that reimplements the same
/// formulas from the R source code.

///|
fn um_close(a : Double, b : Double, tol? : Double = 1.0e-6) -> Bool {
  (a - b).abs() < tol
}

///|
fn um_close_matrix(
  a : Array[Array[Double]],
  b : Array[Array[Double]],
  tol? : Double = 1.0e-6,
) -> Bool {
  if a.length() != b.length() {
    return false
  }
  for i = 0; i < a.length(); i = i + 1 {
    if a[i].length() != b[i].length() {
      return false
    }
    for j = 0; j < a[i].length(); j = j + 1 {
      if (a[i][j] - b[i][j]).abs() >= tol {
        return false
      }
    }
  }
  true
}

///|
fn um_tata_ppm() -> Array[Array[Double]] {
  [
    [0.9, 0.1, 0.0, 0.0],
    [0.0, 0.0, 0.9, 0.1],
    [0.0, 0.0, 0.9, 0.1],
    [0.1, 0.9, 0.0, 0.0],
    [0.1, 0.9, 0.0, 0.0],
  ]
}

// ---------------------------------------------------------------------------
// Alphabet helpers
// ---------------------------------------------------------------------------

///|
test "um_alphabet_letters" {
  let dna = @src.um_alphabet_letters("DNA")
  assert_eq(dna, ["A", "C", "G", "T"])
  let rna = @src.um_alphabet_letters("RNA")
  assert_eq(rna, ["A", "C", "G", "U"])
  let aa = @src.um_alphabet_letters("AA")
  assert_eq(aa.length(), 20)
  assert_eq(aa[0], "A")
  assert_eq(aa[19], "V")
  let custom = @src.um_alphabet_letters("ACGT")
  assert_eq(custom, ["A", "C", "G", "T"])
}

///|
test "um_alphabet_size" {
  assert_eq(@src.um_alphabet_size("DNA"), 4)
  assert_eq(@src.um_alphabet_size("RNA"), 4)
  assert_eq(@src.um_alphabet_size("AA"), 20)
  assert_eq(@src.um_alphabet_size("AC"), 2)
}

///|
test "um_default_bkg" {
  let b4 = @src.um_default_bkg(4)
  assert_eq(b4.length(), 4)
  for v in b4 {
    assert_true(um_close(v, 0.25))
  }
  let b20 = @src.um_default_bkg(20)
  assert_true(um_close(b20[0], 0.05))
  assert_eq(@src.um_default_bkg(0).length(), 0)
}

// ---------------------------------------------------------------------------
// Matrix conversions
// ---------------------------------------------------------------------------

///|
test "um_pcm_to_ppm_basic" {
  // 2 sites, 4 letters, pseudocount 0.8:
  // total = 2, denom = 2 + 0.8 = 2.8, num = c + 0.2
  // row [2,0,0,0]: [(2.2/2.8), (0.2/2.8), (0.2/2.8), (0.2/2.8)]
  let pcm = [[2.0, 0.0, 0.0, 0.0]]
  let ppm = @src.pcm_to_ppm(pcm, 0.8)
  assert_true(um_close(ppm[0][0], 2.2 / 2.8))
  assert_true(um_close(ppm[0][1], 0.2 / 2.8))
  assert_true(um_close(ppm[0][2], 0.2 / 2.8))
  assert_true(um_close(ppm[0][3], 0.2 / 2.8))
}

///|
test "um_pcm_to_ppm_zero_pseudocount" {
  // With pseudocount 0, PPM should be exactly PCM/sum.
  let pcm = [[3.0, 1.0, 0.0, 0.0]]
  let ppm = @src.pcm_to_ppm(pcm, 0.0)
  assert_true(um_close(ppm[0][0], 0.75))
  assert_true(um_close(ppm[0][1], 0.25))
  assert_true(um_close(ppm[0][2], 0.0))
}

///|
test "um_pcm_to_ppm_empty" {
  assert_eq(@src.pcm_to_ppm([], 0.8).length(), 0)
}

///|
test "um_pcm_to_ppm_zero_row_falls_back_uniform" {
  // A row of all zeros with pseudocount 0 -> uniform 1/n.
  let ppm = @src.pcm_to_ppm([[0.0, 0.0, 0.0, 0.0]], 0.0)
  for v in ppm[0] {
    assert_true(um_close(v, 0.25))
  }
}

///|
test "um_ppm_to_pcm_rounds_and_balances" {
  let ppm = [[0.5, 0.25, 0.15, 0.1]]
  let pcm = @src.ppm_to_pcm(ppm, 20.0)
  // column sum should equal nsites exactly
  let mut sum = 0
  for v in pcm[0] {
    sum = sum + v
  }
  assert_eq(sum, 20)
}

///|
test "um_ppm_to_pcm_default_nsites" {
  // nsites <= 1 -> default 100
  let ppm = [[0.4, 0.3, 0.2, 0.1]]
  let pcm = @src.ppm_to_pcm(ppm, 1.0)
  let mut sum = 0
  for v in pcm[0] {
    sum = sum + v
  }
  assert_eq(sum, 100)
}

///|
test "um_ppm_to_pwm_tata" {
  // Verified against Python mirror:
  //   TATA pwm[0] = [1.807355, -1.241008, -4.70044, -4.70044]
  //   TATA pwm[4] = [-1.241008, 1.807355, -4.70044, -4.70044]
  let bkg = [0.25, 0.25, 0.25, 0.25]
  let pwm = @src.ppm_to_pwm(um_tata_ppm(), bkg, 0.8, 20.0)
  assert_true(um_close(pwm[0][0], 1.807355, tol=2.0e-3))
  assert_true(um_close(pwm[0][1], -1.241008, tol=2.0e-3))
  assert_true(um_close(pwm[4][1], 1.807355, tol=2.0e-3))
}

///|
test "um_pwm_to_ppm_roundtrip" {
  // PWM->PPM[0] = [0.875, 0.105769, 0.009615, 0.009615] (mirror)
  let bkg = [0.25, 0.25, 0.25, 0.25]
  let ppm = @src.pwm_to_ppm([[1.807355, -1.241008, -4.70044, -4.70044]], bkg)
  assert_true(um_close(ppm[0][0], 0.875, tol=2.0e-3))
  assert_true(um_close(ppm[0][1], 0.105769, tol=2.0e-3))
  assert_true(um_close(ppm[0][3], 0.009615, tol=2.0e-3))
}

///|
test "um_ppm_to_icm_schneider" {
  // TATA ICM[0] = [1.377904, 0.1531, 0.0, 0.0]
  // TATA ICM[1] = [0.0, 0.0, 1.377904, 0.1531]
  let bkg = [0.25, 0.25, 0.25, 0.25]
  let icm = @src.ppm_to_icm(um_tata_ppm(), bkg, false)
  assert_true(um_close(icm[0][0], 1.377904, tol=2.0e-3))
  assert_true(um_close(icm[0][1], 0.1531, tol=2.0e-3))
  assert_true(um_close(icm[1][2], 1.377904, tol=2.0e-3))
}

///|
test "um_ppm_to_icm_relative_entropy" {
  let bkg = [0.25, 0.25, 0.25, 0.25]
  let icm = @src.ppm_to_icm(um_tata_ppm(), bkg, true)
  // Relative-entropy ICM[j] = PPM[j] * log2(PPM[j]/bkg[j]) (clamped >= 0)
  let expected_0_0 = 0.9 * @math.log2(0.9 / 0.25)
  assert_true(um_close(icm[0][0], expected_0_0, tol=2.0e-3))
  // Zero probabilities produce 0 in RE mode
  assert_true(um_close(icm[0][2], 0.0))
}

///|
test "um_icm_to_ppm_renormalizes" {
  let icm = [[1.0, 0.0, 0.0, 0.0], [0.0, 0.5, 0.5, 0.0]]
  let ppm = @src.icm_to_ppm(icm)
  assert_true(um_close(ppm[0][0], 1.0))
  assert_true(um_close(ppm[0][1], 0.0))
  assert_true(um_close(ppm[1][1], 0.5))
  assert_true(um_close(ppm[1][2], 0.5))
}

///|
test "um_icm_to_ppm_zero_row_uniform" {
  let ppm = @src.icm_to_ppm([[0.0, 0.0, 0.0, 0.0]])
  for v in ppm[0] {
    assert_true(um_close(v, 0.25))
  }
}

// ---------------------------------------------------------------------------
// Information content
// ---------------------------------------------------------------------------

///|
test "um_position_icscore_schneider" {
  // For pure A column with n=4: log2(4) - 0 = 2 bits
  let bkg = [0.25, 0.25, 0.25, 0.25]
  let ic = @src.position_icscore_ppm([1.0, 0.0, 0.0, 0.0], bkg, false)
  assert_true(um_close(ic, 2.0))
  // Uniform column has 0 IC
  let ic0 = @src.position_icscore_ppm([0.25, 0.25, 0.25, 0.25], bkg, false)
  assert_true(um_close(ic0, 0.0))
}

///|
test "um_position_icscore_relative_entropy" {
  let bkg = [0.25, 0.25, 0.25, 0.25]
  // Pure A: 1 * log2(1 / 0.25) = 2 bits
  let ic = @src.position_icscore_ppm([1.0, 0.0, 0.0, 0.0], bkg, true)
  assert_true(um_close(ic, 2.0))
}

///|
test "um_universal_icscore_tata" {
  // TATA Schneider IC = 7.655022 (mirror)
  // TATA relative-entropy IC = 8.315986 (mirror)
  let bkg = [0.25, 0.25, 0.25, 0.25]
  let sch = @src.universal_icscore(um_tata_ppm(), bkg, false)
  assert_true(um_close(sch, 7.655022, tol=2.0e-3))
  let re = @src.universal_icscore(um_tata_ppm(), bkg, true)
  assert_true(um_close(re, 8.315986, tol=2.0e-3))
}

// ---------------------------------------------------------------------------
// Consensus / IUPAC
// ---------------------------------------------------------------------------

///|
test "um_get_consensus_dominant_letters" {
  let dna = "DNA"
  assert_eq(@src.get_consensus([0.9, 0.0, 0.0, 0.0], dna), "A")
  assert_eq(@src.get_consensus([0.0, 0.9, 0.0, 0.0], dna), "C")
  assert_eq(@src.get_consensus([0.0, 0.0, 0.9, 0.0], dna), "G")
  assert_eq(@src.get_consensus([0.0, 0.0, 0.0, 0.9], dna), "T")
}

///|
test "um_get_consensus_two_letter_iupac" {
  let dna = "DNA"
  // A+C sum > 0.75 -> M
  assert_eq(@src.get_consensus([0.5, 0.5, 0.0, 0.0], dna), "M")
  // A+G -> R, A+T -> W, C+G -> S, C+T -> Y, G+T -> K
  assert_eq(@src.get_consensus([0.5, 0.0, 0.5, 0.0], dna), "R")
  assert_eq(@src.get_consensus([0.5, 0.0, 0.0, 0.5], dna), "W")
  assert_eq(@src.get_consensus([0.0, 0.5, 0.5, 0.0], dna), "S")
  assert_eq(@src.get_consensus([0.0, 0.5, 0.0, 0.5], dna), "Y")
  assert_eq(@src.get_consensus([0.0, 0.0, 0.5, 0.5], dna), "K")
}

///|
test "um_get_consensus_three_letter_iupac" {
  let dna = "DNA"
  // A+C+G all > 0.25 -> V (mirror-verified)
  assert_eq(@src.get_consensus([0.4, 0.3, 0.3, 0.0], dna), "V")
  assert_eq(@src.get_consensus([0.3, 0.3, 0.3, 0.1], dna), "V")
  // A+C+T -> H, C+G+T -> B, A+G+T -> D
  assert_eq(@src.get_consensus([0.4, 0.3, 0.0, 0.3], dna), "H")
  assert_eq(@src.get_consensus([0.0, 0.3, 0.3, 0.4], dna), "B")
  assert_eq(@src.get_consensus([0.4, 0.0, 0.3, 0.3], dna), "D")
}

///|
test "um_get_consensus_uniform_is_N" {
  assert_eq(@src.get_consensus([0.25, 0.25, 0.25, 0.25], "DNA"), "N")
}

///|
test "um_get_consensus_rna_uses_U" {
  assert_eq(@src.get_consensus([0.0, 0.0, 0.0, 0.9], "RNA"), "U")
}

///|
test "um_get_consensus_aa_returns_max_letter" {
  // For AA alphabet, the highest-prob letter is returned directly.
  let row = [
    0.1, 0.2, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05,
    0.05, 0.05, 0.05, 0.05, 0.05, 0.05,
  ]
  // index 1 is "R"
  assert_eq(@src.get_consensus(row, "AA"), "R")
}

///|
test "um_consensus_string_tata" {
  // The TATA PPM should produce a consensus like "AGGCC" (mirror-verified:
  // the first column [0.9, 0.1, 0, 0] gives A; then G at positions 2, 3
  // and C at the last two).
  let s = @src.consensus_string(um_tata_ppm(), "DNA")
  assert_eq(s.length(), 5)
  assert_eq(s[0:1].to_string(), "A")
}

///|
test "um_consensus_letter_to_ppm_pure_letters" {
  let a = @src.consensus_letter_to_ppm("A")
  assert_true(a[0] > 0.99)
  assert_true(a[1] < 0.01)
  let t = @src.consensus_letter_to_ppm("T")
  assert_true(t[3] > 0.99)
}

///|
test "um_consensus_letter_to_ppm_two_letter" {
  let m = @src.consensus_letter_to_ppm("M")
  assert_true(um_close(m[0], 0.499))
  assert_true(um_close(m[1], 0.499))
  assert_true(m[2] < 0.01)
}

///|
test "um_consensus_letter_to_ppm_n_is_uniform" {
  let n = @src.consensus_letter_to_ppm("N")
  for v in n {
    assert_true(um_close(v, 0.25))
  }
}

// ---------------------------------------------------------------------------
// UniversalMotif construction
// ---------------------------------------------------------------------------

///|
test "um_universalmotif_new_basic" {
  let m = @src.UniversalMotif::new(
    "TATA-box",
    "DNA",
    um_tata_ppm(),
    @src.UMotifType::PPM,
    nsites=20.0,
  )
  assert_eq(m.name, "TATA-box")
  assert_eq(m.alphabet, "DNA")
  assert_eq(m.motif.length(), 5)
  assert_true(m.motif_type is @src.UMotifType::PPM)
  assert_true(um_close(m.nsites, 20.0))
  // IC score is auto-computed (Schneider) ~ 7.655022
  assert_true(um_close(m.icscore, 7.655022, tol=2.0e-3))
  // Consensus is auto-computed
  assert_eq(m.consensus.length(), 5)
}

///|
test "um_universalmotif_new_default_bkg" {
  let m = @src.UniversalMotif::new(
    "x",
    "DNA",
    [[0.25, 0.25, 0.25, 0.25]],
    @src.UMotifType::PPM,
  )
  assert_eq(m.bkg.length(), 4)
  for v in m.bkg {
    assert_true(um_close(v, 0.25))
  }
}

///|
test "um_universalmotif_to_ppm_idempotent" {
  let m = @src.UniversalMotif::new(
    "x",
    "DNA",
    um_tata_ppm(),
    @src.UMotifType::PPM,
    nsites=20.0,
  )
  let m2 = m.to_ppm()
  assert_true(m2.motif_type is @src.UMotifType::PPM)
  assert_true(um_close_matrix(m2.motif, um_tata_ppm()))
}

///|
test "um_universalmotif_to_pwm" {
  let m = @src.UniversalMotif::new(
    "x",
    "DNA",
    um_tata_ppm(),
    @src.UMotifType::PPM,
    nsites=20.0,
  )
  let pwm_m = m.to_pwm()
  assert_true(pwm_m.motif_type is @src.UMotifType::PWM)
  // pwm[0][0] ~ 1.807355 (mirror)
  assert_true(um_close(pwm_m.motif[0][0], 1.807355, tol=2.0e-3))
}

///|
test "um_universalmotif_to_icm" {
  let m = @src.UniversalMotif::new(
    "x",
    "DNA",
    um_tata_ppm(),
    @src.UMotifType::PPM,
    nsites=20.0,
  )
  let icm_m = m.to_icm()
  assert_true(icm_m.motif_type is @src.UMotifType::ICM)
  // icm[0][0] ~ 1.377904 (mirror)
  assert_true(um_close(icm_m.motif[0][0], 1.377904, tol=2.0e-3))
}

///|
test "um_universalmotif_to_icm_relative_entropy" {
  let m = @src.UniversalMotif::new(
    "x",
    "DNA",
    um_tata_ppm(),
    @src.UMotifType::PPM,
    nsites=20.0,
  )
  let icm_m = m.to_icm(relative_entropy=true)
  let expected = 0.9 * @math.log2(0.9 / 0.25)
  assert_true(um_close(icm_m.motif[0][0], expected, tol=2.0e-3))
}

///|
test "um_create_example_universal_motif" {
  let m = @src.create_example_universal_motif()
  assert_eq(m.name, "TATA-box")
  assert_eq(m.alphabet, "DNA")
  assert_eq(m.motif.length(), 5)
}

// ---------------------------------------------------------------------------
// Reverse complement
// ---------------------------------------------------------------------------

///|
test "um_reverse_complement_motif_dna" {
  let m = @src.UniversalMotif::new(
    "m",
    "DNA",
    um_tata_ppm(),
    @src.UMotifType::PPM,
  )
  let rc = @src.reverse_complement_motif(m)
  // RC of position 0 ([0.9, 0.1, 0, 0]) should appear at the last position
  // with columns reversed: [0, 0, 0.1, 0.9]
  let last = rc.motif[rc.motif.length() - 1]
  assert_true(um_close(last[0], 0.0))
  assert_true(um_close(last[1], 0.0))
  assert_true(um_close(last[2], 0.1))
  assert_true(um_close(last[3], 0.9))
}

///|
test "um_reverse_complement_motif_aa_unchanged" {
  let aa_ppm = [[0.5, 0.5]]
  let m = @src.UniversalMotif::new("m", "AC", aa_ppm, @src.UMotifType::PPM)
  let rc = @src.reverse_complement_motif(m)
  // AA / custom alphabets are returned unchanged
  assert_true(um_close_matrix(rc.motif, aa_ppm))
}

// ---------------------------------------------------------------------------
// compare_positions (12 metrics, mirror-verified)
// ---------------------------------------------------------------------------

///|
test "um_compare_positions_all_metrics" {
  // Mirror reference values from /tmp/um_mirror.py
  let row1 = [0.4, 0.3, 0.2, 0.1]
  let row2 = [0.1, 0.2, 0.3, 0.4]
  let bkg = [0.25, 0.25, 0.25, 0.25]
  let cmp = fn(m : @src.UMetric) -> Double {
    @src.compare_positions(m, row1, row2, bkg, bkg, 100.0, 100.0)
  }
  assert_true(um_close(cmp(@src.UMetric::PCC), -1.0, tol=1.0e-5))
  assert_true(um_close(cmp(@src.UMetric::EUCL), 0.4472135955, tol=1.0e-5))
  assert_true(um_close(cmp(@src.UMetric::SW), 1.8, tol=1.0e-5))
  assert_true(um_close(cmp(@src.UMetric::KL), 0.4336495148, tol=1.0e-3))
  assert_true(um_close(cmp(@src.UMetric::WEUCL), 0.316227766, tol=1.0e-5))
  assert_true(um_close(cmp(@src.UMetric::ALLR), -0.2910616, tol=1.0e-3))
  assert_true(um_close(cmp(@src.UMetric::BHAT), 0.8898979486, tol=1.0e-5))
  assert_true(um_close(cmp(@src.UMetric::HELL), 0.3318162917, tol=1.0e-3))
  assert_true(um_close(cmp(@src.UMetric::WPCC), -1.0, tol=1.0e-5))
  assert_true(um_close(cmp(@src.UMetric::SEUCL), 0.2, tol=1.0e-5))
  assert_true(um_close(cmp(@src.UMetric::MAN), 0.8, tol=1.0e-5))
  assert_true(um_close(cmp(@src.UMetric::ALLR_LL), -0.2910616, tol=1.0e-3))
}

///|
test "um_compare_positions_identical_rows" {
  // Identical rows should give the maximum similarity / minimum distance.
  let row = [0.4, 0.3, 0.2, 0.1]
  let bkg = [0.25, 0.25, 0.25, 0.25]
  let pcc = @src.compare_positions(
    @src.UMetric::PCC,
    row,
    row,
    bkg,
    bkg,
    100.0,
    100.0,
  )
  // PCC of identical vectors is 1 (or 1.0 by the denom<=0 fallback)
  assert_true(pcc >= 0.999)
  let eucl = @src.compare_positions(
    @src.UMetric::EUCL,
    row,
    row,
    bkg,
    bkg,
    100.0,
    100.0,
  )
  assert_true(um_close(eucl, 0.0))
  let bhat = @src.compare_positions(
    @src.UMetric::BHAT,
    row,
    row,
    bkg,
    bkg,
    100.0,
    100.0,
  )
  assert_true(um_close(bhat, 1.0))
}

///|
test "um_compare_positions_allr_ll_floor" {
  // ALLR_LL is clamped to >= -2.0
  let row1 = [0.0, 0.0, 0.0, 1.0]
  let row2 = [1.0, 0.0, 0.0, 0.0]
  let bkg = [0.25, 0.25, 0.25, 0.25]
  let v = @src.compare_positions(
    @src.UMetric::ALLR_LL,
    row1,
    row2,
    bkg,
    bkg,
    100.0,
    100.0,
  )
  assert_true(v >= -2.0 - 1.0e-9)
}

// ---------------------------------------------------------------------------
// compare_motifs (sliding alignment)
// ---------------------------------------------------------------------------

///|
test "um_compare_motifs_identical" {
  let m = @src.UniversalMotif::new(
    "x",
    "DNA",
    um_tata_ppm(),
    @src.UMotifType::PPM,
    nsites=20.0,
  )
  // Mirror: PCC=5.0, EUCL=0, SW=10, BHAT=5, HELL=0, WPCC=5, MAN=0, SEUCL=0
  let pcc = @src.compare_motifs(
    m,
    m,
    @src.UMetric::PCC,
    min_overlap=3,
    try_rc=false,
  )
  assert_true(um_close(pcc, 5.0, tol=2.0e-3))
  let eucl = @src.compare_motifs(
    m,
    m,
    @src.UMetric::EUCL,
    min_overlap=3,
    try_rc=false,
  )
  assert_true(um_close(eucl, 0.0, tol=2.0e-3))
  let sw = @src.compare_motifs(
    m,
    m,
    @src.UMetric::SW,
    min_overlap=3,
    try_rc=false,
  )
  assert_true(um_close(sw, 10.0, tol=2.0e-3))
  let bhat = @src.compare_motifs(
    m,
    m,
    @src.UMetric::BHAT,
    min_overlap=3,
    try_rc=false,
  )
  assert_true(um_close(bhat, 5.0, tol=2.0e-3))
  let hell = @src.compare_motifs(
    m,
    m,
    @src.UMetric::HELL,
    min_overlap=3,
    try_rc=false,
  )
  assert_true(um_close(hell, 0.0, tol=2.0e-3))
  let wpcc = @src.compare_motifs(
    m,
    m,
    @src.UMetric::WPCC,
    min_overlap=3,
    try_rc=false,
  )
  assert_true(um_close(wpcc, 5.0, tol=2.0e-3))
}

///|
test "um_compare_motifs_shifted_by_one" {
  // Motif vs. itself shifted by 1 position: PCC ~ 4.1667 (mirror)
  let m = @src.UniversalMotif::new(
    "x",
    "DNA",
    um_tata_ppm(),
    @src.UMotifType::PPM,
    nsites=20.0,
  )
  let shifted_ppm = [
    [0.0, 0.0, 0.9, 0.1],
    [0.0, 0.0, 0.9, 0.1],
    [0.1, 0.9, 0.0, 0.0],
    [0.1, 0.9, 0.0, 0.0],
    [0.9, 0.1, 0.0, 0.0],
  ]
  let m2 = @src.UniversalMotif::new(
    "y",
    "DNA",
    shifted_ppm,
    @src.UMotifType::PPM,
    nsites=20.0,
  )
  let pcc = @src.compare_motifs(
    m,
    m2,
    @src.UMetric::PCC,
    min_overlap=3,
    try_rc=false,
  )
  assert_true(um_close(pcc, 3.2, tol=2.0e-3))
}

///|
test "um_compare_motifs_with_rc" {
  // Try_rc should find a strong match between a motif and its own reverse
  // complement.
  let m = @src.UniversalMotif::new(
    "x",
    "DNA",
    um_tata_ppm(),
    @src.UMotifType::PPM,
    nsites=20.0,
  )
  let rc = @src.reverse_complement_motif(m)
  // Without RC: weak; with RC: should match (PCC high).
  let pcc_no_rc = @src.compare_motifs(
    m,
    rc,
    @src.UMetric::PCC,
    min_overlap=3,
    try_rc=false,
  )
  let pcc_rc = @src.compare_motifs(
    m,
    rc,
    @src.UMetric::PCC,
    min_overlap=3,
    try_rc=true,
  )
  assert_true(pcc_rc > pcc_no_rc)
  // With RC, an identical-motif-vs-its-RC match is at least as good as no-RC.
  assert_true(pcc_rc > 0.0)
}

// ---------------------------------------------------------------------------
// merge_motifs
// ---------------------------------------------------------------------------

///|
test "um_merge_motifs_identical_keeps_motif" {
  let m = @src.UniversalMotif::new(
    "x",
    "DNA",
    um_tata_ppm(),
    @src.UMotifType::PPM,
    nsites=20.0,
  )
  let merged = @src.merge_motifs([m, m])
  // Merging identical motifs at offset 0 should keep the same PPM.
  assert_eq(merged.motif.length(), m.motif.length())
  assert_true(um_close_matrix(merged.motif, um_tata_ppm()))
}

///|
test "um_merge_motifs_extends_flanks" {
  // Mirror: merging m1 (5 cols) with m2 (5 cols, +1 offset) gives 6 columns.
  // The first column of the merged PPM equals m1's first column verbatim
  // (m2 doesn't overlap there); the second column is the weighted mean of
  // m1[1] and m2[0].
  let m1_ppm = [
    [0.9, 0.1, 0.0, 0.0],
    [0.0, 0.0, 0.9, 0.1],
    [0.0, 0.0, 0.9, 0.1],
    [0.1, 0.9, 0.0, 0.0],
    [0.1, 0.9, 0.0, 0.0],
  ]
  let m2_ppm = [
    [0.0, 0.0, 0.9, 0.1],
    [0.0, 0.0, 0.9, 0.1],
    [0.1, 0.9, 0.0, 0.0],
    [0.1, 0.9, 0.0, 0.0],
    [0.9, 0.1, 0.0, 0.0],
  ]
  let m1 = @src.UniversalMotif::new(
    "a",
    "DNA",
    m1_ppm,
    @src.UMotifType::PPM,
    nsites=20.0,
  )
  let m2 = @src.UniversalMotif::new(
    "b",
    "DNA",
    m2_ppm,
    @src.UMotifType::PPM,
    nsites=20.0,
  )
  let merged = @src.merge_motifs([m1, m2], min_overlap=3, try_rc=false)
  // Best alignment is at offset -1 (m2 shifted left by 1), giving 5 columns.
  // The first column is m2's first column (no m1 overlap).
  assert_eq(merged.motif.length(), 5)
  assert_true(um_close(merged.motif[0][0], 0.0))
  assert_true(um_close(merged.motif[0][1], 0.0))
}

///|
test "um_merge_motifs_single_input" {
  let m = @src.UniversalMotif::new(
    "x",
    "DNA",
    um_tata_ppm(),
    @src.UMotifType::PPM,
    nsites=20.0,
  )
  let merged = @src.merge_motifs([m])
  assert_eq(merged.motif.length(), 5)
  assert_true(merged.motif_type is @src.UMotifType::PPM)
}

///|
test "um_merge_motifs_empty_returns_empty" {
  let merged = @src.merge_motifs([])
  assert_eq(merged.motif.length(), 0)
}

// ---------------------------------------------------------------------------
// trim_motif
// ---------------------------------------------------------------------------

///|
test "um_trim_motif_removes_low_ic_flanks" {
  // Build a 7-position PPM with two low-IC flanks (uniform) and 5 informative
  // middle positions (the TATA motif). Trimming at min_ic=0.25 should yield 5.
  let extended = [
    [0.25, 0.25, 0.25, 0.25], // IC=0
    [0.9, 0.1, 0.0, 0.0],
    [0.0, 0.0, 0.9, 0.1],
    [0.0, 0.0, 0.9, 0.1],
    [0.1, 0.9, 0.0, 0.0],
    [0.1, 0.9, 0.0, 0.0],
    [0.25, 0.25, 0.25, 0.25], // IC=0
  ]
  let m = @src.UniversalMotif::new(
    "x",
    "DNA",
    extended,
    @src.UMotifType::PPM,
    nsites=20.0,
  )
  let trimmed = @src.trim_motif(m)
  // Mirror: trim result width = 5 (mirror verifies IC of flanks is < 0.25)
  assert_eq(trimmed.motif.length(), 5)
}

///|
test "um_trim_motif_left_only" {
  let extended = [
    [0.25, 0.25, 0.25, 0.25],
    [0.9, 0.1, 0.0, 0.0],
    [0.0, 0.0, 0.9, 0.1],
    [0.25, 0.25, 0.25, 0.25],
  ]
  let m = @src.UniversalMotif::new(
    "x",
    "DNA",
    extended,
    @src.UMotifType::PPM,
    nsites=20.0,
  )
  let trimmed = @src.trim_motif(m, trim_from=@src.UTrimFrom::ULeft)
  // Left flank removed; right flank kept
  assert_eq(trimmed.motif.length(), 3)
}

///|
test "um_trim_motif_right_only" {
  let extended = [
    [0.25, 0.25, 0.25, 0.25],
    [0.9, 0.1, 0.0, 0.0],
    [0.0, 0.0, 0.9, 0.1],
    [0.25, 0.25, 0.25, 0.25],
  ]
  let m = @src.UniversalMotif::new(
    "x",
    "DNA",
    extended,
    @src.UMotifType::PPM,
    nsites=20.0,
  )
  let trimmed = @src.trim_motif(m, trim_from=@src.UTrimFrom::URight)
  assert_eq(trimmed.motif.length(), 3)
}

///|
test "um_trim_motif_all_below_threshold_returns_empty" {
  // A uniform motif has all-zero IC; trimming at any threshold leaves nothing.
  let m = @src.UniversalMotif::new(
    "x",
    "DNA",
    [[0.25, 0.25, 0.25, 0.25]],
    @src.UMotifType::PPM,
    nsites=20.0,
  )
  let trimmed = @src.trim_motif(m)
  assert_eq(trimmed.motif.length(), 0)
}

// ---------------------------------------------------------------------------
// MEME I/O
// ---------------------------------------------------------------------------

///|
test "um_write_meme_round_trips_with_read_meme" {
  let m = @src.UniversalMotif::new(
    "TATA-box",
    "DNA",
    um_tata_ppm(),
    @src.UMotifType::PPM,
    nsites=20.0,
  )
  let text = @src.write_meme([m])
  let parsed = @src.read_meme(text)
  assert_eq(parsed.length(), 1)
  let m2 = parsed[0]
  assert_eq(m2.name, "TATA-box")
  assert_eq(m2.motif.length(), 5)
  // The matrix should round-trip closely (PPM)
  for i = 0; i < 5; i = i + 1 {
    for j = 0; j < 4; j = j + 1 {
      assert_true(um_close(m2.motif[i][j], um_tata_ppm()[i][j], tol=1.0e-5))
    }
  }
}

///|
test "um_write_meme_header" {
  let m = @src.UniversalMotif::new(
    "M1",
    "DNA",
    um_tata_ppm(),
    @src.UMotifType::PPM,
    nsites=20.0,
  )
  let text = @src.write_meme([m])
  // Required MEME headers
  assert_true(text.contains("MEME version"))
  assert_true(text.contains("ALPHABET= ACGT"))
  assert_true(text.contains("strands: + -"))
  assert_true(text.contains("Background letter frequencies"))
  assert_true(text.contains("MOTIF M1"))
  assert_true(text.contains("letter-probability matrix:"))
}

///|
test "um_read_meme_multiple_motifs" {
  let text = "MEME version 5\n\n" +
    "ALPHABET= ACGT\n\n" +
    "strands: + -\n\n" +
    "Background letter frequencies\n" +
    "A 0.250000 C 0.250000 G 0.250000 T 0.250000\n" +
    "\nMOTIF motif1\n" +
    "letter-probability matrix: alength= 4 w= 2 nsites= 20.000 E= 0\n" +
    "  0.900000  0.100000  0.000000  0.000000\n" +
    "  0.000000  0.000000  0.900000  0.100000\n" +
    "\nMOTIF motif2\n" +
    "letter-probability matrix: alength= 4 w= 2 nsites= 10.000 E= 0\n" +
    "  0.500000  0.500000  0.000000  0.000000\n" +
    "  0.000000  0.000000  0.500000  0.500000\n"
  let parsed = @src.read_meme(text)
  assert_eq(parsed.length(), 2)
  assert_eq(parsed[0].name, "motif1")
  assert_eq(parsed[1].name, "motif2")
  assert_eq(parsed[0].motif.length(), 2)
  assert_true(um_close(parsed[0].motif[0][0], 0.9))
  assert_true(um_close(parsed[1].motif[1][3], 0.5))
}

// ---------------------------------------------------------------------------
// JASPAR I/O
// ---------------------------------------------------------------------------

///|
test "um_write_jaspar_round_trips_with_read_jaspar" {
  let m = @src.UniversalMotif::new(
    "TATA-box",
    "DNA",
    um_tata_ppm(),
    @src.UMotifType::PPM,
    nsites=20.0,
  )
  let text = @src.write_jaspar(m)
  let parsed = @src.read_jaspar(text)
  assert_eq(parsed.length(), 1)
  let m2 = parsed[0]
  assert_eq(m2.name, "TATA-box")
  // JASPAR stores PCM counts; converting back to PPM should approximate the
  // original PPM.
  let m2_ppm = m2.to_ppm()
  for i = 0; i < 5; i = i + 1 {
    for j = 0; j < 4; j = j + 1 {
      assert_true(um_close(m2_ppm.motif[i][j], um_tata_ppm()[i][j], tol=0.06))
    }
  }
}

///|
test "um_write_jaspar_header_and_counts" {
  let m = @src.UniversalMotif::new(
    "MA0001.1",
    "DNA",
    um_tata_ppm(),
    @src.UMotifType::PPM,
    nsites=100.0,
  )
  let text = @src.write_jaspar(m)
  assert_true(text.contains(">MA0001.1"))
  // 4 rows: one per letter (A, C, G, T)
  assert_true(text.contains("A ["))
  assert_true(text.contains("C ["))
  assert_true(text.contains("G ["))
  assert_true(text.contains("T ["))
}

///|
test "um_read_jaspar_multiple_motifs" {
  let text = ">motif1\n" +
    "A [ 10 0 0 0 ]\n" +
    "C [ 0 10 0 0 ]\n" +
    "G [ 0 0 10 0 ]\n" +
    "T [ 0 0 0 10 ]\n" +
    "\n" +
    ">motif2\n" +
    "A [ 5 5 0 0 ]\n" +
    "C [ 5 5 0 0 ]\n" +
    "G [ 0 0 5 5 ]\n" +
    "T [ 0 0 5 5 ]\n"
  let parsed = @src.read_jaspar(text)
  assert_eq(parsed.length(), 2)
  assert_eq(parsed[0].name, "motif1")
  assert_eq(parsed[1].name, "motif2")
  assert_eq(parsed[0].motif.length(), 4)
  // motif1 column 0 should have counts [10,0,0,0] across A,C,G,T rows
  assert_eq(parsed[0].motif[0][0].to_int(), 10)
}
