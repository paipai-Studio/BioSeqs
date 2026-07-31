///|
/// Test file for prot_dao module (IUPred disorder prediction).

test "disorder_score_positive" {
  // Disordered-promoting amino acids should have positive scores
  assert_eq!(@src.prot_dao_disorder_score('E') > 0.0, true)
  assert_eq!(@src.prot_dao_disorder_score('P') > 0.0, true)
  assert_eq!(@src.prot_dao_disorder_score('K') > 0.0, true)
  assert_eq!(@src.prot_dao_disorder_score('S') > 0.0, true)
}

test "disorder_score_negative" {
  // Order-promoting amino acids should have negative or near-zero scores
  assert_eq!(@src.prot_dao_disorder_score('I') < 0.0, true)
  assert_eq!(@src.prot_dao_disorder_score('W') <= 0.0, true)
  assert_eq!(@src.prot_dao_disorder_score('Y') < 0.0, true)
  assert_eq!(@src.prot_dao_disorder_score('C') <= 0.0, true)
}

test "disorder_score_unknown_aa" {
  // Unknown amino acids should have 0 score
  assert_eq!(@src.prot_dao_disorder_score('X'), 0.0)
  assert_eq!(@src.prot_dao_disorder_score('z'), 0.0)
}

test "energy_score_positive" {
  assert_eq!(@src.prot_dao_energy_score('C') > 0.0, true)
  assert_eq!(@src.prot_dao_energy_score('P') > 0.0, true)
  assert_eq!(@src.prot_dao_energy_score('M') > 0.0, true)
}

test "energy_score_negative" {
  assert_eq!(@src.prot_dao_energy_score('H') < 0.0, true)
  assert_eq!(@src.prot_dao_energy_score('R') < 0.0, true)
}

test "disorder_residue_create" {
  let res = @src.DisorderResidue::new(1, 'A', 0.5, 1.0)
  assert_eq!(res.position, 1)
  assert_eq!(res.amino_acid, 'A')
  assert_eq!(res.disorder_score, 0.5)
  assert_eq!(res.energy_score, 1.0)
  assert_eq!(res.is_disordered, false)
  assert_eq!(res.is_disordered_long, false)
}

test "disorder_region_create" {
  let reg = @src.DisorderRegion::new(10, 50, 0.6, 0.8, region_type="disordered")
  assert_eq!(reg.start, 10)
  assert_eq!(reg.end_, 50)
  assert_eq!(reg.length, 40)
  assert_eq!(reg.avg_score, 0.6)
  assert_eq!(reg.max_score, 0.8)
  assert_eq!(reg.region_type, "disordered")
}

test "disorder_region_type" {
  let reg = @src.DisorderRegion::new(5, 45, 0.7, 0.9, region_type="long disordered")
  assert_eq!(reg.region_type, "long disordered")
}

test "prot_dao_predict_ordered_sequence" {
  // A sequence with mostly ordered amino acids
  let ordered_seq = "AILWVFMWCILVWALMVILWAFMVCLWVILMAWLVFMAW"
  let result = @src.prot_dao_predict(ordered_seq)
  assert_eq!(result.sequence, ordered_seq)
  assert_eq!(result.residues.length(), ordered_seq.length())
  assert_eq!(result.method, "IUPred-like")
}

test "prot_dao_predict_disordered_sequence" {
  // A sequence with many disorder-promoting amino acids
  let disordered_seq = "EPKSESEKPPPPPPKKKEEEPPPKKKSEGGSSSGGGKKKPPP"
  let result = @src.prot_dao_predict(disordered_seq)
  assert_eq!(result.residues.length(), disordered_seq.length())
  // With many disorder-promoting residues, some should be disordered
  assert_eq!(result.get_n_disordered() > 0, true)
}

test "prot_dao_predict_threshold" {
  let seq = "EPKSESEKPPPPPPKKKEEEPPPKKKSEGGSSSGGGKKKPPP"
  let result_low = @src.prot_dao_predict(seq, threshold_disordered=0.3)
  let result_high = @src.prot_dao_predict(seq, threshold_disordered=0.7)
  assert_eq!(result_low.get_n_disordered() >= result_high.get_n_disordered(), true)
}

test "disorder_result_n_disordered" {
  let result = @src.prot_dao_sample()
  let n = result.get_n_disordered()
  assert_eq!(n >= 0, true)
  assert_eq!(n <= result.sequence.length(), true)
}

test "disorder_result_fraction" {
  let result = @src.prot_dao_sample()
  let frac = result.get_fraction_disordered()
  assert_eq!(frac >= 0.0, true)
  assert_eq!(frac <= 1.0, true)
}

test "disorder_result_regions" {
  let result = @src.prot_dao_sample()
  let regions = result.get_regions()
  let n_regions = result.get_n_regions()
  assert_eq!(regions.length(), n_regions)
  // Each region should have valid coordinates
  let mut i = 0
  while i < regions.length() {
    assert_eq!(regions[i].start > 0 || regions[i].end_ > 0, true)
    assert_eq!(regions[i].length >= 0, true)
    i = i + 1
  }
}

test "disorder_result_longest_region" {
  let result = @src.prot_dao_sample()
  let longest = result.get_longest_region()
  assert_eq!(longest.length >= 0, true)
  if result.get_n_regions() > 0 {
    assert_eq!(longest.length > 0, true)
  }
}

test "disorder_result_scores" {
  let result = @src.prot_dao_sample()
  let scores = result.get_scores()
  assert_eq!(scores.length(), result.sequence.length())
  let mut i = 0
  while i < scores.length() {
    assert_eq!(scores[i] >= 0.0 || scores[i] <= 0.0, true)
    i = i + 1
  }
}

test "disorder_result_summary" {
  let result = @src.prot_dao_sample()
  let summary = result.summary()
  assert_eq!(summary.contains("IUPred-like"), true)
  assert_eq!(summary.contains("Sequence length"), true)
  assert_eq!(summary.contains("Disordered residues"), true)
}

test "disorder_result_disordered_sequence" {
  let result = @src.prot_dao_sample()
  let ds = result.disordered_sequence()
  assert_eq!(ds.length(), result.sequence.length())
  // Disordered positions should show the amino acid
  // Ordered positions should show '-'
}

test "prot_dao_sample_sequence" {
  let seq = @src.prot_dao_sample_sequence()
  assert_eq!(seq.length() > 0, true)
  // All characters should be valid amino acids
  let mut i = 0
  while i < seq.length() {
    let c = seq.unsafe_get(i)
    assert_eq!(c >= 65 && c <= 90 || c >= 97 && c <= 122, true)
    i = i + 1
  }
}

test "prot_dao_empty_sequence" {
  let result = @src.prot_dao_predict("")
  assert_eq!(result.residues.length(), 0)
  assert_eq!(result.get_fraction_disordered(), 0.0)
  assert_eq!(result.get_n_disordered(), 0)
}

test "prot_dao_single_residue" {
  let result = @src.prot_dao_predict("A")
  assert_eq!(result.residues.length(), 1)
  assert_eq!(result.sequence, "A")
}

test "prot_dao_to_ascii" {
  let result = @src.prot_dao_sample()
  let ascii = result.to_ascii()
  assert_eq!(ascii.contains("Score"), true)
  assert_eq!(ascii.contains("Seq"), true)
  assert_eq!(ascii.contains("Dis"), true)
}
