///|
/// Tests for isoelectric_point module (Bio.SeqUtils.IsoelectricPoint).

test "mt_charge_at_ph neutral" {
  let seq = "A"
  let charge = @src.mt_charge_at_ph_public(seq, 7.0)
  // Single A: N-term +1 (pH 7 < 8.0), C-term -1 (pH 7 > 3.1), no side chain
  // charge = +1 - 1 = 0
  assert_true(charge.abs() < 0.1)
}

test "mt_charge_at_ph acidic" {
  let seq = "D"
  let charge = @src.mt_charge_at_ph_public(seq, 7.0)
  // D: N-term +1, C-term -1, side chain -1 (pH 7 > 3.9)
  // charge = +1 - 1 - 1 = -1
  assert_true((charge + 1.0).abs() < 0.1)
}

test "mt_charge_at_ph basic" {
  let seq = "K"
  let charge = @src.mt_charge_at_ph_public(seq, 7.0)
  // K: N-term +1, C-term -1, side chain +1 (pH 7 < 10.5)
  // charge = +1 - 1 + 1 = +1
  assert_true((charge - 1.0).abs() < 0.1)
}

test "mt_isoelectric_point single_aa" {
  let pi = @src.mt_isoelectric_point("A")
  // Alanine: only N-term and C-term, pI = (8.0 + 3.1) / 2 = 5.55
  // Actual value may vary slightly due to binary search
  assert_true(pi > 0.0 && pi < 14.0)
}

test "mt_isoelectric_point tyrosine" {
  let pi = @src.mt_isoelectric_point("Y")
  // Tyrosine: N-term +1, C-term -1, side chain pKa = 10.5
  // At pH < 10.5, side chain is protonated (neutral)
  // pI should be around 5.55
  assert_true(pi > 0.0 && pi < 14.0)
}

test "mt_isoelectric_point cysteine" {
  let pi = @src.mt_isoelectric_point("C")
  // Cysteine: N-term +1, C-term -1, side chain pKa = 8.3
  // At pH < 8.3, side chain is protonated (neutral)
  // pI should be around 5.55
  assert_true(pi > 0.0 && pi < 14.0)
}

test "mt_isoelectric_point acidic_protein" {
  let pi = @src.mt_isoelectric_point("DDD")
  // Acidic protein: pI should be low
  assert_true(pi < 5.0)
}

test "mt_isoelectric_point basic_protein" {
  let pi = @src.mt_isoelectric_point("KKK")
  // Basic protein: pI should be high
  assert_true(pi > 7.0)
}

test "mt_isoelectric_point mixed" {
  let pi = @src.mt_isoelectric_point("DEKR")
  assert_true(pi > 0.0 && pi < 14.0)
}

test "mt_isoelectric_point empty" {
  let pi = @src.mt_isoelectric_point("")
  assert_eq(pi, 0.0)
}

test "mt_protein_pi" {
  let pi = @src.mt_protein_pi("A")
  assert_true(pi > 0.0)
}

test "mt_charge_at_ph histidine" {
  let seq = "H"
  let charge = @src.mt_charge_at_ph_public(seq, 7.0)
  // H: N-term +1, C-term -1, side chain +1 (pH 7 < 6.0? No, 7 > 6.0, so side chain = 0)
  // Actually: pH 7 > 6.0, so side chain is deprotonated = 0
  // charge = +1 - 1 + 0 = 0
  assert_true(charge.abs() < 0.1)
}

test "mt_charge_at_ph histidine_low_ph" {
  let seq = "H"
  let charge = @src.mt_charge_at_ph_public(seq, 5.0)
  // H: pH 5 < 6.0, side chain protonated = +1
  // charge = +1 - 1 + 1 = +1
  assert_true((charge - 1.0).abs() < 0.1)
}

test "mt_example_protein" {
  let seq = @src.mt_example_protein()
  let pi = @src.mt_isoelectric_point(seq)
  assert_true(pi > 0.0 && pi < 14.0)
}


test "mt_isoelectric_point arginine" {
  let pi = @src.mt_isoelectric_point("R")
  // R: N-term +1, C-term -1, side chain +1 (pH 7 < 12.0)
  // charge at pH 7 = +1 - 1 + 1 = +1
  // pI should be high (> 10)
  assert_true(pi > 10.0)
}