///|
test "universalmotif_motif_from_pwm" {
  let pwm = [
    [0.9, 0.1, 0.0, 0.0],
    [0.0, 0.0, 0.9, 0.1],
    [0.0, 0.0, 0.9, 0.1],
    [0.1, 0.9, 0.0, 0.0],
  ]
  
  let motif = @src.S4Motif::from_pwm("TATA-box", "DNA", pwm)
  
  assert_eq(motif.name, "TATA-box")
  assert_eq(motif.alphabet, "DNA")
  assert_eq(motif.pwm.length(), 4)
}

///|
test "universalmotif_calculate_consensus" {
  let pwm = [
    [0.9, 0.1, 0.0, 0.0],
    [0.0, 0.0, 0.9, 0.1],
    [0.1, 0.9, 0.0, 0.0],
  ]
  
  let motif = @src.S4Motif::from_pwm("test", "DNA", pwm)
  
  assert_eq(motif.consensus.length(), 3)
}

///|
test "universalmotif_create_example_motif" {
  let motif = @src.create_example_s4motif()
  
  assert_eq(motif.name, "TATA-box")
  assert_eq(motif.pwm.length(), 5)
}