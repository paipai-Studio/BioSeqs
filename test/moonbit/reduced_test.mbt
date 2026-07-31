///|
/// Tests for reduced module (reduced amino acid alphabets).

///|
test "rad_alphabet_creation" {
  let rad = @src.rad_alphabet()
  assert_eq(rad.get_name(), "RAD")
  assert_eq(rad.get_n_groups(), 6)
}

///|
test "rad_alphabet_reduce_a" {
  let rad = @src.rad_alphabet()
  let g = rad.reduce_aa("A")
  assert_eq(g, "0") // (A,G) group
}

///|
test "rad_alphabet_reduce_c" {
  let rad = @src.rad_alphabet()
  let g = rad.reduce_aa("C")
  assert_eq(g, "1") // (C) group
}

///|
test "rad_alphabet_reduce_d" {
  let rad = @src.rad_alphabet()
  let g = rad.reduce_aa("D")
  assert_eq(g, "2") // (D,E,N,Q) group
}

///|
test "rad_alphabet_reduce_i" {
  let rad = @src.rad_alphabet()
  let g = rad.reduce_aa("I")
  assert_eq(g, "3") // (I,L,M,V) group
}

///|
test "rad_alphabet_reduce_f" {
  let rad = @src.rad_alphabet()
  let g = rad.reduce_aa("F")
  assert_eq(g, "4") // (F,Y,W) group
}

///|
test "rad_alphabet_reduce_h" {
  let rad = @src.rad_alphabet()
  let g = rad.reduce_aa("H")
  assert_eq(g, "5") // (H,K,R,S,T) group
}

///|
test "rad_alphabet_reduce_unknown" {
  let rad = @src.rad_alphabet()
  let g = rad.reduce_aa("Z")
  assert_eq(g, "X") // Unknown
}

///|
test "rad_alphabet_reduce_sequence" {
  let rad = @src.rad_alphabet()
  let reduced = rad.reduce_sequence("ACDFH")
  assert_eq(reduced, "01245")
}

///|
test "dayhoff_alphabet_creation" {
  let dh = @src.dayhoff_alphabet()
  assert_eq(dh.get_name(), "Dayhoff")
  assert_eq(dh.get_n_groups(), 6)
}

///|
test "dayhoff_alphabet_reduce_a" {
  let dh = @src.dayhoff_alphabet()
  let g = dh.reduce_aa("A")
  assert_eq(g, "0") // (A,G,P,S,T) group
}

///|
test "dayhoff_alphabet_reduce_p" {
  let dh = @src.dayhoff_alphabet()
  let g = dh.reduce_aa("P")
  assert_eq(g, "0") // (A,G,P,S,T) group
}

///|
test "dayhoff_alphabet_reduce_c" {
  let dh = @src.dayhoff_alphabet()
  let g = dh.reduce_aa("C")
  assert_eq(g, "1") // (C) group
}

///|
test "dayhoff_alphabet_reduce_d" {
  let dh = @src.dayhoff_alphabet()
  let g = dh.reduce_aa("D")
  assert_eq(g, "2") // (D,E,N,Q) group
}

///|
test "dayhoff_alphabet_reduce_i" {
  let dh = @src.dayhoff_alphabet()
  let g = dh.reduce_aa("I")
  assert_eq(g, "3") // (I,L,M,V) group
}

///|
test "dayhoff_alphabet_reduce_f" {
  let dh = @src.dayhoff_alphabet()
  let g = dh.reduce_aa("F")
  assert_eq(g, "4") // (F,W,Y) group
}

///|
test "dayhoff_alphabet_reduce_h" {
  let dh = @src.dayhoff_alphabet()
  let g = dh.reduce_aa("H")
  assert_eq(g, "5") // (H,K,R) group
}

///|
test "dayhoff_alphabet_reduce_sequence" {
  let dh = @src.dayhoff_alphabet()
  let reduced = dh.reduce_sequence("APDCIFH")
  assert_eq(reduced, "0021345")
}

///|
test "charm_alphabet_creation" {
  let charm = @src.charm_alphabet()
  assert_eq(charm.get_name(), "CHARM")
  assert_eq(charm.get_n_groups(), 4)
}

///|
test "charm_alphabet_reduce_a" {
  let charm = @src.charm_alphabet()
  let g = charm.reduce_aa("A")
  assert_eq(g, "0") // (A,C,I,L,M,F,V,W,Y) group
}

///|
test "charm_alphabet_reduce_s" {
  let charm = @src.charm_alphabet()
  let g = charm.reduce_aa("S")
  assert_eq(g, "1") // (S,T,N,Q,G,P) group
}

///|
test "charm_alphabet_reduce_d" {
  let charm = @src.charm_alphabet()
  let g = charm.reduce_aa("D")
  assert_eq(g, "2") // (D,E) group
}

///|
test "charm_alphabet_reduce_h" {
  let charm = @src.charm_alphabet()
  let g = charm.reduce_aa("H")
  assert_eq(g, "3") // (H,K,R) group
}

///|
test "sdm12_alphabet_creation" {
  let sdm = @src.sdm12_alphabet()
  assert_eq(sdm.get_name(), "SDM12")
  assert_eq(sdm.get_n_groups(), 12)
}

///|
test "sdm12_alphabet_reduce_a" {
  let sdm = @src.sdm12_alphabet()
  let g = sdm.reduce_aa("A")
  assert_eq(g, "0")
}

///|
test "sdm12_alphabet_reduce_r" {
  let sdm = @src.sdm12_alphabet()
  let g = sdm.reduce_aa("R")
  assert_eq(g, "1")
}

///|
test "sdm12_alphabet_reduce_n" {
  let sdm = @src.sdm12_alphabet()
  let g = sdm.reduce_aa("N")
  assert_eq(g, "2")
}

///|
test "sdm12_alphabet_reduce_d" {
  let sdm = @src.sdm12_alphabet()
  let g = sdm.reduce_aa("D")
  assert_eq(g, "3")
}

///|
test "sdm12_alphabet_reduce_c" {
  let sdm = @src.sdm12_alphabet()
  let g = sdm.reduce_aa("C")
  assert_eq(g, "4")
}

///|
test "sdm12_alphabet_reduce_e" {
  let sdm = @src.sdm12_alphabet()
  let g = sdm.reduce_aa("E")
  assert_eq(g, "6")
}

///|
test "sdm12_alphabet_reduce_g" {
  let sdm = @src.sdm12_alphabet()
  let g = sdm.reduce_aa("G")
  assert_eq(g, "7")
}

///|
test "sdm12_alphabet_reduce_h" {
  let sdm = @src.sdm12_alphabet()
  let g = sdm.reduce_aa("H")
  assert_eq(g, "8")
}

///|
test "sdm12_alphabet_reduce_i" {
  let sdm = @src.sdm12_alphabet()
  let g = sdm.reduce_aa("I")
  assert_eq(g, "9")
}

///|
test "sdm12_alphabet_reduce_l" {
  let sdm = @src.sdm12_alphabet()
  let g = sdm.reduce_aa("L")
  assert_eq(g, "10")
}

///|
test "sdm12_alphabet_reduce_k" {
  let sdm = @src.sdm12_alphabet()
  let g = sdm.reduce_aa("K")
  assert_eq(g, "11")
}

///|
test "all_reduced_alphabets" {
  let alphabets = @src.all_reduced_alphabets()
  assert_eq(alphabets.length(), 4)
}

///|
test "get_reduced_alphabet_by_name" {
  match @src.get_reduced_alphabet("RAD") {
    Some(a) => assert_eq(a.get_name(), "RAD")
    None => assert_true(false)
  }
  match @src.get_reduced_alphabet("dayhoff") {
    Some(a) => assert_eq(a.get_name(), "Dayhoff")
    None => assert_true(false)
  }
  match @src.get_reduced_alphabet("unknown") {
    Some(_) => assert_true(false)
    None => assert_true(true)
  }
}

///|
test "reduced_identity_same" {
  let rad = @src.rad_alphabet()
  let identity = @src.reduced_identity("ACGT", "ACGT", rad)
  assert_true(identity > 0.99)
}

///|
test "reduced_identity_different" {
  let rad = @src.rad_alphabet()
  let identity = @src.reduced_identity("AAAA", "CCCC", rad)
  assert_true(identity < 0.5)
}

///|
test "reduced_identity_empty" {
  let rad = @src.rad_alphabet()
  let identity = @src.reduced_identity("", "", rad)
  assert_eq(identity, 0.0)
}

///|
test "reduced_identity_same_group" {
  let rad = @src.rad_alphabet()
  // A and G are in the same RAD group (0)
  let identity = @src.reduced_identity("AG", "GA", rad)
  assert_true(identity > 0.99)
}

///|
test "reduced_alphabet_summary" {
  let rad = @src.rad_alphabet()
  let s = rad.summary()
  assert_true(s.length() > 0)
}

///|
test "reduced_sample_sequence" {
  let seq = @src.reduced_sample_sequence()
  assert_true(seq.length() > 0)
}

///|
test "reduced_alphabet_get_groups" {
  let rad = @src.rad_alphabet()
  let groups = rad.get_groups()
  assert_eq(groups.length(), 6)
}
