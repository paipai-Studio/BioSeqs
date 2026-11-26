///|
/// Tests for biotite_sequence_rules.mbt.
///
/// Synthetic operon on the forward strand:
///   promoter "TATAAT" at 0
///   RBS      "AGGAGGT" at 12
///   start    "ATG"     at 25 (gap 6 after RBS end, within [5,13])
///   31 lysine codons "AAG" (93 bp, 31 aa >= min 30)
///   stop     "TAA"     at 121
///   terminator "AAAAAA" at 124
fn bsr_forward_operon() -> String {
  let mut seq = "TATAAT" + "GCGCGC" + "AGGAGGT" + "AAAAAA" + "ATG"
  for _i = 0; _i < 31; _i = _i + 1 {
    seq = seq + "AAG"
  }
  seq + "TAA" + "AAAAAA"
}

///|
test "find genes forward strand" {
  let seq = bsr_forward_operon()
  let genes = @src.biotite_find_genes(
    seq,
    "TATAAT",
    "AAAAAA",
    "AGGAGGT",
    ["ATG"],
    ["TAA", "TAG", "TGA"],
  )
  assert_eq(genes.length(), 1)
  let g = genes[0]
  assert_eq(g.direction, 1)
  assert_eq(g.start, 25)
  assert_eq(g.stop, 124)
  assert_eq(g.codon_start, 28)
  assert_false(g.stop_far_before_start)
}

///|
test "find genes reverse strand" {
  let seq = @src.biostrings_reverse_complement(bsr_forward_operon())
  let genes = @src.biotite_find_genes(
    seq,
    "TATAAT",
    "AAAAAA",
    "AGGAGGT",
    ["ATG"],
    ["TAA", "TAG", "TGA"],
  )
  assert_eq(genes.length(), 1)
  let g = genes[0]
  assert_eq(g.direction, -1)
  let n = 130
  // the gene maps back to [n - 124, n - 25) = [6, 105)
  assert_eq(g.start, n - 124)
  assert_eq(g.stop, n - 25)
  assert_eq(g.codon_start, n - 28)
  assert_false(g.stop_far_before_start)
}

///|
test "find genes min protein length filter" {
  let seq = bsr_forward_operon()
  let genes = @src.biotite_find_genes(
    seq,
    "TATAAT",
    "AAAAAA",
    "AGGAGGT",
    ["ATG"],
    ["TAA", "TAG", "TGA"],
    min_protein_len=40,
  )
  assert_eq(genes.length(), 0)
}

///|
test "find genes max protein length filter" {
  let seq = bsr_forward_operon()
  let genes = @src.biotite_find_genes(
    seq,
    "TATAAT",
    "AAAAAA",
    "AGGAGGT",
    ["ATG"],
    ["TAA", "TAG", "TGA"],
    max_protein_len=Some(20),
  )
  assert_eq(genes.length(), 0)
}

///|
test "find genes requires terminator" {
  // operon without the downstream terminator
  let mut seq = "TATAAT" + "GCGCGC" + "AGGAGGT" + "AAAAAA" + "ATG"
  for _i = 0; _i < 31; _i = _i + 1 {
    seq = seq + "AAG"
  }
  seq = seq + "TAA" + "CGCGCG"
  let genes = @src.biotite_find_genes(
    seq,
    "TATAAT",
    "AAAAAA",
    "AGGAGGT",
    ["ATG"],
    ["TAA", "TAG", "TGA"],
  )
  assert_eq(genes.length(), 0)
}

///|
test "find genes no promoter no genes" {
  let seq = bsr_forward_operon()
  let genes = @src.biotite_find_genes(
    seq,
    "GGGGGG",
    "AAAAAA",
    "AGGAGGT",
    ["ATG"],
    ["TAA", "TAG", "TGA"],
  )
  assert_eq(genes.length(), 0)
}

///|
test "find genes incomplete without stop codon" {
  // construct a candidate that never hits an in-frame stop: the gene
  // runs to the sequence end and is flagged stop_far_before_start
  let mut seq = "TATAAT" + "GCGCGC" + "AGGAGGT" + "AAAAAA" + "ATG"
  for _i = 0; _i < 31; _i = _i + 1 {
    seq = seq + "AAG"
  }
  // no stop codon, but terminator present at the very end
  seq = seq + "AAAAAA"
  let genes = @src.biotite_find_genes(
    seq,
    "TATAAT",
    "AAAAAA",
    "AGGAGGT",
    ["ATG"],
    ["TAA", "TAG", "TGA"],
  )
  assert_eq(genes.length(), 1)
  assert_true(genes[0].stop_far_before_start)
}

///|
test "find genes alternate start codon" {
  let mut seq = "TATAAT" + "GCGCGC" + "AGGAGGT" + "AAAAAA" + "GTG"
  for _i = 0; _i < 31; _i = _i + 1 {
    seq = seq + "AAG"
  }
  seq = seq + "TGA" + "AAAAAA"
  let genes = @src.biotite_find_genes(
    seq,
    "TATAAT",
    "AAAAAA",
    "AGGAGGT",
    ["ATG", "GTG", "TTG"],
    ["TAA", "TAG", "TGA"],
  )
  assert_eq(genes.length(), 1)
  assert_eq(genes[0].direction, 1)
}

///|
test "find genes gap filter" {
  // RBS-to-start gap of 20 exceeds rbs_max_gap=13
  let mut seq = "TATAAT" + "GCGCGC" + "AGGAGGT" + "AAAAAAAAAAAAAAAAAAAA" + "ATG"
  for _i = 0; _i < 31; _i = _i + 1 {
    seq = seq + "AAG"
  }
  seq = seq + "TAA" + "AAAAAA"
  let genes = @src.biotite_find_genes(
    seq,
    "TATAAT",
    "AAAAAA",
    "AGGAGGT",
    ["ATG"],
    ["TAA", "TAG", "TGA"],
  )
  assert_eq(genes.length(), 0)
}

///|
test "find genes terminator must directly follow stop codon" {
  // Terminator appears 10 bp downstream of the stop codon instead of
  // directly after it -> no gene is reported.
  let mut seq = "TATAAT" + "GCGCGC" + "AGGAGGT" + "AAAAAA" + "ATG"
  for _i = 0; _i < 31; _i = _i + 1 {
    seq = seq + "AAG"
  }
  seq = seq + "TAA" + "CGCGCGCGC" + "AAAAAA"
  let genes = @src.biotite_find_genes(
    seq,
    "TATAAT",
    "AAAAAA",
    "AGGAGGT",
    ["ATG"],
    ["TAA", "TAG", "TGA"],
  )
  assert_eq(genes.length(), 0)
}
