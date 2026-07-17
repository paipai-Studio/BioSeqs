///|
test "codon_usage_count_codons" {
  let seq = "ATGGCCATG"
  let counts = @src.count_codons(seq)

  assert_eq(counts.get("ATG").unwrap_or(0), 2)
  assert_eq(counts.get("GCC").unwrap_or(0), 1)
}

///|
test "codon_usage_count_codons_empty" {
  let counts = @src.count_codons("")
  assert_eq(counts.length(), 0)
}

///|
test "codon_usage_count_codons_short" {
  let counts = @src.count_codons("AT")
  assert_eq(counts.length(), 0)
}

///|
test "codon_usage_rscu" {
  let counts = Map([], capacity=6)
  counts.set("GCA", 2)
  counts.set("GCC", 2)
  counts.set("GCG", 2)
  counts.set("GCT", 2)
  counts.set("AAA", 3)
  counts.set("AAG", 1)
  let rscu_map = @src.codon_rscu(counts)

  assert_true(rscu_map.get("GCA").unwrap_or(0.0) == 1.0)
  assert_true(rscu_map.get("GCC").unwrap_or(0.0) == 1.0)
  assert_true(rscu_map.get("AAA").unwrap_or(0.0) > 1.0)
  assert_true(rscu_map.get("AAG").unwrap_or(0.0) < 1.0)
}

///|
test "codon_usage_cai" {
  let seq = "ATGGCCATG"
  let optimal = Map([], capacity=2)
  optimal.set("ATG", 0.9)
  optimal.set("GCC", 0.8)
  let cai_val = @src.codon_cai(seq, optimal)

  assert_true(cai_val > 0.0 && cai_val <= 1.0)
}

///|
test "codon_usage_cai_empty" {
  let optimal = Map([], capacity=1)
  optimal.set("ATG", 1.5)
  let cai_val = @src.codon_cai("", optimal)

  assert_eq(cai_val, 0.0)
}

///|
test "codon_usage_enc" {
  let counts = Map([], capacity=8)
  counts.set("GCA", 1)
  counts.set("GCC", 1)
  counts.set("GCG", 1)
  counts.set("GCT", 1)
  counts.set("AAA", 1)
  counts.set("AAG", 1)
  counts.set("ATG", 1)
  counts.set("TGG", 1)
  let enc_val = @src.codon_enc(counts)

  assert_true(enc_val > 20.0 && enc_val <= 61.0)
}

///|
test "codon_usage_enc_high_bias" {
  let counts = Map([], capacity=8)
  counts.set("GCA", 10)
  counts.set("GCC", 0)
  counts.set("GCG", 0)
  counts.set("GCT", 0)
  counts.set("AAA", 10)
  counts.set("AAG", 0)
  counts.set("ATG", 1)
  counts.set("TGG", 1)
  let enc_val = @src.codon_enc(counts)

  assert_true(enc_val < 61.0)
}

///|
test "codon_usage_gc3" {
  let seq = "ATGCCCAAA"
  let gc3_val = @src.codon_gc3(seq)

  assert_eq(gc3_val, 2.0 / 3.0)
}

///|
test "codon_usage_gc3_empty" {
  let gc3_val = @src.codon_gc3("")

  assert_eq(gc3_val, 0.0)
}

///|
test "codon_usage_translate_codon" {
  assert_eq(@src.codon_translate("ATG"), "M")
  assert_eq(@src.codon_translate("TAA"), "*")
  assert_eq(@src.codon_translate("GCA"), "A")
  assert_eq(@src.codon_translate("XXX"), "X")
}

///|
test "codon_usage_reverse_complement_codon" {
  assert_eq(@src.codon_reverse_complement("ATG"), "CAT")
  assert_eq(@src.codon_reverse_complement("GCA"), "TGC")
}

///|
test "codon_usage_optimal_codons" {
  let rscu_map = Map([], capacity=6)
  rscu_map.set("GCA", 1.5)
  rscu_map.set("GCC", 1.0)
  rscu_map.set("GCG", 0.8)
  rscu_map.set("GCT", 0.7)
  rscu_map.set("AAA", 1.2)
  rscu_map.set("AAG", 0.8)
  let optimal = @src.codon_optimal(rscu_map)

  assert_true(optimal.contains("GCA"))
  assert_true(optimal.contains("AAA"))
}

///|
test "codon_usage_cbi" {
  let seq = "ATGGCCATG"
  let optimal = ["ATG", "GCC"]
  let cbi_val = @src.codon_cbi(seq, optimal)

  assert_true(cbi_val >= -1.0 && cbi_val <= 1.0)
}

///|
test "codon_usage_fop" {
  let seq = "ATGGCCATGGCA"
  let optimal = ["ATG", "GCC"]
  let fop_val = @src.codon_fop(seq, optimal)

  assert_true(fop_val >= 0.0 && fop_val <= 1.0)
}

///|
test "codon_usage_fop_all_optimal" {
  let seq = "ATGATG"
  let optimal = ["ATG"]
  let fop_val = @src.codon_fop(seq, optimal)

  assert_eq(fop_val, 1.0)
}

///|
test "codon_usage_fop_no_optimal" {
  let seq = "GCA"
  let optimal = ["ATG"]
  let fop_val = @src.codon_fop(seq, optimal)

  assert_eq(fop_val, 0.0)
}

///|
test "codon_usage_full_analysis" {
  let seq = "ATGGCCATGAAATGGTGA"
  let result = @src.codon_usage_analysis(seq)

  assert_true(result.codon_counts.length() > 0)
  assert_true(result.rscu.length() > 0)
  assert_true(result.cai > 0.0)
  assert_true(result.enc > 20.0)
  assert_true(result.gc3 >= 0.0)
  assert_true(result.optimal_codons.length() > 0)
}

///|
test "codon_usage_with_reference" {
  let seq = "ATGGCCATG"
  let ref_seq = "ATGATGATGGCCGCCGCC"
  let result = @src.codon_usage_with_ref(seq, ref_seq)

  assert_true(result.cai > 0.0)
}
