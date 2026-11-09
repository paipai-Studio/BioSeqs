///|
/// Tests for seq_filter_catalog.mbt and seq_variant_library.mbt.

///|
fn seq_ft_seq(s : String) -> @src.Seq {
  @src.Seq::new(s)
}

///|
test "seq_filter_poly_a_and_clean" {
  let hits = @src.seq_apply_filter_catalog(seq_ft_seq("ACGTAAAAAAAAACGT"))
  let mut found_polya = false
  for h in hits {
    if h.rule_id == "poly_a_tail" {
      assert_eq(h.start, 4)
      found_polya = true
    }
  }
  assert_true(found_polya)
  // clean sequence: no hits at all
  let clean = @src.seq_apply_filter_catalog(seq_ft_seq("ACGTACGTACGT"))
  assert_eq(clean.length(), 0)
}

///|
test "seq_filter_adapter_contamination" {
  let hits = @src.seq_apply_filter_catalog(seq_ft_seq("AGATCGGAAGAGCACGT"))
  let mut found = false
  for h in hits {
    if h.rule_id == "illumina_adapter" {
      assert_eq(h.category, "contamination")
      assert_eq(h.severity, "fail")
      found = true
    }
  }
  assert_true(found)
}

///|
test "seq_filter_n_run_is_fatal" {
  assert_true(@src.seq_filter_has_fatal(seq_ft_seq("ACNNNNNACGT")))
  // plain sequence without N runs is not fatal
  assert_true(!@src.seq_filter_has_fatal(seq_ft_seq("ACGTACGTNN")))
}

///|
test "seq_filter_restriction_site" {
  let hits = @src.seq_apply_filter_catalog(seq_ft_seq("TTTGAATTCTTT"))
  let mut found_ecori = false
  for h in hits {
    if h.rule_id == "ecori_site" {
      assert_eq(h.start, 3)
      assert_eq(h.reason, "EcoRI restriction site GAATTC")
      found_ecori = true
    }
  }
  assert_true(found_ecori)
}

///|
test "seq_saturate_site_excludes_wildtype" {
  let variants = @src.seq_saturate_site(@src.Seq::new("ACGT"), 1)
  assert_eq(variants.length(), 3)
  for v in variants {
    assert_eq(v.to_string().length(), 4)
    assert_true(v.to_string() != "ACGT")
  }
}

///|
test "seq_degenerate_sequence_writes_scheme" {
  let deg = @src.seq_degenerate_sequence(
    @src.Seq::new("AAACCC"),
    [0],
    scheme="NNK",
  )
  assert_eq(deg.to_string(), "NNKCCC")
  // two codon sites
  let deg2 = @src.seq_degenerate_sequence(
    @src.Seq::new("AAACCC"),
    [0, 1],
    scheme="NDT",
  )
  assert_eq(deg2.to_string(), "NDTNDT")
}

///|
test "seq_degenerate_variant_count_and_expansion" {
  // NNK covers 4 * 4 * 2 = 32 codons per site
  assert_eq(@src.seq_degenerate_variant_count([0], scheme="NNK"), 32)
  assert_eq(@src.seq_degenerate_variant_count([0, 1], scheme="NNK"), 1024)
  // small scheme expands into concrete sequences
  let variants = @src.seq_degenerate_variants(
    @src.Seq::new("AAACCC"),
    [0],
    scheme="KK",
    max_variants=16,
  )
  // K = {G, T} on positions 0,1 -> 4 concrete variants
  assert_eq(variants.length(), 4)
  let mut seen = 0
  for v in variants {
    let s = v.to_string()
    assert_eq(s.length(), 6)
    assert_true(
      s == "GGACCC" || s == "GTACCC" || s == "TGACCC" || s == "TTACCC",
    )
    seen = seen + 1
  }
  assert_eq(seen, 4)
}

///|
test "seq_combo_library_cartesian" {
  let variants : Map[Int32, Array[String]] = Map([], capacity=2)
  variants[0] = ["C", "G"]
  variants[2] = ["T"]
  let combos = @src.seq_combo_library(@src.Seq::new("AAA"), variants)
  assert_eq(combos.length(), 2)
  let first = combos[0].to_string()
  let second = combos[1].to_string()
  assert_true(first == "CAT" || first == "GAT")
  assert_true(second == "CAT" || second == "GAT")
  assert_true(first != second)
}

///|
test "seq_reverse_translate_degenerate_codons" {
  // M is encoded by ATG only
  match @src.seq_reverse_translate_degenerate(@src.Seq::new("M")) {
    Some(dna) => assert_eq(dna, "ATG")
    None => fail("M should reverse-translate")
  }
  // K = AAA/AAG -> second position is R (A/G)
  match @src.seq_reverse_translate_degenerate(@src.Seq::new("K")) {
    Some(dna) => assert_eq(dna, "AAR")
    None => fail("K should reverse-translate")
  }
  // X has no codons
  assert_true(@src.seq_reverse_translate_degenerate(@src.Seq::new("X")) is None)
}
