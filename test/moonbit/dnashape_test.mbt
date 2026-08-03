///|
/// Tests for the dnashape module (Bioconductor DNAshapeR).
///
/// Covers DNA shape feature prediction, encoding/decoding, lookup tables,
/// and shape result statistics for MGW, Roll, ProT, HelT, and EP features.

// ===========================================================================
// 1. DnaShapeFeature enum
// ===========================================================================

///|
test "ds_feature_to_string" {
  assert_eq(@src.dnashape_make_feature("MGW").to_string(), "MGW")
  assert_eq(@src.dnashape_make_feature("Roll").to_string(), "Roll")
  assert_eq(@src.dnashape_make_feature("ProT").to_string(), "ProT")
  assert_eq(@src.dnashape_make_feature("HelT").to_string(), "HelT")
  assert_eq(@src.dnashape_make_feature("EP").to_string(), "EP")
}

///|
test "ds_feature_kmer_size" {
  let mgw = @src.dnashape_make_feature("MGW")
  let roll = @src.dnashape_make_feature("Roll")
  let prot = @src.dnashape_make_feature("ProT")
  let helt = @src.dnashape_make_feature("HelT")
  let ep = @src.dnashape_make_feature("EP")
  assert_eq(mgw.kmer_size(), 5)
  assert_eq(roll.kmer_size(), 4)
  assert_eq(prot.kmer_size(), 4)
  assert_eq(helt.kmer_size(), 5)
  assert_eq(ep.kmer_size(), 5)
}

///|
test "ds_make_feature" {
  assert_eq(@src.dnashape_make_feature("MGW").to_string(), "MGW")
  assert_eq(@src.dnashape_make_feature("Roll").to_string(), "Roll")
  assert_eq(@src.dnashape_make_feature("ProT").to_string(), "ProT")
  assert_eq(@src.dnashape_make_feature("HelT").to_string(), "HelT")
  assert_eq(@src.dnashape_make_feature("EP").to_string(), "EP")
  // Unknown name defaults to MGW
  assert_eq(@src.dnashape_make_feature("XYZ").to_string(), "MGW")
}

// ===========================================================================
// 2. Base and k-mer encoding
// ===========================================================================

///|
test "ds_encode_base" {
  assert_eq(@src.dnashape_encode_base("A"), 0)
  assert_eq(@src.dnashape_encode_base("C"), 1)
  assert_eq(@src.dnashape_encode_base("G"), 2)
  assert_eq(@src.dnashape_encode_base("T"), 3)
  assert_eq(@src.dnashape_encode_base("a"), 0)
  assert_eq(@src.dnashape_encode_base("c"), 1)
  assert_eq(@src.dnashape_encode_base("g"), 2)
  assert_eq(@src.dnashape_encode_base("t"), 3)
  assert_eq(@src.dnashape_encode_base("N"), -1)
}

///|
test "ds_encode_kmer" {
  // AA = 0*4+0 = 0
  assert_eq(@src.dnashape_encode_kmer("AA"), 0)
  // AC = 0*4+1 = 1
  assert_eq(@src.dnashape_encode_kmer("AC"), 1)
  // TT = 3*4+3 = 15
  assert_eq(@src.dnashape_encode_kmer("TT"), 15)
  // AAAA = 0
  assert_eq(@src.dnashape_encode_kmer("AAAA"), 0)
  // ACGT = 0*64+1*16+2*4+3 = 27
  assert_eq(@src.dnashape_encode_kmer("ACGT"), 27)
  // Contains N -> -1
  assert_eq(@src.dnashape_encode_kmer("AN"), -1)
}

///|
test "ds_decode_kmer" {
  assert_eq(@src.dnashape_decode_kmer(0, 2), "AA")
  assert_eq(@src.dnashape_decode_kmer(15, 2), "TT")
  assert_eq(@src.dnashape_decode_kmer(1, 2), "AC")
  assert_eq(@src.dnashape_decode_kmer(0, 4), "AAAA")
  assert_eq(@src.dnashape_decode_kmer(27, 4), "ACGT")
}

///|
test "ds_encode_decode_roundtrip" {
  let mut i = 0
  while i < 256 {
    let kmer = @src.dnashape_decode_kmer(i, 4)
    assert_eq(@src.dnashape_encode_kmer(kmer), i)
    i = i + 1
  }
}

///|
test "ds_revcomp" {
  assert_eq(@src.dnashape_revcomp("ATCG"), "CGAT")
  assert_eq(@src.dnashape_revcomp("AAAA"), "TTTT")
  assert_eq(@src.dnashape_revcomp("ACGT"), "ACGT")
  assert_eq(@src.dnashape_revcomp(""), "")
  assert_eq(@src.dnashape_revcomp("a"), "T")
}

// ===========================================================================
// 3. Lookup tables
// ===========================================================================

///|
test "ds_roll_table" {
  let t = @src.dnashape_get_roll_table()
  assert_eq(t.length(), 16)
  // AA (index 0) = -6.0
  let v0 = t.get(0).unwrap_or(999.0)
  assert_true((v0 - (-6.0)).abs() < 0.001)
  // TT (index 15) = -6.0
  let v15 = t.get(15).unwrap_or(999.0)
  assert_true((v15 - (-6.0)).abs() < 0.001)
  // CG (index 6) = 6.0
  let v6 = t.get(6).unwrap_or(999.0)
  assert_true((v6 - 6.0).abs() < 0.001)
}

///|
test "ds_prot_table" {
  let t = @src.dnashape_get_prot_table()
  assert_eq(t.length(), 16)
  // AA = 15.0
  let v0 = t.get(0).unwrap_or(0.0)
  assert_true((v0 - 15.0).abs() < 0.001)
  // TA = -2.0
  let v12 = t.get(12).unwrap_or(0.0)
  assert_true((v12 - (-2.0)).abs() < 0.001)
}

///|
test "ds_mgw_table" {
  let t = @src.dnashape_get_mgw_table()
  assert_eq(t.length(), 16)
  // AA = 3.5 (A-tract narrow groove)
  let v0 = t.get(0).unwrap_or(0.0)
  assert_true((v0 - 3.5).abs() < 0.001)
  // CG = 5.0
  let v6 = t.get(6).unwrap_or(0.0)
  assert_true((v6 - 5.0).abs() < 0.001)
}

///|
test "ds_helt_table" {
  let t = @src.dnashape_get_helt_table()
  assert_eq(t.length(), 16)
  // AA = 35.0
  let v0 = t.get(0).unwrap_or(0.0)
  assert_true((v0 - 35.0).abs() < 0.001)
  // AT = 31.0
  let v3 = t.get(3).unwrap_or(0.0)
  assert_true((v3 - 31.0).abs() < 0.001)
}

///|
test "ds_ep_table" {
  let t = @src.dnashape_get_ep_table()
  assert_eq(t.length(), 16)
  // AA = -1.5
  let v0 = t.get(0).unwrap_or(0.0)
  assert_true((v0 - (-1.5)).abs() < 0.001)
  // CG = -0.6
  let v6 = t.get(6).unwrap_or(0.0)
  assert_true((v6 - (-0.6)).abs() < 0.001)
}

// ===========================================================================
// 4. Shape prediction
// ===========================================================================

///|
test "ds_predict_mgw_short_seq" {
  // Sequence shorter than k-mer (5) returns all NaN values
  let mgw = @src.dnashape_make_feature("MGW")
  let result = @src.dnashape_predict("ACG", mgw)
  assert_eq(result.values.length(), 3)
  for v in result.values {
    assert_true(v.is_nan())
  }
}

///|
test "ds_predict_roll_basic" {
  // k=4, sequence length 8 produces some values
  let roll = @src.dnashape_make_feature("Roll")
  let result = @src.dnashape_predict("ATCGATCG", roll)
  assert_eq(result.values.length(), 8)
  // At least one value should be non-NaN
  let mut has_value = false
  for v in result.values {
    if !v.is_nan() {
      has_value = true
    }
  }
  assert_true(has_value)
}

///|
test "ds_predict_mgw_basic" {
  let mgw = @src.dnashape_make_feature("MGW")
  let seq = "ATCGATCGATCGATC"
  let result = @src.dnashape_predict(seq, mgw)
  assert_eq(result.values.length(), seq.length())
  // Mean should be in a reasonable range for MGW (3-6 Angstroms)
  let m = result.mean()
  assert_true(m > 0.0)
  assert_true(m < 10.0)
}

///|
test "ds_predict_helt_basic" {
  let helt = @src.dnashape_make_feature("HelT")
  let seq = "GCGCGCGCGCGC"
  let result = @src.dnashape_predict(seq, helt)
  assert_eq(result.values.length(), seq.length())
  // HelT should be near 34 degrees for GC-rich
  let m = result.mean()
  assert_true(m > 25.0 && m < 45.0)
}

///|
test "ds_predict_prot_basic" {
  let prot = @src.dnashape_make_feature("ProT")
  let seq = "ATCGATCGATCG"
  let result = @src.dnashape_predict(seq, prot)
  assert_eq(result.values.length(), seq.length())
  // ProT should be in range
  let m = result.mean()
  assert_true(m > 0.0)
}

///|
test "ds_predict_all_features" {
  let seq = "ATCGATCGATCGATCGATCG"
  let pred = @src.dnashape_predict_all(seq)
  assert_eq(pred.sequence, seq)
  assert_eq(pred.mgw.values.length(), seq.length())
  assert_eq(pred.roll.values.length(), seq.length())
  assert_eq(pred.prot.values.length(), seq.length())
  assert_eq(pred.helt.values.length(), seq.length())
}

///|
test "ds_predict_atract_narrow_mgw" {
  // A-tracts are known to have narrow minor groove
  let mgw_feat = @src.dnashape_make_feature("MGW")
  let atract = @src.dnashape_sample_atract()
  let gc_rich = @src.dnashape_sample_gc_rich()
  let mgw_atract = @src.dnashape_predict(atract, mgw_feat)
  let mgw_gc = @src.dnashape_predict(gc_rich, mgw_feat)
  // A-tract should have lower (narrower) MGW than GC-rich
  assert_true(mgw_atract.mean() < mgw_gc.mean())
}

// ===========================================================================
// 5. DnaShapeResult statistics
// ===========================================================================

///|
test "ds_result_mean_min_max" {
  let roll = @src.dnashape_make_feature("Roll")
  let result = @src.dnashape_predict("ATCGATCGATCGATCG", roll)
  let m = result.mean()
  let lo = result.min()
  let hi = result.max()
  assert_true(lo <= m)
  assert_true(m <= hi)
}

///|
test "ds_result_get_value" {
  let roll = @src.dnashape_make_feature("Roll")
  let result = @src.dnashape_predict("ATCGATCGATCGATCG", roll)
  // All values should be accessible
  let mut i = 0
  while i < result.values.length() {
    let _ = result.get_value(i)
    i = i + 1
  }
}

///|
test "ds_result_length" {
  let roll = @src.dnashape_make_feature("Roll")
  let result = @src.dnashape_predict("ATCGATCG", roll)
  assert_eq(result.length(), 8)
}

///|
test "ds_result_feature_name" {
  let roll = @src.dnashape_make_feature("Roll")
  let result = @src.dnashape_predict("ATCGATCG", roll)
  assert_eq(result.feature_name(), "Roll")
}

///|
test "ds_result_to_string" {
  let roll = @src.dnashape_make_feature("Roll")
  let result = @src.dnashape_predict("ATCGATCGATCG", roll)
  let s = @src.dnashape_result_to_string(result)
  assert_true(s.contains("Roll"))
  assert_true(s.contains("["))
}

// ===========================================================================
// 6. Composition and utility
// ===========================================================================

///|
test "ds_dinuc_composition" {
  let counts = @src.dnashape_dinuc_composition("ATCG")
  // AT, TC, CG => indices 3, 13, 6
  assert_eq(counts.length(), 16)
  assert_eq(counts[3], 1)
  assert_eq(counts[13], 1)
  assert_eq(counts[6], 1)
}

///|
test "ds_gc_content" {
  assert_true((@src.dnashape_gc_content("GCGC") - 1.0).abs() < 0.001)
  assert_true((@src.dnashape_gc_content("ATAT") - 0.0).abs() < 0.001)
  assert_true((@src.dnashape_gc_content("ATGC") - 0.5).abs() < 0.001)
  assert_eq(@src.dnashape_gc_content(""), 0.0)
}

///|
test "ds_sample_sequences" {
  let seq = @src.dnashape_sample_sequence()
  assert_true(seq.length() >= 49)
  let atract = @src.dnashape_sample_atract()
  assert_true(atract.length() >= 49)
  let gc = @src.dnashape_sample_gc_rich()
  assert_true(gc.length() >= 49)
  let mixed = @src.dnashape_sample_mixed()
  assert_true(mixed.length() >= 49)
}

///|
test "ds_heatmap" {
  let roll = @src.dnashape_make_feature("Roll")
  let result = @src.dnashape_predict("ATCGATCGATCGATCGATCG", roll)
  let hm = @src.dnashape_heatmap(result)
  // Heatmap string should have same length as sequence
  assert_eq(hm.length(), 20)
}
