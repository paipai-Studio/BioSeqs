///|
test "kf_lgamma known values" {
  // lgamma(1) = 0
  assert_true(@src.kf_lgamma(1.0).abs() < 1.0e-9)
  // lgamma(0.5) = ln(sqrt(pi)) = 0.5*ln(pi)
  let expected = 0.5 * @math.ln(@math.PI)
  assert_true((@src.kf_lgamma(0.5) - expected).abs() < 1.0e-9)
  // lgamma(5) = ln(24) = ln(4!)
  assert_true((@src.kf_lgamma(5.0) - @math.ln(24.0)).abs() < 1.0e-9)
}

///|
test "kf_erfc known values" {
  // erfc(0) = 1
  assert_true((@src.kf_erfc(0.0) - 1.0).abs() < 1.0e-9)
  // erfc(inf) = 0, erfc(large) ≈ 0
  assert_true(@src.kf_erfc(10.0) < 1.0e-40)
  // erfc(-x) = 2 - erfc(x)
  let e5 = @src.kf_erfc(5.0)
  let em5 = @src.kf_erfc(-5.0)
  assert_true((em5 - (2.0 - e5)).abs() < 1.0e-9)
  // erfc(1) ≈ 0.157299207050285
  assert_true((@src.kf_erfc(1.0) - 0.1572992070502851).abs() < 1.0e-9)
}

///|
test "kf_gammap gammaq complement" {
  // gammap(s,z) + gammaq(s,z) = 1
  for s in [1.0, 2.0, 3.5, 0.5, 10.0] {
    for z in [0.5, 1.0, 2.0, 5.0, 10.0] {
      let p = @src.kf_gammap(s, z)
      let q = @src.kf_gammaq(s, z)
      assert_true((p + q - 1.0).abs() < 1.0e-9)
    }
  }
}

///|
test "kf_gammap known values" {
  // gammap(1, z) = 1 - exp(-z)
  for z in [0.5, 1.0, 2.0, 5.0] {
    let expected = 1.0 - @math.exp(-z)
    assert_true((@src.kf_gammap(1.0, z) - expected).abs() < 1.0e-9)
  }
  // gammap(0.5, 1) ≈ 0.8427 (erf(1) = 1 - erfc(1))
  assert_true((@src.kf_gammap(0.5, 1.0) - 0.8427007929497149).abs() < 1.0e-9)
}

///|
test "kf_gammaq known values" {
  // gammaq(1, z) = exp(-z)
  for z in [0.5, 1.0, 2.0, 5.0] {
    let expected = @math.exp(-z)
    assert_true((@src.kf_gammaq(1.0, z) - expected).abs() < 1.0e-9)
  }
}

///|
test "kf_betai known values" {
  // betai(1, 1, x) = x (uniform distribution CDF)
  for x in [0.1, 0.25, 0.5, 0.75, 0.9] {
    assert_true((@src.kf_betai(1.0, 1.0, x) - x).abs() < 1.0e-9)
  }
  // betai(a, b, 0) = 0
  assert_true(@src.kf_betai(2.0, 3.0, 0.0).abs() < 1.0e-9)
  // betai(a, b, 1) = 1
  assert_true((@src.kf_betai(2.0, 3.0, 1.0) - 1.0).abs() < 1.0e-9)
  // betai(2, 2, 0.5) = 0.5 (symmetric beta)
  assert_true((@src.kf_betai(2.0, 2.0, 0.5) - 0.5).abs() < 1.0e-9)
}

///|
test "hts_parse_region_str full" {
  let r = @src.hts_parse_region_str("chr1:100-200")
  assert_true(r.is_some())
  let rr = r.unwrap()
  assert_eq(rr.chrom, "chr1")
  assert_eq(rr.beg, 99) // 1-based → 0-based
  assert_eq(rr.end, 200)
}

///|
test "hts_parse_region_str single_pos" {
  let r = @src.hts_parse_region_str("chr2:500")
  let rr = r.unwrap()
  assert_eq(rr.chrom, "chr2")
  assert_eq(rr.beg, 499)
  assert_eq(rr.end, 500)
}

///|
test "hts_parse_region_str chrom_only" {
  let r = @src.hts_parse_region_str("chr3")
  let rr = r.unwrap()
  assert_eq(rr.chrom, "chr3")
  assert_eq(rr.beg, 0)
  assert_eq(rr.end, 0x7FFFFFFF)
}

///|
test "hts_parse_region_list multi" {
  let regions = @src.hts_parse_region_list("chr1:100-200,chr2:300-400,chr3")
  assert_eq(regions.length(), 3)
  assert_eq(regions[0].chrom, "chr1")
  assert_eq(regions[0].beg, 99)
  assert_eq(regions[0].end, 200)
  assert_eq(regions[1].chrom, "chr2")
  assert_eq(regions[2].chrom, "chr3")
}

///|
test "hts_regidx basic" {
  let idx = @src.HtsRegIdx::new()
  idx.insert("chr1", 100, 200)
  idx.insert("chr1", 300, 400)
  idx.insert("chr2", 0, 1000)
  idx.compact_all()
  assert_eq(idx.nregs(), 3)
  assert_eq(idx.seq_nregs("chr1"), 2)
  assert_eq(idx.seq_nregs("chr2"), 1)
  assert_eq(idx.seq_nregs("chr3"), 0)
}

///|
test "hts_regidx overlap" {
  let idx = @src.HtsRegIdx::new()
  idx.insert("chr1", 100, 200)
  idx.insert("chr1", 300, 400)
  idx.compact_all()
  // Overlapping query
  assert_true(idx.overlap("chr1", 150, 160))
  assert_true(idx.overlap("chr1", 190, 210))
  assert_true(idx.overlap("chr1", 350, 360))
  // Non-overlapping
  assert_true(not(idx.overlap("chr1", 200, 300)))
  assert_true(not(idx.overlap("chr1", 250, 290)))
  // Different chrom
  assert_true(not(idx.overlap("chr2", 150, 160)))
}

///|
test "hts_regidx merge" {
  let idx = @src.HtsRegIdx::new()
  idx.insert("chr1", 100, 200)
  idx.insert("chr1", 150, 250) // overlaps with first
  idx.insert("chr1", 300, 400)
  idx.insert("chr1", 350, 450) // overlaps with third
  idx.compact_all()
  // After merge: [100-250, 300-450]
  assert_eq(idx.seq_nregs("chr1"), 2)
  let overlaps = idx.overlapping("chr1", 0, 1000)
  assert_eq(overlaps.length(), 2)
  let (o0b, o0e) = overlaps[0]
  assert_eq(o0b, 100)
  assert_eq(o0e, 250)
  let (o1b, o1e) = overlaps[1]
  assert_eq(o1b, 300)
  assert_eq(o1e, 450)
}

///|
test "hts_regidx from_string" {
  let idx = @src.hts_regidx_from_string("chr1:100-200,chr2:300-400")
  assert_eq(idx.nregs(), 2)
  assert_true(idx.overlap("chr1", 150, 160))
  assert_true(idx.overlap("chr2", 350, 360))
  assert_true(not(idx.overlap("chr1", 250, 290)))
}

///|
test "hts_regidx from_bed" {
  let bed = "chr1\t100\t200\nchr1\t300\t400\nchr2\t0\t1000\n"
  let idx = @src.hts_regidx_from_bed(bed)
  assert_eq(idx.nregs(), 3)
  assert_true(idx.overlap("chr1", 150, 160))
  assert_true(idx.overlap("chr2", 500, 600))
}

///|
test "hts_regidx seq_names" {
  let idx = @src.hts_regidx_from_string("chr1:100-200,chr2:300-400,chr3:500-600")
  let names = idx.seq_names()
  assert_eq(names.length(), 3)
}
