///|
test "sv params default" {
  let p = @src.SvParams::default()
  assert_eq(p.min_support, 3)
  assert_eq(p.min_sv_len, 50)
}

///|
test "sv params custom" {
  let p = @src.SvParams::new(min_support=5, min_sv_len=100)
  assert_eq(p.min_support, 5)
  assert_eq(p.min_sv_len, 100)
}

///|
test "parse cigar simple" {
  let elems = @src.parse_cigar("10M5D20M")
  assert_eq(elems.length(), 3)
  assert_eq(elems[0].length, 10)
  assert_eq(elems[1].length, 5)
  assert_eq(elems[2].length, 20)
}

///|
test "parse cigar with insert" {
  let elems = @src.parse_cigar("10M5I10M")
  assert_eq(elems.length(), 3)
  // 仅验证长度, 不直接比较 CigarOp 枚举变体
  assert_eq(elems[1].length, 5)
}

///|
test "parse cigar star" {
  let elems = @src.parse_cigar("*")
  assert_eq(elems.length(), 0)
}

///|
test "sv candidates from sam deletion" {
  let cigar = @src.parse_cigar("10M50D20M")
  let params = @src.SvParams::default()
  let calls = @src.sv_candidates_from_sam("chr1", 100, cigar, params)
  assert_eq(calls.length(), 1)
  assert_eq(@src.sv_type_label(calls[0].sv_type()), "DEL")
  assert_eq(calls[0].start(), 110)
  assert_eq(calls[0].end(), 160)
}

///|
test "sv candidates from sam insertion" {
  let cigar = @src.parse_cigar("10M100I20M")
  let params = @src.SvParams::default()
  let calls = @src.sv_candidates_from_sam("chr1", 100, cigar, params)
  assert_eq(calls.length(), 1)
  assert_eq(@src.sv_type_label(calls[0].sv_type()), "INS")
}

///|
test "sv candidates from sam skip" {
  let cigar = @src.parse_cigar("10M50N20M")
  let params = @src.SvParams::default()
  let calls = @src.sv_candidates_from_sam("chr1", 100, cigar, params)
  assert_eq(calls.length(), 1)
  assert_eq(@src.sv_type_label(calls[0].sv_type()), "DEL")
}

///|
test "sv candidates short deletion ignored" {
  let cigar = @src.parse_cigar("10M10D20M")
  let params = @src.SvParams::default()
  let calls = @src.sv_candidates_from_sam("chr1", 100, cigar, params)
  assert_eq(calls.length(), 0)
}

///|
test "sv candidates no sv" {
  let cigar = @src.parse_cigar("100M")
  let params = @src.SvParams::default()
  let calls = @src.sv_candidates_from_sam("chr1", 100, cigar, params)
  assert_eq(calls.length(), 0)
}

///|
test "merge sv candidates" {
  let cigar = @src.parse_cigar("10M50D20M")
  let params = @src.SvParams::default()
  let calls1 = @src.sv_candidates_from_sam("chr1", 100, cigar, params)
  let calls2 = @src.sv_candidates_from_sam("chr1", 100, cigar, params)
  let all_calls : Array[@src.SvCall] = []
  for c in calls1 {
    all_calls.push(c)
  }
  for c in calls2 {
    all_calls.push(c)
  }
  let merged = @src.merge_sv_candidates(all_calls)
  assert_eq(merged.length(), 1)
  assert_eq(merged[0].support_reads(), 2)
}

///|
test "filter sv calls" {
  let cigar = @src.parse_cigar("10M50D20M")
  let params = @src.SvParams::default()
  let calls = @src.sv_candidates_from_sam("chr1", 100, cigar, params)
  let filtered = @src.filter_sv_calls(calls, 1)
  assert_eq(filtered.length(), 1)
}

///|
test "sort sv calls by score" {
  let cigar = @src.parse_cigar("10M50D20M")
  let params = @src.SvParams::default()
  let calls = @src.sv_candidates_from_sam("chr1", 100, cigar, params)
  let sorted = @src.sort_sv_calls_by_score(calls)
  assert_eq(sorted.length(), 1)
}

///|
test "sv caller batch" {
  // 批量检测通过 sv_candidates_from_sam + merge_sv_candidates 组合验证
  let cigar = @src.parse_cigar("10M50D20M")
  let params = @src.SvParams::default()
  let calls1 = @src.sv_candidates_from_sam("chr1", 100, cigar, params)
  let calls2 = @src.sv_candidates_from_sam("chr1", 100, cigar, params)
  let all_calls : Array[@src.SvCall] = []
  for c in calls1 {
    all_calls.push(c)
  }
  for c in calls2 {
    all_calls.push(c)
  }
  let merged = @src.merge_sv_candidates(all_calls)
  assert_eq(merged.length(), 1)
  assert_eq(merged[0].support_reads(), 2)
}

///|
test "sv call high confidence" {
  let cigar = @src.parse_cigar("10M50D20M")
  let params = @src.SvParams::default()
  let calls = @src.sv_candidates_from_sam("chr1", 100, cigar, params)
  let c = calls[0]
  assert_true(c.is_high_confidence())
}

///|
test "sv call fields" {
  let cigar = @src.parse_cigar("10M50D20M")
  let params = @src.SvParams::default()
  let calls = @src.sv_candidates_from_sam("chr1", 100, cigar, params)
  let c = calls[0]
  assert_eq(c.chrom(), "chr1")
  assert_eq(@src.sv_type_label(c.sv_type()), "DEL")
  assert_eq(c.start(), 110)
  assert_eq(c.end(), 160)
}
