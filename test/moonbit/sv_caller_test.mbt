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
  // 3 条同位置读数达到默认 min_support=3, 应报告 1 个 DEL
  let records : Array[(String, Int32, Array[@src.CigarElement])] = [
    ("chr1", 100, @src.parse_cigar("10M50D20M")),
    ("chr1", 100, @src.parse_cigar("10M50D20M")),
    ("chr1", 100, @src.parse_cigar("10M50D20M")),
  ]
  let params = @src.SvParams::default()
  let calls = @src.sv_caller_batch(records, params)
  assert_eq(calls.length(), 1)
  assert_eq(calls[0].support_reads(), 3)
  assert_eq(@src.sv_type_label(calls[0].sv_type()), "DEL")
}

///|
test "sv caller batch filters low support" {
  // 仅 2 条读数低于默认 min_support=3, 应被过滤为空
  let records : Array[(String, Int32, Array[@src.CigarElement])] = [
    ("chr1", 100, @src.parse_cigar("10M50D20M")),
    ("chr1", 100, @src.parse_cigar("10M50D20M")),
  ]
  let params = @src.SvParams::default()
  let calls = @src.sv_caller_batch(records, params)
  assert_eq(calls.length(), 0)
}

///|
test "sv caller batch with custom min support" {
  let records : Array[(String, Int32, Array[@src.CigarElement])] = [
    ("chr1", 100, @src.parse_cigar("10M50D20M")),
    ("chr1", 100, @src.parse_cigar("10M50D20M")),
  ]
  let params = @src.SvParams::new(min_support=1)
  let calls = @src.sv_caller_batch(records, params)
  assert_eq(calls.length(), 1)
  assert_eq(calls[0].support_reads(), 2)
}

///|
test "sv caller batch with depth" {
  let records : Array[(String, Int32, Array[@src.CigarElement])] = [
    ("chr1", 100, @src.parse_cigar("10M50D20M")),
    ("chr1", 100, @src.parse_cigar("10M50D20M")),
    ("chr1", 100, @src.parse_cigar("10M50D20M")),
  ]
  let params : @src.SvParams = @src.SvParams::default()
  let depth_by_chrom : Map[String, Int32] = Map([("chr1", 30)])
  let calls = @src.sv_caller_batch(records, params, depth_by_chrom~)
  assert_eq(calls.length(), 1)
  assert_eq(calls[0].total_reads(), 30)
  assert_true(calls[0].score() == 0.1)
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

///|
test "merge sv candidates with depth" {
  let cigar = @src.parse_cigar("10M50D20M")
  let params = @src.SvParams::default()
  let calls = @src.sv_candidates_from_sam("chr1", 100, cigar, params)
  let all_calls : Array[@src.SvCall] = []
  for _ in 0..<3 {
    for c in calls {
      all_calls.push(c)
    }
  }
  let depth_by_chrom : Map[String, Int32] = Map([("chr1", 12)])
  let merged = @src.merge_sv_candidates(all_calls, depth_by_chrom~)
  assert_eq(merged.length(), 1)
  assert_eq(merged[0].support_reads(), 3)
  assert_eq(merged[0].total_reads(), 12)
  assert_true(merged[0].score() == 0.25)
  assert_true(!merged[0].is_high_confidence())
}

///|
test "merge sv candidates depth below support" {
  let cigar = @src.parse_cigar("10M50D20M")
  let params = @src.SvParams::default()
  let calls = @src.sv_candidates_from_sam("chr1", 100, cigar, params)
  let all_calls : Array[@src.SvCall] = []
  for _ in 0..<2 {
    for c in calls {
      all_calls.push(c)
    }
  }
  let depth_by_chrom : Map[String, Int32] = Map([("chr1", 1)])
  let merged = @src.merge_sv_candidates(all_calls, depth_by_chrom~)
  assert_eq(merged[0].support_reads(), 2)
  assert_eq(merged[0].total_reads(), 2)
  assert_true(merged[0].score() == 1.0)
}
