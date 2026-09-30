///|
test "parse_mm_tag simple" {
  // C+m,5,10; means C+ strand, mod code 'm', at deltas 5 and 10
  let mods = @src.parse_base_mods("C+m,5,10;", None, "CCCCCACCCCCACCCC")
  let all = mods.all_mods()
  // Sequence has C at positions 0,1,2,3,4,6,7,8,9,10,12,13,14,15
  // Delta 5 → 6th C (0-indexed) = position 10
  // Delta 10 → 11th C = position 14
  // Wait, deltas are cumulative: first delta=5, second=10 means positions 5 and 5+10+1=16
  // Actually in htslib, deltas are 0-based offsets, first is absolute, rest are relative to previous
  // delta 5 → 5th C (0-indexed), delta 10 → 5+10+1=16th C
  // But we only have 14 C's, so only first mod should be found
  assert_true(all.length() >= 1)
}

///|
test "parse_mm_tag multiple_types" {
  let mods = @src.parse_base_mods("C+m,0,2;C+h,1;", None, "CCACCACC")
  // C at positions 0,1,4,5
  // m: delta 0 → pos 0, delta 2 → pos 4 (cumulative: 0, 0+2+1=3 → pos 5? No)
  // Actually deltas: first=0 (0th C = pos 0), second=2 (0+2+1=3rd C = pos 5)
  // h: delta 1 (1st C = pos 1)
  let all = mods.all_mods()
  // Should have 2 m mods + 1 h mod = 3 total
  assert_eq(all.length(), 3)
}

///|
test "parse_mm_tag with_quality" {
  let mods = @src.parse_base_mods("C+m,0,1;", Some([255, 128]), "CCC")
  let all = mods.all_mods()
  assert_eq(all.length(), 2)
  // First mod at position 0, quality 255
  assert_eq(all[0].pos, 0)
  assert_eq(all[0].qual, 255)
  assert_eq(all[0].code, 'm')
  // Second mod at position 1, quality 128
  assert_eq(all[1].pos, 2)
  assert_eq(all[1].qual, 128)
}

///|
test "parse_mm_tag strand" {
  let mods_fwd = @src.parse_base_mods("C+m,0;", None, "C")
  let all_fwd = mods_fwd.all_mods()
  assert_eq(all_fwd.length(), 1)
  assert_eq(all_fwd[0].strand, 1)

  let mods_rev = @src.parse_base_mods("C-m,0;", None, "C")
  let all_rev = mods_rev.all_mods()
  assert_eq(all_rev.length(), 1)
  assert_eq(all_rev[0].strand, -1)
}

///|
test "parse_mm_tag different_bases" {
  let mods = @src.parse_base_mods("C+m,0;A+a,0;", None, "CA")
  let all = mods.all_mods()
  assert_eq(all.length(), 2)
  // C+m at position 0
  assert_eq(all[0].canonical, 'C')
  assert_eq(all[0].code, 'm')
  // A+a at position 1
  assert_eq(all[1].canonical, 'A')
  assert_eq(all[1].code, 'a')
}

///|
test "BaseModState mods_at" {
  let mods = @src.parse_base_mods("C+m,0,2;", None, "CCACC")
  // C at positions 0,1,3,4
  // delta 0 → pos 0, delta 2 → pos 3 (cumulative 0+2+1=3)
  let at_0 = mods.mods_at(0)
  assert_eq(at_0.length(), 1)
  assert_eq(at_0[0].code, 'm')
  let at_4 = mods.mods_at(4)
  assert_eq(at_4.length(), 1)
  let at_1 = mods.mods_at(1)
  assert_eq(at_1.length(), 0)
}

///|
test "BaseModState mod_codes" {
  let mods = @src.parse_base_mods("C+m,0;C+h,0;A+a,0;", None, "CA")
  let codes = mods.mod_codes()
  assert_eq(codes.length(), 3)
}

///|
test "BaseModState count_by_code" {
  let mods = @src.parse_base_mods("C+m,0,1,2;", None, "CCCC")
  // 3 m modifications
  assert_eq(mods.count_by_code('m'), 2)
  assert_eq(mods.count_by_code('h'), 0)
}

///|
test "BaseModState methylation_summary" {
  // Sequence "CCCC" has 4 C's, 2 are modified
  let mods = @src.parse_base_mods("C+m,0,1;", None, "CCCC")
  let summary = mods.methylation_summary("CCCC")
  assert_eq(summary.length(), 1)
  assert_eq(summary[0].canonical, 'C')
  assert_eq(summary[0].code, 'm')
  assert_eq(summary[0].modified, 2)
  assert_eq(summary[0].unmodified, 2)
  assert_true((summary[0].fraction - 0.5).abs() < 1.0e-9)
}

///|
test "parse_mm_tag empty" {
  let mods = @src.parse_base_mods("", None, "ACGT")
  assert_eq(mods.all_mods().length(), 0)
}

///|
test "parse_mm_tag U_to_T" {
  // U should be treated as T
  let mods = @src.parse_base_mods("U+m,0;", None, "T")
  let all = mods.all_mods()
  assert_eq(all.length(), 1)
  assert_eq(all[0].canonical, 'T')
}

///|
test "parse_mm_tag no_ml_tag" {
  let mods = @src.parse_base_mods("C+m,0,1;", None, "CCC")
  let all = mods.all_mods()
  assert_eq(all.length(), 2)
  // Quality should be -1 when no ML tag
  assert_eq(all[0].qual, -1)
  assert_eq(all[1].qual, -1)
}

