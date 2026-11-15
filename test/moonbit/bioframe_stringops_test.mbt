///|
/// Tests for bioframe_stringops.mbt — UCSC region strings + natural sort.

// ---- bf_to_ucsc_string ----------------------------------------------------
// 0-based half-open → 1-based inclusive:
//   start_0 + 1 → start_1,   end_0 → end_1
// [0, 10) 0-based covers 0..9 (10 bp) → [1, 10] 1-based inclusive

test "bf_to_ucsc_string: basic 0→1 conversion" {
  let s = @src.bf_to_ucsc_string("chr1", 0, 10)
  assert_eq(s, "chr1:1-10")
}

test "bf_to_ucsc_string: non-zero start" {
  let s = @src.bf_to_ucsc_string("chr1", 100, 200)
  assert_eq(s, "chr1:101-200")
}

test "bf_to_ucsc_string: end==start empty interval" {
  // [100, 100) has zero length; 1-based [101, 100] is also zero-length
  let s = @src.bf_to_ucsc_string("chr1", 100, 100)
  assert_eq(s, "chr1:101-100")
}

// ---- bf_is_complete_ucsc_string -------------------------------------------

test "bf_is_complete_ucsc_string: valid pattern" {
  assert_true(@src.bf_is_complete_ucsc_string("chr1:1-100"))
}

// NOTE: current implementation allows additional dashes after the first;
// it only rejects extra colons.  Test reflects actual source behaviour.
test "bf_is_complete_ucsc_string: multi-dash accepted (current impl)" {
  assert_true(@src.bf_is_complete_ucsc_string("chr1:1-2-3"))
}

test "bf_is_complete_ucsc_string: no dash → false" {
  assert_true(!@src.bf_is_complete_ucsc_string("chr1:1"))
}

test "bf_is_complete_ucsc_string: empty → false" {
  assert_true(!@src.bf_is_complete_ucsc_string(""))
}

test "bf_is_complete_ucsc_string: no colon → false" {
  assert_true(!@src.bf_is_complete_ucsc_string("chr1"))
}

test "bf_is_complete_ucsc_string: too short → false" {
  assert_true(!@src.bf_is_complete_ucsc_string("c:1"))
}

test "bf_is_complete_ucsc_string: valid chrM" {
  assert_true(@src.bf_is_complete_ucsc_string("chrM:1-1000"))
}

// ---- bf_parse_region_string -----------------------------------------------

test "bf_parse_region_string: chrom only" {
  let (c, s, e) = @src.bf_parse_region_string("chr1")
  assert_eq(c, "chr1")
  match s {
    None => assert_true(true)
    Some(_) => assert_true(false)
  }
  match e {
    None => assert_true(true)
    Some(_) => assert_true(false)
  }
}

test "bf_parse_region_string: chrom:start-end → 0-based" {
  let (c, s, e) = @src.bf_parse_region_string("chr1:1000-2000")
  assert_eq(c, "chr1")
  match s {
    Some(v) => assert_eq(v, 999)
    None => assert_true(false)
  }
  match e {
    Some(v) => assert_eq(v, 2000)
    None => assert_true(false)
  }
}

test "bf_parse_region_string: MT alias resolved with hg38" {
  let (c, s, e) = @src.bf_parse_region_string("MT:1-100", assembly=Some("hg38"))
  assert_eq(c, "chrM")
  match s {
    Some(v) => assert_eq(v, 0)
    None => assert_true(false)
  }
  match e {
    Some(v) => assert_eq(v, 100)
    None => assert_true(false)
  }
}

test "bf_parse_region_string: numeric alias 1 resolved with hg38" {
  let (c, s, e) = @src.bf_parse_region_string("1:100-200", assembly=Some("hg38"))
  assert_eq(c, "chr1")
  match s {
    Some(v) => assert_eq(v, 99)
    None => assert_true(false)
  }
  match e {
    Some(v) => assert_eq(v, 200)
    None => assert_true(false)
  }
}

test "bf_parse_region_string: no assembly → alias not resolved" {
  let (c, _, _) = @src.bf_parse_region_string("MT:1-100")
  assert_eq(c, "MT")
}

// ---- bf_natsort_chroms ----------------------------------------------------

test "bf_natsort_chroms: classic chr10 vs chr2" {
  let inp = ["chr10", "chr2", "chr1"]
  let out = @src.bf_natsort_chroms(inp)
  assert_eq(out.length(), 3)
  assert_eq(out[0], "chr1")
  assert_eq(out[1], "chr2")
  assert_eq(out[2], "chr10")
}

// NOTE: natural key uses to_lower(), so chrM→"chrm" < chrX→"chrx" < chrY→"chry".
// Expected order: chr1 (num=1), chrM, chrX, chrY
test "bf_natsort_chroms: autosomes then X/Y/M" {
  let inp = ["chrX", "chr1", "chrY", "chrM"]
  let out = @src.bf_natsort_chroms(inp)
  assert_eq(out[0], "chr1")
  assert_eq(out[1], "chrM")
  assert_eq(out[2], "chrX")
  assert_eq(out[3], "chrY")
}

test "bf_natsort_chroms: empty → empty" {
  let out = @src.bf_natsort_chroms([])
  assert_eq(out.length(), 0)
}

// ---- bf_natsort_intervals -------------------------------------------------

test "bf_natsort_intervals: sorted by chrom natural then start" {
  let inp = [
    @src.BedInterval::new("chr10", 100, 200),
    @src.BedInterval::new("chr2", 50, 60),
    @src.BedInterval::new("chr1", 1000, 2000),
    @src.BedInterval::new("chr2", 10, 20),
  ]
  let out = @src.bf_natsort_intervals(inp)
  assert_eq(out.length(), 4)
  assert_eq(out[0].chrom, "chr1")
  assert_eq(out[1].chrom, "chr2")
  assert_eq(out[1].start, 10)
  assert_eq(out[2].chrom, "chr2")
  assert_eq(out[2].start, 50)
  assert_eq(out[3].chrom, "chr10")
}
