///|
/// Tests for bioframe_ucsc.mbt — UCSC chrom.sizes/centromeres fetch with fallback.

// ---- bf_ucsc_chromsizes_url ----------------------------------------------

test "bf_ucsc_chromsizes_url: builds expected UCSC URL" {
  let u = @src.bf_ucsc_chromsizes_url("hg38")
  assert_true(u.contains("hgdownload.soe.ucsc.edu"))
  assert_true(u.contains("hg38"))
  assert_true(u.has_suffix(".chrom.sizes"))
}

// ---- bf_ucsc_centromeres_url --------------------------------------------

///|
test "bf_ucsc_centromeres_url: builds expected UCSC URL" {
  let u = @src.bf_ucsc_centromeres_url("mm10")
  assert_true(u.has_suffix(".centromeres"))
}

// ---- bf_fetch_chromsizes_from_ucsc (no fetch_url → fallback) ------------

///|
test "bf_fetch_chromsizes_from_ucsc: no fetch → built-in fallback" {
  let cs = @src.bf_fetch_chromsizes_from_ucsc("hg38")
  assert_true(cs.length() >= 25)
  let mut found_chr1 = false
  for pair in cs {
    if pair.0 == "chr1" {
      assert_eq(pair.1.to_int64(), 248956422L)
      found_chr1 = true
      break
    }
  }
  assert_true(found_chr1)
}

///|
test "bf_fetch_chromsizes_from_ucsc: unknown assembly → empty" {
  let cs = @src.bf_fetch_chromsizes_from_ucsc("fake")
  assert_eq(cs.length(), 0)
}

// ---- bf_fetch_chromsizes_from_ucsc with fetch_url callback ---------------

///|
test "bf_fetch_chromsizes_from_ucsc: fetch returns None → fallback" {
  let fetcher = fn(_url : String) -> String? { None }
  let cs = @src.bf_fetch_chromsizes_from_ucsc("hg38", fetch_url=Some(fetcher))
  assert_true(cs.length() >= 25)
}

///|
test "bf_fetch_chromsizes_from_ucsc: fetch returns valid TSV → used" {
  let fetcher = fn(_url : String) -> String? {
    Some("chrFake\t1000\nchrOther\t2000\n")
  }
  let cs = @src.bf_fetch_chromsizes_from_ucsc("hg38", fetch_url=Some(fetcher))
  assert_eq(cs.length(), 2)
  assert_eq(cs[0].0, "chrFake")
  assert_eq(cs[1].0, "chrOther")
}

// NOTE: source's bf_parse_chromsizes_tsv accepts any two-tab-line because
// parse_int never fails — so "not\tproper\ttsv" would actually parse as
// chrom="not", length="tsv" → valid length=0 (or parse_int fallback).
// To trigger fallback we need an empty body that yields 0 rows.

///|
test "bf_fetch_chromsizes_from_ucsc: fetch returns empty → fallback" {
  let fetcher = fn(_url : String) -> String? { Some("") }
  let cs = @src.bf_fetch_chromsizes_from_ucsc("hg38", fetch_url=Some(fetcher))
  assert_true(cs.length() >= 25)
}

// ---- bf_fetch_centromeres_from_ucsc (no fetch_url → fallback) ------------

///|
test "bf_fetch_centromeres_from_ucsc: hg19 fallback → 24 rows" {
  let cents = @src.bf_fetch_centromeres_from_ucsc("hg19")
  assert_eq(cents.length(), 24)
  assert_eq(cents[0].0, "chr1")
}

///|
test "bf_fetch_centromeres_from_ucsc: fetch returns None → fallback" {
  let fetcher = fn(_url : String) -> String? { None }
  let cents = @src.bf_fetch_centromeres_from_ucsc(
    "hg19",
    fetch_url=Some(fetcher),
  )
  assert_eq(cents.length(), 24)
}

///|
test "bf_fetch_centromeres_from_ucsc: fetch returns valid → used" {
  let fetcher = fn(_url : String) -> String? {
    Some("chrX\t100\t200\nchrY\t50\t150\n")
  }
  let cents = @src.bf_fetch_centromeres_from_ucsc(
    "hg19",
    fetch_url=Some(fetcher),
  )
  assert_eq(cents.length(), 2)
  assert_eq(cents[0].0, "chrX")
}
