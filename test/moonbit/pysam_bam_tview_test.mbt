///|
/// Tests for `pysam_bam_tview` module — samtools tview style read alignment browser.
/// Reuses the same base64 BAM fixture as bam_test.mbt: read1@100(8M), read2@200(4M1I), read3@300(6M).

let tview_base64 = "QkFNAR8AAABASEQJVk46MS42CkBTUQlTTjpjaHIxCUxOOjEwMDAKAQAAAAUAAABjaHIxAOgDAABNAAAAAAAAAGQAAAAGPEkSAQADAAgAAAD//////////wAAAAByZWFkMQCAAAAAEkgSSAAAAAAAAAAATk1pAAAAAE1EWjgAQVNpCAAAAFhTQStXAAAAAAAAAMgAAAAGFEkSAgAABgUAAAD//////////wAAAAByZWFkMgBAAAAAEQAAABJIEAAAAAAATk1pAQAAAEFTZgAAkEBSR1pncnAxAFhBQkMDAAAAAQIDRQAAAAAAAAAsAQAABihJEgEAEAAGAAAA//////////8AAAAAcmVhZDMAYAAAAP///wAAAAAAAE1DYwFTQ0P6U0lz/f9VSVNg6g=="

let tview_ref : String = "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT"

///|
test "tview_basic_output_nonempty" {
  let bam = @src.parse_bam_base64(tview_base64)
  let out = @src.bam_tview(bam, "chr1", 90, 320, None, 10)
  assert_true(out.length() > 0)
}

///|
test "tview_contains_ruler_prefix" {
  let bam = @src.parse_bam_base64(tview_base64)
  let out = @src.bam_tview(bam, "chr1", 90, 320, None, 10)
  assert_true(out.has_prefix("    "))
}

///|
test "tview_wide_window_sees_all_three_reads" {
  let bam = @src.parse_bam_base64(tview_base64)
  let out = @src.bam_tview(bam, "chr1", 90, 320, None, 10)
  assert_true(out.contains("read1"))
  assert_true(out.contains("read2"))
  assert_true(out.contains("read3"))
}

///|
test "tview_narrow_window_sees_only_read1" {
  let bam = @src.parse_bam_base64(tview_base64)
  // read1 at pos 100, read2 at 200 — window 90-150 only covers read1
  let out = @src.bam_tview(bam, "chr1", 90, 150, None, 10)
  assert_true(out.contains("read1"))
  assert_true(!out.contains("read2"))
}

///|
test "tview_with_refseq_adds_ref_row" {
  let bam = @src.parse_bam_base64(tview_base64)
  let out_no_ref = @src.bam_tview(bam, "chr1", 100, 140, None, 10)
  assert_true(!out_no_ref.contains("REF"))
  let out_with_ref = @src.bam_tview(bam, "chr1", 100, 140, Some(tview_ref), 10)
  assert_true(out_with_ref.contains("REF"))
}

///|
test "tview_max_reads_caps" {
  let bam = @src.parse_bam_base64(tview_base64)
  // Wide window sees all 3, cap to 1 → only read1
  let out = @src.bam_tview(bam, "chr1", 90, 320, None, 1)
  assert_true(out.contains("read1"))
  assert_true(!out.contains("read2"))
  assert_true(!out.contains("read3"))
}

///|
test "tview_max_reads_caps_to_two" {
  let bam = @src.parse_bam_base64(tview_base64)
  let out = @src.bam_tview(bam, "chr1", 90, 320, None, 2)
  assert_true(out.contains("read1"))
  assert_true(out.contains("read2"))
  assert_true(!out.contains("read3"))
}

///|
test "tview_no_overlap_graceful" {
  let bam = @src.parse_bam_base64(tview_base64)
  let out = @src.bam_tview(bam, "chr1", 9999, 10039, None, 10)
  assert_true(out.contains("no reads overlap"))
}

///|
test "tview_autowiden_narrow_window" {
  let bam = @src.parse_bam_base64(tview_base64)
  // 5 bp window auto-widens to 40 — still covers read1 at 100
  let out = @src.bam_tview(bam, "chr1", 100, 105, None, 10)
  assert_true(out.contains("read1"))
}

///|
test "tview_contig_filter_excludes_other_chrom" {
  let bam = @src.parse_bam_base64(tview_base64)
  let out = @src.bam_tview(bam, "chr2", 100, 140, None, 10)
  assert_true(out.contains("no reads overlap"))
}

///|
test "tview_empty_contig_ignores_filter" {
  let bam = @src.parse_bam_base64(tview_base64)
  let out = @src.bam_tview(bam, "", 90, 320, None, 10)
  assert_true(out.contains("read1"))
  assert_true(out.contains("read2"))
  assert_true(out.contains("read3"))
}

///|
test "tview_read2_not_in_read1_window" {
  let bam = @src.parse_bam_base64(tview_base64)
  // read1: [100,108), read2: [200,204) — they don't overlap
  let out = @src.bam_tview(bam, "chr1", 100, 108, None, 10)
  assert_true(out.contains("read1"))
  assert_true(!out.contains("read2"))
}
