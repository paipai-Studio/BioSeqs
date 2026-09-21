///|
/// Test mpileup tabular output (`BamFile::mpileup`) against a real
/// pysam-generated BAM fixture, plus the `to_uppercase` helper used to
/// normalize the reference base column.

///|
/// Base64-encoded raw BAM (pysam 0.24.1), three reads on chr1 (LN:1000):
///   read1: pos=100 8M,     seq ACGTACGT, mapq 60, forward
///   read2: pos=200 4M1I,   seq ACGTA,    mapq 20, forward (secondary/qcfail)
///   read3: pos=300 6M,     seq NNNNNN,   mapq 40, reverse
let mp_bam_base64 = "QkFNAR8AAABASEQJVk46MS42CkBTUQlTTjpjaHIxCUxOOjEwMDAKAQAAAAUAAABjaHIxAOgDAABNAAAAAAAAAGQAAAAGPEkSAQADAAgAAAD//////////wAAAAByZWFkMQCAAAAAEkgSSAAAAAAAAAAATk1pAAAAAE1EWjgAQVNpCAAAAFhTQStXAAAAAAAAAMgAAAAGFEkSAgAABgUAAAD//////////wAAAAByZWFkMgBAAAAAEQAAABJIEAAAAAAATk1pAQAAAEFTZgAAkEBSR1pncnAxAFhBQkMDAAAAAQIDRQAAAAAAAAAsAQAABihJEgEAEAAGAAAA//////////8AAAAAcmVhZDMAYAAAAP///wAAAAAAAE1DYwFTQ0P6U0lz/f9VSVNg6g=="

///|
/// Reference sequence covering the three read regions: read1/read2 match,
/// read3 (NNNNNN) mismatches against ACGTAC.
fn mp_ref() -> String {
  let mut r = ""
  let mut i = 0
  while i < 306 {
    if i >= 100 && i < 108 {
      r = r + "ACGTACGT"[i - 100:i - 100 + 1].to_owned()
    } else if i >= 200 && i < 204 {
      r = r + "ACGT"[i - 200:i - 200 + 1].to_owned()
    } else if i >= 300 && i < 306 {
      r = r + "ACGTAC"[i - 300:i - 300 + 1].to_owned()
    } else {
      r = r + "A"
    }
    i = i + 1
  }
  r
}

///|
test "mpileup_read1_head_match" {
  let bam = @src.parse_bam_base64(mp_bam_base64)
  let lines = bam.mpileup("chr1", 0, 306, mp_ref())
  assert_eq(lines.length(), 18)
  let l = lines[0]
  assert_eq(l.chrom, "chr1")
  assert_eq(l.pos, 101)
  assert_eq(l.ref_base, "A")
  assert_eq(l.depth, 1)
  assert_eq(l.bases, "^].")
  assert_eq(l.quals, "!")
}

///|
test "mpileup_read1_tail" {
  let bam = @src.parse_bam_base64(mp_bam_base64)
  let lines = bam.mpileup("chr1", 0, 306, mp_ref())
  let l = lines[7]
  assert_eq(l.pos, 108)
  assert_eq(l.bases, ".$")
  assert_eq(l.quals, "!")
}

///|
test "mpileup_read2_insertion_anchor" {
  let bam = @src.parse_bam_base64(mp_bam_base64)
  let lines = bam.mpileup("chr1", 0, 306, mp_ref())
  let l = lines[11]
  assert_eq(l.pos, 204)
  assert_eq(l.ref_base, "T")
  assert_eq(l.bases, ".+1A$")
  assert_eq(l.quals, "!")
}

///|
test "mpileup_read3_reverse_mismatch" {
  let bam = @src.parse_bam_base64(mp_bam_base64)
  let lines = bam.mpileup("chr1", 0, 306, mp_ref())
  let head = lines[12]
  assert_eq(head.pos, 301)
  assert_eq(head.bases, "^In")
  assert_eq(head.quals, "!")
  let tail = lines[17]
  assert_eq(tail.pos, 306)
  assert_eq(tail.bases, "n$")
  assert_eq(tail.quals, "!")
}

///|
test "to_uppercase_basic" {
  assert_eq(@src.to_uppercase("acgtacgtn"), "ACGTACGTN")
  assert_eq(@src.to_uppercase("ACGTT"), "ACGTT")
  assert_eq(@src.to_uppercase(""), "")
}