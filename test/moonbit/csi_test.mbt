///|
/// Tests for the CSI index parser and reg2bins, cross-validated against the
/// hts-specs CSIv1 format and htslib.

///|
let csi_base64 = "Q1NJAQ4AAAAGAAAAJgAAAAIAAAABAAAAAgAAAAAAAAAjAAAAAAAAAAoAAABjaHIxAGNocjIAAgAAAAIAAABJkgAAeQAAAAAAAAABAAAAeQAAAAAAAACjAAAAAAAAAEqSBAAAAAAAAAAAAAIAAAB5AAAAAAAAAKMAAAAAAAAAAgAAAAAAAAAAAAAAAAAAAAIAAABJkgAAowAAAAAAAAABAAAAowAAAAAAAAAAAJ4AAAAAAEqSBAAAAAAAAAAAAAIAAACjAAAAAAAAAAAAngAAAAAAAQAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"

///|
fn csi_ints(a : Array[Int], b : Array[Int]) -> Bool {
  if a.length() != b.length() {
    return false
  }
  for i = 0; i < a.length(); i = i + 1 {
    if a[i] != b[i] {
      return false
    }
  }
  true
}

///|
test "parse_csi_structure" {
  let idx = @src.parse_csi(@src.base64_decode(csi_base64))
  assert_eq(idx.min_shift, 14)
  assert_eq(idx.n_lvls, 6)
  assert_eq(idx.n_no_coor, 0L)
  assert_eq(idx.references.length(), 2)
  // ref 0 (chr1): two bins
  let r0 = idx.references[0]
  assert_eq(r0.bins.length(), 2)
  assert_eq(r0.bins[0].bin, 37449)
  assert_eq(r0.bins[0].loffset, 121L)
  assert_eq(r0.bins[0].chunks.length(), 1)
  assert_eq(r0.bins[0].chunks[0].chunk_beg, 121L)
  assert_eq(r0.bins[0].chunks[0].chunk_end, 163L)
  assert_eq(r0.bins[1].bin, 299594) // meta bin (linear index)
  assert_eq(r0.bins[1].chunks.length(), 2)
  // ref 1 (chr2)
  let r1 = idx.references[1]
  assert_eq(r1.bins[0].bin, 37449)
  assert_eq(r1.bins[0].loffset, 163L)
  assert_eq(r1.bins[0].chunks[0].chunk_beg, 163L)
  assert_eq(r1.bins[0].chunks[0].chunk_end, 10354688L)
}

///|
test "csi_reg2bins_matches_htslib" {
  assert_true(csi_ints(@src.csi_reg2bins(100, 101, 14, 6), [0, 1, 9, 73, 585, 4681, 37449]))
  assert_true(csi_ints(@src.csi_reg2bins(100, 20000, 14, 6), [0, 1, 9, 73, 585, 4681, 37449, 37450]))
  assert_true(csi_ints(@src.csi_reg2bins(0, 1000, 14, 6), [0, 1, 9, 73, 585, 4681, 37449]))
  assert_true(csi_ints(@src.csi_reg2bins(16384, 16385, 14, 6), [0, 1, 9, 73, 585, 4681, 37450]))
  // empty region
  assert_eq(@src.csi_reg2bins(5, 5, 14, 6).length(), 0)
}