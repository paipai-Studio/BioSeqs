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
///|
test "csi_query" {
  let idx = @src.parse_csi(@src.base64_decode(csi_base64))
  // ref 0 (chr1): bin 37449 loffset=121 chunk(121,163)
  let q0 = @src.csi_query(idx, 0, 0, 1000)
  assert_eq(q0.min_off, 121L)
  assert_eq(q0.chunks.length(), 1)
  assert_eq(q0.chunks[0].chunk_beg, 121L)
  assert_eq(q0.chunks[0].chunk_end, 163L)
  // sub-region -> same bin
  let q1 = @src.csi_query(idx, 0, 100, 200)
  assert_eq(q1.min_off, 121L)
  assert_eq(q1.chunks.length(), 1)
  // region beyond data: bin walk 37450 -> 37449, empty chunks
  let q2 = @src.csi_query(idx, 0, 16384, 20000)
  assert_eq(q2.min_off, 121L)
  assert_eq(q2.chunks.length(), 0)
  // ref 1 (chr2): bin 37449 loffset=163 chunk(163, 10354688)
  let q3 = @src.csi_query(idx, 1, 0, 500)
  assert_eq(q3.min_off, 163L)
  assert_eq(q3.chunks.length(), 1)
  assert_eq(q3.chunks[0].chunk_beg, 163L)
  assert_eq(q3.chunks[0].chunk_end, 10354688L)
  // out-of-range ref
  let raised = try {
    ignore(@src.csi_query(idx, 5, 0, 100))
    false
  } catch {
    CsiError(_) => true
  }
  assert_true(raised)
}
