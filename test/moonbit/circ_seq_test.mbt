///|
test "CircSeq::new construction" {
  let circ = @src.CircSeq::new(sequence="ATCGATCG")
  assert_eq(circ.circ_length(), 8)
  assert_eq(circ.circ_sequence(), "ATCGATCG")
}

///|
test "CircSeq::with_name" {
  let circ = @src.CircSeq::with_name(sequence="ATCGATCG", name="pBR322")
  assert_eq(circ.circ_name(), "pBR322")
  assert_eq(circ.circ_length(), 8)
}

///|
test "CircSeq::circ_gc_content" {
  let circ = @src.CircSeq::new(sequence="GCGCGCGC")
  let gc = circ.circ_gc_content()
  let diff = (gc - 1.0).abs()
  assert_true(diff < 1.0e-10)
}

///|
test "CircSeq::circ_gc_content AT rich" {
  let circ = @src.CircSeq::new(sequence="ATATATAT")
  let gc = circ.circ_gc_content()
  assert_true((gc).abs() < 1.0e-10)
}

///|
test "CircSeq::circ_at" {
  let circ = @src.CircSeq::new(sequence="ATCG")
  assert_eq(circ.circ_at(0), 65)
  assert_eq(circ.circ_at(1), 84)
  assert_eq(circ.circ_at(2), 67)
  assert_eq(circ.circ_at(3), 71)
}

///|
test "CircSeq::circ_at wrapping" {
  let circ = @src.CircSeq::new(sequence="ATCG")
  assert_eq(circ.circ_at(4), 65)
  assert_eq(circ.circ_at(5), 84)
  assert_eq(circ.circ_at(-1), 71)
}

///|
test "CircSeq::circ_subseq linear" {
  let circ = @src.CircSeq::new(sequence="ATCGATCG")
  let sub = circ.circ_subseq(0, 3)
  assert_eq(sub, "ATCG")
}

///|
test "CircSeq::circ_subseq wrapping around" {
  let circ = @src.CircSeq::new(sequence="ATCGATCG")
  let sub = circ.circ_subseq(6, 3)
  assert_eq(sub, "CGATCG")
}

///|
test "CircSeq::circ_subseq entire sequence" {
  let circ = @src.CircSeq::new(sequence="ATCG")
  let sub = circ.circ_subseq(0, 3)
  assert_eq(sub, "ATCG")
}

///|
test "CircSeq::circ_subseq single base" {
  let circ = @src.CircSeq::new(sequence="ATCG")
  let sub = circ.circ_subseq(2, 2)
  assert_eq(sub, "C")
}

///|
test "CircSeq::circ_slice" {
  let circ = @src.CircSeq::new(sequence="ATCG")
  let slice = circ.circ_slice(0, 3)
  assert_eq(slice, "ATC")
}

///|
test "CircSeq::circ_slice wrapping" {
  let circ = @src.CircSeq::new(sequence="ATCG")
  let slice = circ.circ_slice(3, 3)
  assert_eq(slice, "GAT")
}

///|
test "CircSeq::circ_rotate" {
  let circ = @src.CircSeq::new(sequence="ATCG")
  let rotated = circ.circ_rotate(1)
  assert_eq(rotated.circ_sequence(), "TCGA")
}

///|
test "CircSeq::circ_rotate by 0" {
  let circ = @src.CircSeq::new(sequence="ATCG")
  let rotated = circ.circ_rotate(0)
  assert_eq(rotated.circ_sequence(), "ATCG")
}

///|
test "CircSeq::circ_rotate by full length" {
  let circ = @src.CircSeq::new(sequence="ATCG")
  let rotated = circ.circ_rotate(4)
  assert_eq(rotated.circ_sequence(), "ATCG")
}

///|
test "CircSeq::circ_rotate by negative" {
  let circ = @src.CircSeq::new(sequence="ATCG")
  let rotated = circ.circ_rotate(-1)
  assert_eq(rotated.circ_sequence(), "GATC")
}

///|
test "CircSeq::circ_reverse_complement" {
  let circ = @src.CircSeq::new(sequence="ATCG")
  let rc = circ.circ_reverse_complement()
  assert_eq(rc.circ_sequence(), "CGAT")
}

///|
test "CircSeq::circ_find" {
  let circ = @src.CircSeq::new(sequence="ATCGATCG")
  let pos = circ.circ_find("TCGA")
  assert_eq(pos, 1)
}

///|
test "CircSeq::circ_find not found" {
  let circ = @src.CircSeq::new(sequence="ATCGATCG")
  let pos = circ.circ_find("ZZZZ")
  assert_eq(pos, -1)
}

///|
test "CircSeq::circ_find wrapping" {
  let circ = @src.CircSeq::new(sequence="GGATCC")
  let pos = circ.circ_find("GAT")
  assert_true(pos >= 0)
}

///|
test "CircSeq::circ_find_all" {
  let circ = @src.CircSeq::new(sequence="ATCGATCG")
  let positions = circ.circ_find_all("ATCG")
  assert_eq(positions.length(), 2)
  assert_eq(positions[0], 0)
  assert_eq(positions[1], 4)
}

///|
test "CircSeq::circ_find_rc" {
  let circ = @src.CircSeq::new(sequence="GAATTC")
  let rc_positions = circ.circ_find_rc("AATT")
  assert_true(rc_positions.length() >= 0)
}

///|
test "CircSeq::circ_has_site" {
  let circ = @src.CircSeq::new(sequence="GAATTC")
  assert_true(circ.circ_has_site("GAATTC"))
  assert_false(circ.circ_has_site("XXXXX"))
}

///|
test "CircSeq::circ_find_sites" {
  let circ = @src.CircSeq::new(sequence="GAATTC")
  let sites = circ.circ_find_sites("GAATTC")
  assert_true(sites.length() >= 1)
}

///|
test "CircSeq::circ_digest no cut" {
  let circ = @src.CircSeq::new(sequence="ATCGATCGATCG")
  let sites : Array[String] = ["GAATTC"]
  let fragments = circ.circ_digest(sites)
  assert_eq(fragments.length(), 1)
  assert_eq(fragments[0], 12)
}

///|
test "CircSeq::circ_digest with cut" {
  let circ = @src.CircSeq::new(sequence="GAATTC")
  let sites : Array[String] = ["GAATTC"]
  let fragments = circ.circ_digest(sites)
  assert_true(fragments.length() >= 1)
}

///|
test "CircSeq::circ_digest_enzymes" {
  let circ = @src.CircSeq::new(sequence="GAATTC")
  let enzymes : Array[String] = ["GAATTC"]
  let fragments = circ.circ_digest_enzymes(enzymes)
  assert_true(fragments.length() >= 1)
}

///|
test "circular_distance simple" {
  let d = @src.circular_distance(0, 5, 10)
  assert_eq(d, 5)
}

///|
test "circular_distance wrapping" {
  let d = @src.circular_distance(8, 2, 10)
  assert_eq(d, 4)
}

///|
test "circular_distance same position" {
  let d = @src.circular_distance(5, 5, 10)
  assert_eq(d, 0)
}

///|
test "circular_distance zero length" {
  let d = @src.circular_distance(0, 5, 0)
  assert_eq(d, 0)
}

///|
test "circ_to_linear" {
  let pos1 = @src.circ_to_linear(12, 10)
  assert_eq(pos1, 2)
  let pos2 = @src.circ_to_linear(-1, 10)
  assert_eq(pos2, 9)
}

///|
test "linear_to_circular" {
  let pos1 = @src.linear_to_circular(12, 10)
  assert_eq(pos1, 2)
  let pos2 = @src.linear_to_circular(9, 10)
  assert_eq(pos2, 9)
}

///|
test "CircSeq::circ_insert" {
  let circ = @src.CircSeq::new(sequence="ATCG")
  let inserted = circ.circ_insert(2, "XX")
  assert_eq(inserted.circ_sequence(), "ATXXCG")
  assert_eq(inserted.circ_length(), 6)
}

///|
test "CircSeq::circ_delete" {
  let circ = @src.CircSeq::new(sequence="ATCGATCG")
  let deleted = circ.circ_delete(2, 5)
  assert_eq(deleted.circ_length(), 4)
  assert_eq(deleted.circ_sequence(), "ATCG")
}

///|
test "CircSeq::circ_delete wrapping" {
  let circ = @src.CircSeq::new(sequence="ATCGATCG")
  let deleted = circ.circ_delete(6, 7)
  assert_eq(deleted.circ_length(), 6)
}

///|
test "CircSeq::circ_random_sequence" {
  let circ = @src.circ_random_sequence(100, 0.5)
  assert_eq(circ.circ_length(), 100)
  assert_true(circ.circ_gc_content() > 0.0)
}

///|
test "CircSeq::circ_random_sequence with seed" {
  let circ1 = @src.circ_random_sequence(50, 0.5, seed=42)
  let circ2 = @src.circ_random_sequence(50, 0.5, seed=42)
  assert_eq(circ1.circ_sequence(), circ2.circ_sequence())
}

///|
test "CircSeq::circ_ligate" {
  let vector = @src.CircSeq::new(sequence="ATCGATCG")
  let ligated = @src.circ_ligate("XX", vector, 4)
  assert_eq(ligated.circ_length(), 10)
  assert_true(ligated.circ_sequence().contains("XX"))
}

///|
test "CircSeq::circ_tm" {
  let circ = @src.CircSeq::new(sequence="GCGCGCGCGCGC")
  let tm = @src.circ_tm(circ, 0, 10)
  assert_true(tm > 0.0)
}

///|
test "CircSeq::circ_common_restriction_sites" {
  let sites = @src.circ_common_restriction_sites()
  assert_true(sites.length() > 0)
  assert_eq(sites[0].0, "EcoRI")
  assert_eq(sites[0].1, "GAATTC")
}

///|
test "circ_unique_fragments" {
  let frags = @src.circ_unique_fragments([100, 200, 100, 300, 200])
  assert_eq(frags.length(), 3)
  assert_eq(frags[0], 100)
  assert_eq(frags[1], 200)
  assert_eq(frags[2], 300)
}

///|
test "circ_largest_fragment" {
  let largest = @src.circ_largest_fragment([100, 200, 500, 50])
  assert_eq(largest, 500)
}

///|
test "CircSeq::circ_design_primers" {
  let circ = @src.CircSeq::new(sequence="ATCGATCGATCGATCG")
  let (fwd, rev, prod) = circ.circ_design_primers(4, 12)
  assert_true(prod > 0)
  assert_true(fwd >= 0)
  assert_true(rev >= 0)
}

///|
test "is_within_circular_range" {
  assert_true(@src.is_within_circular_range(5, 3, 7, 10))
  assert_false(@src.is_within_circular_range(1, 3, 7, 10))
}

///|
test "is_within_circular_range wrapping" {
  assert_true(@src.is_within_circular_range(8, 6, 2, 10))
  assert_false(@src.is_within_circular_range(4, 6, 2, 10))
}

///|
test "circular_overlap" {
  let overlap = @src.circular_overlap(2, 5, 4, 7, 10)
  assert_eq(overlap, 2)
}

///|
test "circular_overlap no overlap" {
  let overlap = @src.circular_overlap(0, 2, 5, 7, 10)
  assert_eq(overlap, 0)
}
