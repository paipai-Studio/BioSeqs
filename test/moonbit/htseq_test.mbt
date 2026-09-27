///|
/// Tests for htseq-count (htseq.mbt).

let htseq_test_bam_base64 = "QkFNAR8AAABASEQJVk46MS42CkBTUQlTTjpjaHIxCUxOOjEwMDAKAQAAAAUAAABjaHIxAOgDAABNAAAAAAAAAGQAAAAGPEkSAQADAAgAAAD//////////wAAAAByZWFkMQCAAAAAEkgSSAAAAAAAAAAATk1pAAAAAE1EWjgAQVNpCAAAAFhTQStXAAAAAAAAAMgAAAAGFEkSAgAABgUAAAD//////////wAAAAByZWFkMgBAAAAAEQAAABJIEAAAAAAATk1pAQAAAEFTZgAAkEBSR1pncnAxAFhBQkMDAAAAAQIDRQAAAAAAAAAsAQAABihJEgEAEAAGAAAA//////////8AAAAAcmVhZDMAYAAAAP///wAAAAAAAE1DYwFTQ0P6U0lz/f9VSVNg6g=="

fn ht_gene(id : String, start : Int, end_ : Int) -> @src.GFFFeature {
  @src.GFFFeature::new("chr1", "src", "gene", start, end_, "+").add_attribute("ID", id)
}

///|
/// 三个 read 分别精确落在三个基因上 → 各 1 个 count，且无未计数类别。
test "htseq_count_basic" {
  let bam = @src.parse_bam_base64(htseq_test_bam_base64)
  let g1 = ht_gene("g1", 101, 108)
  let g2 = ht_gene("g2", 201, 204)
  let g3 = ht_gene("g3", 301, 306)
  let r = @src.htseq_count(bam, [g1, g2, g3])
  assert_eq(r.count_of("g1"), 1)
  assert_eq(r.count_of("g2"), 1)
  assert_eq(r.count_of("g3"), 1)
  assert_eq(r.get_no_feature(), 0)
  assert_eq(r.get_ambiguous(), 0)
  assert_eq(r.get_not_aligned(), 0)
}

///|
/// read1 同时落在 g1 和 g4（重叠）→ ambiguous；read3 不落任何基因 → no_feature。
test "htseq_count_ambiguous_and_no_feature" {
  let bam = @src.parse_bam_base64(htseq_test_bam_base64)
  let g1 = ht_gene("g1", 101, 108)
  let g4 = ht_gene("g4", 100, 110)  // 覆盖 read1 的整个区间
  let g2 = ht_gene("g2", 201, 204)
  let r = @src.htseq_count(bam, [g1, g4, g2])
  assert_eq(r.get_ambiguous(), 1)  // read1
  assert_eq(r.count_of("g2"), 1)   // read2
  assert_eq(r.get_no_feature(), 1) // read3
  assert_eq(r.get_not_aligned(), 0)
}

///|
/// intersection-strict：read 必须完全落在 feature 内。
test "htseq_count_intersection_strict" {
  let bam = @src.parse_bam_base64(htseq_test_bam_base64)
  // g1 只覆盖 read1 的一半 [101,104]：union 算重叠、strict 不算
  let g1 = ht_gene("g1", 101, 104)
  let r = @src.htseq_count(bam, [g1], mode="intersection-strict")
  assert_eq(r.count_of("g1"), 0)
  assert_eq(r.get_no_feature(), 3)

  let r2 = @src.htseq_count(bam, [g1], mode="union")
  assert_eq(r2.count_of("g1"), 1)  // read1 与 [101,104] 重叠
  assert_eq(r2.get_no_feature(), 2)
}
