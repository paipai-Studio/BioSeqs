///|
/// Tests for the SeqKit record operations + faidx/subseq ports.

///|
fn sq_records() -> Array[@src.SeqRecord] {
  @src.parse_fasta(">s1 d1\nACGT\n>s2\nTT\n>s3 d3\nACGT\n>s4\nGGGGGG\n")
}

///|
test "seqkit_head_range" {
  let r = sq_records()
  assert_eq(@src.seqkit_head(r, 2).length(), 2)
  assert_eq(@src.seqkit_head(r, 2)[0].id, "s1")
  let rg = @src.seqkit_range(r, 2, 3)
  assert_eq(rg.length(), 2)
  assert_eq(rg[0].id, "s2")
  assert_eq(rg[1].id, "s3")
}

///|
test "seqkit_grep" {
  let r = sq_records()
  // full name "s3 d3" contains "d3"
  let g = @src.seqkit_grep(r, "d3", false)
  assert_eq(g.length(), 1)
  assert_eq(g[0].id, "s3")
  // invert keeps the other three
  let inv = @src.seqkit_grep(r, "d3", true)
  assert_eq(inv.length(), 3)
}

///|
test "seqkit_sort" {
  let r = sq_records()
  // length descending: 6,4,4,2 (stable keeps s1 before s3)
  let by_len = @src.seqkit_sort_length(r, true)
  assert_eq(by_len[0].id, "s4")
  assert_eq(by_len[3].id, "s2")
  let by_name = @src.seqkit_sort_name(r)
  assert_eq((by_name[0].id, by_name[1].id, by_name[2].id, by_name[3].id), ("s1", "s2", "s3", "s4"))
}

///|
test "seqkit_renumber_dup_rmdup_common" {
  let r = sq_records()
  // renumber
  let rn = @src.seqkit_renumber(r)
  assert_eq((rn[0].id, rn[3].id), ("1", "4"))
  // duplicate 2x
  assert_eq(@src.seqkit_duplicate(r, 2).length(), 8)
  // rmdup by sequence: "ACGT" duplicated (s1, s3) → keep s1
  let rm = @src.seqkit_rmdup(r)
  assert_eq(rm.length(), 3)
  assert_eq((rm[0].id, rm[1].id, rm[2].id), ("s1", "s2", "s4"))
  // common: a=[s1,s2], b=[s2,s4] → s2
  let a = @src.parse_fasta(">s1\nAC\n>s2\nGG\n")
  let b = @src.parse_fasta(">s2\nGG\n>s4\nTTT\n")
  let c = @src.seqkit_common(a, b)
  assert_eq(c.length(), 1)
  assert_eq(c[0].id, "s2")
}

///|
test "seqkit_faidx_subseq" {
  let content = ">chr1\nACGTACGT\n"
  let faidx = @src.seqkit_faidx(content)
  assert_true(@src.starts_with(faidx, 0, "chr1\t8"))
  let fasta = @src.Fasta::from_content(content)
  assert_eq(@src.seqkit_subseq(fasta, "chr1", 1, 4), "ACGT")
  assert_eq(@src.seqkit_subseq(fasta, "chr1", 5, 8), "ACGT")
  assert_eq(@src.seqkit_subseq(fasta, "chr1", 1, 8), "ACGTACGT")
  assert_eq(@src.seqkit_subseq(fasta, "chr1", 1, 100), "")
  assert_eq(@src.seqkit_subseq(fasta, "chrX", 1, 4), "")
}///|
test "seqkit_sliding_translate_mutate" {
  // sliding: "ACGTACGT", window 4, step 3 -> "ACGT","TACG"
  let win = @src.seqkit_sliding("ACGTACGT", 4, 3)
  assert_eq(win.length(), 2)
  assert_eq(win[0], "ACGT")
  assert_eq(win[1], "TACG")
  // translate: ATG=M, GCC=A
  assert_eq(@src.seqkit_translate("ATGGCC", 1), "MA")
  // frame 2 skips the first base: TGG=W
  assert_eq(@src.seqkit_translate("ATGGCC", 2), "W")
  // mutate point / delete / insert
  assert_eq(@src.seqkit_mutate_point("ACGTACGT", 2, "X"), "AXGTACGT")
  assert_eq(@src.seqkit_mutate_delete("ACGTACGT", 2, 4), "AACGT")
  assert_eq(@src.seqkit_mutate_insert("ACGT", 0, "TT"), "TTACGT")
  assert_eq(@src.seqkit_mutate_insert("ACGT", 2, "NN"), "ACNNGT")
}

///|
test "seqkit_shuffle_sample" {
  let r = @src.seqkit_head(@src.seqkit_renumber(sq_records()), 100)
  assert_eq(@src.seqkit_shuffle(r, 11).length(), 4)
  let s = @src.seqkit_sample(r, 2, 11)
  assert_eq(s.length(), 2)
  // sampling is without replacement (distinct ids)
  assert_true(s[0].id != s[1].id)
}