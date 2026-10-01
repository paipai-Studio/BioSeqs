///|
test "minimap_minimizers basic extraction" {
  let seq = "ACGTACGTACGT"
  let mins = @src.minimap_minimizers(seq, 4, 5)
  assert_true(mins.length() > 0)
  for m in mins {
    assert_true(m.pos >= 0)
    assert_true(m.pos < seq.length())
    assert_true(m.strand == 0 || m.strand == 1)
  }
}

///|
test "minimap_minimizers empty sequence" {
  let mins = @src.minimap_minimizers("", 4, 5)
  assert_eq(mins.length(), 0)
}

///|
test "minimap_minimizers short sequence" {
  let mins = @src.minimap_minimizers("AC", 4, 5)
  assert_eq(mins.length(), 0)
}

///|
test "minimap_minimizers with N characters" {
  let seq = "ACGTNNNNACGT"
  let mins = @src.minimap_minimizers(seq, 4, 5)
  for m in mins {
    assert_true(m.hash != -1)
  }
}

///|
test "minimap_index build and seed" {
  let ref_seq = "ACGTACGTACGTACGTACGTACGTACGTACGT"
  let refs = [@src.MinimapRefSeq::new("ref1", ref_seq)]
  let index = @src.minimap_index_build(refs, 4, 5)
  let query = "ACGTACGTACGT"
  let seeds = @src.minimap_seed(index, query)
  assert_true(seeds.length() > 0)
}

///|
test "minimap_chain basic chaining" {
  let ref_seq = "ACGTACGTACGTACGTACGTACGTACGTACGT"
  let refs = [@src.MinimapRefSeq::new("ref1", ref_seq)]
  let index = @src.minimap_index_build(refs, 4, 5)
  let query = "ACGTACGTACGT"
  let seeds = @src.minimap_seed(index, query)
  let chains = @src.minimap_chain(seeds, 80, 1)
  assert_true(chains.length() > 0)
  for c in chains {
    assert_true(c.q_start >= 0)
    assert_true(c.q_end > c.q_start)
    assert_true(c.r_start >= 0)
    assert_true(c.r_end > c.r_start)
  }
}

///|
test "minimap_align exact match" {
  let ref_seq = "ACGTACGTACGTACGTACGTACGTACGTACGT"
  let refs = [@src.MinimapRefSeq::new("ref1", ref_seq)]
  let index = @src.minimap_index_build(refs, 4, 5)
  let query = "ACGTACGTACGT"
  let results = @src.minimap_align(index, query, "query1")
  assert_true(results.length() > 0)
  let r = results[0]
  assert_eq(r.ref_name, "ref1")
  assert_true(r.ref_start >= 0)
}

///|
test "minimap_align with mismatch" {
  let ref_seq = "ACGTACGTACGTACGTACGTACGTACGTACGT"
  let refs = [@src.MinimapRefSeq::new("ref1", ref_seq)]
  let index = @src.minimap_index_build(refs, 4, 5)
  let query = "ACGTACGTACCT"
  let results = @src.minimap_align(index, query, "query1")
  assert_true(results.length() > 0)
}

///|
test "minimap_to_sam produces valid SAM" {
  let ref_seq = "ACGTACGTACGTACGTACGTACGTACGTACGT"
  let refs = [@src.MinimapRefSeq::new("ref1", ref_seq)]
  let index = @src.minimap_index_build(refs, 4, 5)
  let query = "ACGTACGTACGT"
  let results = @src.minimap_align(index, query, "query1")
  assert_true(results.length() > 0)
  let sam = @src.minimap_to_sam(results[0], query)
  assert_true(sam.length() > 0)
  assert_true(sam.contains("\t"))
}

///|
test "minimap_align_batch multiple queries" {
  let ref_seq = "ACGTACGTACGTACGTACGTACGTACGTACGT"
  let refs = [@src.MinimapRefSeq::new("ref1", ref_seq)]
  let index = @src.minimap_index_build(refs, 4, 5)
  let queries = [("q1", "ACGTACGTACGT"), ("q2", "ACGTACGTACGT"), ("q3", "TTACGTACGTAC")]
  let sam = @src.minimap_align_batch(index, queries)
  assert_true(sam.length() > 0)
}


