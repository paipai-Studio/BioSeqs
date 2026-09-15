///|
/// Tests for primersearch module (Bio.Emboss.PrimerSearch).

test "primersearch_read_primers" {
  let content = "# comment line\nGeneA ATGAAACGTT TTACAGCTAA\nGeneB GGCCTTAACC GGAATTCCGG\n"
  let pairs = @src.primersearch_read_primers(content)
  assert_eq(pairs.length(), 2)
  assert_eq(pairs[0].name, "GeneA")
  assert_eq(pairs[0].forward, "ATGAAACGTT")
  assert_eq(pairs[0].reverse, "TTACAGCTAA")
}

test "primersearch_revcomp" {
  assert_eq(@src.primersearch_revcomp("ACGT"), "ACGT")
  assert_eq(@src.primersearch_revcomp("AAAA"), "TTTT")
  assert_eq(@src.primersearch_revcomp("ATGC"), "GCAT")
}

test "primersearch_find_hits exact" {
  let seq = "GGGATGAAACGTTCCCTTTACAGCTAAAGGG"
  let hits = @src.primersearch_find_hits(seq, "ATGAAACGTT")
  // Forward match at position 3
  assert_true(hits.length() >= 1)
  assert_eq(hits[0].position, 3)
  assert_eq(hits[0].strand, "+")
  assert_eq(hits[0].mismatches, 0)
}

test "primersearch_find_hits with mismatches" {
  let seq = "GGGATGAAACGTTCCCTTTACAGCTAAAGGG"
  // "ATGAAACGTA" differs at last base from target
  let hits0 = @src.primersearch_find_hits(seq, "ATGAAACGTA", max_mismatches=0)
  let hits1 = @src.primersearch_find_hits(seq, "ATGAAACGTA", max_mismatches=1)
  assert_eq(hits0.length(), 0)
  assert_true(hits1.length() >= 1)
  assert_eq(hits1[0].mismatches, 1)
}

test "primersearch_run product" {
  let primers = @src.primersearch_example_primers()
  let seqs = @src.primersearch_example_sequences()
  let result = @src.primersearch_run(primers, seqs)
  assert_eq(result.n_primers, 1)
  assert_eq(result.n_sequences, 1)
  assert_true(result.products.length() >= 1)
  // Product: fwd at 3, rev-rc region ends at 28 -> size 26
  let prod = result.products[0]
  assert_eq(prod.primer_name, "GeneA_F")
  assert_eq(prod.sequence_id, "chr1")
  assert_true(prod.product_size > 0)
}

test "primersearch_no_hit" {
  let primers = [@src.PrimerPair::new("NoHit", "AAAAAAAAAA", "CCCCCCCCCC")]
  let seqs = @src.primersearch_example_sequences()
  let result = @src.primersearch_run(primers, seqs)
  assert_eq(result.products.length(), 0)
}


test "primersearch_to_string" {
  let primers = @src.primersearch_example_primers()
  let seqs = @src.primersearch_example_sequences()
  let result = @src.primersearch_run(primers, seqs)
  let s = @src.primersearch_to_string(result)
  assert_true(s.contains("Primer searching"))
  assert_true(s.contains("GeneA_F"))
  assert_true(s.contains("chr1"))
}

test "primersearch_filter_by_size" {
  let primers = @src.primersearch_example_primers()
  let seqs = @src.primersearch_example_sequences()
  let result = @src.primersearch_run(primers, seqs)
  let prod = result.products[0]
  let in_range = @src.primersearch_filter_by_size(result, prod.product_size, prod.product_size)
  assert_eq(in_range.length(), 1)
  let out_range = @src.primersearch_filter_by_size(result, prod.product_size + 1000, prod.product_size + 2000)
  assert_eq(out_range.length(), 0)
}

test "primersearch no hit" {
  let primers = [@src.PrimerPair::new("NoHit", "AAAAAAAAAA", "CCCCCCCCCC")]
  let seqs = @src.primersearch_example_sequences()
  let result = @src.primersearch_run(primers, seqs)
  assert_eq(result.products.length(), 0)
}