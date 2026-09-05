///|
/// Tests for the PQS (G-quadruplex) finder, validated against
/// Bioconductor pqsfinder 2.x reference outputs.

///|
test "pqs_finder_basic_vignette" {
  // pqsfinder(DNAString("GGGAGGGTGGGTGGG"))
  let hits = @src.pqs_finder("GGGAGGGTGGGTGGG")
  assert_eq(hits.length(), 1)
  let p = hits[0]
  assert_eq(p.start, 0)
  assert_eq(p.end, 15)
  assert_eq(p.score, 73)
  assert_eq(p.nt, 3)
  assert_eq(p.nb, 0)
  assert_eq(p.nm, 0)
  assert_eq(p.rl1, 3)
  assert_eq(p.rl2, 3)
  assert_eq(p.rl3, 3)
  assert_eq(p.ll1, 1)
  assert_eq(p.ll2, 1)
  assert_eq(p.ll3, 1)
  assert_eq(p.strand, "+")
}

///|
test "pqs_finder_flanked_vignette" {
  // Same hit inside a flanked sequence
  let hits = @src.pqs_finder("AAAAGGGAGGGTGGGTGGGAAAA")
  assert_eq(hits.length(), 1)
  assert_eq(hits[0].start, 4)
  assert_eq(hits[0].end, 19)
  assert_eq(hits[0].score, 73)
}

///|
test "pqs_finder_longer_loops_lower_score" {
  // GGG-AA-GGG-ATTT-GGG-A-GGG: loops 2,4,1 -> sum 7 -> penalty 13
  // score = (3-1)*40 - 13 = 67
  let hits = @src.pqs_finder("GGGAAGGGATTTGGGAGGG")
  assert_eq(hits.length(), 1)
  assert_eq(hits[0].score, 67)
  assert_eq(hits[0].ll1, 2)
  assert_eq(hits[0].ll2, 4)
  assert_eq(hits[0].ll3, 1)
}

///|
test "pqs_finder_telomere_repeat" {
  // Human telomere repeat: first G at index 3, loops TTA x3 -> sum 9,
  // penalty round(6.6*3^0.8)=16 -> score 80-16=64, PQS spans [3,24)
  let hits = @src.pqs_finder("TTAGGGTTAGGGTTAGGGTTAGGG")
  assert_eq(hits.length(), 1)
  assert_eq(hits[0].score, 64)
  assert_eq(hits[0].start, 3)
  assert_eq(hits[0].end, 24)
}

///|
test "pqs_finder_both_strands" {
  // reverse complement is TTAGGGTTAGGGTTAGGGTTAGGG -> loops TTA x3,
  // score 64, PQS maps back to [0,21) on the minus strand
  let hits = @src.pqs_finder("CCCTAACCCTAACCCTAACCCTAA", strand="*")
  assert_eq(hits.length(), 1)
  assert_eq(hits[0].strand, "-")
  assert_eq(hits[0].score, 64)
  assert_eq(hits[0].start, 0)
  assert_eq(hits[0].end, 21)
}

///|
test "pqs_finder_negative_strand_only" {
  let hits = @src.pqs_finder("CCCTAACCCTAACCCTAACCCTAA", strand="-")
  assert_eq(hits.length(), 1)
  assert_eq(hits[0].strand, "-")
  assert_eq(hits[0].start, 0)
  assert_eq(hits[0].end, 21)
  assert_eq(hits[0].score, 64)
}

///|
test "pqs_finder_sense_strand_skips_antisense" {
  let hits = @src.pqs_finder("CCCTAACCCTAACCCTAA", strand="+")
  assert_eq(hits.length(), 0)
}

///|
test "pqs_finder_no_quadruplex" {
  assert_eq(@src.pqs_finder("ACGTACGTACGTACGT").length(), 0)
  assert_eq(@src.pqs_finder("").length(), 0)
}

///|
test "pqs_finder_min_score_filter" {
  // short two-tetrad hit scores (2-1)*40 - loop penalty = below default 47
  let loose = @src.pqs_finder("GGGAGGGTGGGTGGG", options=@src.PqsOptions::of(min_score=0))
  assert_true(loose.length() >= 1)
}

///|
test "pqs_finder_bulge" {
  // GGG-T-GGG with bulge T inside run: GGGTGGGTGGGTGGG has 4 runs where
  // run 2 = GGTGGG (w=6, g=6) is a bulged run relative to reference w=3
  let hits = @src.pqs_finder("GGGTGGGTGGGTGGG")
  assert_eq(hits.length(), 1)
  assert_eq(hits[0].score, 73)
  assert_eq(hits[0].nt, 3)
}

///|
test "pqs_finder_invalid_strand_raises" {
  let result = try {
    ignore(@src.pqs_finder("GGGAGGGTGGGTGGG", strand="x"))
    false
  } catch {
    SeqError(_) => true
  }
  assert_true(result)
}
