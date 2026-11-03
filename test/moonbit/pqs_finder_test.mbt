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
  let loose = @src.pqs_finder(
    "GGGAGGGTGGGTGGG",
    options=@src.PqsOptions::of(min_score=0),
  )
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

///|
/// pqs_search returns per-position density / max_scores vectors (deep
/// search, mirroring pqsfinder density()/maxScores()).
test "pqs_search_density_vignette" {
  let r = @src.pqs_search("GGGAGGGTGGGTGGG")
  assert_eq(r.hits.length(), 1)
  // exactly one scoreable conformation covers all 15 positions
  assert_eq(r.density, [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1])
  assert_eq(r.max_scores, [
    73, 73, 73, 73, 73, 73, 73, 73, 73, 73, 73, 73, 73, 73, 73,
  ])
}

///|
test "pqs_search_density_flanked" {
  let r = @src.pqs_search("AAAAGGGAGGGTGGGTGGGAAAA")
  assert_eq(r.density, [
    0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0,
  ])
  assert_eq(r.max_scores, [
    0, 0, 0, 0, 73, 73, 73, 73, 73, 73, 73, 73, 73, 73, 73, 73, 73, 73, 73, 0, 0,
    0, 0,
  ])
}

///|
test "pqs_search_density_telomere" {
  // hit [3,24) score 64 -> three leading zeros, 21 covered positions
  let r = @src.pqs_search("TTAGGGTTAGGGTTAGGGTTAGGG")
  assert_eq(r.density[0], 0)
  assert_eq(r.density[2], 0)
  for i = 3; i < 24; i = i + 1 {
    assert_eq(r.density[i], 1)
    assert_eq(r.max_scores[i], 64)
  }
  assert_eq(r.max_scores[0], 0)
}

///|
test "pqs_search_density_minus_strand_mirrored" {
  // C-rich sequence: the hit is on the minus strand at [0,21); density
  // vectors are reported on original coordinates
  let r = @src.pqs_search("CCCTAACCCTAACCCTAACCCTAA", strand="*")
  assert_eq(r.hits[0].strand, "-")
  for i = 0; i < 21; i = i + 1 {
    assert_eq(r.density[i], 1)
    assert_eq(r.max_scores[i], 64)
  }
  for i = 21; i < 24; i = i + 1 {
    assert_eq(r.density[i], 0)
    assert_eq(r.max_scores[i], 0)
  }
}

///|
test "pqs_search_strand_plus_excludes_antisense" {
  let r = @src.pqs_search("CCCTAACCCTAACCCTAACCCTAA", strand="+")
  assert_eq(r.hits.length(), 0)
  for x in r.density {
    assert_eq(x, 0)
  }
  for x in r.max_scores {
    assert_eq(x, 0)
  }
}

///|
test "pqs_search_density_spanning_conformations" {
  // Two separated G4s share scoreable spanning conformations. Reference
  // values from a Python mirror of the exhaustive search.
  let seq = "GGGAGGGTGGGTGGGAAAAAGGGAGGGTGGGTGGG"
  let r = @src.pqs_search(seq)
  assert_eq(r.hits.length(), 2)
  // density peaks at the cores of each quadruplex
  assert_eq(r.density[0], 28)
  assert_eq(r.density[14], 92)
  assert_eq(r.density[20], 92)
  assert_eq(r.density[34], 28)
  // positions in the AAAAA linker are covered by spanning conformations
  // scoring 67 (loops 1,1,5 etc.)
  assert_eq(r.max_scores[0], 73)
  for i = 15; i < 20; i = i + 1 {
    assert_eq(r.max_scores[i], 67)
  }
  assert_eq(r.max_scores[20], 73)
}

///|
test "pqs_search_overlapping_reports_all_conformations" {
  let seq = "GGGAGGGTGGGTGGGAAAAAGGGAGGGTGGGTGGG"
  let nonoverlap = @src.pqs_search(seq, overlapping=false)
  let overlap = @src.pqs_search(seq, overlapping=true)
  assert_eq(nonoverlap.hits.length(), 2)
  // mirror reference: 24 distinct (start,end) conformations
  assert_eq(overlap.hits.length(), 24)
  // every reported hit passes the score threshold and is start-sorted
  let mut prev_start = -1
  for p in overlap.hits {
    assert_true(p.score >= 47)
    assert_true(p.start >= prev_start)
    prev_start = p.start
  }
  // the two canonical non-overlapping hits are present with score 73
  let mut found_first = false
  let mut found_second = false
  let mut found_spanning = false
  for p in overlap.hits {
    if p.start == 0 && p.end == 15 && p.score == 73 {
      found_first = true
    }
    if p.start == 20 && p.end == 35 && p.score == 73 {
      found_second = true
    }
    // spanning conformation across the AAAAA linker, loops (1,1,5)
    if p.start == 4 && p.end == 23 && p.score == 67 {
      found_spanning = true
    }
  }
  assert_true(found_first)
  assert_true(found_second)
  assert_true(found_spanning)
  // overlapping export never reports fewer hits than non-overlapping
  assert_true(overlap.hits.length() >= nonoverlap.hits.length())
}

///|
test "pqs_finder_overlapping_flag" {
  let seq = "GGGAGGGTGGGTGGGAAAAAGGGAGGGTGGGTGGG"
  let few = @src.pqs_finder(seq)
  let many = @src.pqs_finder(seq, overlapping=true)
  assert_eq(few.length(), 2)
  assert_eq(many.length(), 24)
}

///|
test "pqs_sequence_extraction" {
  let seq = "AAAAGGGAGGGTGGGTGGGAAAA"
  let hits = @src.pqs_finder(seq)
  assert_eq(hits[0].sequence(seq), "GGGAGGGTGGGTGGG")
  // minus-strand hit: returns the C-rich substring in original coordinates
  let cseq = "CCCTAACCCTAACCCTAACCCTAA"
  let mhits = @src.pqs_finder(cseq, strand="-")
  assert_eq(mhits.length(), 1)
  assert_eq(mhits[0].sequence(cseq), "CCCTAACCCTAACCCTAACCC")
}

///|
test "pqs_search_empty_sequence" {
  let r = @src.pqs_search("")
  assert_eq(r.hits.length(), 0)
  assert_eq(r.density.length(), 0)
  assert_eq(r.max_scores.length(), 0)
}
