///|
fn alignment_counts_pair(
  target : String,
  query : String,
  target_coordinates : Array[Int],
  query_coordinates : Array[Int],
) -> @src.CoordinatePairwiseAlignment {
  @src.coordinate_pairwise_alignment(
    "target", target, "query", query, target_coordinates, query_coordinates,
  ) catch {
    _ => abort("valid alignment should build")
  }
}

///|
fn alignment_counts_basic(
  target : String,
  query : String,
  target_coordinates : Array[Int],
  query_coordinates : Array[Int],
) -> @src.AlignmentCounts {
  alignment_counts_pair(target, query, target_coordinates, query_coordinates).alignment_counts() catch {
    _ => abort("valid alignment should count")
  }
}

///|
fn alignment_counts_affine(
  open_score : Double,
  extend_score : Double,
) -> @src.AlignmentGapScores {
  @src.AlignmentGapScores::affine(open_score, extend_score) catch {
    _ => abort("finite affine scores should build")
  }
}

///|
test "alignment counts: ungapped identities" {
  let counts = alignment_counts_basic("ACGT", "ACGT", [0, 4], [0, 4])
  assert_eq(counts.aligned, 4)
  assert_eq(counts.identities, 4)
  assert_eq(counts.mismatches, 0)
  assert_eq(counts.gaps(), 0)
}

///|
test "alignment counts: ungapped mismatches" {
  let counts = alignment_counts_basic("ACGT", "AGCT", [0, 4], [0, 4])
  assert_eq(counts.aligned, 4)
  assert_eq(counts.identities, 2)
  assert_eq(counts.mismatches, 2)
}

///|
test "alignment counts: default scores are unknown" {
  let counts = alignment_counts_basic("AC", "AC", [0, 2], [0, 2])
  assert_eq(counts.positives, None)
  assert_eq(counts.substitution_score, None)
  assert_eq(counts.gap_score, None)
  assert_eq(counts.score, None)
}

///|
test "alignment counts: wildcard excludes mismatch" {
  let alignment = alignment_counts_pair("AC?T", "ACGT", [0, 4], [0, 4])
  let config = @src.AlignmentCountsConfig::with_wildcard('?')
  let counts = alignment.alignment_counts(config~) catch {
    _ => abort("wildcard counting should succeed")
  }
  assert_eq(counts.aligned, 4)
  assert_eq(counts.identities, 3)
  assert_eq(counts.mismatches, 0)
}

///|
test "alignment counts: wildcard excludes identity" {
  let alignment = alignment_counts_pair("AC?T", "AC?T", [0, 4], [0, 4])
  let config = @src.AlignmentCountsConfig::with_wildcard('?')
  let counts = alignment.alignment_counts(config~) catch {
    _ => abort("wildcard counting should succeed")
  }
  assert_eq(counts.identities, 3)
  assert_eq(counts.mismatches, 0)
}

///|
test "alignment counts: wildcard remains aligned" {
  let alignment = alignment_counts_pair("NN", "AA", [0, 2], [0, 2])
  let counts = alignment.alignment_counts(
    config=@src.AlignmentCountsConfig::with_wildcard('N'),
  ) catch {
    _ => abort("wildcard counting should succeed")
  }
  assert_eq(counts.aligned, 2)
  assert_eq(counts.identities + counts.mismatches, 0)
}

///|
test "alignment counts: BLOSUM62 official protein example" {
  let alignment = alignment_counts_pair(
    "EPQSDPSVEPPLSQETFSDLWKLLPE",
    "EPSSETGMDPPLSQETFEDLWSLLPD",
    [0, 26],
    [0, 26],
  )
  let config = @src.AlignmentCountsConfig::with_matrix(@src.blosum62())
  let counts = alignment.alignment_counts(config~) catch {
    _ => abort("BLOSUM62 counting should succeed")
  }
  assert_eq(counts.identities, 17)
  assert_eq(counts.mismatches, 9)
  assert_eq(counts.positives, Some(21))
  assert_eq(counts.substitution_score, Some(101.0))
  assert_eq(counts.score, None)
}

///|
test "alignment counts: BLOSUM45 official protein example" {
  let alignment = alignment_counts_pair(
    "EPQSDPSVEPPLSQETFSDLWKLLPE",
    "EPSSETGMDPPLSQETFEDLWSLLPD",
    [0, 26],
    [0, 26],
  )
  let config = @src.AlignmentCountsConfig::with_matrix(@src.blosum45())
  let counts = alignment.alignment_counts(config~) catch {
    _ => abort("BLOSUM45 counting should succeed")
  }
  assert_eq(counts.positives, Some(21))
  assert_eq(counts.substitution_score, Some(122.0))
  let matrix = @src.blosum45()
  assert_eq(matrix.get_score("C", "C"), 12)
  assert_eq(matrix.get_score("P", "P"), 9)
  assert_eq(matrix.get_score("W", "W"), 15)
  assert_eq(matrix.get_score("K", "S"), -1)
  assert_eq(matrix.get_score("S", "K"), -1)
}

///|
test "alignment counts: matrix accepts lower case residues" {
  let alignment = alignment_counts_pair("ar", "AR", [0, 2], [0, 2])
  let counts = alignment.alignment_counts(
    config=@src.AlignmentCountsConfig::with_matrix(@src.blosum62()),
  ) catch {
    _ => abort("case-insensitive matrix scoring should succeed")
  }
  assert_eq(counts.identities, 0)
  assert_eq(counts.mismatches, 2)
  assert_eq(counts.positives, Some(2))
  assert_eq(counts.substitution_score, Some(9.0))
}

///|
test "alignment counts: matrix rejects unknown residue" {
  let alignment = alignment_counts_pair("AX", "AA", [0, 2], [0, 2])
  let failed = try {
    ignore(
      alignment.alignment_counts(
        config=@src.AlignmentCountsConfig::with_matrix(@src.blosum62()),
      ),
    )
    false
  } catch {
    AlignmentCountsError(_) => true
  }
  assert_true(failed)
}

///|
test "alignment counts: matrix wildcard skips unknown residue" {
  let alignment = alignment_counts_pair("AX", "AA", [0, 2], [0, 2])
  let config = @src.AlignmentCountsConfig::with_matrix(
    @src.blosum62(),
    wildcard=Some('X'),
  )
  let counts = alignment.alignment_counts(config~) catch {
    _ => abort("matrix wildcard should skip X")
  }
  assert_eq(counts.aligned, 2)
  assert_eq(counts.identities, 1)
  assert_eq(counts.positives, Some(1))
  assert_eq(counts.substitution_score, Some(4.0))
}

///|
test "alignment counts: match mismatch substitution score" {
  let alignment = alignment_counts_pair("ACGT", "AGGT", [0, 4], [0, 4])
  let config = @src.AlignmentCountsConfig::create(
    match_score=Some(2.0),
    mismatch_score=Some(-1.0),
  ) catch {
    _ => abort("valid score config should build")
  }
  let counts = alignment.alignment_counts(config~) catch {
    _ => abort("match scoring should succeed")
  }
  assert_eq(counts.substitution_score, Some(5.0))
  assert_eq(counts.gap_score, None)
  assert_eq(counts.score, None)
}

///|
test "alignment counts: matrix takes precedence over match score" {
  let alignment = alignment_counts_pair("AA", "AA", [0, 2], [0, 2])
  let config = @src.AlignmentCountsConfig::create(
    substitution_matrix=Some(@src.blosum62()),
    match_score=Some(100.0),
    mismatch_score=Some(-100.0),
  ) catch {
    _ => abort("combined config should build")
  }
  let counts = alignment.alignment_counts(config~) catch {
    _ => abort("matrix precedence should count")
  }
  assert_eq(counts.substitution_score, Some(8.0))
}

///|
test "alignment counts: left insertion" {
  let counts = alignment_counts_basic("AAAAA", "AAAAA", [2, 2, 5], [0, 2, 5])
  assert_eq(counts.open_left_insertions, 1)
  assert_eq(counts.extend_left_insertions, 1)
  assert_eq(counts.left_insertions(), 2)
  assert_eq(counts.aligned, 3)
}

///|
test "alignment counts: left deletion" {
  let counts = alignment_counts_basic("AAAAA", "AAAAA", [0, 2, 5], [2, 2, 5])
  assert_eq(counts.open_left_deletions, 1)
  assert_eq(counts.extend_left_deletions, 1)
  assert_eq(counts.left_deletions(), 2)
}

///|
test "alignment counts: internal insertion" {
  let counts = alignment_counts_basic("AAAA", "AAAAAA", [0, 2, 2, 4], [
    0, 2, 4, 6,
  ])
  assert_eq(counts.open_internal_insertions, 1)
  assert_eq(counts.extend_internal_insertions, 1)
  assert_eq(counts.internal_insertions(), 2)
}

///|
test "alignment counts: internal deletion" {
  let counts = alignment_counts_basic("AAAAAA", "AAAA", [0, 2, 4, 6], [
    0, 2, 2, 4,
  ])
  assert_eq(counts.open_internal_deletions, 1)
  assert_eq(counts.extend_internal_deletions, 1)
  assert_eq(counts.internal_deletions(), 2)
}

///|
test "alignment counts: right insertion" {
  let counts = alignment_counts_basic("AAA", "AAAAA", [0, 3, 3], [0, 3, 5])
  assert_eq(counts.open_right_insertions, 1)
  assert_eq(counts.extend_right_insertions, 1)
  assert_eq(counts.right_insertions(), 2)
}

///|
test "alignment counts: right deletion" {
  let counts = alignment_counts_basic("AAAAA", "AAA", [0, 3, 5], [0, 3, 3])
  assert_eq(counts.open_right_deletions, 1)
  assert_eq(counts.extend_right_deletions, 1)
  assert_eq(counts.right_deletions(), 2)
}

///|
test "alignment counts: one-letter gap has no extension" {
  let counts = alignment_counts_basic("AAA", "AAAA", [0, 2, 2, 3], [0, 2, 3, 4])
  assert_eq(counts.open_internal_insertions, 1)
  assert_eq(counts.extend_internal_insertions, 0)
}

///|
test "alignment counts: four-letter gap has three extensions" {
  let counts = alignment_counts_basic("AAAA", "AAAAAAAA", [0, 2, 2, 4], [
    0, 2, 6, 8,
  ])
  assert_eq(counts.open_internal_insertions, 1)
  assert_eq(counts.extend_internal_insertions, 3)
}

///|
test "alignment counts: split consecutive insertion remains one gap" {
  let counts = alignment_counts_basic("AA", "AAAAA", [0, 0, 0, 2], [0, 1, 3, 5])
  assert_eq(counts.open_left_insertions, 1)
  assert_eq(counts.extend_left_insertions, 2)
}

///|
test "alignment counts: split consecutive deletion remains one gap" {
  let counts = alignment_counts_basic("AAAAA", "AA", [0, 1, 3, 5], [0, 0, 0, 2])
  assert_eq(counts.open_left_deletions, 1)
  assert_eq(counts.extend_left_deletions, 2)
}

///|
test "alignment counts: diagonal resets insertion path" {
  let counts = alignment_counts_basic("AAAA", "AAAAAA", [0, 1, 1, 3, 3, 4], [
    0, 1, 2, 4, 5, 6,
  ])
  assert_eq(counts.open_internal_insertions, 2)
  assert_eq(counts.extend_internal_insertions, 0)
}

///|
test "alignment counts: changing gap direction opens both gaps" {
  let counts = alignment_counts_basic("AAAA", "AAAA", [0, 0, 2, 4], [0, 2, 2, 4])
  assert_eq(counts.open_left_insertions, 1)
  assert_eq(counts.open_internal_deletions, 1)
  assert_eq(counts.open_gaps(), 2)
}

///|
test "alignment counts: left classification precedes right" {
  let counts = alignment_counts_basic("", "", [], [])
  assert_eq(counts.gaps(), 0)
  let alignment = @src.coordinate_pairwise_alignment_with_lengths(
    "target",
    0,
    "query",
    2,
    [0, 0],
    [0, 2],
  ) catch {
    _ => abort("gap-only alignment should build")
  }
  let gap_counts = alignment.alignment_counts() catch {
    _ => abort("gap-only alignment should count")
  }
  assert_eq(gap_counts.left_insertions(), 2)
  assert_eq(gap_counts.right_insertions(), 0)
}

///|
test "alignment counts: aggregate insertion properties" {
  let counts = alignment_counts_basic("AAAA", "AAAAAAAA", [0, 0, 2, 2, 4, 4], [
    0, 1, 3, 5, 7, 8,
  ])
  assert_eq(counts.left_insertions(), 1)
  assert_eq(counts.internal_insertions(), 2)
  assert_eq(counts.right_insertions(), 1)
  assert_eq(counts.insertions(), 4)
  assert_eq(counts.open_insertions(), 3)
  assert_eq(counts.extend_insertions(), 1)
}

///|
test "alignment counts: aggregate deletion properties" {
  let counts = alignment_counts_basic("AAAAAAAA", "AAAA", [0, 1, 3, 5, 7, 8], [
    0, 0, 2, 2, 4, 4,
  ])
  assert_eq(counts.left_deletions(), 1)
  assert_eq(counts.internal_deletions(), 2)
  assert_eq(counts.right_deletions(), 1)
  assert_eq(counts.deletions(), 4)
  assert_eq(counts.open_deletions(), 3)
  assert_eq(counts.extend_deletions(), 1)
}

///|
test "alignment counts: aggregate side properties" {
  let counts = alignment_counts_basic("AAAA", "AAAA", [0, 0, 2, 4, 4], [
    0, 1, 3, 3, 4,
  ])
  assert_eq(counts.left_gaps(), 1)
  assert_eq(counts.internal_gaps(), 2)
  assert_eq(counts.right_gaps(), 1)
  assert_eq(counts.gaps(), 4)
  assert_eq(counts.open_gaps(), 3)
  assert_eq(counts.extend_gaps(), 1)
}

///|
test "alignment counts: affine complete score" {
  let alignment = alignment_counts_pair("GAACT", "GAT", [0, 2, 4, 5], [
    0, 2, 2, 3,
  ])
  let config = @src.AlignmentCountsConfig::scored(
    match_score=1.0,
    mismatch_score=0.0,
    gap_scores=alignment_counts_affine(-1.0, -1.0),
  ) catch {
    _ => abort("scored config should build")
  }
  let counts = alignment.alignment_counts(config~) catch {
    _ => abort("scored counts should succeed")
  }
  assert_eq(counts.substitution_score, Some(3.0))
  assert_eq(counts.gap_score, Some(-2.0))
  assert_eq(counts.score, Some(1.0))
}

///|
test "alignment counts: affine open extend score" {
  let alignment = alignment_counts_pair("GAACT", "GAT", [0, 2, 4, 5], [
    0, 2, 2, 3,
  ])
  let config = @src.AlignmentCountsConfig::scored(
    match_score=2.0,
    mismatch_score=-1.0,
    gap_scores=alignment_counts_affine(-5.0, -1.0),
  ) catch {
    _ => abort("scored config should build")
  }
  let counts = alignment.alignment_counts(config~) catch {
    _ => abort("scored counts should succeed")
  }
  assert_eq(counts.substitution_score, Some(6.0))
  assert_eq(counts.gap_score, Some(-6.0))
  assert_eq(counts.score, Some(0.0))
}

///|
test "alignment counts: directional insertion score" {
  let scores = @src.AlignmentGapScores::directional(
    insertion_open_score=-3.0,
    insertion_extend_score=-0.5,
    deletion_open_score=-8.0,
    deletion_extend_score=-2.0,
  ) catch {
    _ => abort("directional scores should build")
  }
  let config = @src.AlignmentCountsConfig::scored(
    match_score=0.0,
    mismatch_score=0.0,
    gap_scores=scores,
  ) catch {
    _ => abort("scored config should build")
  }
  let counts = alignment_counts_pair("AAA", "AAAAA", [0, 3, 3], [0, 3, 5]).alignment_counts(
    config~,
  ) catch {
    _ => abort("directional score should count")
  }
  assert_eq(counts.gap_score, Some(-3.5))
}

///|
test "alignment counts: directional deletion score" {
  let scores = @src.AlignmentGapScores::directional(
    insertion_open_score=-3.0,
    insertion_extend_score=-0.5,
    deletion_open_score=-8.0,
    deletion_extend_score=-2.0,
  ) catch {
    _ => abort("directional scores should build")
  }
  let config = @src.AlignmentCountsConfig::scored(
    match_score=0.0,
    mismatch_score=0.0,
    gap_scores=scores,
  ) catch {
    _ => abort("scored config should build")
  }
  let counts = alignment_counts_pair("AAAAA", "AAA", [0, 3, 5], [0, 3, 3]).alignment_counts(
    config~,
  ) catch {
    _ => abort("directional score should count")
  }
  assert_eq(counts.gap_score, Some(-10.0))
}

///|
test "alignment counts: all twelve gap scores are applied" {
  let scores = @src.AlignmentGapScores::create(
    open_left_insertion_score=-1.0,
    extend_left_insertion_score=-2.0,
    open_left_deletion_score=-3.0,
    extend_left_deletion_score=-4.0,
    open_internal_insertion_score=-5.0,
    extend_internal_insertion_score=-6.0,
    open_internal_deletion_score=-7.0,
    extend_internal_deletion_score=-8.0,
    open_right_insertion_score=-9.0,
    extend_right_insertion_score=-10.0,
    open_right_deletion_score=-11.0,
    extend_right_deletion_score=-12.0,
  ) catch {
    _ => abort("position-specific scores should build")
  }
  let config = @src.AlignmentCountsConfig::scored(
    match_score=0.0,
    mismatch_score=0.0,
    gap_scores=scores,
  ) catch {
    _ => abort("scored config should build")
  }
  let alignment = alignment_counts_pair("AAAA", "AAAAA", [0, 0, 2, 4, 4], [
    0, 2, 4, 4, 5,
  ])
  let counts = alignment.alignment_counts(config~) catch {
    _ => abort("position-specific scoring should succeed")
  }
  assert_eq(counts.open_left_insertions, 1)
  assert_eq(counts.extend_left_insertions, 1)
  assert_eq(counts.open_internal_deletions, 1)
  assert_eq(counts.extend_internal_deletions, 1)
  assert_eq(counts.open_right_insertions, 1)
  assert_eq(counts.gap_score, Some(-27.0))
}

///|
test "alignment counts: matrix and gap score total" {
  let config = @src.AlignmentCountsConfig::scored_matrix(
    matrix=@src.blosum62(),
    gap_scores=alignment_counts_affine(-5.0, -1.0),
  )
  let alignment = alignment_counts_pair("AAAA", "AAA", [0, 2, 3, 4], [
    0, 2, 2, 3,
  ])
  let counts = alignment.alignment_counts(config~) catch {
    _ => abort("matrix gap scoring should succeed")
  }
  assert_eq(counts.substitution_score, Some(12.0))
  assert_eq(counts.gap_score, Some(-5.0))
  assert_eq(counts.score, Some(7.0))
  assert_eq(counts.positives, Some(3))
}

///|
test "alignment counts: reverse query identities" {
  let counts = alignment_counts_basic("GAACT", "ATC", [0, 2, 4, 5], [3, 1, 1, 0])
  assert_eq(counts.aligned, 3)
  assert_eq(counts.identities, 3)
  assert_eq(counts.mismatches, 0)
  assert_eq(counts.internal_deletions(), 2)
}

///|
test "alignment counts: reverse query mismatch" {
  let counts = alignment_counts_basic("GACT", "AATC", [0, 4], [4, 0])
  assert_eq(counts.aligned, 4)
  assert_eq(counts.identities, 3)
  assert_eq(counts.mismatches, 1)
}

///|
test "alignment counts: reverse matrix score" {
  let alignment = alignment_counts_pair("AAAA", "TTTT", [0, 4], [4, 0])
  let counts = alignment.alignment_counts(
    config=@src.AlignmentCountsConfig::with_matrix(@src.blosum62()),
  ) catch {
    _ => abort("reverse matrix scoring should succeed")
  }
  assert_eq(counts.identities, 4)
  assert_eq(counts.substitution_score, Some(16.0))
  assert_eq(counts.positives, Some(4))
}

///|
test "alignment counts: reverse right gap classification" {
  let counts = alignment_counts_basic("AAAAA", "TTT", [0, 3, 5], [3, 0, 0])
  assert_eq(counts.right_deletions(), 2)
  assert_eq(counts.open_right_deletions, 1)
}

///|
test "alignment counts: reverse internal gap classification" {
  let counts = alignment_counts_basic("AAAAAA", "TTTT", [0, 2, 4, 6], [
    4, 2, 2, 0,
  ])
  assert_eq(counts.internal_deletions(), 2)
  assert_eq(counts.open_internal_deletions, 1)
}

///|
test "alignment counts: empty path" {
  let alignment = @src.coordinate_pairwise_alignment_with_lengths(
    "target",
    10,
    "query",
    20,
    [],
    [],
  ) catch {
    _ => abort("empty alignment should build")
  }
  let counts = alignment.alignment_counts() catch {
    _ => abort("empty alignment should count")
  }
  assert_eq(counts.aligned, 0)
  assert_eq(counts.gaps(), 0)
  assert_eq(
    counts.summary(),
    "AlignmentCounts(aligned=0, identities=0, mismatches=0, gaps=0)",
  )
}

///|
test "alignment counts: unknown sequence still counts coordinates" {
  let alignment = @src.coordinate_pairwise_alignment_with_lengths(
    "target",
    5,
    "query",
    3,
    [0, 2, 4, 5],
    [0, 2, 2, 3],
  ) catch {
    _ => abort("length-only alignment should build")
  }
  let counts = alignment.alignment_counts() catch {
    _ => abort("length-only alignment should count")
  }
  assert_eq(counts.aligned, 3)
  assert_eq(counts.identities, 0)
  assert_eq(counts.mismatches, 0)
  assert_eq(counts.internal_deletions(), 2)
}

///|
test "alignment counts: unknown sequence scored composition is zero" {
  let alignment = @src.coordinate_pairwise_alignment_with_lengths(
    "target",
    4,
    "query",
    4,
    [0, 4],
    [0, 4],
  ) catch {
    _ => abort("length-only alignment should build")
  }
  let config = @src.AlignmentCountsConfig::create(
    match_score=Some(2.0),
    mismatch_score=Some(-1.0),
  ) catch {
    _ => abort("score config should build")
  }
  let counts = alignment.alignment_counts(config~) catch {
    _ => abort("length-only scoring should count")
  }
  assert_eq(counts.aligned, 4)
  assert_eq(counts.substitution_score, Some(0.0))
}

///|
test "alignment counts: rejects unequal aligned steps" {
  let alignment = alignment_counts_pair("AAA", "AA", [0, 3], [0, 2])
  let failed = try {
    ignore(alignment.alignment_counts())
    false
  } catch {
    AlignmentCountsError(_) => true
  }
  assert_true(failed)
}

///|
test "alignment counts: config rejects lone match score" {
  let failed = try {
    ignore(@src.AlignmentCountsConfig::create(match_score=Some(1.0)))
    false
  } catch {
    AlignmentCountsError(_) => true
  }
  assert_true(failed)
}

///|
test "alignment counts: config rejects lone mismatch score" {
  let failed = try {
    ignore(@src.AlignmentCountsConfig::create(mismatch_score=Some(-1.0)))
    false
  } catch {
    AlignmentCountsError(_) => true
  }
  assert_true(failed)
}

///|
test "alignment counts: affine rejects non-finite open score" {
  let failed = try {
    ignore(@src.AlignmentGapScores::affine(1.0e301, -1.0))
    false
  } catch {
    AlignmentCountsError(_) => true
  }
  assert_true(failed)
}

///|
test "alignment counts: directional rejects non-finite extend score" {
  let failed = try {
    ignore(
      @src.AlignmentGapScores::directional(
        insertion_open_score=-1.0,
        insertion_extend_score=-1.0,
        deletion_open_score=-1.0,
        deletion_extend_score=1.0e301,
      ),
    )
    false
  } catch {
    AlignmentCountsError(_) => true
  }
  assert_true(failed)
}

///|
test "alignment counts: legacy basic counts remain unchanged" {
  let alignment = alignment_counts_pair("AAAAAA", "AAAA", [0, 2, 4, 6], [
    0, 2, 2, 4,
  ])
  let legacy = alignment.counts()
  let detailed = alignment.alignment_counts() catch {
    _ => abort("detailed counts should succeed")
  }
  assert_eq(legacy.aligned, detailed.aligned)
  assert_eq(legacy.query_gap_bases, detailed.deletions())
  assert_eq(legacy.query_gap_events, detailed.open_deletions())
}

///|
test "alignment counts: summary reports core statistics" {
  let counts = alignment_counts_basic("GAACT", "GAT", [0, 2, 4, 5], [0, 2, 2, 3])
  assert_eq(
    counts.summary(),
    "AlignmentCounts(aligned=3, identities=3, mismatches=0, gaps=2)",
  )
}

///|
test "alignment counts: MSA official pairwise sums" {
  let alignment = @src.coordinate_multiple_alignment(
    ["target", "query1", "query2"],
    ["CGGTTTTT", "AG-TTT--", "AGGTTT--"],
  ) catch {
    _ => abort("official MSA should build")
  }
  let counts = alignment.alignment_counts() catch {
    _ => abort("MSA counts should succeed")
  }
  assert_eq(counts.aligned, 16)
  assert_eq(counts.identities, 14)
  assert_eq(counts.mismatches, 2)
  assert_eq(counts.insertions(), 1)
  assert_eq(counts.deletions(), 5)
  assert_eq(counts.left_gaps(), 0)
  assert_eq(counts.internal_gaps(), 2)
  assert_eq(counts.right_gaps(), 4)
}

///|
test "alignment counts: MSA shared gap columns are ignored" {
  let alignment = @src.coordinate_multiple_alignment(["a", "b", "c"], [
    "AA--CC", "AA--CC", "AATTCC",
  ]) catch {
    _ => abort("MSA with shared gaps should build")
  }
  let counts = alignment.alignment_counts() catch {
    _ => abort("shared gap columns should count")
  }
  assert_eq(counts.aligned, 12)
  assert_eq(counts.identities, 12)
  assert_eq(counts.gaps(), 4)
}

///|
test "alignment counts: MSA matrix sums all unordered pairs" {
  let alignment = @src.coordinate_multiple_alignment(["a", "b", "c"], [
    "AA", "AR", "RR",
  ]) catch {
    _ => abort("protein MSA should build")
  }
  let counts = alignment.alignment_counts(
    config=@src.AlignmentCountsConfig::with_matrix(@src.blosum62()),
  ) catch {
    _ => abort("MSA matrix scoring should succeed")
  }
  assert_eq(counts.aligned, 6)
  assert_eq(counts.identities, 2)
  assert_eq(counts.mismatches, 4)
  assert_eq(counts.positives, Some(2))
  assert_eq(counts.substitution_score, Some(5.0))
}

///|
test "alignment counts: MSA wildcard sums all unordered pairs" {
  let alignment = @src.coordinate_multiple_alignment(["a", "b", "c"], [
    "AN", "AA", "AT",
  ]) catch {
    _ => abort("wildcard MSA should build")
  }
  let counts = alignment.alignment_counts(
    config=@src.AlignmentCountsConfig::with_wildcard('N'),
  ) catch {
    _ => abort("MSA wildcard counting should succeed")
  }
  assert_eq(counts.aligned, 6)
  assert_eq(counts.identities, 3)
  assert_eq(counts.mismatches, 1)
}

///|
test "alignment counts: MSA gap score sums pair gaps" {
  let alignment = @src.coordinate_multiple_alignment(["a", "b", "c"], [
    "A-AA", "AAAA", "A--A",
  ]) catch {
    _ => abort("gapped MSA should build")
  }
  let config = @src.AlignmentCountsConfig::scored(
    match_score=1.0,
    mismatch_score=-1.0,
    gap_scores=alignment_counts_affine(-2.0, -1.0),
  ) catch {
    _ => abort("scored config should build")
  }
  let counts = alignment.alignment_counts(config~) catch {
    _ => abort("MSA gap scoring should succeed")
  }
  assert_eq(counts.gaps(), 4)
  assert_eq(counts.open_gaps(), 3)
  assert_eq(counts.extend_gaps(), 1)
  assert_eq(counts.gap_score, Some(-7.0))
  assert_eq(counts.score, Some(0.0))
}

///|
test "alignment counts: one-row MSA has zero pair counts" {
  let alignment = @src.coordinate_multiple_alignment(["only"], ["ACGT"]) catch {
    _ => abort("one-row MSA should build")
  }
  let counts = alignment.alignment_counts() catch {
    _ => abort("one-row MSA should count")
  }
  assert_eq(counts.aligned, 0)
  assert_eq(counts.gaps(), 0)
}

///|
test "alignment counts: scored one-row MSA has zero scores" {
  let alignment = @src.coordinate_multiple_alignment(["only"], ["ACGT"]) catch {
    _ => abort("one-row MSA should build")
  }
  let config = @src.AlignmentCountsConfig::scored(
    match_score=1.0,
    mismatch_score=-1.0,
    gap_scores=alignment_counts_affine(-2.0, -1.0),
  ) catch {
    _ => abort("scored config should build")
  }
  let counts = alignment.alignment_counts(config~) catch {
    _ => abort("one-row scored MSA should count")
  }
  assert_eq(counts.substitution_score, Some(0.0))
  assert_eq(counts.gap_score, Some(0.0))
  assert_eq(counts.score, Some(0.0))
}
