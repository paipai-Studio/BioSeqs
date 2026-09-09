///|
/// Tests for GenomicScores. Expected values were computed with a Python
/// mirror of the step/linear interpolation and range summaries.

///|
fn gs_test_scores() -> @src.GScores raise @src.SeqIOError {
  let gss = @src.GScores::new("phastCons7way", "hg38")
  gss.set_scores(
    "chr1",
    "default",
    [100, 200, 300, 400],
    [0.1, 0.5, 0.9, 0.4],
  )
  gss
}

///|
fn gs_close(a : Double, b : Double, tol : Double) -> Bool {
  if a != a {
    return b != b
  }
  (a - b).abs() <= tol
}

///|
test "gscores_provider_metadata" {
  let gss = gs_test_scores()
  assert_eq(gss.provider_name(), "phastCons7way")
  assert_eq(gss.organism_name(), "hg38")
  assert_eq(gss.n_stored_positions("chr1"), 4)
  assert_eq(gss.n_stored_positions("chrX"), 0)
}

///|
test "gscores_score_at_constant" {
  let gss = gs_test_scores()
  // exact stored positions
  assert_true(gs_close(gss.score_at("chr1", 100), 0.1, 1.0e-12))
  assert_true(gs_close(gss.score_at("chr1", 400), 0.4, 1.0e-12))
  // between stored positions: step function keeps the previous score
  assert_true(gs_close(gss.score_at("chr1", 150), 0.1, 1.0e-12))
  assert_true(gs_close(gss.score_at("chr1", 350), 0.9, 1.0e-12))
  // before the first stored position: NaN
  let missing = gss.score_at("chr1", 50)
  assert_true(missing != missing)
  // unknown chromosome: NaN
  let absent = gss.score_at("chrX", 100)
  assert_true(absent != absent)
}

///|
test "gscores_score_at_linear" {
  let gss = gs_test_scores()
  gss.set_interpolation(@src.GsInterpolation::GsLinear)
  // mirror values for linear interpolation
  assert_true(gs_close(gss.score_at("chr1", 150), 0.3, 1.0e-12))
  assert_true(gs_close(gss.score_at("chr1", 250), 0.7, 1.0e-12))
  assert_true(gs_close(gss.score_at("chr1", 350), 0.65, 1.0e-12))
  // exact positions unchanged
  assert_true(gs_close(gss.score_at("chr1", 300), 0.9, 1.0e-12))
  // beyond the last stored position: clamps to the last score
  assert_true(gs_close(gss.score_at("chr1", 450), 0.4, 1.0e-12))
  // before the first stored position: NaN
  let missing = gss.score_at("chr1", 50)
  assert_true(missing != missing)
}

///|
test "gscores_range_summaries" {
  let gss = gs_test_scores()
  // closed range 100-400 covers all four stored scores
  assert_true(gs_close(gss.range_score("chr1", 100, 400), 0.475, 1.0e-12))
  assert_true(gs_close(gss.range_score("chr1", 100, 400, summary="min"), 0.1, 1.0e-12))
  assert_true(gs_close(gss.range_score("chr1", 100, 400, summary="max"), 0.9, 1.0e-12))
  assert_true(gs_close(gss.range_score("chr1", 100, 400, summary="median"), 0.45, 1.0e-12))
  // range 100-200 covers only the two stored bases (150 has no data)
  assert_true(gs_close(gss.range_score("chr1", 100, 200), 0.3, 1.0e-12))
  // range without data: NaN
  let empty = gss.range_score("chr1", 1, 10)
  assert_true(empty != empty)
}

///|
test "gscores_multiple_populations" {
  let gss = gs_test_scores()
  gss.set_scores(
    "chr1",
    "phyloP",
    [100, 200, 300, 400],
    [0.9, 0.7, 0.2, 0.2],
  )
  let pops = gss.populations()
  assert_eq(pops.length(), 2)
  assert_true(gs_close(gss.range_score("chr1", 100, 400, population="default"), 0.475, 1.0e-12))
  assert_true(gs_close(gss.range_score("chr1", 100, 400, population="phyloP"), 0.5, 1.0e-12))
  gss.set_default_population("phyloP")
  assert_true(gs_close(gss.score_at("chr1", 200), 0.7, 1.0e-12))
}

///|
test "gscores_multiple_chromosomes_and_unsorted_input" {
  let gss = gs_test_scores()
  gss.set_scores(
    "chr2",
    "default",
    [500, 100, 300],
    [0.3, 0.8, 0.2],
  )
  assert_eq(gss.n_stored_positions("chr2"), 3)
  // unsorted input is sorted: position 100 now holds 0.8
  assert_true(gs_close(gss.score_at("chr2", 100), 0.8, 1.0e-12))
  assert_true(gs_close(gss.score_at("chr2", 300), 0.2, 1.0e-12))
  assert_true(gs_close(gss.score_at("chr2", 400), 0.2, 1.0e-12))
  // extending chr1 with new positions keeps earlier series aligned
  gss.set_scores("chr1", "default", [150], [0.25])
  assert_eq(gss.n_stored_positions("chr1"), 5)
  assert_true(gs_close(gss.score_at("chr1", 150), 0.25, 1.0e-12))
  assert_true(gs_close(gss.score_at("chr1", 200), 0.5, 1.0e-12))
}

///|
test "gscores_length_mismatch_raises" {
  let gss = gs_test_scores()
  let failed = try {
    gss.set_scores("chr1", "bad", [1, 2, 3], [0.1, 0.2])
    false
  } catch {
    _ => true
  }
  assert_true(failed)
}

///|
test "gscores_fetch_batch" {
  let gss = gs_test_scores()
  let fetched = @src.gscores_fetch(
    gss,
    ["chr1", "chr1", "chr1"],
    [100, 100, 250],
    [200, 400, 260],
  )
  assert_eq(fetched.length(), 3)
  assert_true(gs_close(fetched[0], 0.3, 1.0e-12))
  assert_true(gs_close(fetched[1], 0.475, 1.0e-12))
  // 250 and 260 both carry the step score 0.5
  assert_true(gs_close(fetched[2], 0.5, 1.0e-12))
  // mismatched arrays raise
  let failed = try {
    let _ = @src.gscores_fetch(gss, ["chr1"], [1], [2, 3])
    false
  } catch {
    _ => true
  }
  assert_true(failed)
}
