///|
/// Tests for biotite_kmer_selector.mbt — minimizer, syncmer, mincode.
test "biotite_minimizer_select" {
  let selector = @src.BiotiteMinimizerSelector::new(k=3, window_size=2)
  let positions = selector.select("AAAA")
  assert_eq(positions.length(), 1)
  assert_eq(positions[0], 0)
  assert_eq(selector.select("AC").length(), 0)
}

///|
test "biotite_syncmer_select" {
  let selector = @src.BiotiteSyncmerSelector::new(k=3, s=1, t=0)
  let positions = selector.select("AAAA")
  assert_eq(positions.length(), 2)
  assert_eq(positions[0], 0)
  assert_eq(positions[1], 1)
}

///|
test "biotite_mincode_density_and_select" {
  let selector = @src.BiotiteMincodeSelector::new(k=3, density=0.1)
  let density = selector.density()
  assert_true(density >= 0.0)
  assert_true(density <= 1.0)
  let positions = selector.select("ACGTACGT")
  for p in positions {
    assert_true(p >= 0)
    assert_true(p <= 5)
  }
}

///|
test "biotite_selectors_short_sequence" {
  let minimizer = @src.BiotiteMinimizerSelector::new(k=5, window_size=2)
  assert_eq(minimizer.select("ACG").length(), 0)
  let syncmer = @src.BiotiteSyncmerSelector::new(k=5, s=2, t=0)
  assert_eq(syncmer.select("ACG").length(), 0)
  let mincode = @src.BiotiteMincodeSelector::new(k=5, density=0.1)
  assert_eq(mincode.select("ACG").length(), 0)
}
