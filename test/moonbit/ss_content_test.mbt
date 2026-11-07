///|
/// Tests for ss_content.mbt (secondary structure content summary).

///|
test "calculate_ss_content_mixed" {
  // 10 residues: 4 helix (H,G,I count as helix), 3 sheet (E,B),
  // 2 turn (T,S), 1 coil (C).
  let codes = ['H', 'G', 'I', 'H', 'E', 'B', 'E', 'T', 'S', 'C']
  let ss = @src.calculate_ss_content(codes)
  assert_eq(ss.total_residues, 10)
  assert_eq(ss.helix_count, 4)
  assert_eq(ss.sheet_count, 3)
  assert_eq(ss.turn_count, 2)
  assert_eq(ss.coil_count, 1)
  assert_true(
    @moonbit.assert_num_close(
      @src.SSContent::get_helix_fraction(ss),
      0.4,
      1.0e-12,
    ),
  )
  assert_true(
    @moonbit.assert_num_close(
      @src.SSContent::get_sheet_fraction(ss),
      0.3,
      1.0e-12,
    ),
  )
  assert_true(@moonbit.assert_num_close(ss.turn_fraction, 0.2, 1.0e-12))
  assert_true(
    @moonbit.assert_num_close(
      @src.SSContent::get_coil_fraction(ss),
      0.1,
      1.0e-12,
    ),
  )
}

///|
test "calculate_ss_content_empty" {
  // Empty input: all counts 0, fractions 0 (division guarded).
  let ss = @src.calculate_ss_content([])
  assert_eq(ss.total_residues, 0)
  assert_true(@moonbit.assert_num_close(ss.helix_fraction, 0.0, 1.0e-12))
  assert_true(@moonbit.assert_num_close(ss.coil_fraction, 0.0, 1.0e-12))
}

///|
test "calculate_ss_content_all_helix" {
  let ss = @src.calculate_ss_content(['H', 'H', 'G'])
  assert_eq(ss.helix_count, 3)
  assert_true(@moonbit.assert_num_close(ss.helix_fraction, 1.0, 1.0e-12))
  assert_true(@moonbit.assert_num_close(ss.sheet_fraction, 0.0, 1.0e-12))
}

///|
test "count_helix_segments" {
  // 'HH---GGH': G and H are both helix codes, so GGH is ONE contiguous
  // run -> 2 segments.
  assert_eq(
    @src.count_helix_segments(['H', 'H', 'C', 'C', 'C', 'G', 'G', 'H']),
    2,
  )
  // Fully separated runs -> 3 segments.
  assert_eq(
    @src.count_helix_segments(['H', 'H', 'C', 'C', 'C', 'G', 'G', 'C', 'I']),
    3,
  )
  // No helix at all.
  assert_eq(@src.count_helix_segments(['E', 'C', 'T']), 0)
  // Empty input.
  assert_eq(@src.count_helix_segments([]), 0)
}

///|
test "count_strand_segments" {
  // Strand runs: EE-C-BBB -> 2 segments (E and B both count as sheet).
  assert_eq(@src.count_strand_segments(['E', 'E', 'C', 'T', 'B', 'B', 'B']), 2)
  assert_eq(@src.count_strand_segments(['H', 'H']), 0)
}

///|
test "helix_boundaries_runs" {
  // 'C H H C G C' -> one helix run at [1, 2] and a single-residue run [4, 4].
  let b = @src.helix_boundaries(['C', 'H', 'H', 'C', 'G', 'C'])
  assert_eq(b.length(), 2)
  assert_eq(b[0].0, 1)
  assert_eq(b[0].1, 2)
  assert_eq(b[1].0, 4)
  assert_eq(b[1].1, 4)
  // Helix running to the end of the chain closes at the last index.
  let b2 = @src.helix_boundaries(['C', 'H', 'H'])
  assert_eq(b2.length(), 1)
  assert_eq(b2[0].0, 1)
  assert_eq(b2[0].1, 2)
  // No helices -> no boundaries.
  assert_eq(@src.helix_boundaries(['C', 'E']).length(), 0)
}

///|
test "strand_boundaries_runs" {
  // 'E C E E C' -> runs [0,0] and [2,3].
  let b = @src.strand_boundaries(['E', 'C', 'E', 'E', 'C'])
  assert_eq(b.length(), 2)
  assert_eq(b[0].0, 0)
  assert_eq(b[0].1, 0)
  assert_eq(b[1].0, 2)
  assert_eq(b[1].1, 3)
  // Empty input.
  assert_eq(@src.strand_boundaries([]).length(), 0)
}
