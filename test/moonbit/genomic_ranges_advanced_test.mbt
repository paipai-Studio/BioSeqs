///|
/// GenomicRanges Advanced tests - Tiling, windowing, and range operations.

///|
test "tile_by_width_correct_number_of_tiles" {
  let gr = @src.granges(["chr1"], [(1, 10)], [@src.strand_plus()])
  let tiled = @src.granges_tile(gr, width=5)
  assert_eq(tiled.seqnames.length(), 2)
  assert_eq(tiled.starts[0], 1)
  assert_eq(tiled.ends[0], 5)
  assert_eq(tiled.starts[1], 6)
  assert_eq(tiled.ends[1], 10)
}

///|
test "tile_by_n_correct_number_of_tiles" {
  let gr = @src.granges(["chr1"], [(1, 10)], [@src.strand_plus()])
  let tiled = @src.granges_tile(gr, n=2)
  assert_eq(tiled.seqnames.length(), 2)
  assert_eq(tiled.starts[0], 1)
  assert_eq(tiled.ends[0], 5)
  assert_eq(tiled.starts[1], 6)
  assert_eq(tiled.ends[1], 10)
}

///|
test "tile_width_larger_than_range_single_tile" {
  let gr = @src.granges(["chr1"], [(1, 10)], [@src.strand_plus()])
  let tiled = @src.granges_tile(gr, width=20)
  assert_eq(tiled.seqnames.length(), 1)
  assert_eq(tiled.starts[0], 1)
  assert_eq(tiled.ends[0], 10)
}

///|
test "tile_default_returns_original" {
  let gr = @src.granges(
    ["chr1", "chr1"], [(1, 10), (20, 30)], [@src.strand_plus(), @src.strand_minus()],
  )
  let tiled = @src.granges_tile(gr)
  assert_eq(tiled.seqnames.length(), 2)
  assert_eq(tiled.starts[0], 1)
  assert_eq(tiled.ends[0], 10)
  assert_eq(tiled.starts[1], 20)
  assert_eq(tiled.ends[1], 30)
}

///|
test "sliding_windows_correct_overlap" {
  let gr = @src.granges(["chr1"], [(1, 10)], [@src.strand_plus()])
  let windows = @src.granges_sliding_windows(gr, 5, 2)
  assert_true(windows.seqnames.length() >= 2)
  // Consecutive windows overlap (end of first >= start of second)
  assert_true(windows.ends[0] >= windows.starts[1])
}

///|
test "sliding_windows_step_greater_than_width_no_overlap" {
  let gr = @src.granges(["chr1"], [(1, 20)], [@src.strand_plus()])
  let windows = @src.granges_sliding_windows(gr, 5, 10)
  assert_true(windows.seqnames.length() >= 2)
  // No overlap: end of first < start of second
  assert_true(windows.ends[0] < windows.starts[1])
}

///|
test "tile_genome_covers_entire_genome" {
  let tiled = @src.tile_genome([("chr1", 10)], 3)
  assert_eq(tiled.seqnames.length(), 4)
  assert_eq(tiled.starts[0], 1)
  assert_eq(tiled.ends[0], 3)
  assert_eq(tiled.starts[3], 10)
  assert_eq(tiled.ends[3], 10)
}

///|
test "tile_genome_multiple_chromosomes" {
  let tiled = @src.tile_genome([("chr1", 5), ("chr2", 5)], 2)
  assert_eq(tiled.seqnames.length(), 6)
  assert_eq(tiled.seqnames[0], "chr1")
  assert_eq(tiled.seqnames[3], "chr2")
}

///|
test "tile_by_ranges_correct_number_per_range" {
  let gr = @src.granges(["chr1"], [(1, 10)], [@src.strand_plus()])
  let tiled = @src.tile_by_ranges(gr, 2)
  assert_eq(tiled.seqnames.length(), 2)
  assert_eq(tiled.starts[0], 1)
  assert_eq(tiled.ends[0], 5)
  assert_eq(tiled.starts[1], 6)
  assert_eq(tiled.ends[1], 10)
}

///|
test "tile_by_ranges_distributes_remainder" {
  let gr = @src.granges(["chr1"], [(1, 11)], [@src.strand_plus()])
  let tiled = @src.tile_by_ranges(gr, 3)
  assert_eq(tiled.seqnames.length(), 3)
  // range_width=11, tile_width=3, remainder=2
  // First two tiles get +1 (width 4), last gets width 3
  assert_eq(tiled.starts[0], 1)
  assert_eq(tiled.ends[0], 4)
  assert_eq(tiled.starts[1], 5)
  assert_eq(tiled.ends[1], 8)
  assert_eq(tiled.starts[2], 9)
  assert_eq(tiled.ends[2], 11)
}

///|
test "coverage_by_window_correct_counts" {
  let gr = @src.granges(
    ["chr1", "chr1"], [(1, 5), (3, 8)], [@src.strand_plus(), @src.strand_minus()],
  )
  let cov = @src.coverage_by_window(gr, "chr1", 10, 5, 5)
  assert_eq(cov.length(), 2)
  assert_eq(cov[0].0, 1)
  assert_eq(cov[0].1, 5)
  assert_eq(cov[0].2, 2)
  assert_eq(cov[1].0, 6)
  assert_eq(cov[1].1, 10)
  assert_eq(cov[1].2, 1)
}

///|
test "coverage_by_window_no_overlapping_ranges" {
  let gr = @src.granges(["chr1"], [(1, 5)], [@src.strand_plus()])
  let cov = @src.coverage_by_window(gr, "chr2", 10, 5, 5)
  assert_eq(cov.length(), 2)
  assert_eq(cov[0].2, 0)
  assert_eq(cov[1].2, 0)
}

///|
test "bin_genome_correct_bin_count" {
  let gr = @src.granges(["chr1"], [(1, 5)], [@src.strand_plus()])
  let bins = @src.bin_genome(gr, [("chr1", 10)], 5)
  assert_eq(bins.length(), 2)
  assert_eq(bins[0].0, "chr1")
  assert_eq(bins[0].1, 1)
  assert_eq(bins[0].2, 5)
  assert_eq(bins[0].3, 1)
  assert_eq(bins[1].3, 0)
}

///|
test "bin_genome_density_values" {
  let gr = @src.granges(
    ["chr1", "chr1"], [(1, 5), (3, 8)], [@src.strand_plus(), @src.strand_minus()],
  )
  let bins = @src.bin_genome(gr, [("chr1", 10)], 5)
  assert_eq(bins.length(), 2)
  // First bin (1-5) overlaps 2 ranges, density = 2/2 = 1.0
  assert_true(bins[0].4 == 1.0)
  // Second bin (6-10) overlaps 1 range, density = 1/2 = 0.5
  assert_true(bins[1].4 == 0.5)
}

///|
test "granges_resize_advanced_fix_start" {
  let gr = @src.granges(["chr1"], [(10, 20)], [@src.strand_plus()])
  let resized = @src.granges_resize_advanced(gr, 15, "start")
  assert_eq(resized.starts[0], 10)
  assert_eq(resized.ends[0], 24)
  assert_eq(resized.widths[0], 15)
}

///|
test "granges_resize_advanced_fix_end" {
  let gr = @src.granges(["chr1"], [(10, 20)], [@src.strand_plus()])
  let resized = @src.granges_resize_advanced(gr, 15, "end")
  assert_eq(resized.starts[0], 6)
  assert_eq(resized.ends[0], 20)
  assert_eq(resized.widths[0], 15)
}

///|
test "granges_resize_advanced_fix_center" {
  let gr = @src.granges(["chr1"], [(10, 20)], [@src.strand_plus()])
  let resized = @src.granges_resize_advanced(gr, 15, "center")
  // center = (10+20)/2 = 15, new_start = 15 - 15/2 = 15 - 7 = 8
  assert_eq(resized.starts[0], 8)
  assert_eq(resized.ends[0], 22)
  assert_eq(resized.widths[0], 15)
}

///|
test "granges_promoters_on_plus_strand" {
  let gr = @src.granges_single("chr1", 100, 200, @src.strand_plus())
  let promoters = @src.granges_promoters(gr, 50, 20)
  assert_eq(promoters.starts[0], 50)
  assert_eq(promoters.ends[0], 119)
}

///|
test "granges_promoters_on_minus_strand" {
  let gr = @src.granges_single("chr1", 100, 200, @src.strand_minus())
  let promoters = @src.granges_promoters(gr, 50, 20)
  assert_eq(promoters.starts[0], 181)
  assert_eq(promoters.ends[0], 250)
}

///|
test "granges_gaps_finds_correct_gaps" {
  let gr = @src.granges(
    ["chr1", "chr1"], [(1, 10), (20, 30)], [@src.strand_plus(), @src.strand_minus()],
  )
  let gaps = @src.granges_gaps(gr, end_limit=100)
  assert_eq(gaps.seqnames.length(), 2)
  assert_eq(gaps.starts[0], 11)
  assert_eq(gaps.ends[0], 19)
  assert_eq(gaps.starts[1], 31)
  assert_eq(gaps.ends[1], 100)
}

///|
test "granges_gaps_with_end_limit" {
  let gr = @src.granges(
    ["chr1", "chr1"], [(1, 10), (20, 30)], [@src.strand_plus(), @src.strand_minus()],
  )
  let gaps = @src.granges_gaps(gr, end_limit=50)
  assert_eq(gaps.seqnames.length(), 2)
  assert_eq(gaps.starts[0], 11)
  assert_eq(gaps.ends[0], 19)
  assert_eq(gaps.starts[1], 31)
  assert_eq(gaps.ends[1], 50)
}

///|
test "granges_subtract_removes_overlapping_region" {
  let gr1 = @src.granges(["chr1"], [(1, 30)], [@src.strand_plus()])
  let gr2 = @src.granges(["chr1"], [(10, 20)], [@src.strand_minus()])
  let result = @src.granges_subtract(gr1, gr2)
  assert_eq(result.seqnames.length(), 2)
  assert_eq(result.starts[0], 1)
  assert_eq(result.ends[0], 9)
  assert_eq(result.starts[1], 21)
  assert_eq(result.ends[1], 30)
}

///|
test "granges_subtract_no_overlap_returns_original" {
  let gr1 = @src.granges(["chr1"], [(1, 10)], [@src.strand_plus()])
  let gr2 = @src.granges(["chr1"], [(20, 30)], [@src.strand_minus()])
  let result = @src.granges_subtract(gr1, gr2)
  assert_eq(result.seqnames.length(), 1)
  assert_eq(result.starts[0], 1)
  assert_eq(result.ends[0], 10)
}

///|
test "granges_subtract_splits_range" {
  let gr1 = @src.granges(["chr1"], [(1, 100)], [@src.strand_plus()])
  let gr2 = @src.granges(["chr1"], [(30, 50)], [@src.strand_minus()])
  let result = @src.granges_subtract(gr1, gr2)
  assert_eq(result.seqnames.length(), 2)
  assert_eq(result.starts[0], 1)
  assert_eq(result.ends[0], 29)
  assert_eq(result.starts[1], 51)
  assert_eq(result.ends[1], 100)
}
