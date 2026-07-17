///|
test "genomic_ranges_granges_construction" {
  let seqnames = ["chr1", "chr1", "chr2"]
  let ranges = [(1, 10), (5, 15), (20, 30)]
  let strands = [@src.strand_plus(), @src.strand_minus(), @src.strand_star()]
  let gr = @src.granges(seqnames, ranges, strands)

  assert_eq(gr.seqnames[0], "chr1")
  assert_eq(gr.starts[0], 1)
  assert_eq(gr.ends[0], 10)
  assert_eq(@src.granges_length(gr), 3)
}

///|
test "genomic_ranges_granges_single" {
  let gr = @src.granges_single("chr1", 10, 20, @src.strand_plus())

  assert_eq(gr.seqnames[0], "chr1")
  assert_eq(gr.starts[0], 10)
  assert_eq(gr.ends[0], 20)
  assert_eq(gr.widths[0], 11)
}

///|
test "genomic_ranges_shift" {
  let seqnames = ["chr1", "chr2"]
  let ranges = [(1, 10), (20, 30)]
  let strands = [@src.strand_plus(), @src.strand_star()]
  let gr = @src.granges(seqnames, ranges, strands)

  let shifted = @src.granges_shift(gr, 5)

  assert_eq(shifted.starts[0], 6)
  assert_eq(shifted.ends[0], 15)
  assert_eq(shifted.starts[1], 25)
  assert_eq(shifted.ends[1], 35)
}

///|
test "genomic_ranges_shift_negative" {
  let seqnames = ["chr1"]
  let ranges = [(10, 20)]
  let strands = [@src.strand_minus()]
  let gr = @src.granges(seqnames, ranges, strands)

  let shifted = @src.granges_shift(gr, -3)

  assert_eq(shifted.starts[0], 7)
  assert_eq(shifted.ends[0], 17)
}

///|
test "genomic_ranges_narrow" {
  let seqnames = ["chr1"]
  let ranges = [(1, 20)]
  let strands = [@src.strand_plus()]
  let gr = @src.granges(seqnames, ranges, strands)

  let narrowed = @src.granges_narrow(gr, 5, 15)

  assert_eq(narrowed.starts[0], 5)
  assert_eq(narrowed.ends[0], 15)
  assert_eq(narrowed.widths[0], 11)
}

///|
test "genomic_ranges_resize_start" {
  let seqnames = ["chr1"]
  let ranges = [(10, 20)]
  let strands = [@src.strand_plus()]
  let gr = @src.granges(seqnames, ranges, strands)

  let resized = @src.granges_resize(gr, 15, "start")

  assert_eq(resized.starts[0], 10)
  assert_eq(resized.ends[0], 24)
  assert_eq(resized.widths[0], 15)
}

///|
test "genomic_ranges_resize_end" {
  let seqnames = ["chr1"]
  let ranges = [(10, 20)]
  let strands = [@src.strand_plus()]
  let gr = @src.granges(seqnames, ranges, strands)

  let resized = @src.granges_resize(gr, 15, "end")

  assert_eq(resized.starts[0], 6)
  assert_eq(resized.ends[0], 20)
  assert_eq(resized.widths[0], 15)
}

///|
test "genomic_ranges_resize_center" {
  let seqnames = ["chr1"]
  let ranges = [(10, 20)]
  let strands = [@src.strand_plus()]
  let gr = @src.granges(seqnames, ranges, strands)

  let resized = @src.granges_resize(gr, 15, "center")

  assert_eq(resized.starts[0], 8)
  assert_eq(resized.ends[0], 22)
  assert_eq(resized.widths[0], 15)
}

///|
test "genomic_ranges_flank_start" {
  let seqnames = ["chr1"]
  let ranges = [(10, 20)]
  let strands = [@src.strand_plus()]
  let gr = @src.granges(seqnames, ranges, strands)

  let flanked = @src.granges_flank(gr, 5, true, false)

  assert_eq(flanked.starts[0], 5)
  assert_eq(flanked.ends[0], 9)
  assert_eq(flanked.widths[0], 5)
}

///|
test "genomic_ranges_flank_end" {
  let seqnames = ["chr1"]
  let ranges = [(10, 20)]
  let strands = [@src.strand_plus()]
  let gr = @src.granges(seqnames, ranges, strands)

  let flanked = @src.granges_flank(gr, 5, false, false)

  assert_eq(flanked.starts[0], 21)
  assert_eq(flanked.ends[0], 25)
  assert_eq(flanked.widths[0], 5)
}

///|
test "genomic_ranges_flank_both" {
  let seqnames = ["chr1"]
  let ranges = [(10, 20)]
  let strands = [@src.strand_plus()]
  let gr = @src.granges(seqnames, ranges, strands)

  let flanked = @src.granges_flank(gr, 5, true, true)

  assert_eq(flanked.starts[0], 5)
  assert_eq(flanked.ends[0], 25)
  assert_eq(flanked.widths[0], 21)
}

///|
test "genomic_ranges_reduce" {
  let seqnames = ["chr1", "chr1", "chr1", "chr2"]
  let ranges = [(1, 10), (5, 15), (20, 30), (1, 5)]
  let strands = [
    @src.strand_plus(),
    @src.strand_minus(),
    @src.strand_star(),
    @src.strand_plus(),
  ]
  let gr = @src.granges(seqnames, ranges, strands)

  let reduced = @src.granges_reduce(gr, 0)

  assert_eq(reduced.seqnames.length(), 3)
  assert_eq(reduced.starts[0], 1)
  assert_eq(reduced.ends[0], 15)
  assert_eq(reduced.starts[1], 20)
  assert_eq(reduced.ends[1], 30)
  assert_eq(reduced.starts[2], 1)
  assert_eq(reduced.ends[2], 5)
}

///|
test "genomic_ranges_disjoin" {
  let seqnames = ["chr1", "chr1"]
  let ranges = [(1, 5), (3, 10)]
  let strands = [@src.strand_plus(), @src.strand_minus()]
  let gr = @src.granges(seqnames, ranges, strands)

  let disjoined = @src.granges_disjoin(gr)

  assert_eq(disjoined.seqnames.length(), 3)
  assert_eq(disjoined.starts[0], 1)
  assert_eq(disjoined.ends[0], 2)
  assert_eq(disjoined.starts[1], 3)
  assert_eq(disjoined.ends[1], 5)
  assert_eq(disjoined.starts[2], 6)
  assert_eq(disjoined.ends[2], 10)
}

///|
test "genomic_ranges_union" {
  let gr1 = @src.granges(["chr1", "chr1"], [(1, 10), (20, 30)], [
    @src.strand_plus(),
    @src.strand_minus(),
  ])
  let gr2 = @src.granges(["chr1"], [(5, 25)], [@src.strand_star()])

  let union_result = @src.granges_union(gr1, gr2)

  assert_eq(union_result.seqnames.length(), 1)
  assert_eq(union_result.starts[0], 1)
  assert_eq(union_result.ends[0], 30)
}

///|
test "genomic_ranges_intersect" {
  let gr1 = @src.granges(["chr1", "chr2"], [(1, 20), (10, 30)], [
    @src.strand_plus(),
    @src.strand_star(),
  ])
  let gr2 = @src.granges(["chr1", "chr2"], [(10, 15), (20, 25)], [
    @src.strand_minus(),
    @src.strand_plus(),
  ])

  let intersect_result = @src.granges_intersect(gr1, gr2)

  assert_eq(intersect_result.seqnames.length(), 2)
  assert_eq(intersect_result.starts[0], 10)
  assert_eq(intersect_result.ends[0], 15)
  assert_eq(intersect_result.starts[1], 20)
  assert_eq(intersect_result.ends[1], 25)
}

///|
test "genomic_ranges_setdiff" {
  let gr1 = @src.granges(["chr1"], [(1, 30)], [@src.strand_plus()])
  let gr2 = @src.granges(["chr1"], [(10, 20)], [@src.strand_minus()])

  let diff_result = @src.granges_setdiff(gr1, gr2)

  assert_eq(diff_result.seqnames.length(), 2)
  assert_eq(diff_result.starts[0], 1)
  assert_eq(diff_result.ends[0], 9)
  assert_eq(diff_result.starts[1], 21)
  assert_eq(diff_result.ends[1], 30)
}

///|
test "genomic_ranges_count_overlaps" {
  let query = @src.granges(["chr1", "chr1"], [(1, 5), (3, 8)], [
    @src.strand_plus(),
    @src.strand_minus(),
  ])
  let subject = @src.granges(
    ["chr1", "chr1", "chr2"],
    [(2, 4), (6, 10), (1, 5)],
    [@src.strand_star(), @src.strand_plus(), @src.strand_star()],
  )

  let counts = @src.count_overlaps(query, subject)

  assert_eq(counts.length(), 2)
  assert_eq(counts[0], 1)
  assert_eq(counts[1], 2)
}

///|
test "genomic_ranges_find_overlaps" {
  let query = @src.granges(["chr1", "chr1"], [(1, 5), (3, 8)], [
    @src.strand_plus(),
    @src.strand_minus(),
  ])
  let subject = @src.granges(
    ["chr1", "chr1", "chr2"],
    [(2, 4), (6, 10), (1, 5)],
    [@src.strand_star(), @src.strand_plus(), @src.strand_star()],
  )

  let overlaps = @src.find_overlaps(query, subject)

  assert_eq(overlaps.length(), 3)
}

///|
test "genomic_ranges_distance" {
  let gr1 = @src.granges(
    ["chr1", "chr1", "chr2"],
    [(1, 10), (20, 30), (1, 5)],
    [@src.strand_plus(), @src.strand_minus(), @src.strand_star()],
  )
  let gr2 = @src.granges(["chr1", "chr1"], [(12, 15), (25, 28)], [
    @src.strand_star(),
    @src.strand_plus(),
  ])

  let distances = @src.granges_distance(gr1, gr2)

  assert_eq(distances.length(), 3)
  assert_eq(distances[0].to_int(), 1)
  assert_eq(distances[1].to_int(), 0)
  assert_true(distances[2].is_nan())
}

///|
test "genomic_ranges_nearest" {
  let gr1 = @src.granges(["chr1", "chr1"], [(1, 5), (15, 20)], [
    @src.strand_plus(),
    @src.strand_minus(),
  ])
  let gr2 = @src.granges(
    ["chr1", "chr1", "chr2"],
    [(10, 12), (25, 30), (1, 5)],
    [@src.strand_star(), @src.strand_plus(), @src.strand_star()],
  )

  let indices = @src.nearest(gr1, gr2)

  assert_eq(indices.length(), 2)
  assert_eq(indices[0], 0)
  assert_eq(indices[1], 0)
}

///|
test "genomic_ranges_strand_to_char" {
  assert_eq(@src.strand_to_char(@src.strand_plus()), "+")
  assert_eq(@src.strand_to_char(@src.strand_minus()), "-")
  assert_eq(@src.strand_to_char(@src.strand_star()), "*")
}

///|
test "genomic_ranges_empty" {
  let gr = @src.granges([], [], [])

  assert_eq(@src.granges_length(gr), 0)

  let shifted = @src.granges_shift(gr, 5)
  assert_eq(@src.granges_length(shifted), 0)

  let reduced = @src.granges_reduce(gr, 0)
  assert_eq(@src.granges_length(reduced), 0)

  let disjoined = @src.granges_disjoin(gr)
  assert_eq(@src.granges_length(disjoined), 0)
}
