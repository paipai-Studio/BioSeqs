///|
fn main {
  println("=== GenomicRanges Demo ===")
  println("")

  println("1. GRanges Construction")
  let seqnames = ["chr1", "chr1", "chr2"]
  let ranges = [(1, 10), (5, 15), (20, 30)]
  let strands = [@bio.strand_plus(), @bio.strand_minus(), @bio.strand_star()]
  let gr = @bio.granges(seqnames, ranges, strands)

  let mut i = 0
  while i < @bio.granges_length(gr) {
    let strand_char = @bio.strand_to_char(gr.strands[i])
    println(
      "  \{gr.seqnames[i]}:\{gr.starts[i]}-\{gr.ends[i]}(\{strand_char}) w=\{gr.widths[i]}",
    )
    i = i + 1
  }
  println("")

  println("2. Interval Operations")
  let shifted = @bio.granges_shift(gr, 5)
  println("After shift(5):")
  i = 0
  while i < @bio.granges_length(shifted) {
    let strand_char = @bio.strand_to_char(shifted.strands[i])
    println(
      "  \{shifted.seqnames[i]}:\{shifted.starts[i]}-\{shifted.ends[i]}(\{strand_char})",
    )
    i = i + 1
  }

  let narrowed = @bio.granges_narrow(gr, 3, 8)
  println("After narrow(3,8):")
  i = 0
  while i < @bio.granges_length(narrowed) {
    let strand_char = @bio.strand_to_char(narrowed.strands[i])
    println(
      "  \{narrowed.seqnames[i]}:\{narrowed.starts[i]}-\{narrowed.ends[i]}(\{strand_char})",
    )
    i = i + 1
  }

  let resized = @bio.granges_resize(gr, 15, "center")
  println("After resize(15, center):")
  i = 0
  while i < @bio.granges_length(resized) {
    let strand_char = @bio.strand_to_char(resized.strands[i])
    println(
      "  \{resized.seqnames[i]}:\{resized.starts[i]}-\{resized.ends[i]}(\{strand_char})",
    )
    i = i + 1
  }

  let flanked = @bio.granges_flank(gr, 3, true, false)
  println("After flank(3, start):")
  i = 0
  while i < @bio.granges_length(flanked) {
    let strand_char = @bio.strand_to_char(flanked.strands[i])
    println(
      "  \{flanked.seqnames[i]}:\{flanked.starts[i]}-\{flanked.ends[i]}(\{strand_char})",
    )
    i = i + 1
  }
  println("")

  println("3. Reduce & Disjoin")
  let reduced = @bio.granges_reduce(gr, 0)
  println("After reduce:")
  i = 0
  while i < @bio.granges_length(reduced) {
    let strand_char = @bio.strand_to_char(reduced.strands[i])
    println(
      "  \{reduced.seqnames[i]}:\{reduced.starts[i]}-\{reduced.ends[i]}(\{strand_char})",
    )
    i = i + 1
  }

  let disjoined = @bio.granges_disjoin(gr)
  println("After disjoin:")
  i = 0
  while i < @bio.granges_length(disjoined) {
    let strand_char = @bio.strand_to_char(disjoined.strands[i])
    println(
      "  \{disjoined.seqnames[i]}:\{disjoined.starts[i]}-\{disjoined.ends[i]}(\{strand_char})",
    )
    i = i + 1
  }
  println("")

  println("4. Set Operations")
  let gr1 = @bio.granges(["chr1", "chr1"], [(1, 10), (20, 30)], [
    @bio.strand_plus(),
    @bio.strand_minus(),
  ])
  let gr2 = @bio.granges(["chr1"], [(5, 25)], [@bio.strand_star()])

  let union_result = @bio.granges_union(gr1, gr2)
  println("A union B:")
  i = 0
  while i < @bio.granges_length(union_result) {
    let strand_char = @bio.strand_to_char(union_result.strands[i])
    println(
      "  \{union_result.seqnames[i]}:\{union_result.starts[i]}-\{union_result.ends[i]}(\{strand_char})",
    )
    i = i + 1
  }

  let intersect_result = @bio.granges_intersect(gr1, gr2)
  println("A intersect B:")
  i = 0
  while i < @bio.granges_length(intersect_result) {
    let strand_char = @bio.strand_to_char(intersect_result.strands[i])
    println(
      "  \{intersect_result.seqnames[i]}:\{intersect_result.starts[i]}-\{intersect_result.ends[i]}(\{strand_char})",
    )
    i = i + 1
  }

  let diff_result = @bio.granges_setdiff(gr1, gr2)
  println("A setdiff B:")
  i = 0
  while i < @bio.granges_length(diff_result) {
    let strand_char = @bio.strand_to_char(diff_result.strands[i])
    println(
      "  \{diff_result.seqnames[i]}:\{diff_result.starts[i]}-\{diff_result.ends[i]}(\{strand_char})",
    )
    i = i + 1
  }
  println("")

  println("5. Overlap Detection")
  let query = @bio.granges(["chr1", "chr1"], [(1, 5), (3, 8)], [
    @bio.strand_plus(),
    @bio.strand_minus(),
  ])
  let subject = @bio.granges(["chr1", "chr1"], [(2, 4), (6, 10)], [
    @bio.strand_star(),
    @bio.strand_plus(),
  ])

  let counts = @bio.count_overlaps(query, subject)
  println("count_overlaps: \{counts}")

  let overlaps = @bio.find_overlaps(query, subject)
  println("find_overlaps: \{overlaps}")
  println("")

  println("6. Distance & Nearest")
  let gr_a = @bio.granges(["chr1", "chr1"], [(1, 10), (20, 30)], [
    @bio.strand_plus(),
    @bio.strand_minus(),
  ])
  let gr_b = @bio.granges(["chr1", "chr1"], [(12, 15), (25, 28)], [
    @bio.strand_star(),
    @bio.strand_plus(),
  ])

  let distances = @bio.granges_distance(gr_a, gr_b)
  println("distances: \{distances}")

  let indices = @bio.nearest(gr_a, gr_b)
  println("nearest indices: \{indices}")
}
