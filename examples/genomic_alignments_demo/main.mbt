///|
fn main {
  println("=== GenomicAlignments Demo ===")
  println("Inspired by Bioconductor GenomicAlignments package")
  println("")

  println("1. Create GAlignments object")
  let galn = @bio.GAlignments::new(
    ["chr1", "chr1", "chr1", "chr2"],
    [100, 105, 150, 200],
    [110, 115, 160, 210],
    [
      @bio.strand_plus(),
      @bio.strand_plus(),
      @bio.strand_minus(),
      @bio.strand_plus(),
    ],
    ["read1", "read2", "read3", "read4"],
    [30, 30, 60, 45],
    [0, 0, 16, 0],
  )
  println("Created GAlignments with \{galn.length()} alignments")
  println("")

  println("2. Filter by mapping quality (MAPQ >= 40)")
  let filtered = galn.filter_by_mapq(40)
  println("After filtering: \{filtered.length()} alignments")
  println("")

  println("3. Filter by strand (+)")
  let plus_strand = galn.filter_by_strand(@bio.strand_plus())
  println("Plus strand alignments: \{plus_strand.length()}")
  println("")

  println("4. Convert to GRanges")
  let gr = galn.to_granges()
  println("Converted to GRanges with \{@bio.granges_length(gr)} ranges")
  println("")

  println("5. Calculate coverage")
  let cov = @bio.coverage(galn)
  let mut i = 0
  while i < cov.length() {
    println(
      "  \{cov[i].seqname}: positions=\{cov[i].positions.length()}, depths=\{cov[i].depths.length()}",
    )
    i = i + 1
  }
  println("")

  println("6. Summarize overlaps with genes")
  let features = @bio.granges(["chr1", "chr1"], [(100, 180), (190, 260)], [
    @bio.strand_plus(),
    @bio.strand_plus(),
  ])
  let summaries = @bio.summarize_overlaps(galn, features, ["gene1", "gene2"])
  i = 0
  while i < summaries.length() {
    println(
      "  \{summaries[i].feature_id}: count=\{summaries[i].count}, coverage=\{summaries[i].coverage_fraction}",
    )
    i = i + 1
  }
  println("")

  println("7. Pileup analysis")
  let galn_pileup = @bio.GAlignments::new(
    ["chr1", "chr1"],
    [1, 1],
    [3, 3],
    [@bio.strand_plus(), @bio.strand_plus()],
    ["ACG", "ATG"],
    [30, 30],
    [0, 0],
  )
  let pileup = @bio.pileup(galn_pileup, @hashmap.HashMap([]))
  i = 0
  while i < pileup.length() {
    println(
      "  Position \{pileup[i].position}: ref=\{pileup[i].ref_base}, depth=\{pileup[i].depth}",
    )
    i = i + 1
  }
}
