///|
test "genomic_alignments_galignments_basic" {
  let galn = @src.GAlignments::new(
    ["chr1", "chr1", "chr2"],
    [100, 150, 200],
    [200, 250, 300],
    [@src.strand_plus(), @src.strand_minus(), @src.strand_plus()],
    ["read1", "read2", "read3"],
    [30, 60, 45],
    [0, 16, 0],
  )
  assert_eq(galn.length(), 3)
}

///|
test "genomic_alignments_coverage" {
  let galn = @src.GAlignments::new(
    ["chr1", "chr1"],
    [100, 105],
    [110, 115],
    [@src.strand_plus(), @src.strand_plus()],
    ["read1", "read2"],
    [30, 30],
    [0, 0],
  )
  let cov = @src.coverage(galn)
  assert_eq(cov.length(), 1)
  assert_eq(cov[0].seqname, "chr1")
  assert_eq(cov[0].positions.length(), 16)
}

///|
test "genomic_alignments_summarize_overlaps" {
  let galn = @src.GAlignments::new(
    ["chr1", "chr1", "chr1"],
    [100, 120, 200],
    [150, 180, 250],
    [@src.strand_plus(), @src.strand_plus(), @src.strand_plus()],
    ["read1", "read2", "read3"],
    [30, 30, 30],
    [0, 0, 0],
  )
  let features = @src.granges(["chr1", "chr1"], [(100, 180), (190, 260)], [
    @src.strand_plus(),
    @src.strand_plus(),
  ])
  let summaries = @src.summarize_overlaps(galn, features, ["gene1", "gene2"])
  assert_eq(summaries.length(), 2)
  assert_eq(summaries[0].count, 2)
  assert_eq(summaries[1].count, 1)
}

///|
test "genomic_alignments_filter_by_mapq" {
  let galn = @src.GAlignments::new(
    ["chr1", "chr1", "chr1"],
    [100, 150, 200],
    [200, 250, 300],
    [@src.strand_plus(), @src.strand_plus(), @src.strand_plus()],
    ["read1", "read2", "read3"],
    [10, 30, 50],
    [0, 0, 0],
  )
  let filtered = galn.filter_by_mapq(30)
  assert_eq(filtered.length(), 2)
}

///|
test "genomic_alignments_filter_by_strand" {
  let galn = @src.GAlignments::new(
    ["chr1", "chr1", "chr1"],
    [100, 150, 200],
    [200, 250, 300],
    [@src.strand_plus(), @src.strand_minus(), @src.strand_plus()],
    ["read1", "read2", "read3"],
    [30, 30, 30],
    [0, 16, 0],
  )
  let plus = galn.filter_by_strand(@src.strand_plus())
  assert_eq(plus.length(), 2)
  let minus = galn.filter_by_strand(@src.strand_minus())
  assert_eq(minus.length(), 1)
}

///|
test "genomic_alignments_to_granges" {
  let galn = @src.GAlignments::new(
    ["chr1", "chr2"],
    [100, 200],
    [200, 300],
    [@src.strand_plus(), @src.strand_minus()],
    ["read1", "read2"],
    [30, 30],
    [0, 0],
  )
  let gr = galn.to_granges()
  assert_eq(@src.granges_length(gr), 2)
}

///|
test "genomic_alignments_pileup" {
  let galn = @src.GAlignments::new(
    ["chr1", "chr1"],
    [1, 1],
    [3, 3],
    [@src.strand_plus(), @src.strand_plus()],
    ["ACG", "ATG"],
    [30, 30],
    [0, 0],
  )
  let pileup = @src.pileup(galn, @hashmap.HashMap([]))
  assert_eq(pileup.length(), 3)
}
