///|
test "bedpe parse full 10 columns" {
  let content = "chr1\t100\t200\tchr1\t500\t600\tloop1\t500\t+\t-\nchr2\t50\t150\tchr3\t250\t350\tsv1\t900\t.\t.\n"
  let pairs = @src.parse_bedpe(content)
  assert_eq(pairs.length(), 2)
  assert_eq(pairs[0].chrom1, "chr1")
  assert_eq(pairs[0].start1, 100)
  assert_eq(pairs[0].end1, 200)
  assert_eq(pairs[0].chrom2, "chr1")
  assert_eq(pairs[0].start2, 500)
  assert_eq(pairs[0].end2, 600)
  assert_eq(pairs[0].name, "loop1")
  assert_eq(pairs[0].score, 500)
  assert_eq(pairs[0].strand1, "+")
  assert_eq(pairs[0].strand2, "-")
  assert_eq(pairs[0].fields.length(), 0)
  assert_eq(pairs[1].chrom1, "chr2")
  assert_eq(pairs[1].chrom2, "chr3")
  assert_eq(pairs[1].strand1, ".")
}

///|
test "bedpe parse minimal 6 columns and extra fields" {
  let content = "chr1\t100\t200\tchr1\t500\t600\nchr4\t0\t10\tchr5\t20\t30\te\t100\t+\t-\ta1\ta2\n"
  let pairs = @src.parse_bedpe(content)
  assert_eq(pairs.length(), 2)
  assert_eq(pairs[0].name, ".")
  assert_eq(pairs[0].score, 0)
  assert_eq(pairs[0].strand1, ".")
  assert_eq(pairs[1].name, "e")
  assert_eq(pairs[1].fields, ["a1", "a2"])
}

///|
test "bedpe parse skips comments and blank lines" {
  let content = "# comment\ntrack name=loops\nbrowser position chr1:1-1000\n\nchr1\t10\t20\tchr1\t30\t40\tL\t100\t+\t+\n"
  let pairs = @src.parse_bedpe(content)
  assert_eq(pairs.length(), 1)
  assert_eq(pairs[0].name, "L")
}

///|
test "bedpe parse space separated" {
  let content = "chr1 100 200 chr1 500 600 loop1 500 + -\n"
  let pairs = @src.parse_bedpe(content)
  assert_eq(pairs.length(), 1)
  assert_eq(pairs[0].name, "loop1")
  assert_eq(pairs[0].score, 500)
}

///|
test "bedpe parse crlf" {
  let content = "chr1\t100\t200\tchr1\t500\t600\r\n"
  let pairs = @src.parse_bedpe(content)
  assert_eq(pairs.length(), 1)
  assert_eq(pairs[0].end2, 600)
}

///|
test "bedpe roundtrip write" {
  let content = "chr1\t100\t200\tchr1\t500\t600\tloop1\t500\t+\t-\tX\tY\n"
  let pairs = @src.parse_bedpe(content)
  let out = @src.write_bedpe(pairs)
  let reparsed = @src.parse_bedpe(out)
  assert_eq(reparsed.length(), 1)
  assert_eq(reparsed[0].name, "loop1")
  assert_eq(reparsed[0].score, 500)
  assert_eq(reparsed[0].strand1, "+")
  assert_eq(reparsed[0].strand2, "-")
  assert_eq(reparsed[0].fields, ["X", "Y"])
}

///|
fn bedpe_raises(content : String) -> String {
  try {
    ignore(@src.parse_bedpe(content))
    "NO_ERROR"
  } catch {
    BedpeError(msg) => msg
  }
}

///|
test "bedpe parse too few columns raises" {
  let msg = bedpe_raises("chr1\t100\t200\tchr1\n")
  assert_true(msg.contains("at least 6 columns"))
}

///|
test "bedpe parse invalid coordinate raises" {
  let msg = bedpe_raises("chr1\tabc\t200\tchr1\t500\t600\n")
  assert_true(msg.contains("Invalid start1"))
}

///|
test "bedpe parse end before start raises" {
  let msg = bedpe_raises("chr1\t300\t200\tchr1\t500\t600\n")
  assert_true(msg.contains("end before start"))
}

///|
test "bedpe parse invalid strand raises" {
  let msg = bedpe_raises("chr1\t100\t200\tchr1\t500\t600\tn\t100\tx\ty\n")
  assert_true(msg.contains("strand1"))
}

///|
test "bedpe parse score out of range raises" {
  let msg = bedpe_raises("chr1\t100\t200\tchr1\t500\t600\tn\t2000\t+\t-\n")
  assert_true(msg.contains("out of range"))
}

///|
test "bedpe loop_length" {
  let same = @src.BedpePair::new("chr1", 100, 200, "chr1", 500, 600)
  assert_eq(same.loop_length(), 300)
  let overlap = @src.BedpePair::new("chr1", 100, 250, "chr1", 200, 600)
  assert_eq(overlap.loop_length(), 0)
  let inter = @src.BedpePair::new("chr1", 100, 200, "chr2", 500, 600)
  assert_eq(inter.loop_length(), -1)
}

///|
test "bedpe anchor length and same chrom" {
  let p = @src.BedpePair::new("chr1", 100, 200, "chr1", 500, 560)
  assert_eq(p.anchor_length1(), 100)
  assert_eq(p.anchor_length2(), 60)
  assert_true(p.is_same_chrom())
  assert_false(
    @src.BedpePair::new("chr1", 100, 200, "chr2", 500, 600).is_same_chrom(),
  )
}

///|
test "bedpe to_bed" {
  let p = @src.BedpePair::new("chr1", 100, 200, "chr2", 500, 600)
  let (a, b) = p.to_bed()
  assert_eq(a.chrom, "chr1")
  assert_eq(a.start, 100)
  assert_eq(a.end, 200)
  assert_eq(b.chrom, "chr2")
  assert_eq(b.start, 500)
  assert_eq(b.end, 600)
}

///|
test "bedpe overlaps_region" {
  let p = @src.BedpePair::new("chr1", 100, 200, "chr2", 500, 600)
  assert_true(p.overlaps_region("chr1", 150, 250))
  assert_true(p.overlaps_region("chr2", 550, 700))
  assert_false(p.overlaps_region("chr3", 0, 1000))
  assert_false(p.overlaps_region("chr1", 200, 300))
  assert_true(p.overlaps_region("chr1", 199, 201))
}

///|
test "bedpe filter same chrom" {
  let pairs = [
    @src.BedpePair::new("chr1", 100, 200, "chr1", 500, 600),
    @src.BedpePair::new("chr1", 100, 200, "chr2", 500, 600),
    @src.BedpePair::new("chr2", 0, 10, "chr2", 20, 30),
  ]
  let same = @src.bedpe_filter_same_chrom(pairs)
  assert_eq(same.length(), 2)
}

///|
test "bedpe count chrom pairs" {
  let pairs = [
    @src.BedpePair::new("chr1", 0, 10, "chr2", 0, 10),
    @src.BedpePair::new("chr1", 0, 10, "chr2", 20, 30),
    @src.BedpePair::new("chr1", 0, 10, "chr1", 20, 30),
  ]
  let counts = @src.bedpe_count_chrom_pairs(pairs)
  assert_eq(counts.get("chr1|chr2"), Some(2))
  assert_eq(counts.get("chr1|chr1"), Some(1))
}

///|
test "bedpe distance histogram" {
  let pairs = [
    @src.BedpePair::new("chr1", 0, 10, "chr1", 110, 120),
    @src.BedpePair::new("chr1", 0, 10, "chr1", 210, 220),
    @src.BedpePair::new("chr1", 0, 10, "chr1", 310, 320),
    @src.BedpePair::new("chr1", 0, 10, "chr2", 500, 600),
  ]
  let hist = @src.bedpe_distance_histogram(pairs, 0, 100, 3)
  assert_eq(hist, [0, 1, 2])
  let hist2 = @src.bedpe_distance_histogram(pairs, 0, 1000, 1)
  assert_eq(hist2, [3])
}

///|
test "bedpe pair overlaps" {
  let a = @src.BedpePair::new("chr1", 100, 200, "chr2", 500, 600)
  let b = @src.BedpePair::new("chr1", 150, 250, "chr3", 0, 10)
  let c = @src.BedpePair::new("chr9", 0, 10, "chr8", 0, 10)
  assert_true(@src.bedpe_pair_overlaps(a, b))
  assert_false(@src.bedpe_pair_overlaps(a, c))
}
