///|
/// Tests for plyranges - tidy genomic data manipulation

///|
test "plyranges_create_tidified_granges" {
  let ranges = @src.tgr_create(
    ["chr1", "chr1", "chr2"],
    [100, 200, 300],
    [200, 300, 400],
    ["+", "-", "+"],
  )
  assert_eq(ranges.seqnames.length(), 3)
  assert_eq(ranges.starts.length(), 3)
  assert_eq(ranges.ends.length(), 3)
}

///|
test "plyranges_length" {
  let ranges = @src.tgr_create(
    ["chr1", "chr1", "chr2", "chr2"],
    [100, 200, 300, 400],
    [200, 300, 400, 500],
    ["+", "-", "+", "-"],
  )
  assert_eq(@src.tgr_length(ranges), 4)
}

///|
test "plyranges_filter_seqname" {
  let ranges = @src.tgr_create(
    ["chr1", "chr1", "chr2", "chr2"],
    [100, 200, 300, 400],
    [200, 300, 400, 500],
    ["+", "-", "+", "-"],
  )
  let filtered = @src.tgr_filter_seqname(ranges, "chr1")
  assert_eq(filtered.seqnames.length(), 2)
  assert_eq(filtered.seqnames[0], "chr1")
}

///|
test "plyranges_filter_strand" {
  let ranges = @src.tgr_create(
    ["chr1", "chr1", "chr2", "chr2"],
    [100, 200, 300, 400],
    [200, 300, 400, 500],
    ["+", "-", "+", "-"],
  )
  let filtered = @src.tgr_filter_strand(ranges, "+")
  assert_eq(filtered.seqnames.length(), 2)
}

///|
test "plyranges_mutate" {
  let ranges = @src.tgr_create(
    ["chr1", "chr2"],
    [100, 200],
    [200, 300],
    ["+", "-"],
  )
  let mutated = @src.tgr_mutate(ranges, "gc_content", [0.45, 0.62])
  assert_eq(mutated.metadata["gc_content"].length(), 2)
  assert_eq(mutated.metadata["gc_content"][0], 0.45)
}

///|
test "plyranges_mutate_str" {
  let ranges = @src.tgr_create(
    ["chr1", "chr2"],
    [100, 200],
    [200, 300],
    ["+", "-"],
  )
  let mutated = @src.tgr_mutate_str(ranges, "gene_name", ["BRCA1", "TP53"])
  assert_eq(mutated.metadata_str["gene_name"].length(), 2)
  assert_eq(mutated.metadata_str["gene_name"][0], "BRCA1")
}

///|
test "plyranges_select" {
  let ranges = @src.tgr_create(
    ["chr1", "chr2"],
    [100, 200],
    [200, 300],
    ["+", "-"],
  )
  let mutated = @src.tgr_mutate(ranges, "score", [1.0, 2.0])
  let selected = @src.tgr_select(mutated, ["score"])
  assert_eq(selected.metadata["score"].length(), 2)
}

///|
test "plyranges_arrange" {
  let ranges = @src.tgr_create(
    ["chr1", "chr1", "chr2"],
    [300, 100, 200],
    [400, 200, 300],
    ["+", "-", "+"],
  )
  let arranged = @src.tgr_arrange(ranges)
  // Sorted by seqname then start: chr1/100, chr1/300, chr2/200
  assert_eq(arranged.starts[0], 100)
  assert_eq(arranged.starts[1], 300)
  assert_eq(arranged.starts[2], 200)
}

///|
test "plyranges_rename" {
  let ranges = @src.tgr_create(
    ["chr1", "chr2"],
    [100, 200],
    [200, 300],
    ["+", "-"],
  )
  let mutated = @src.tgr_mutate(ranges, "score", [1.0, 2.0])
  let renamed = @src.tgr_rename(mutated, "score", "alignment_score")
  assert_eq(renamed.metadata["alignment_score"].length(), 2)
}

///|
test "plyranges_width" {
  let ranges = @src.tgr_create(
    ["chr1", "chr2"],
    [100, 200],
    [200, 350],
    ["+", "-"],
  )
  assert_eq(ranges.widths[0], 101)
  assert_eq(ranges.widths[1], 151)
}

///|
test "plyranges_empty" {
  let ranges = @src.tgr_create([], [], [], [])
  assert_eq(ranges.seqnames.length(), 0)
  assert_eq(@src.tgr_length(ranges), 0)
}

///|
test "plyranges_filter_start" {
  let ranges = @src.tgr_create(
    ["chr1", "chr1", "chr2"],
    [100, 500, 200],
    [200, 600, 300],
    ["+", "-", "+"],
  )
  let filtered = @src.tgr_filter_start(ranges, 300)
  assert_eq(filtered.seqnames.length(), 1)
}

///|
test "plyranges_filter_end" {
  let ranges = @src.tgr_create(
    ["chr1", "chr1", "chr2"],
    [100, 500, 200],
    [200, 600, 300],
    ["+", "-", "+"],
  )
  let filtered = @src.tgr_filter_end(ranges, 250)
  assert_eq(filtered.seqnames.length(), 1)
}

///|
test "plyranges_multiple_metadata" {
  let ranges = @src.tgr_create(
    ["chr1", "chr2"],
    [100, 200],
    [200, 300],
    ["+", "-"],
  )
  let mutated = @src.tgr_mutate(ranges, "score", [1.5, 2.5])
  let mutated2 = @src.tgr_mutate(mutated, "p_value", [0.01, 0.05])
  assert_eq(mutated2.metadata.length(), 2)
}

///|
test "plyranges_summary" {
  let ranges = @src.tgr_create(
    ["chr1", "chr1", "chr2", "chr2"],
    [100, 200, 300, 400],
    [200, 300, 400, 500],
    ["+", "-", "+", "-"],
  )
  let summary = @src.tgr_summary(ranges)
  assert_eq(summary.n_ranges, 4)
}
