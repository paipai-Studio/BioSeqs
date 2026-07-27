///|
/// GenomicRanges Gaps 模块测试。
/// 测试 @src.grg_* 系列函数：GRanges 高级操作。

///|
test "grg_granges_constructor" {
  let gr = @src.grg_granges(["chr1", "chr1"], [1, 10], [10, 20])
  assert_eq(gr.seqnames.length(), 2)
  assert_eq(gr.starts[0], 1)
  assert_eq(gr.ends[1], 20)
  assert_eq(gr.strands[0], "+")
  assert_eq(gr.strands[1], "+")
}

///|
test "grg_granges_with_custom_strands" {
  let gr = @src.grg_granges(
    ["chr1", "chr2"],
    [5, 15],
    [10, 25],
    strands=["+", "-"],
  )
  assert_eq(gr.seqnames.length(), 2)
  assert_eq(gr.strands[0], "+")
  assert_eq(gr.strands[1], "-")
}

///|
test "grg_length_basic" {
  let gr = @src.grg_granges(
    ["chr1", "chr1", "chr2"],
    [1, 10, 20],
    [5, 15, 25],
  )
  assert_eq(@src.grg_length(gr), 3)
}

///|
test "grg_empty_creates_empty_granges" {
  let gr = @src.grg_empty()
  assert_eq(@src.grg_length(gr), 0)
  assert_eq(gr.seqnames.length(), 0)
  assert_eq(gr.starts.length(), 0)
  assert_eq(gr.ends.length(), 0)
  assert_eq(gr.strands.length(), 0)
}

///|
test "grg_gaps_finds_gaps_between_ranges" {
  let gr = @src.grg_granges(
    ["chr1", "chr1"],
    [1, 30],
    [10, 50],
  )
  let sl : Map[String, Int] = Map([], capacity=1)
  sl["chr1"] = 100
  let gaps = @src.grg_gaps(gr, sl)
  assert_eq(gaps.seqnames.length(), 2)
  assert_eq(gaps.seqnames[0], "chr1")
  assert_eq(gaps.starts[0], 11)
  assert_eq(gaps.ends[0], 29)
  assert_eq(gaps.gap_lengths[0], 19)
  assert_eq(gaps.starts[1], 51)
  assert_eq(gaps.ends[1], 100)
  assert_eq(gaps.gap_lengths[1], 50)
}

///|
test "grg_gaps_no_gaps_when_adjacent" {
  let gr = @src.grg_granges(
    ["chr1", "chr1"],
    [1, 11],
    [10, 20],
  )
  let sl : Map[String, Int] = Map([], capacity=1)
  sl["chr1"] = 20
  let gaps = @src.grg_gaps(gr, sl)
  assert_eq(gaps.seqnames.length(), 0)
}

///|
test "grg_gaps_empty_input" {
  let gr = @src.grg_empty()
  let sl : Map[String, Int] = Map([], capacity=1)
  sl["chr1"] = 100
  let gaps = @src.grg_gaps(gr, sl)
  assert_eq(gaps.seqnames.length(), 0)
  assert_eq(gaps.starts.length(), 0)
  assert_eq(gaps.ends.length(), 0)
  assert_eq(gaps.gap_lengths.length(), 0)
}

///|
test "grg_gaps_multiple_chromosomes" {
  let gr = @src.grg_granges(
    ["chr1", "chr2"],
    [1, 1],
    [10, 10],
  )
  let sl : Map[String, Int] = Map([], capacity=2)
  sl["chr1"] = 20
  sl["chr2"] = 30
  let gaps = @src.grg_gaps(gr, sl)
  assert_eq(gaps.seqnames.length(), 2)
  assert_eq(gaps.seqnames[0], "chr1")
  assert_eq(gaps.starts[0], 11)
  assert_eq(gaps.ends[0], 20)
  assert_eq(gaps.seqnames[1], "chr2")
  assert_eq(gaps.starts[1], 11)
  assert_eq(gaps.ends[1], 30)
}

///|
test "grg_gaps_overlapping_ranges_merged" {
  let gr = @src.grg_granges(
    ["chr1", "chr1", "chr1"],
    [1, 5, 20],
    [10, 15, 30],
  )
  let sl : Map[String, Int] = Map([], capacity=1)
  sl["chr1"] = 50
  let gaps = @src.grg_gaps(gr, sl)
  assert_eq(gaps.seqnames.length(), 2)
  assert_eq(gaps.starts[0], 16)
  assert_eq(gaps.ends[0], 19)
  assert_eq(gaps.starts[1], 31)
  assert_eq(gaps.ends[1], 50)
}

///|
test "grg_reduce_merges_overlapping_ranges" {
  let gr = @src.grg_granges(
    ["chr1", "chr1"],
    [1, 5],
    [10, 15],
  )
  let result = @src.grg_reduce(gr)
  assert_eq(result.seqnames.length(), 1)
  assert_eq(result.starts[0], 1)
  assert_eq(result.ends[0], 15)
  assert_eq(result.merged_count[0], 2)
}

///|
test "grg_reduce_merges_adjacent_ranges" {
  let gr = @src.grg_granges(
    ["chr1", "chr1"],
    [1, 11],
    [10, 20],
  )
  let result = @src.grg_reduce(gr)
  assert_eq(result.seqnames.length(), 1)
  assert_eq(result.starts[0], 1)
  assert_eq(result.ends[0], 20)
}

///|
test "grg_reduce_no_overlapping_ranges" {
  let gr = @src.grg_granges(
    ["chr1", "chr1"],
    [1, 30],
    [10, 50],
  )
  let result = @src.grg_reduce(gr)
  assert_eq(result.seqnames.length(), 2)
  assert_eq(result.starts[0], 1)
  assert_eq(result.ends[0], 10)
  assert_eq(result.starts[1], 30)
  assert_eq(result.ends[1], 50)
}

///|
test "grg_reduce_with_min_gap_parameter" {
  let gr = @src.grg_granges(
    ["chr1", "chr1"],
    [1, 15],
    [10, 25],
  )
  let result = @src.grg_reduce(gr, min_gap=5)
  assert_eq(result.seqnames.length(), 1)
  assert_eq(result.starts[0], 1)
  assert_eq(result.ends[0], 25)
}

///|
test "grg_reduce_empty_input" {
  let gr = @src.grg_empty()
  let result = @src.grg_reduce(gr)
  assert_eq(result.seqnames.length(), 0)
  assert_eq(result.merged_count.length(), 0)
}

///|
test "grg_disjoin_splits_overlapping_ranges" {
  let gr = @src.grg_granges(
    ["chr1", "chr1"],
    [1, 5],
    [10, 15],
  )
  let result = @src.grg_disjoin(gr)
  assert_eq(result.seqnames.length(), 3)
  assert_eq(result.starts[0], 1)
  assert_eq(result.ends[0], 4)
  assert_eq(result.starts[1], 5)
  assert_eq(result.ends[1], 10)
  assert_eq(result.starts[2], 11)
  assert_eq(result.ends[2], 15)
}

///|
test "grg_disjoin_no_overlapping_ranges" {
  let gr = @src.grg_granges(
    ["chr1", "chr1"],
    [1, 20],
    [10, 30],
  )
  let result = @src.grg_disjoin(gr)
  assert_eq(result.seqnames.length(), 3)
  assert_eq(result.starts[0], 1)
  assert_eq(result.ends[0], 10)
  assert_eq(result.starts[1], 11)
  assert_eq(result.ends[1], 19)
  assert_eq(result.starts[2], 20)
  assert_eq(result.ends[2], 30)
}

///|
test "grg_disjoin_empty_input" {
  let gr = @src.grg_empty()
  let result = @src.grg_disjoin(gr)
  assert_eq(result.seqnames.length(), 0)
  assert_eq(result.segment_ids.length(), 0)
}

///|
test "grg_setdiff_removes_overlap" {
  let gr1 = @src.grg_granges(["chr1"], [1], [30])
  let gr2 = @src.grg_granges(["chr1"], [10], [20])
  let result = @src.grg_setdiff(gr1, gr2)
  assert_eq(result.seqnames.length(), 2)
  assert_eq(result.starts[0], 1)
  assert_eq(result.ends[0], 9)
  assert_eq(result.starts[1], 21)
  assert_eq(result.ends[1], 30)
}

///|
test "grg_setdiff_no_overlap" {
  let gr1 = @src.grg_granges(["chr1"], [1], [10])
  let gr2 = @src.grg_granges(["chr1"], [20], [30])
  let result = @src.grg_setdiff(gr1, gr2)
  assert_eq(result.seqnames.length(), 1)
  assert_eq(result.starts[0], 1)
  assert_eq(result.ends[0], 10)
}

///|
test "grg_setdiff_full_containment" {
  let gr1 = @src.grg_granges(["chr1"], [1], [10])
  let gr2 = @src.grg_granges(["chr1"], [1], [20])
  let result = @src.grg_setdiff(gr1, gr2)
  assert_eq(result.seqnames.length(), 0)
}

///|
test "grg_intersect_finds_overlap" {
  let gr1 = @src.grg_granges(["chr1"], [1], [30])
  let gr2 = @src.grg_granges(["chr1"], [10], [20])
  let result = @src.grg_intersect(gr1, gr2)
  assert_eq(result.seqnames.length(), 1)
  assert_eq(result.starts[0], 10)
  assert_eq(result.ends[0], 20)
}

///|
test "grg_intersect_no_intersection" {
  let gr1 = @src.grg_granges(["chr1"], [1], [10])
  let gr2 = @src.grg_granges(["chr1"], [20], [30])
  let result = @src.grg_intersect(gr1, gr2)
  assert_eq(result.seqnames.length(), 0)
}

///|
test "grg_intersect_different_chromosomes" {
  let gr1 = @src.grg_granges(["chr1"], [1], [10])
  let gr2 = @src.grg_granges(["chr2"], [1], [10])
  let result = @src.grg_intersect(gr1, gr2)
  assert_eq(result.seqnames.length(), 0)
}

///|
test "grg_union_merges_ranges" {
  let gr1 = @src.grg_granges(["chr1"], [1], [10])
  let gr2 = @src.grg_granges(["chr1"], [5], [20])
  let result = @src.grg_union(gr1, gr2)
  assert_eq(result.seqnames.length(), 1)
  assert_eq(result.starts[0], 1)
  assert_eq(result.ends[0], 20)
}

///|
test "grg_union_disjoint_ranges" {
  let gr1 = @src.grg_granges(["chr1"], [1], [10])
  let gr2 = @src.grg_granges(["chr1"], [20], [30])
  let result = @src.grg_union(gr1, gr2)
  assert_eq(result.seqnames.length(), 2)
}

///|
test "grg_compute_coverage_basic" {
  let gr = @src.grg_granges(
    ["chr1", "chr1"],
    [1, 3],
    [5, 8],
  )
  let sl : Map[String, Int] = Map([], capacity=1)
  sl["chr1"] = 10
  let cov = @src.grg_compute_coverage(gr, sl)
  let arr = cov["chr1"]
  assert_eq(arr.length(), 10)
  assert_eq(arr[0], 1)
  assert_eq(arr[1], 1)
  assert_eq(arr[2], 2)
  assert_eq(arr[3], 2)
  assert_eq(arr[4], 2)
  assert_eq(arr[5], 1)
  assert_eq(arr[6], 1)
  assert_eq(arr[7], 1)
  assert_eq(arr[8], 0)
  assert_eq(arr[9], 0)
}

///|
test "grg_compute_coverage_no_ranges" {
  let gr = @src.grg_empty()
  let sl : Map[String, Int] = Map([], capacity=1)
  sl["chr1"] = 5
  let cov = @src.grg_compute_coverage(gr, sl)
  let arr = cov["chr1"]
  assert_eq(arr.length(), 5)
  assert_eq(arr[0], 0)
  assert_eq(arr[4], 0)
}

///|
test "grg_trim_clips_to_chromosome_bounds" {
  let gr = @src.grg_granges(
    ["chr1", "chr1"],
    [-5, 8],
    [10, 30],
  )
  let sl : Map[String, Int] = Map([], capacity=1)
  sl["chr1"] = 20
  let result = @src.grg_trim(gr, sl)
  assert_eq(result.seqnames.length(), 2)
  assert_eq(result.starts[0], 1)
  assert_eq(result.ends[0], 10)
  assert_eq(result.starts[1], 8)
  assert_eq(result.ends[1], 20)
}

///|
test "grg_trim_removes_completely_outside_ranges" {
  let gr = @src.grg_granges(
    ["chr1", "chr1"],
    [1, 25],
    [10, 30],
  )
  let sl : Map[String, Int] = Map([], capacity=1)
  sl["chr1"] = 20
  let result = @src.grg_trim(gr, sl)
  assert_eq(result.seqnames.length(), 1)
  assert_eq(result.starts[0], 1)
  assert_eq(result.ends[0], 10)
}

///|
test "grg_promoters_plus_strand" {
  let gr = @src.grg_granges(
    ["chr1"],
    [100],
    [200],
    strands=["+"],
  )
  let result = @src.grg_promoters(gr, upstream=50, downstream=20)
  assert_eq(result.seqnames.length(), 1)
  assert_eq(result.starts[0], 50)
  assert_eq(result.ends[0], 119)
}

///|
test "grg_promoters_minus_strand" {
  let gr = @src.grg_granges(
    ["chr1"],
    [100],
    [200],
    strands=["-"],
  )
  let result = @src.grg_promoters(gr, upstream=50, downstream=20)
  assert_eq(result.seqnames.length(), 1)
  assert_eq(result.starts[0], 181)
  assert_eq(result.ends[0], 250)
}

///|
test "grg_promoters_default_parameters" {
  let gr = @src.grg_granges(
    ["chr1"],
    [100],
    [200],
    strands=["+"],
  )
  let result = @src.grg_promoters(gr)
  assert_eq(result.seqnames.length(), 1)
  assert_eq(result.starts[0], 100 - 2000)
  assert_eq(result.ends[0], 100 + 200 - 1)
}

///|
test "grg_edge_case_single_nucleotide_range" {
  let gr = @src.grg_granges(["chr1"], [50], [50])
  assert_eq(@src.grg_length(gr), 1)
  let sl : Map[String, Int] = Map([], capacity=1)
  sl["chr1"] = 100
  let gaps = @src.grg_gaps(gr, sl)
  assert_eq(gaps.seqnames.length(), 2)
  assert_eq(gaps.starts[0], 1)
  assert_eq(gaps.ends[0], 49)
  assert_eq(gaps.starts[1], 51)
  assert_eq(gaps.ends[1], 100)
}

///|
test "grg_edge_case_multiple_chromosomes_reduce" {
  let gr = @src.grg_granges(
    ["chr1", "chr1", "chr2", "chr2"],
    [1, 5, 1, 10],
    [10, 15, 5, 20],
  )
  let result = @src.grg_reduce(gr)
  assert_eq(result.seqnames.length(), 3)
  assert_eq(result.seqnames[0], "chr1")
  assert_eq(result.starts[0], 1)
  assert_eq(result.ends[0], 15)
  assert_eq(result.seqnames[1], "chr2")
  assert_eq(result.starts[1], 1)
  assert_eq(result.ends[1], 5)
  assert_eq(result.seqnames[2], "chr2")
  assert_eq(result.starts[2], 10)
  assert_eq(result.ends[2], 20)
}

///|
test "grg_edge_case_setdiff_different_chromosomes" {
  let gr1 = @src.grg_granges(["chr1", "chr2"], [1, 1], [10, 10])
  let gr2 = @src.grg_granges(["chr1"], [5], [15])
  let result = @src.grg_setdiff(gr1, gr2)
  assert_eq(result.seqnames.length(), 2)
  assert_eq(result.seqnames[0], "chr1")
  assert_eq(result.starts[0], 1)
  assert_eq(result.ends[0], 4)
  assert_eq(result.seqnames[1], "chr2")
  assert_eq(result.starts[1], 1)
  assert_eq(result.ends[1], 10)
}

///|
test "grg_edge_case_intersect_empty_grg1" {
  let gr1 = @src.grg_empty()
  let gr2 = @src.grg_granges(["chr1"], [1], [10])
  let result = @src.grg_intersect(gr1, gr2)
  assert_eq(result.seqnames.length(), 0)
}

///|
test "grg_edge_case_coverage_multiple_chromosomes" {
  let gr = @src.grg_granges(
    ["chr1", "chr2"],
    [1, 1],
    [5, 3],
  )
  let sl : Map[String, Int] = Map([], capacity=2)
  sl["chr1"] = 5
  sl["chr2"] = 5
  let cov = @src.grg_compute_coverage(gr, sl)
  let arr1 = cov["chr1"]
  let arr2 = cov["chr2"]
  assert_eq(arr1.length(), 5)
  assert_eq(arr2.length(), 5)
  assert_eq(arr1[0], 1)
  assert_eq(arr1[4], 1)
  assert_eq(arr2[0], 1)
  assert_eq(arr2[3], 0)
  assert_eq(arr2[4], 0)
}

///|
test "grg_edge_case_trim_no_clipping_needed" {
  let gr = @src.grg_granges(
    ["chr1"],
    [5],
    [10],
  )
  let sl : Map[String, Int] = Map([], capacity=1)
  sl["chr1"] = 100
  let result = @src.grg_trim(gr, sl)
  assert_eq(result.seqnames.length(), 1)
  assert_eq(result.starts[0], 5)
  assert_eq(result.ends[0], 10)
}