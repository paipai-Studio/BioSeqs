///|
fn grl_test_object() -> @src.GRangesList {
  @src.granges_list(
    [
      @src.granges(["chr1", "chr1"], [(100, 109), (200, 209)], [
        @src.strand_plus(),
        @src.strand_plus(),
      ]),
      @src.granges(["chr1", "chr1"], [(150, 159), (300, 309)], [
        @src.strand_minus(),
        @src.strand_minus(),
      ]),
      @src.granges_single("chr2", 50, 59, @src.strand_star()),
      @src.granges([], [], []),
    ],
    names=["txA", "txB", "txC", "empty"],
    element_metadata=[
      Map([("gene", "geneA")]),
      Map([("gene", "geneB")]),
      Map([("gene", "geneC")]),
      Map([("gene", "none")]),
    ],
  ) catch {
    _ => abort("failed to construct GRangesList fixture")
  }
}

///|
fn grl_test_subject() -> @src.GRanges {
  @src.granges(
    ["chr1", "chr1", "chr1", "chr2", "chr3"],
    [(105, 205), (155, 156), (205, 206), (55, 55), (1, 10)],
    [
      @src.strand_plus(),
      @src.strand_minus(),
      @src.strand_plus(),
      @src.strand_star(),
      @src.strand_plus(),
    ],
  )
}

///|
fn grl_test_experiment() -> @src.RangedSummarizedExperiment {
  @src.RangedSummarizedExperiment::new_with_range_groups(
    assays=Map([
      ("counts", [[10.0, 11.0], [20.0, 21.0], [30.0, 31.0], [40.0, 41.0]]),
    ]),
    row_range_groups=grl_test_object(),
    col_data=[Map([("sample", "S1")]), Map([("sample", "S2")])],
    row_names=["txA", "txB", "txC", "empty"],
    metadata=Map([("organism", "human")]),
  ) catch {
    _ => abort("failed to construct grouped RangedSummarizedExperiment")
  }
}

///|
test "granges_list: construct and access" {
  let ranges = grl_test_object()
  assert_true(ranges.is_valid())
  assert_eq(ranges.length(), 4)
  assert_eq(ranges.names(), ["txA", "txB", "txC", "empty"])
  assert_eq(ranges.element_lengths(), [2, 2, 1, 0])
  assert_eq(ranges.element_metadata()[1]["gene"], "geneB")
  assert_false(ranges.is_empty())
  assert_eq(ranges.summary(), "GRangesList(4 elements, 5 ranges)")
}

///|
test "granges_list: get by index and name returns ranges" {
  let ranges = grl_test_object()
  match ranges.get(1) {
    Some(element) => {
      assert_eq(element.starts, [150, 300])
      assert_true(element.strands == [@src.strand_minus(), @src.strand_minus()])
    }
    None => assert_true(false)
  }
  match ranges.get_by_name("txC") {
    Some(element) => assert_eq(element.seqnames, ["chr2"])
    None => assert_true(false)
  }
  assert_true(ranges.get(-1) is None)
  assert_true(ranges.get_by_name("missing") is None)
}

///|
test "granges_list: empty semantics include all-empty elements" {
  let no_elements = @src.granges_list([]) catch {
    _ => abort("empty list should be valid")
  }
  let empty_elements = @src.granges_list(
    [@src.granges([], [], []), @src.granges([], [], [])],
    names=["a", "b"],
  ) catch {
    _ => abort("all-empty list should be valid")
  }
  assert_true(no_elements.is_empty())
  assert_true(empty_elements.is_empty())
}

///|
test "granges_list: rejects outer dimension mismatch" {
  let names_raised = try {
    ignore(
      @src.granges_list(
        [@src.granges_single("chr1", 1, 2, @src.strand_plus())],
        names=["a", "b"],
      ),
    )
    false
  } catch {
    GRangesListError(_) => true
  }
  let metadata_raised = try {
    ignore(
      @src.granges_list(
        [@src.granges_single("chr1", 1, 2, @src.strand_plus())],
        element_metadata=[Map([]), Map([])],
      ),
    )
    false
  } catch {
    GRangesListError(_) => true
  }
  assert_true(names_raised)
  assert_true(metadata_raised)
}

///|
test "granges_list: rejects malformed inner GRanges" {
  let raised = try {
    ignore(
      @src.granges_list([
        @src.granges(["chr1", "chr2"], [(1, 10), (20, 30)], [@src.strand_plus()]),
      ]),
    )
    false
  } catch {
    GRangesListError(_) => true
  }
  assert_true(raised)
}

///|
test "granges_list: split uses first-appearance group order" {
  let flat = @src.granges(
    ["chr1", "chr2", "chr1", "chr2", "chr1"],
    [(1, 5), (10, 15), (20, 25), (30, 35), (40, 45)],
    [
      @src.strand_plus(),
      @src.strand_minus(),
      @src.strand_plus(),
      @src.strand_minus(),
      @src.strand_plus(),
    ],
  )
  let split = @src.granges_split_as_list(flat, [
    "txB", "txA", "txB", "txA", "txC",
  ])
  assert_eq(split.names(), ["txB", "txA", "txC"])
  assert_eq(split.element_lengths(), [2, 2, 1])
  match split.get_by_name("txB") {
    Some(element) => assert_eq(element.starts, [1, 20])
    None => assert_true(false)
  }
}

///|
test "granges_list: split rejects group length mismatch" {
  let raised = try {
    ignore(
      @src.granges_split_as_list(
        @src.granges_single("chr1", 1, 10, @src.strand_plus()),
        [],
      ),
    )
    false
  } catch {
    GRangesListError(_) => true
  }
  assert_true(raised)
}

///|
test "granges_list: partition supports empty features" {
  let flat = @src.granges(
    ["chr1", "chr1", "chr2"],
    [(1, 5), (10, 15), (20, 25)],
    [@src.strand_plus(), @src.strand_plus(), @src.strand_minus()],
  )
  let partitioned = @src.granges_list_from_partition(flat, [2, 0, 1], names=[
    "tx1", "empty", "tx2",
  ])
  assert_eq(partitioned.element_lengths(), [2, 0, 1])
  match partitioned.get(1) {
    Some(element) => assert_eq(@src.granges_length(element), 0)
    None => assert_true(false)
  }
}

///|
test "granges_list: partition validates lengths" {
  let flat = @src.granges_single("chr1", 1, 10, @src.strand_plus())
  let sum_raised = try {
    ignore(@src.granges_list_from_partition(flat, [2]))
    false
  } catch {
    GRangesListError(_) => true
  }
  let negative_raised = try {
    ignore(@src.granges_list_from_partition(flat, [-1, 2]))
    false
  } catch {
    GRangesListError(_) => true
  }
  assert_true(sum_raised)
  assert_true(negative_raised)
}

///|
test "granges_list: unlist preserves member order" {
  let flat = grl_test_object().unlist()
  assert_eq(flat.seqnames, ["chr1", "chr1", "chr1", "chr1", "chr2"])
  assert_eq(flat.starts, [100, 200, 150, 300, 50])
  assert_eq(flat.ends, [109, 209, 159, 309, 59])
}

///|
test "granges_list: unlist and relist round trip" {
  let original = grl_test_object()
  let rebuilt = original.relist(original.unlist())
  assert_eq(rebuilt.names(), original.names())
  assert_eq(rebuilt.element_lengths(), original.element_lengths())
  assert_eq(rebuilt.element_metadata()[2]["gene"], "geneC")
  for index = 0; index < original.length(); index = index + 1 {
    match (original.get(index), rebuilt.get(index)) {
      (Some(left), Some(right)) => {
        assert_eq(left.seqnames, right.seqnames)
        assert_eq(left.starts, right.starts)
        assert_eq(left.ends, right.ends)
        assert_true(left.strands == right.strands)
      }
      _ => assert_true(false)
    }
  }
}

///|
test "granges_list: relist rejects incompatible flat length" {
  let raised = try {
    ignore(
      grl_test_object().relist(
        @src.granges_single("chr1", 1, 10, @src.strand_plus()),
      ),
    )
    false
  } catch {
    GRangesListError(_) => true
  }
  assert_true(raised)
}

///|
test "granges_list: subset coordinates names and metadata" {
  let subset = grl_test_object().subset([2, 0, 2])
  assert_eq(subset.names(), ["txC", "txA", "txC"])
  assert_eq(subset.element_lengths(), [1, 2, 1])
  assert_eq(subset.element_metadata()[0]["gene"], "geneC")
  assert_eq(subset.element_metadata()[1]["gene"], "geneA")
}

///|
test "granges_list: invalid subset index raises" {
  let raised = try {
    ignore(grl_test_object().subset([4]))
    false
  } catch {
    GRangesListError(_) => true
  }
  assert_true(raised)
}

///|
test "granges_list: feature bounds summarize compound features" {
  let bounds = grl_test_object().feature_bounds()
  assert_eq(bounds.seqnames, ["chr1", "chr1", "chr2", ""])
  assert_eq(bounds.starts, [100, 150, 50, 0])
  assert_eq(bounds.ends, [209, 309, 59, -1])
  assert_eq(bounds.widths, [110, 160, 10, 0])
  assert_true(bounds.strands[0] == @src.strand_plus())
  assert_true(bounds.strands[1] == @src.strand_minus())
}

///|
test "granges_list: mixed feature bounds use wildcard identifiers" {
  let mixed = @src.granges_list([
    @src.granges(["chr1", "chr2"], [(10, 20), (30, 40)], [
      @src.strand_plus(),
      @src.strand_minus(),
    ]),
  ]) catch {
    _ => abort("mixed feature should be valid")
  }
  let bounds = mixed.feature_bounds()
  assert_eq(bounds.seqnames, ["*"])
  assert_true(bounds.strands == [@src.strand_star()])
  assert_eq(bounds.starts, [10])
  assert_eq(bounds.ends, [40])
}

///|
test "granges_list: concat appends outer elements" {
  let left = grl_test_object().subset([0])
  let right = grl_test_object().subset([1, 2])
  let combined = left.concat(right)
  assert_eq(combined.names(), ["txA", "txB", "txC"])
  assert_eq(combined.element_lengths(), [2, 2, 1])
  assert_eq(combined.element_metadata()[2]["gene"], "geneC")
}

///|
test "granges_list: parallel concat combines corresponding features" {
  let shifted = grl_test_object().shift(1000)
  let combined = grl_test_object().parallel_concat(shifted)
  assert_eq(combined.element_lengths(), [4, 4, 2, 0])
  match combined.get(0) {
    Some(element) => assert_eq(element.starts, [100, 200, 1100, 1200])
    None => assert_true(false)
  }
}

///|
test "granges_list: parallel operations require equal element counts" {
  let raised = try {
    ignore(grl_test_object().parallel_concat(grl_test_object().subset([0, 1])))
    false
  } catch {
    GRangesListError(_) => true
  }
  assert_true(raised)
}

///|
test "granges_list: intra-range operations preserve groups" {
  let shifted = grl_test_object().shift(10)
  let narrowed = grl_test_object().narrow(2, 5)
  let resized = grl_test_object().resize(5, "start")
  assert_eq(shifted.element_lengths(), [2, 2, 1, 0])
  match shifted.get(0) {
    Some(element) => assert_eq(element.starts, [110, 210])
    None => assert_true(false)
  }
  match narrowed.get(0) {
    Some(element) => {
      assert_eq(element.starts, [101, 201])
      assert_eq(element.ends, [104, 204])
    }
    None => assert_true(false)
  }
  match resized.get(1) {
    Some(element) => assert_eq(element.ends, [154, 304])
    None => assert_true(false)
  }
}

///|
test "granges_list: flank and promoters preserve outer metadata" {
  let flanked = grl_test_object().flank(5, true, false)
  let promoted = grl_test_object().promoters(upstream=10, downstream=3)
  assert_eq(flanked.names(), grl_test_object().names())
  assert_eq(promoted.element_metadata()[0]["gene"], "geneA")
  match flanked.get(0) {
    Some(element) => assert_eq(element.starts, [95, 195])
    None => assert_true(false)
  }
  match promoted.get(1) {
    Some(element) => {
      assert_eq(element.starts, [157, 307])
      assert_eq(element.ends, [169, 319])
    }
    None => assert_true(false)
  }
}

///|
test "granges_list: promoters reject negative widths" {
  let raised = try {
    ignore(grl_test_object().promoters(upstream=-1))
    false
  } catch {
    GRangesListError(_) => true
  }
  assert_true(raised)
}

///|
test "granges_list: reduce disjoin and sort operate per feature" {
  let ranges = @src.granges_list(
    [
      @src.granges(["chr1", "chr1", "chr1"], [(20, 30), (1, 10), (8, 15)], [
        @src.strand_plus(),
        @src.strand_plus(),
        @src.strand_plus(),
      ]),
    ],
    names=["tx"],
  ) catch {
    _ => abort("failed to construct interval fixture")
  }
  match ranges.sort_ranges().get(0) {
    Some(element) => assert_eq(element.starts, [1, 8, 20])
    None => assert_true(false)
  }
  match ranges.reduce().get(0) {
    Some(element) => {
      assert_eq(element.starts, [1, 20])
      assert_eq(element.ends, [15, 30])
    }
    None => assert_true(false)
  }
  match ranges.disjoin().get(0) {
    Some(element) => assert_eq(element.starts, [1, 8, 11, 16, 20])
    None => assert_true(false)
  }
}

///|
test "granges_list: width summaries exclude introns" {
  assert_eq(grl_test_object().total_widths(), [20, 20, 10, 0])
}

///|
test "granges_list: parallel set operations preserve grouping" {
  let left = @src.granges_list([
    @src.granges_single("chr1", 1, 10, @src.strand_star()),
    @src.granges_single("chr2", 20, 30, @src.strand_star()),
  ]) catch {
    _ => abort("failed to build left set")
  }
  let right = @src.granges_list([
    @src.granges_single("chr1", 5, 15, @src.strand_star()),
    @src.granges_single("chr2", 25, 35, @src.strand_star()),
  ]) catch {
    _ => abort("failed to build right set")
  }
  match left.parallel_union(right).get(0) {
    Some(element) => {
      assert_eq(element.starts, [1])
      assert_eq(element.ends, [15])
    }
    None => assert_true(false)
  }
  match left.parallel_intersect(right).get(1) {
    Some(element) => {
      assert_eq(element.starts, [25])
      assert_eq(element.ends, [30])
    }
    None => assert_true(false)
  }
  match left.parallel_setdiff(right).get(0) {
    Some(element) => {
      assert_eq(element.starts, [1])
      assert_eq(element.ends, [4])
    }
    None => assert_true(false)
  }
}

///|
test "granges_list: feature overlaps are deduplicated and strand-aware" {
  let ranges = grl_test_object()
  let subject = grl_test_subject()
  assert_eq(ranges.find_overlaps(subject), [(0, 0), (0, 2), (1, 1), (2, 3)])
  assert_eq(ranges.count_overlaps(subject), [2, 1, 1, 0])
  assert_eq(ranges.overlaps_any(subject), [true, true, true, false])
  assert_eq(ranges.find_overlaps(subject, ignore_strand=true), [
    (0, 0),
    (0, 2),
    (1, 0),
    (1, 1),
    (2, 3),
  ])
}

///|
test "granges_list: subset by overlap keeps each feature once" {
  let subset = grl_test_object().subset_by_overlaps(grl_test_subject())
  assert_eq(subset.names(), ["txA", "txB", "txC"])
  assert_eq(subset.element_lengths(), [2, 2, 1])
}

///|
test "granges_list: grouped subjects use feature indices" {
  let subject = @src.granges_list(
    [
      @src.granges(["chr1", "chr1"], [(105, 106), (205, 206)], [
        @src.strand_plus(),
        @src.strand_plus(),
      ]),
      @src.granges_single("chr2", 55, 56, @src.strand_plus()),
    ],
    names=["feature1", "feature2"],
  ) catch {
    _ => abort("failed to construct grouped subject")
  }
  let query = grl_test_object()
  assert_eq(query.find_group_overlaps(subject), [(0, 0), (2, 1)])
  assert_eq(query.count_group_overlaps(subject), [1, 0, 1, 0])
  assert_eq(@src.granges_find_overlaps_list(query.unlist(), subject), [
    (0, 0),
    (1, 0),
    (4, 1),
  ])
}

///|
test "granges_list: nearest uses minimum member distance" {
  let subject = @src.granges(
    ["chr1", "chr1", "chr2"],
    [(120, 125), (170, 175), (70, 75)],
    [@src.strand_plus(), @src.strand_minus(), @src.strand_plus()],
  )
  let ranges = grl_test_object()
  assert_eq(ranges.nearest(subject), [0, 1, 2, -1])
  assert_eq(ranges.distance_to_nearest(subject), [10, 10, 10, -1])
}

///|
test "granges_list: coverage uses all member ranges" {
  let coverage = grl_test_object().coverage(Map([("chr1", 320), ("chr2", 80)]))
  assert_eq(coverage["chr1"][99], 1)
  assert_eq(coverage["chr1"][149], 1)
  assert_eq(coverage["chr1"][199], 1)
  assert_eq(coverage["chr1"][299], 1)
  assert_eq(coverage["chr1"][250], 0)
  assert_eq(coverage["chr2"][49], 1)
}

///|
test "grouped RangedSummarizedExperiment: construct and access" {
  let experiment = grl_test_experiment()
  assert_true(experiment.is_valid())
  assert_true(experiment.has_grouped_row_ranges())
  assert_eq(experiment.nrow(), 4)
  assert_eq(experiment.ncol(), 2)
  assert_eq(experiment.row_ranges().starts, [100, 150, 50, 0])
  assert_eq(experiment.row_ranges().ends, [209, 309, 59, -1])
  assert_eq(experiment.row_data()[0]["gene"], "geneA")
  assert_eq(experiment.metadata()["organism"], "human")
  assert_eq(
    experiment.summary(),
    "RangedSummarizedExperiment(4 range groups x 2 samples, assays=[counts])",
  )
  match experiment.row_range_groups() {
    Some(groups) => assert_eq(groups.element_lengths(), [2, 2, 1, 0])
    None => assert_true(false)
  }
}

///|
test "grouped RangedSummarizedExperiment: rejects assay mismatch" {
  let raised = try {
    ignore(
      @src.RangedSummarizedExperiment::new_with_range_groups(
        assays=Map([("counts", [[1.0], [2.0]])]),
        row_range_groups=grl_test_object(),
        col_data=[Map([])],
      ),
    )
    false
  } catch {
    RangedSummarizedExperimentError(_) => true
  }
  assert_true(raised)
}

///|
test "grouped RangedSummarizedExperiment: exact overlaps avoid introns" {
  let query = @src.granges(
    ["chr1", "chr1", "chr1"],
    [(120, 130), (205, 206), (250, 260)],
    [@src.strand_plus(), @src.strand_plus(), @src.strand_minus()],
  )
  let experiment = grl_test_experiment()
  assert_eq(experiment.find_overlaps(query), [(0, 1)])
  assert_eq(experiment.count_overlaps(query), [1, 0, 0, 0])
}

///|
test "grouped RangedSummarizedExperiment: row subset stays coordinated" {
  let subset = grl_test_experiment().subset_rows([2, 0, 2])
  assert_true(subset.is_valid())
  assert_eq(subset.row_names(), ["txC", "txA", "txC"])
  assert_eq(subset.row_data()[0]["gene"], "geneC")
  match subset.row_range_groups() {
    Some(groups) => {
      assert_eq(groups.names(), ["txC", "txA", "txC"])
      assert_eq(groups.element_lengths(), [1, 2, 1])
    }
    None => assert_true(false)
  }
  match subset.assay("counts") {
    Some(assay) => assert_eq(assay, [[30.0, 31.0], [10.0, 11.0], [30.0, 31.0]])
    None => assert_true(false)
  }
}

///|
test "grouped RangedSummarizedExperiment: column subset preserves groups" {
  let subset = grl_test_experiment().subset_cols([1])
  assert_true(subset.has_grouped_row_ranges())
  assert_eq(subset.ncol(), 1)
  assert_eq(subset.col_data()[0]["sample"], "S2")
  match subset.row_range_groups() {
    Some(groups) => assert_eq(groups.element_lengths(), [2, 2, 1, 0])
    None => assert_true(false)
  }
}

///|
test "grouped RangedSummarizedExperiment: transformations preserve groups" {
  let shifted = grl_test_experiment().shift(10)
  let promoted = grl_test_experiment().promoters(upstream=10, downstream=3)
  assert_true(shifted.has_grouped_row_ranges())
  assert_eq(shifted.row_ranges().starts, [110, 160, 60, 0])
  match shifted.row_range_groups() {
    Some(groups) =>
      match groups.get(0) {
        Some(element) => assert_eq(element.starts, [110, 210])
        None => assert_true(false)
      }
    None => assert_true(false)
  }
  match promoted.row_range_groups() {
    Some(groups) =>
      match groups.get(1) {
        Some(element) => assert_eq(element.starts, [157, 307])
        None => assert_true(false)
      }
    None => assert_true(false)
  }
}

///|
test "grouped RangedSummarizedExperiment: nearest and coverage use members" {
  let subject = @src.granges(
    ["chr1", "chr1", "chr2"],
    [(120, 125), (170, 175), (70, 75)],
    [@src.strand_plus(), @src.strand_minus(), @src.strand_plus()],
  )
  let experiment = grl_test_experiment()
  assert_eq(experiment.nearest(subject), [0, 1, 2, -1])
  assert_eq(experiment.distance_to_nearest(subject), [10, 10, 10, -1])
  let coverage = experiment.coverage(Map([("chr1", 320), ("chr2", 80)]))
  assert_eq(coverage["chr1"][120], 0)
  assert_eq(coverage["chr1"][199], 1)
}

///|
test "grouped RangedSummarizedExperiment: flat replacement clears grouping" {
  let flat = grl_test_experiment().with_row_ranges(
    @src.granges(
      ["chr1", "chr1", "chr2", "chr3"],
      [(1, 10), (20, 30), (40, 50), (60, 70)],
      [
        @src.strand_plus(),
        @src.strand_minus(),
        @src.strand_star(),
        @src.strand_plus(),
      ],
    ),
  )
  assert_false(flat.has_grouped_row_ranges())
  assert_true(flat.row_range_groups() is None)
  assert_eq(
    flat.summary(),
    "RangedSummarizedExperiment(4 ranges x 2 samples, assays=[counts])",
  )
}
