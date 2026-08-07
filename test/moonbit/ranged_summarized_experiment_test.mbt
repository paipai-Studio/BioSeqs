///|
fn rse_test_object() -> @src.RangedSummarizedExperiment {
  @src.RangedSummarizedExperiment::new(
    assays=Map([
      (
        "counts",
        [
          [30.0, 31.0, 32.0],
          [10.0, 11.0, 12.0],
          [20.0, 21.0, 22.0],
          [40.0, 41.0, 42.0],
        ],
      ),
      (
        "normalized",
        [[3.0, 3.1, 3.2], [1.0, 1.1, 1.2], [2.0, 2.1, 2.2], [4.0, 4.1, 4.2]],
      ),
    ]),
    row_ranges=@src.granges(
      ["chr2", "chr1", "chr1", "chr3"],
      [(300, 349), (100, 149), (200, 249), (50, 99)],
      [
        @src.strand_plus(),
        @src.strand_plus(),
        @src.strand_minus(),
        @src.strand_star(),
      ],
    ),
    col_data=[
      Map([("sample", "S1")]),
      Map([("sample", "S2")]),
      Map([("sample", "S3")]),
    ],
    row_data=[
      Map([("type", "geneC")]),
      Map([("type", "geneA")]),
      Map([("type", "geneB")]),
      Map([("type", "geneD")]),
    ],
    row_names=["geneC", "geneA", "geneB", "geneD"],
    metadata=Map([("study", "airway-like")]),
  ) catch {
    _ => abort("failed to construct RangedSummarizedExperiment fixture")
  }
}

///|
test "ranged_summarized_experiment: construct and access" {
  let rse = rse_test_object()
  assert_true(rse.is_valid())
  assert_eq(rse.nrow(), 4)
  assert_eq(rse.ncol(), 3)
  assert_eq(rse.assay_names().length(), 2)
  assert_eq(rse.row_ranges().seqnames[0], "chr2")
  assert_eq(rse.row_data()[1]["type"], "geneA")
  assert_eq(rse.row_names()[2], "geneB")
  assert_eq(rse.col_data()[2]["sample"], "S3")
  assert_eq(rse.metadata()["study"], "airway-like")
  assert_eq(
    rse.summary(),
    "RangedSummarizedExperiment(4 ranges x 3 samples, assays=[counts, normalized])",
  )
}

///|
test "ranged_summarized_experiment: rejects assay row mismatch" {
  let raised = try {
    ignore(
      @src.RangedSummarizedExperiment::new(
        assays=Map([("counts", [[1.0, 2.0], [3.0, 4.0]])]),
        row_ranges=@src.granges_single("chr1", 1, 10, @src.strand_plus()),
        col_data=[Map([]), Map([])],
      ),
    )
    false
  } catch {
    RangedSummarizedExperimentError(_) => true
  }
  assert_true(raised)
}

///|
test "ranged_summarized_experiment: rejects parallel annotation mismatch" {
  let raised = try {
    ignore(
      @src.RangedSummarizedExperiment::new(
        assays=Map([("counts", [[1.0], [2.0]])]),
        row_ranges=@src.granges(["chr1", "chr1"], [(1, 10), (20, 30)], [
          @src.strand_plus(),
          @src.strand_plus(),
        ]),
        col_data=[Map([])],
        row_data=[Map([("id", "only-one")])],
      ),
    )
    false
  } catch {
    RangedSummarizedExperimentError(_) => true
  }
  assert_true(raised)
}

///|
test "ranged_summarized_experiment: attach ranges to existing experiment" {
  let experiment = @src.summarized_experiment(
    Map([("counts", [[1.0, 2.0], [3.0, 4.0]])]),
    [("legacy", 0, 0), ("legacy", 0, 0)],
    [Map([("sample", "A")]), Map([("sample", "B")])],
    Map([("source", "existing")]),
  )
  let rse = @src.RangedSummarizedExperiment::from_experiment(
    experiment~,
    row_ranges=@src.granges(["chr1", "chr2"], [(10, 20), (30, 40)], [
      @src.strand_plus(),
      @src.strand_minus(),
    ]),
    row_names=["a", "b"],
  ) catch {
    _ => abort("failed to attach ranges")
  }
  assert_true(rse.is_valid())
  assert_eq(rse.row_ranges().seqnames, ["chr1", "chr2"])
  assert_eq(rse.metadata()["source"], "existing")
}

///|
test "ranged_summarized_experiment: subset rows stays coordinated" {
  let subset = rse_test_object().subset_rows([2, 0, 2])
  assert_true(subset.is_valid())
  assert_eq(subset.row_names(), ["geneB", "geneC", "geneB"])
  assert_eq(subset.row_ranges().seqnames, ["chr1", "chr2", "chr1"])
  assert_eq(subset.row_data()[0]["type"], "geneB")
  match subset.assay("counts") {
    Some(assay) =>
      assert_eq(assay, [
        [20.0, 21.0, 22.0],
        [30.0, 31.0, 32.0],
        [20.0, 21.0, 22.0],
      ])
    None => assert_true(false)
  }
}

///|
test "ranged_summarized_experiment: subset columns stays coordinated" {
  let subset = rse_test_object().subset_cols([2, 0])
  assert_true(subset.is_valid())
  assert_eq(subset.ncol(), 2)
  assert_eq(subset.col_data()[0]["sample"], "S3")
  assert_eq(subset.col_data()[1]["sample"], "S1")
  assert_eq(subset.row_names(), ["geneC", "geneA", "geneB", "geneD"])
  match subset.assay("counts") {
    Some(assay) => {
      assert_eq(assay[0], [32.0, 30.0])
      assert_eq(assay[2], [22.0, 20.0])
    }
    None => assert_true(false)
  }
}

///|
test "ranged_summarized_experiment: invalid subset index raises" {
  let raised = try {
    ignore(rse_test_object().subset_rows([4]))
    false
  } catch {
    RangedSummarizedExperimentError(_) => true
  }
  assert_true(raised)
}

///|
test "ranged_summarized_experiment: strand-aware overlaps" {
  let subject = @src.granges(
    ["chr1", "chr1", "chr3", "chr2"],
    [(120, 130), (220, 230), (80, 120), (340, 360)],
    [
      @src.strand_plus(),
      @src.strand_plus(),
      @src.strand_minus(),
      @src.strand_plus(),
    ],
  )
  let rse = rse_test_object()
  assert_eq(rse.find_overlaps(subject), [(0, 3), (1, 0), (3, 2)])
  assert_eq(rse.count_overlaps(subject), [1, 1, 0, 1])
  assert_eq(rse.overlaps_any(subject), [true, true, false, true])
  assert_eq(rse.find_overlaps(subject, ignore_strand=true), [
    (0, 3),
    (1, 0),
    (2, 1),
    (3, 2),
  ])
}

///|
test "ranged_summarized_experiment: subset by overlaps keeps each row once" {
  let subject = @src.granges(
    ["chr1", "chr1", "chr1"],
    [(105, 110), (120, 125), (225, 230)],
    [@src.strand_plus(), @src.strand_plus(), @src.strand_minus()],
  )
  let subset = rse_test_object().subset_by_overlaps(subject)
  assert_eq(subset.nrow(), 2)
  assert_eq(subset.row_names(), ["geneA", "geneB"])
  match subset.assay("counts") {
    Some(assay) => assert_eq(assay, [[10.0, 11.0, 12.0], [20.0, 21.0, 22.0]])
    None => assert_true(false)
  }
}

///|
test "ranged_summarized_experiment: nearest and distance" {
  let subject = @src.granges(
    ["chr1", "chr1", "chr3"],
    [(160, 170), (260, 270), (1, 10)],
    [@src.strand_plus(), @src.strand_minus(), @src.strand_plus()],
  )
  let rse = rse_test_object()
  assert_eq(rse.nearest(subject), [-1, 0, 1, 2])
  assert_eq(rse.distance_to_nearest(subject), [-1, 10, 10, 39])
}

///|
test "ranged_summarized_experiment: coverage delegates to row ranges" {
  let sequence_lengths = Map([("chr1", 260), ("chr2", 350), ("chr3", 100)])
  let coverage = rse_test_object().coverage(sequence_lengths)
  assert_eq(coverage["chr1"][99], 1)
  assert_eq(coverage["chr1"][149], 0)
  assert_eq(coverage["chr1"][199], 1)
  assert_eq(coverage["chr2"][299], 1)
  assert_eq(coverage["chr3"][98], 1)
}

///|
test "ranged_summarized_experiment: intra-range transformations preserve data" {
  let rse = rse_test_object()
  let shifted = rse.shift(10)
  assert_eq(shifted.row_ranges().starts, [310, 110, 210, 60])

  let narrowed = rse.narrow(2, 10)
  assert_eq(narrowed.row_ranges().widths, [9, 9, 9, 9])

  let resized = rse.resize(5, "start")
  assert_eq(resized.row_ranges().ends, [304, 104, 204, 54])

  let flanked = rse.flank(10, true, false)
  assert_eq(flanked.row_ranges().starts[1], 90)
  assert_eq(flanked.row_ranges().ends[1], 99)

  match shifted.assay("counts") {
    Some(assay) => assert_eq(assay[1], [10.0, 11.0, 12.0])
    None => assert_true(false)
  }
  assert_eq(shifted.row_names(), rse.row_names())
}

///|
test "ranged_summarized_experiment: promoters are strand-aware" {
  let promoters = rse_test_object().promoters(upstream=50, downstream=20)
  assert_eq(promoters.row_ranges().starts, [250, 50, 230, 0])
  assert_eq(promoters.row_ranges().ends, [319, 119, 299, 69])
  assert_eq(promoters.row_ranges().widths, [70, 70, 70, 70])
  assert_true(promoters.is_valid())
}

///|
test "ranged_summarized_experiment: sort reorders all row components" {
  let sorted = rse_test_object().sort()
  assert_eq(sorted.row_ranges().seqnames, ["chr1", "chr1", "chr2", "chr3"])
  assert_eq(sorted.row_names(), ["geneA", "geneB", "geneC", "geneD"])
  assert_eq(sorted.row_data()[0]["type"], "geneA")
  match sorted.assay("counts") {
    Some(assay) => {
      assert_eq(assay[0], [10.0, 11.0, 12.0])
      assert_eq(assay[2], [30.0, 31.0, 32.0])
    }
    None => assert_true(false)
  }

  let descending = rse_test_object().sort(decreasing=true)
  assert_eq(descending.row_names(), ["geneD", "geneC", "geneB", "geneA"])
}

///|
test "ranged_summarized_experiment: replacing ranges validates length" {
  let raised = try {
    ignore(
      rse_test_object().with_row_ranges(
        @src.granges_single("chr1", 1, 10, @src.strand_plus()),
      ),
    )
    false
  } catch {
    RangedSummarizedExperimentError(_) => true
  }
  assert_true(raised)
}

///|
test "ranged_summarized_experiment: empty container" {
  let rse = @src.RangedSummarizedExperiment::new(
    assays=Map([]),
    row_ranges=@src.granges([], [], []),
    col_data=[],
  ) catch {
    _ => abort("failed to construct empty RangedSummarizedExperiment")
  }
  assert_true(rse.is_valid())
  assert_eq(rse.nrow(), 0)
  assert_eq(rse.ncol(), 0)
  assert_eq(rse.subset_by_overlaps(@src.granges([], [], [])).nrow(), 0)
}
