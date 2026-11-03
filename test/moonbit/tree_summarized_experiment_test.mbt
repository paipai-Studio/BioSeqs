///|
fn tse_test_row_tree() -> @src.Tree {
  let a = @src.Clade::new(name=Some("A"))
  let b = @src.Clade::new(name=Some("B"))
  let c = @src.Clade::new(name=Some("C"))
  let d = @src.Clade::new(name=Some("D"))
  let group_a = @src.Clade::new(name=Some("GroupA"), clades=[a, b])
  let group_b = @src.Clade::new(name=Some("GroupB"), clades=[c, d])
  let root = @src.Clade::new(name=Some("All"), clades=[group_a, group_b])
  @src.Tree::new(root, rooted=true, name=Some("taxa"))
}

///|
fn tse_test_col_tree() -> @src.Tree {
  let s1 = @src.Clade::new(name=Some("S1"))
  let s2 = @src.Clade::new(name=Some("S2"))
  let s3 = @src.Clade::new(name=Some("S3"))
  let s4 = @src.Clade::new(name=Some("S4"))
  let control = @src.Clade::new(name=Some("Control"), clades=[s1, s2])
  let treated = @src.Clade::new(name=Some("Treated"), clades=[s3, s4])
  let root = @src.Clade::new(name=Some("Samples"), clades=[control, treated])
  @src.Tree::new(root, rooted=true, name=Some("samples"))
}

///|
fn tse_test_experiment() -> @src.SummarizedExperiment {
  @src.summarized_experiment(
    Map([
      (
        "counts",
        [
          [1.0, 2.0, 3.0, 4.0],
          [3.0, 4.0, 5.0, 6.0],
          [5.0, 6.0, 7.0, 8.0],
          [7.0, 8.0, 9.0, 10.0],
        ],
      ),
      (
        "normalized",
        [
          [0.1, 0.2, 0.3, 0.4],
          [0.3, 0.4, 0.5, 0.6],
          [0.5, 0.6, 0.7, 0.8],
          [0.7, 0.8, 0.9, 1.0],
        ],
      ),
    ]),
    [("A", 0, 0), ("B", 0, 0), ("C", 0, 0), ("D", 0, 0)],
    [
      Map([("sample", "S1")]),
      Map([("sample", "S2")]),
      Map([("sample", "S3")]),
      Map([("sample", "S4")]),
    ],
    Map([("study", "tree-demo")]),
  )
}

///|
fn tse_test_object() -> @src.TreeSummarizedExperiment {
  @src.TreeSummarizedExperiment::new(
    experiment=tse_test_experiment(),
    row_tree=Some(tse_test_row_tree()),
    row_node_labels=["A", "B", "C", "D"],
    col_tree=Some(tse_test_col_tree()),
    col_node_labels=["S1", "S2", "S3", "S4"],
    reference_sequences=["AAAA", "CCCC", "GGGG", "TTTT"],
  )
}

///|
test "tree_summarized_experiment: construct and access links" {
  let tse = tse_test_object()
  assert_true(tse.is_valid())
  assert_eq(tse.nrow(), 4)
  assert_eq(tse.ncol(), 4)
  assert_eq(tse.assay_names().length(), 2)
  assert_eq(tse.row_links().length(), 4)
  assert_eq(tse.col_links().length(), 4)
  assert_eq(tse.row_links()[0].node_label(), "A")
  assert_eq(tse.row_links()[0].node_alias(), "alias_1")
  assert_eq(tse.row_links()[0].node_number(), 1)
  assert_true(tse.row_links()[0].is_leaf())
  assert_eq(tse.row_links()[0].tree_name(), "taxa")
  assert_eq(tse.reference_sequences()[3], "TTTT")
}

///|
test "tree_summarized_experiment: invalid link is detected" {
  let tse = @src.TreeSummarizedExperiment::new(
    experiment=tse_test_experiment(),
    row_tree=Some(tse_test_row_tree()),
    row_node_labels=["A", "B", "missing", "D"],
  )
  assert_false(tse.is_valid())
  assert_eq(tse.row_links().length(), 3)
}

///|
test "tree_summarized_experiment: tree node queries" {
  let tree = tse_test_row_tree()
  assert_eq(@src.tse_find_descendants(tree, "GroupA", leaves_only=true), [
    "A", "B",
  ])
  assert_eq(@src.tse_find_ancestors(tree, "C"), ["All", "GroupB"])
  assert_true(@src.tse_is_leaf(tree, "D"))
  assert_false(@src.tse_is_leaf(tree, "GroupB"))
  assert_eq(@src.tse_find_descendants(tree, "missing").length(), 0)
}

///|
test "tree_summarized_experiment: subset rows keeps links and sequences" {
  let subset = tse_test_object().subset_rows([1, 3])
  assert_true(subset.is_valid())
  assert_eq(subset.nrow(), 2)
  assert_eq(subset.row_links()[0].node_label(), "B")
  assert_eq(subset.row_links()[1].node_label(), "D")
  assert_eq(subset.reference_sequences(), ["CCCC", "TTTT"])
  match subset.assay("counts") {
    Some(assay) =>
      assert_eq(assay, [[3.0, 4.0, 5.0, 6.0], [7.0, 8.0, 9.0, 10.0]])
    None => assert_true(false)
  }
}

///|
test "tree_summarized_experiment: subset by row tree node" {
  match tse_test_object().subset_by_row_nodes(["GroupA"]) {
    None => assert_true(false)
    Some(subset) => {
      assert_true(subset.is_valid())
      assert_eq(subset.nrow(), 2)
      assert_eq(subset.row_links()[0].node_label(), "A")
      assert_eq(subset.row_links()[1].node_label(), "B")
    }
  }
}

///|
test "tree_summarized_experiment: subset by column tree node" {
  match tse_test_object().subset_by_col_nodes(["Treated"]) {
    None => assert_true(false)
    Some(subset) => {
      assert_true(subset.is_valid())
      assert_eq(subset.ncol(), 2)
      assert_eq(subset.col_links()[0].node_label(), "S3")
      assert_eq(subset.col_links()[1].node_label(), "S4")
    }
  }
}

///|
test "tree_summarized_experiment: aggregate rows by sum" {
  match tse_test_object().aggregate_rows(["GroupA", "GroupB"]) {
    None => assert_true(false)
    Some(aggregated) => {
      assert_true(aggregated.is_valid())
      assert_eq(aggregated.nrow(), 2)
      assert_eq(aggregated.row_links()[0].node_label(), "GroupA")
      assert_false(aggregated.row_links()[0].is_leaf())
      match aggregated.assay("counts") {
        Some(assay) =>
          assert_eq(assay, [[4.0, 6.0, 8.0, 10.0], [12.0, 14.0, 16.0, 18.0]])
        None => assert_true(false)
      }
      match aggregated.assay("normalized") {
        Some(assay) => {
          assert_true((assay[0][0] - 0.4).abs() < 1.0e-12)
          assert_true((assay[0][1] - 0.6).abs() < 1.0e-12)
          assert_true((assay[0][2] - 0.8).abs() < 1.0e-12)
          assert_true((assay[0][3] - 1.0).abs() < 1.0e-12)
        }
        None => assert_true(false)
      }
    }
  }
}

///|
test "tree_summarized_experiment: aggregate rows by mean min and max" {
  let tse = tse_test_object()
  match
    tse.aggregate_rows(["GroupA"], aggregation=@src.tse_aggregation_mean()) {
    Some(result) =>
      match result.assay("counts") {
        Some(assay) => assert_eq(assay[0], [2.0, 3.0, 4.0, 5.0])
        None => assert_true(false)
      }
    None => assert_true(false)
  }
  match tse.aggregate_rows(["GroupA"], aggregation=@src.tse_aggregation_min()) {
    Some(result) =>
      match result.assay("counts") {
        Some(assay) => assert_eq(assay[0], [1.0, 2.0, 3.0, 4.0])
        None => assert_true(false)
      }
    None => assert_true(false)
  }
  match tse.aggregate_rows(["GroupA"], aggregation=@src.tse_aggregation_max()) {
    Some(result) =>
      match result.assay("counts") {
        Some(assay) => assert_eq(assay[0], [3.0, 4.0, 5.0, 6.0])
        None => assert_true(false)
      }
    None => assert_true(false)
  }
}

///|
test "tree_summarized_experiment: aggregate columns" {
  match tse_test_object().aggregate_cols(["Control", "Treated"]) {
    None => assert_true(false)
    Some(aggregated) => {
      assert_true(aggregated.is_valid())
      assert_eq(aggregated.ncol(), 2)
      assert_eq(aggregated.col_links()[0].node_label(), "Control")
      match aggregated.assay("counts") {
        Some(assay) => {
          assert_eq(assay[0], [3.0, 7.0])
          assert_eq(assay[3], [15.0, 19.0])
        }
        None => assert_true(false)
      }
    }
  }
}

///|
test "tree_summarized_experiment: duplicate links can aggregate" {
  let tse = @src.TreeSummarizedExperiment::new(
    experiment=tse_test_experiment(),
    row_tree=Some(tse_test_row_tree()),
    row_node_labels=["A", "A", "C", "D"],
  )
  assert_true(tse.is_valid())
  match tse.aggregate_rows(["A"]) {
    Some(aggregated) =>
      match aggregated.assay("counts") {
        Some(assay) => assert_eq(assay[0], [4.0, 6.0, 8.0, 10.0])
        None => assert_true(false)
      }
    None => assert_true(false)
  }
}

///|
test "tree_summarized_experiment: works without trees" {
  let tse = @src.TreeSummarizedExperiment::new(experiment=tse_test_experiment())
  assert_true(tse.is_valid())
  assert_eq(tse.row_links().length(), 0)
  assert_eq(tse.col_links().length(), 0)
  assert_true(tse.aggregate_rows(["All"]) is None)
  assert_true(tse.subset_by_col_nodes(["Samples"]) is None)
}
