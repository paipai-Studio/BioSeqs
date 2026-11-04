// ===========================================================================
// treetime_lite_test.mbt — 分子钟标定与祖先序列重建测试(TreeTime-lite).
// ===========================================================================

///|
/// 构造树并计算 root-to-tip 拟合(辅助函数).
fn treetime_test_fit(
  newick : String,
  dates : Map[String, Double],
) -> @src.RootToTipResult {
  @src.phylo_clock_root_to_tip(@src.skbio_tree_read(newick), dates) catch {
    ClockError(message) => abort("root-to-tip fit failed: " + message)
  }
}

///|
test "treetime_root_to_tip_perfect_linear_clock" {
  let dates : Map[String, Double] = Map([
    ("A", 2020.0),
    ("B", 2021.0),
    ("C", 2022.0),
  ])
  let fit = treetime_test_fit("(A:0.1,B:0.2,C:0.3);", dates)
  assert_true(fit.rate > 0.0999 && fit.rate < 0.1001)
  assert_true(fit.root_date > 2018.99 && fit.root_date < 2019.01)
  assert_true(fit.r_squared > 0.999)
  assert_eq(fit.tips.length(), 3)
  for tip in fit.tips {
    assert_true(tip.residual.abs() < 1.0e-9)
  }
}

///|
test "treetime_root_to_tip_handles_nested_tree" {
  // 树: ((A:0.05,B:0.05):0.05,C:0.15); A/B 的 root-to-tip 都是 0.1, C 为 0.15.
  let dates : Map[String, Double] = Map([
    ("A", 2020.0),
    ("B", 2020.0),
    ("C", 2022.0),
  ])
  let fit = treetime_test_fit("((A:0.05,B:0.05):0.05,C:0.15);", dates)
  assert_eq(fit.tips.length(), 3)
  assert_true(fit.rate > 0.0249 && fit.rate < 0.0251)
  assert_true(fit.root_date > 2015.99 && fit.root_date < 2016.01)
}

///|
test "treetime_root_to_tip_missing_date_raises" {
  let dates : Map[String, Double] = Map([("A", 2020.0), ("B", 2021.0)])
  let raised = try {
    ignore(
      @src.phylo_clock_root_to_tip(
        @src.skbio_tree_read("(A:0.1,B:0.2,C:0.3);"),
        dates,
      ),
    )
    false
  } catch {
    @src.ClockError(_) => true
  }
  assert_true(raised)
}

///|
test "treetime_root_to_tip_non_increasing_raises" {
  let dates : Map[String, Double] = Map([
    ("A", 2022.0),
    ("B", 2021.0),
    ("C", 2020.0),
  ])
  let raised = try {
    ignore(
      @src.phylo_clock_root_to_tip(
        @src.skbio_tree_read("(A:0.1,B:0.2,C:0.3);"),
        dates,
      ),
    )
    false
  } catch {
    @src.ClockError(_) => true
  }
  assert_true(raised)
}

///|
test "treetime_dated_tree_marks_node_times" {
  let dates : Map[String, Double] = Map([
    ("A", 2020.0),
    ("B", 2021.0),
    ("C", 2022.0),
  ])
  let fit = treetime_test_fit("(A:0.1,B:0.2,C:0.3);", dates)
  let tree = @src.phylo_clock_dated_tree(
    @src.skbio_tree_read("(A:0.1,B:0.2,C:0.3);"),
    fit,
  )
  assert_true(tree.time > 2018.99 && tree.time < 2019.01)
  assert_eq(tree.children.length(), 3)
  let expected : Array[Double] = [2020.0, 2021.0, 2022.0]
  for i = 0; i < 3; i = i + 1 {
    assert_true((tree.children[i].time - expected[i]).abs() < 0.01)
  }
}

///|
test "treetime_ancestral_identical_tips" {
  let sequences : Map[String, String] = Map([
    ("A", "ACGT"),
    ("B", "ACGT"),
    ("C", "ACGT"),
  ])
  let result = @src.phylo_clock_ancestral(
    @src.skbio_tree_read("(A:0.01,B:0.01,C:0.01);"),
    sequences,
  ) catch {
    ClockError(message) => abort("ancestral reconstruction failed: " + message)
  }
  // 后序展开: 3 个叶节点在前, 根节点在最后.
  assert_eq(result.taxa.length(), 4)
  let root = result.taxa[result.taxa.length() - 1]
  assert_eq(root.sequence, "ACGT")
  for site in root.sites {
    assert_true(site.probability > 0.99)
  }
  assert_false(result.log_likelihood.is_nan())
}

///|
test "treetime_ancestral_reports_ambiguity" {
  let sequences : Map[String, String] = Map([("A", "A"), ("B", "C")])
  let result = @src.phylo_clock_ancestral(
    @src.skbio_tree_read("(A:0.01,B:0.01);"),
    sequences,
  ) catch {
    ClockError(message) => abort("ancestral reconstruction failed: " + message)
  }
  let root = result.taxa[2]
  // 两个等长叶节点观测到不同的碱基 -> 根节点 A/C 约各 0.5.
  let posterior = root.sites[0].posterior
  assert_true((posterior[0] - 0.5).abs() < 0.01)
  assert_true((posterior[1] - 0.5).abs() < 0.01)
  assert_true(posterior[2] < 0.01)
  assert_true(posterior[3] < 0.01)
}

///|
test "treetime_ancestral_missing_sequence_raises" {
  let sequences : Map[String, String] = Map([("A", "ACGT")])
  let raised = try {
    ignore(
      @src.phylo_clock_ancestral(
        @src.skbio_tree_read("(A:0.01,B:0.01);"),
        sequences,
      ),
    )
    false
  } catch {
    @src.ClockError(_) => true
  }
  assert_true(raised)
}

///|
test "treetime_root_to_tip_json_roundtrip" {
  let dates : Map[String, Double] = Map([
    ("A", 2020.0),
    ("B", 2021.0),
    ("C", 2022.0),
  ])
  let fit = treetime_test_fit("(A:0.1,B:0.2,C:0.3);", dates)
  let text = @src.phylo_clock_root_to_tip_json(fit)
  let parsed = @src.parse_json(text)
  assert_true(parsed is Some(_))
  assert_true(text.contains("\"rate\""))
  assert_true(text.contains("\"tips\""))
  assert_true(parsed.unwrap() is @src.JsonObj(_))
}
