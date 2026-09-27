///|
/// Tests for maximum-likelihood tree inference (ml_tree.mbt, JC69 model).

fn mt_leaf(name : String) -> @src.Clade {
  @src.Clade::new(name=Some(name))
}

fn mt_internal(children : Array[@src.Clade]) -> @src.Clade {
  @src.Clade::new(clades=children)
}

fn mt_tree(root : @src.Clade) -> @src.Tree {
  @src.Tree::new(root, rooted=false)
}

fn mt_map(pairs : Array[(String, String)]) -> Map[String, String] {
  let m : Map[String, String] = Map([])
  for p in pairs {
    m[p.0] = p.1
  }
  m
}

///|
/// 正确分组 ((A,B),(C,D)) 的似然应高于错误分组 ((A,C),(B,D))。
test "ml_tree_log_likelihood_prefers_correct_topology" {
  let m = mt_map([("A", "AAAA"), ("B", "AAAA"), ("C", "CCCC"), ("D", "CCCC")])
  let correct = mt_tree(
    mt_internal([
      mt_internal([mt_leaf("A"), mt_leaf("B")]),
      mt_internal([mt_leaf("C"), mt_leaf("D")]),
    ]),
  )
  let wrong = mt_tree(
    mt_internal([
      mt_internal([mt_leaf("A"), mt_leaf("C")]),
      mt_internal([mt_leaf("B"), mt_leaf("D")]),
    ]),
  )
  let ll_correct = @src.ml_tree_log_likelihood(correct, m)
  let ll_wrong = @src.ml_tree_log_likelihood(wrong, m)
  assert_true(ll_correct > ll_wrong)
}

///|
/// 完全一致的序列：一点替换都没有，似然应为负有限数（非 0、非 NaN）。
test "ml_tree_log_likelihood_identical" {
  let m = mt_map([("A", "ACGT"), ("B", "ACGT"), ("C", "ACGT")])
  let tree = mt_tree(
    mt_internal([mt_leaf("A"), mt_internal([mt_leaf("B"), mt_leaf("C")])]),
  )
  let ll = @src.ml_tree_log_likelihood(tree, m)
  assert_true(ll < 0.0)
}

///|
/// NNI 搜索应至少不恶化似然（从错误拓扑出发）。
test "ml_nni_search_monotonic" {
  let m = mt_map([("A", "AAAA"), ("B", "AAAA"), ("C", "CCCC"), ("D", "CCCC")])
  let wrong = mt_tree(
    mt_internal([
      mt_internal([mt_leaf("A"), mt_leaf("C")]),
      mt_internal([mt_leaf("B"), mt_leaf("D")]),
    ]),
  )
  let before = @src.ml_tree_log_likelihood(wrong, m)
  let result = @src.ml_nni_search(wrong, m)
  let after = @src.ml_tree_log_likelihood(result, m)
  assert_true(after >= before - 1.0e-9)
}
