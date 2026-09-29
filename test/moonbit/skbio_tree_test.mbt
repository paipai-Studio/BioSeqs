///|
/// Tests for the scikit-bio tree + UniFrac port. Expected values are the exact
/// doctest outputs from scikit-bio 0.7.3 `diversity/beta/_unifrac.py`.

///|
test "skbio_tree_read_and_tips" {
  let tree = @src.skbio_tree_read("(a:1,b:2);")
  assert_false(tree.is_tip())
  assert_eq(tree.children.length(), 2)
  assert_true(tree.children[0].is_tip())
  assert_true(tree.children[1].is_tip())
  assert_eq(tree.tips().length(), 2)
  assert_eq(tree.tips()[0], "a")
  assert_eq(tree.tips()[1], "b")
}

///|
test "skbio_unifrac" {
  let tree = @src.skbio_tree_read(
    "(((((U1:0.5,U2:0.5):0.5,U3:1.0):1.0):0.0," +
    "(U4:0.75,(U5:0.5,((U6:0.33,U7:0.62):0.5,U8:0.5):0.5):0.5):1.25):0.0)root;",
  )
  let taxa = ["U1", "U2", "U3", "U4", "U5", "U6", "U7", "U8"]
  let u = [1.0, 0.0, 0.0, 4.0, 1.0, 2.0, 3.0, 0.0]
  let v = [0.0, 1.0, 1.0, 6.0, 0.0, 1.0, 0.0, 0.0]
  // doctest: uu = 0.37, wu = 1.54
  let uu = @src.skbio_unweighted_unifrac(u, v, taxa, tree)
  assert_true(assert_close(uu, 0.37, 0.005))
  let wu = @src.skbio_weighted_unifrac(u, v, taxa, tree)
  assert_true(assert_close(wu, 1.54, 0.005))
  // identical samples -> distance 0
  let z = @src.skbio_unweighted_unifrac(u, u, taxa, tree)
  assert_true(assert_close(z, 0.0, 1.0e-12))
}

// ============================================================
// faith_pd
// ============================================================

test "skbio_faith_pd" {
  // root(0) → node1(1.0) → [A(0.5), B(0.5)] + C(2.0)
  let a = @src.SkbioNode::new(Some("A"), 0.5, [])
  let b = @src.SkbioNode::new(Some("B"), 0.5, [])
  let node1 = @src.SkbioNode::new(None, 1.0, [a, b])
  let c = @src.SkbioNode::new(Some("C"), 2.0, [])
  let root = @src.SkbioNode::new(None, 0.0, [node1, c])
  // A、C observed，B 未 observed → PD = A(0.5) + node1(1.0) + C(2.0) = 3.5
  let pd = @src.skbio_faith_pd([1.0, 0.0, 1.0], ["A", "B", "C"], root)
  assert_true(assert_close(pd, 3.5, 1.0e-9))
}
