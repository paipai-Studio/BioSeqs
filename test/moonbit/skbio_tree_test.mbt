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