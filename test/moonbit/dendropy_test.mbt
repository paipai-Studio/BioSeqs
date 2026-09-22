///|
/// Tests for the DendroPy tree-measure port. Values are hand-computed.

///|
fn tleaf(n : String) -> @src.Clade {
  @src.Clade::new(name=Some(n))
}

///|
fn tleaf_len(n : String, bl : Double) -> @src.Clade {
  @src.Clade::new(name=Some(n), branch_length=Some(bl))
}

///|
fn tnode(children : Array[@src.Clade]) -> @src.Clade {
  @src.Clade::new(clades=children)
}

///|
fn tnode_len(children : Array[@src.Clade], bl : Double) -> @src.Clade {
  @src.Clade::new(clades=children, branch_length=Some(bl))
}

///|
test "dendropy_shape_balanced" {
  let t = tnode([tnode([tleaf("A"), tleaf("B")]), tnode([tleaf("C"), tleaf("D")])])
  assert_eq(@src.dendropy_colless(t), 0)
  assert_eq(@src.dendropy_sackin(t), 8)
  assert_true(assert_close(@src.dendropy_n_bar(t), 2.0, 1.0e-9))
  assert_true(assert_close(@src.dendropy_b1(t), 2.0, 1.0e-9))
}

///|
test "dendropy_shape_caterpillar" {
  let t = tnode([tnode([tnode([tleaf("A"), tleaf("B")]), tleaf("C")]), tleaf("D")])
  assert_eq(@src.dendropy_colless(t), 3)
  assert_eq(@src.dendropy_sackin(t), 9)
  assert_true(assert_close(@src.dendropy_n_bar(t), 2.25, 1.0e-9))
  assert_true(assert_close(@src.dendropy_b1(t), 1.5, 1.0e-9))
}

///|
test "dendropy_tree_height" {
  let ab = tnode_len([tleaf_len("A", 1.0), tleaf_len("B", 1.0)], 1.0)
  let cd = tnode_len([tleaf_len("C", 1.0), tleaf_len("D", 1.0)], 1.0)
  let root = tnode([ab, cd])
  assert_true(assert_close(@src.dendropy_tree_height(root), 2.0, 1.0e-9))
}

///|
test "dendropy_symmetric_difference" {
  let t1 = tnode([tleaf("A"), tnode([tleaf("B"), tleaf("C")])])
  let t2 = tnode([tnode([tleaf("A"), tleaf("C")]), tleaf("B")])
  assert_eq(@src.dendropy_symmetric_difference(t1, t1), 0)
  assert_eq(@src.dendropy_symmetric_difference(t1, t2), 2)
}