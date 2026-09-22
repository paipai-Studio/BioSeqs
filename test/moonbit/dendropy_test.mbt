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
}///|
test "dendropy_weighted_robinson_foulds" {
  // tree1 = (A,(B,C)) with (B,C) branch 2.0 -> split {A}=2.0
  let t1 = tnode([
    tleaf_len("A", 1.0),
    tnode_len([tleaf_len("B", 3.0), tleaf_len("C", 4.0)], 2.0),
  ])
  // tree2 = ((A,C),B) with (A,C) branch 1.0 -> split {B}=1.0
  let t2 = tnode([
    tnode_len([tleaf_len("A", 2.0), tleaf_len("C", 3.0)], 1.0),
    tleaf_len("B", 4.0),
  ])
  assert_true(assert_close(@src.dendropy_weighted_robinson_foulds(t1, t2), 3.0, 1.0e-9))
}

///|
test "dendropy_patristic_distance" {
  // root: A (2.0), X (1.0); X: B (3.0), C (4.0)
  let root = tnode([tleaf_len("A", 2.0), tnode_len([tleaf_len("B", 3.0), tleaf_len("C", 4.0)], 1.0)])
  assert_true(assert_close(@src.dendropy_patristic_distance(root, "A", "B"), 6.0, 1.0e-9))
  assert_true(assert_close(@src.dendropy_patristic_distance(root, "A", "C"), 7.0, 1.0e-9))
  assert_true(assert_close(@src.dendropy_patristic_distance(root, "B", "C"), 7.0, 1.0e-9))
  let (names, matrix) = @src.dendropy_patristic_distance_matrix(root)
  assert_eq(names.length(), 3)
  assert_true(assert_close(matrix[0][1], 6.0, 1.0e-9))
  assert_true(assert_close(matrix[1][2], 7.0, 1.0e-9))
}

///|
test "dendropy_pybus_harvey_gamma" {
  // ultrametric comb tree, speciation heights 3, 2, 1
  let ab = tnode_len([tleaf_len("A", 1.0), tleaf_len("B", 1.0)], 1.0)
  let abc = tnode_len([ab, tleaf_len("C", 2.0)], 1.0)
  let root = tnode([abc, tleaf_len("D", 3.0)])
  // gamma = -sqrt(24)/9
  let expect = -(24.0).sqrt() / 9.0
  assert_true(assert_close(@src.dendropy_pybus_harvey_gamma(root), expect, 1.0e-6))
}