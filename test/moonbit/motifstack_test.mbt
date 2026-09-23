///|
/// Tests for motifStack — motif stacking and clustering.

///|
test "stackmotif construction and len" {
  let m = @src.StackMotif::new(
    name="M1",
    pwm=[[0.5, 0.5], [0.5, 0.5], [0.0, 0.0], [0.0, 0.0]],
  )
  assert_eq(m.name, "M1")
  assert_eq(m.len(), 2)
}

///|
test "stackmotif distance identical is zero" {
  let m1 = @src.StackMotif::new(
    name="A",
    pwm=[[0.5, 0.5], [0.5, 0.5], [0.0, 0.0], [0.0, 0.0]],
  )
  assert_true(@src.motifstack_distance(m1, m1) < 1.0e-12)
}

///|
test "stackmotif distance matrix symmetric" {
  let motifs = @src.motifstack_sample_motifs()
  let dm = @src.motifstack_distance_matrix(motifs)
  assert_eq(dm.length(), 3)
  for i in 0..<3 {
    assert_eq(dm[i][i], 0.0)
    for j in (i + 1)..<3 {
      assert_eq(dm[i][j], dm[j][i])
    }
  }
}

///|
test "stackmotif distance GATA1 vs GATA2 small" {
  let motifs = @src.motifstack_sample_motifs()
  let d = @src.motifstack_distance(motifs[0], motifs[1])
  assert_true(d >= 0.0)
  assert_true(d < 0.5)
}

///|
test "stackmotif distance GATA1 vs CTCF large" {
  let motifs = @src.motifstack_sample_motifs()
  let d = @src.motifstack_distance(motifs[0], motifs[2])
  assert_true(d >= 0.5)
}

///|
test "stackmotif hclust three leaves" {
  let motifs = @src.motifstack_sample_motifs()
  let dm = @src.motifstack_distance_matrix(motifs)
  match @src.motifstack_hclust(dm) {
    Some(root) => {
      let leaves = root.leaves()
      assert_eq(leaves.length(), 3)
      assert_true(root.height() >= 0.0)
    }
    None => assert_true(false)
  }
}

///|
test "stackmotif hclust empty returns none" {
  let dm : Array[Array[Double]] = Array(capacity=0)
  assert_true(@src.motifstack_hclust(dm) is None)
}

///|
test "stackmotif hclust single returns leaf" {
  let dm = [[0.0]]
  match @src.motifstack_hclust(dm) {
    Some(root) => {
      let leaves = root.leaves()
      assert_eq(leaves.length(), 1)
      assert_eq(leaves[0], 0)
    }
    None => assert_true(false)
  }
}

///|
test "stackmotif stack order clusters GATA1 and GATA2 together" {
  let motifs = @src.motifstack_sample_motifs()
  let order = @src.motifstack_stack_order(motifs)
  assert_eq(order.length(), 3)
  // GATA1 (0) and GATA2 (1) should be closer than either to CTCF (2)
  let d_01 = @src.motifstack_distance(motifs[0], motifs[1])
  let d_02 = @src.motifstack_distance(motifs[0], motifs[2])
  let d_12 = @src.motifstack_distance(motifs[1], motifs[2])
  assert_true(d_01 < d_02)
  assert_true(d_01 < d_12)
}

///|
test "stackmotif column ic positive for non-uniform column" {
  let m = @src.StackMotif::new(
    name="info",
    pwm=[[1.0, 0.0], [0.0, 1.0], [0.0, 0.0], [0.0, 0.0]],
  )
  let ics = @src.motifstack_column_ic(m)
  assert_eq(ics.length(), 2)
  assert_true(ics[0] > 0.0)
  assert_true(ics[1] > 0.0)
}

///|
test "stackmotif column ic zero for uniform column" {
  let m = @src.StackMotif::new(
    name="uniform",
    pwm=[[0.25, 0.25], [0.25, 0.25], [0.25, 0.25], [0.25, 0.25]],
  )
  let ics = @src.motifstack_column_ic(m)
  assert_true(ics[0] < 0.01)
  assert_true(ics[1] < 0.01)
}

///|
test "stackmotif ascii logo non-empty" {
  let motifs = @src.motifstack_sample_motifs()
  let logo = @src.motifstack_ascii_logo(motifs[0], rows=4)
  assert_true(logo.length() > 0)
  assert_true(logo.contains("\n"))
}

///|
test "stackmotif stack ascii contains motif names" {
  let motifs = @src.motifstack_sample_motifs()
  let stacked = @src.motifstack_stack_ascii(motifs, rows=4)
  assert_true(stacked.contains("GATA1"))
  assert_true(stacked.contains("GATA2"))
  assert_true(stacked.contains("CTCF"))
}

///|
test "stackmotif flatten length 4 * width" {
  let m = @src.StackMotif::new(
    name="x",
    pwm=[[0.1, 0.2, 0.3], [0.1, 0.2, 0.3], [0.1, 0.2, 0.3], [0.1, 0.2, 0.3]],
  )
  let flat = @src.motifstack_flatten(m)
  assert_eq(flat.length(), 12)
}

///|
test "stackmotif node leaves merge" {
  let leaf1 = @src.MotifStackNode::Leaf(0)
  let leaf2 = @src.MotifStackNode::Leaf(1)
  let merge = @src.MotifStackNode::Merge(leaf1, leaf2, 0.5)
  let leaves = merge.leaves()
  assert_eq(leaves.length(), 2)
  assert_true(leaves.contains(0))
  assert_true(leaves.contains(1))
}

///|
test "stackmotif node height" {
  let leaf = @src.MotifStackNode::Leaf(7)
  assert_eq(leaf.height(), 0.0)
  let merge = @src.MotifStackNode::Merge(
    @src.MotifStackNode::Leaf(0),
    @src.MotifStackNode::Leaf(1),
    0.75,
  )
  assert_true((merge.height() - 0.75).abs() < 1.0e-12)
}