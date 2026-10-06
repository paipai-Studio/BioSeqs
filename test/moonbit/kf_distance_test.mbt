///|
/// Tests for Kuhner-Felsenstein branch-score distance.

///|
fn kf_leaf(name : String, bl : Double) -> @src.Clade {
  @src.Clade::new(name=Some(name), branch_length=Some(bl))
}

///|
fn kf_internal(children : Array[@src.Clade], bl : Double) -> @src.Clade {
  @src.Clade::new(name=Some("internal"), clades=children, branch_length=Some(bl))
}

///|
/// 4-leaf balanced: ((A:1,B:1):2,(C:1,D:1):2):0
fn kf_balanced4() -> @src.Clade {
  kf_internal([
    kf_internal([kf_leaf("A", 1.0), kf_leaf("B", 1.0)], 2.0),
    kf_internal([kf_leaf("C", 1.0), kf_leaf("D", 1.0)], 2.0),
  ], 0.0)
}

///|
/// 4-leaf caterpillar: ((A:1,B:1):2,C:3):0 with D:3
/// Actually: ((A:1,B:1):3,(C:2,D:2):1):0
fn kf_caterpillar4() -> @src.Clade {
  kf_internal([
    kf_internal([kf_leaf("A", 1.0), kf_leaf("B", 1.0)], 3.0),
    kf_internal([kf_leaf("C", 2.0), kf_leaf("D", 2.0)], 1.0),
  ], 0.0)
}

///|
/// 5-leaf tree: ((A:1,B:1):2,(C:1,(D:1,E:1):2):3):0
fn kf_tree5() -> @src.Clade {
  kf_internal([
    kf_internal([kf_leaf("A", 1.0), kf_leaf("B", 1.0)], 2.0),
    kf_internal([
      kf_leaf("C", 1.0),
      kf_internal([kf_leaf("D", 1.0), kf_leaf("E", 1.0)], 2.0),
    ], 3.0),
  ], 0.0)
}

///|
/// 6-leaf tree: ((A:1,B:1):2,((C:1,D:1):2,(E:1,F:1):2):3):0
fn kf_tree6() -> @src.Clade {
  kf_internal([
    kf_internal([kf_leaf("A", 1.0), kf_leaf("B", 1.0)], 2.0),
    kf_internal([
      kf_internal([kf_leaf("C", 1.0), kf_leaf("D", 1.0)], 2.0),
      kf_internal([kf_leaf("E", 1.0), kf_leaf("F", 1.0)], 2.0),
    ], 3.0),
  ], 0.0)
}

///|
/// Same topology as kf_balanced4 but different branch lengths
fn kf_balanced4_diff_bl() -> @src.Clade {
  kf_internal([
    kf_internal([kf_leaf("A", 2.0), kf_leaf("B", 3.0)], 1.0),
    kf_internal([kf_leaf("C", 0.5), kf_leaf("D", 1.5)], 4.0),
  ], 0.0)
}

///|
/// Same topology and branch lengths as kf_balanced4
fn kf_balanced4_same() -> @src.Clade {
  kf_internal([
    kf_internal([kf_leaf("A", 1.0), kf_leaf("B", 1.0)], 2.0),
    kf_internal([kf_leaf("C", 1.0), kf_leaf("D", 1.0)], 2.0),
  ], 0.0)
}

///|
test "kf: identical trees have zero distance" {
  let t1 = kf_balanced4()
  let t2 = kf_balanced4_same()
  let d = @src.dendropy_kuhner_felsenstein_distance(t1, t2)
  assert_true(d < 1.0e-10)
}

///|
test "kf: same tree to itself is zero" {
  let t = kf_tree5()
  let d = @src.dendropy_kuhner_felsenstein_distance(t, t)
  assert_true(d < 1.0e-10)
}

///|
test "kf: different branch lengths nonzero distance" {
  let t1 = kf_balanced4()
  let t2 = kf_balanced4_diff_bl()
  let d = @src.dendropy_kuhner_felsenstein_distance(t1, t2)
  assert_true(d > 0.0)
}

///|
test "kf: different topology nonzero distance" {
  let t1 = kf_balanced4()
  let t2 = kf_caterpillar4()
  let d = @src.dendropy_kuhner_felsenstein_distance(t1, t2)
  assert_true(d > 0.0)
}

///|
test "kf: symmetric property" {
  let t1 = kf_balanced4()
  let t2 = kf_caterpillar4()
  let d12 = @src.dendropy_kuhner_felsenstein_distance(t1, t2)
  let d21 = @src.dendropy_kuhner_felsenstein_distance(t2, t1)
  assert_true((d12 - d21) * (d12 - d21) < 1.0e-20)
}

///|
test "kf: different taxon count returns -1" {
  let t1 = kf_balanced4()
  let t2 = kf_tree5()
  let d = @src.dendropy_kuhner_felsenstein_distance(t1, t2)
  assert_eq(d, -1.0)
}

///|
test "kf: only internal branches matter" {
  // Two trees with same internal branch lengths but different pendant lengths
  // KF distance should be 0 (pendant branches excluded)
  let t1 = kf_internal([
    kf_internal([kf_leaf("A", 1.0), kf_leaf("B", 1.0)], 2.0),
    kf_internal([kf_leaf("C", 1.0), kf_leaf("D", 1.0)], 2.0),
  ], 0.0)
  let t2 = kf_internal([
    kf_internal([kf_leaf("A", 99.0), kf_leaf("B", 88.0)], 2.0),
    kf_internal([kf_leaf("C", 77.0), kf_leaf("D", 66.0)], 2.0),
  ], 0.0)
  let d = @src.dendropy_kuhner_felsenstein_distance(t1, t2)
  assert_true(d < 1.0e-10)
}

///|
test "kf: root branch excluded" {
  // Root branch length should not affect KF distance
  let t1 = kf_internal([
    kf_internal([kf_leaf("A", 1.0), kf_leaf("B", 1.0)], 2.0),
    kf_internal([kf_leaf("C", 1.0), kf_leaf("D", 1.0)], 2.0),
  ], 5.0)
  let t2 = kf_internal([
    kf_internal([kf_leaf("A", 1.0), kf_leaf("B", 1.0)], 2.0),
    kf_internal([kf_leaf("C", 1.0), kf_leaf("D", 1.0)], 2.0),
  ], 99.0)
  let d = @src.dendropy_kuhner_felsenstein_distance(t1, t2)
  assert_true(d < 1.0e-10)
}

///|
test "kf: 5-leaf tree distance to itself" {
  let t = kf_tree5()
  let d = @src.dendropy_kuhner_felsenstein_distance(t, t)
  assert_true(d < 1.0e-10)
}

///|
test "kf: 6-leaf tree distance to itself" {
  let t = kf_tree6()
  let d = @src.dendropy_kuhner_felsenstein_distance(t, t)
  assert_true(d < 1.0e-10)
}

///|
test "kf: 6-leaf different branch lengths" {
  let t1 = kf_tree6()
  let t2 = kf_internal([
    kf_internal([kf_leaf("A", 2.0), kf_leaf("B", 3.0)], 1.0),
    kf_internal([
      kf_internal([kf_leaf("C", 0.5), kf_leaf("D", 1.5)], 4.0),
      kf_internal([kf_leaf("E", 2.0), kf_leaf("F", 1.0)], 0.5),
    ], 5.0),
  ], 0.0)
  let d = @src.dendropy_kuhner_felsenstein_distance(t1, t2)
  assert_true(d > 0.0)
}

///|
test "kf: triangle inequality with 3 trees" {
  let t1 = kf_balanced4()
  let t2 = kf_caterpillar4()
  let t3 = kf_balanced4_diff_bl()
  let d12 = @src.dendropy_kuhner_felsenstein_distance(t1, t2)
  let d23 = @src.dendropy_kuhner_felsenstein_distance(t2, t3)
  let d13 = @src.dendropy_kuhner_felsenstein_distance(t1, t3)
  // d(t1,t3) <= d(t1,t2) + d(t2,t3)
  assert_true(d13 <= d12 + d23 + 1.0e-10)
}

///|
test "kf: completely different topology 4 leaves" {
  // ((A,C),(B,D)) vs ((A,B),(C,D))
  let t1 = kf_internal([
    kf_internal([kf_leaf("A", 1.0), kf_leaf("B", 1.0)], 2.0),
    kf_internal([kf_leaf("C", 1.0), kf_leaf("D", 1.0)], 2.0),
  ], 0.0)
  let t2 = kf_internal([
    kf_internal([kf_leaf("A", 1.0), kf_leaf("C", 1.0)], 2.0),
    kf_internal([kf_leaf("B", 1.0), kf_leaf("D", 1.0)], 2.0),
  ], 0.0)
  let d = @src.dendropy_kuhner_felsenstein_distance(t1, t2)
  // Both internal branches are 2.0, completely different splits
  // Each tree has 1 non-trivial internal split (the two cherries)
  // Actually for 4 taxa: each tree has 2 internal nodes (root + 1 internal)
  // Non-trivial: k=2, n-k=2, so the cherry nodes count
  // t1 splits: {A,B} and {C,D}, each with bl=2
  // t2 splits: {A,C} and {B,D}, each with bl=2
  // No overlap, so d = sqrt(2^2 + 2^2 + 2^2 + 2^2) = sqrt(16) = 4
  assert_true(d > 0.0)
  // Check exact value: 4 non-overlapping splits each with bl=2
  // d = sqrt(4 * 4) = sqrt(16) = 4.0
  assert_true((d - 4.0) * (d - 4.0) < 1.0e-10)
}

///|
test "kf: non-negativity" {
  let t1 = kf_tree5()
  let t2 = kf_tree6()
  // Different taxon counts return -1, so test with same count
  let t3 = kf_internal([
    kf_internal([kf_leaf("A", 2.0), kf_leaf("B", 0.5)], 1.0),
    kf_internal([
      kf_leaf("C", 3.0),
      kf_internal([kf_leaf("D", 0.5), kf_leaf("E", 2.0)], 1.0),
    ], 2.0),
  ], 0.0)
  let d = @src.dendropy_kuhner_felsenstein_distance(kf_tree5(), t3)
  assert_true(d >= 0.0)
}

///|
test "kf: 3-leaf tree has no non-trivial splits" {
  // ((A,B),C) — only internal node is the cherry, but n=3, k=2, n-k=1 < 2
  // So no non-trivial splits, distance is always 0
  let t1 = kf_internal([
    kf_internal([kf_leaf("A", 1.0), kf_leaf("B", 1.0)], 2.0),
    kf_leaf("C", 3.0),
  ], 0.0)
  let t2 = kf_internal([
    kf_internal([kf_leaf("A", 5.0), kf_leaf("B", 3.0)], 1.0),
    kf_leaf("C", 2.0),
  ], 0.0)
  let d = @src.dendropy_kuhner_felsenstein_distance(t1, t2)
  assert_true(d < 1.0e-10)
}
