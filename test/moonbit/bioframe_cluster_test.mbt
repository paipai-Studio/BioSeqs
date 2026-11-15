///|
/// Tests for bioframe_cluster.mbt — interval clustering.

// ---- Empty input -----------------------------------------------------------

test "bf_cluster: empty input → empty output" {
  let (rows, n) = @src.bf_cluster([])
  assert_eq(rows.length(), 0)
  assert_eq(n, 0)
}

// ---- No overlap → separate clusters ----------------------------------------

test "bf_cluster: no overlap → each interval independent" {
  let a = [
    @src.BedInterval::new("chr1", 0, 100),
    @src.BedInterval::new("chr1", 200, 300),
    @src.BedInterval::new("chr1", 500, 600),
  ]
  let (rows, n) = @src.bf_cluster(a)
  assert_eq(rows.length(), 3)
  assert_eq(n, 3)
  assert_eq(rows[0].cluster_id, 0)
  assert_eq(rows[1].cluster_id, 1)
  assert_eq(rows[2].cluster_id, 2)
}

// ---- Fully overlapping → one cluster ---------------------------------------

test "bf_cluster: overlapping → one cluster" {
  let a = [
    @src.BedInterval::new("chr1", 0, 100),
    @src.BedInterval::new("chr1", 50, 150),
    @src.BedInterval::new("chr1", 120, 200),
  ]
  let (rows, n) = @src.bf_cluster(a)
  assert_eq(n, 1)
  for r in rows {
    assert_eq(r.cluster_id, 0)
  }
}

// ---- Book-ended intervals do NOT cluster with min_dist=None ----------------

test "bf_cluster: book-ended do NOT cluster by default (min_dist=None)" {
  let a = [
    @src.BedInterval::new("chr1", 0, 100),
    @src.BedInterval::new("chr1", 100, 200),
  ]
  let (_, n) = @src.bf_cluster(a)
  assert_eq(n, 2)
}

// ---- min_dist parameter controls gap tolerance ------------------------------

test "bf_cluster: min_dist=5 clusters intervals within 5bp gap" {
  let a = [
    @src.BedInterval::new("chr1", 0, 100),
    @src.BedInterval::new("chr1", 104, 200),
  ]
  let (_, n) = @src.bf_cluster(a, min_dist=Some(5))
  assert_eq(n, 1)
}

test "bf_cluster: min_dist=5 does NOT cluster 6bp gap" {
  let a = [
    @src.BedInterval::new("chr1", 0, 100),
    @src.BedInterval::new("chr1", 106, 200),
  ]
  let (_, n) = @src.bf_cluster(a, min_dist=Some(5))
  assert_eq(n, 2)
}

// ---- Book-ended intervals cluster with min_dist=0 --------------------------

test "bf_cluster: min_dist=0 clusters book-ended intervals" {
  let a = [
    @src.BedInterval::new("chr1", 0, 100),
    @src.BedInterval::new("chr1", 100, 200),
  ]
  let (_, n) = @src.bf_cluster(a, min_dist=Some(0))
  assert_eq(n, 1)
}

// ---- Different chromosomes always separate ---------------------------------

test "bf_cluster: different chroms → separate clusters" {
  let a = [
    @src.BedInterval::new("chr1", 0, 100),
    @src.BedInterval::new("chr2", 0, 100),
  ]
  let (rows, n) = @src.bf_cluster(a)
  assert_eq(n, 2)
  assert_eq(rows[0].cluster_id, 0)
  assert_eq(rows[1].cluster_id, 1)
}

// ---- Mixed chroms + overlap ------------------------------------------------

test "bf_cluster: mixed chroms with multiple groups" {
  let a = [
    @src.BedInterval::new("chr1", 0, 100),
    @src.BedInterval::new("chr1", 50, 150),
    @src.BedInterval::new("chr2", 0, 100),
    @src.BedInterval::new("chr2", 200, 300),
  ]
  let (rows, n) = @src.bf_cluster(a)
  assert_eq(n, 3)
  assert_eq(rows[0].cluster_id, 0)
  assert_eq(rows[1].cluster_id, 0)
  assert_eq(rows[2].cluster_id, 1)
  assert_eq(rows[3].cluster_id, 2)
}

// ---- bf_cluster_with_min_dist wrapper --------------------------------------

test "bf_cluster_with_min_dist: convenience wrapper" {
  let a = [
    @src.BedInterval::new("chr1", 0, 100),
    @src.BedInterval::new("chr1", 105, 200),
    @src.BedInterval::new("chr1", 300, 400),
  ]
  let (rows, n) = @src.bf_cluster_with_min_dist(a, 10)
  assert_eq(n, 2)
  assert_eq(rows[0].cluster_id, 0)
  assert_eq(rows[1].cluster_id, 0)
  assert_eq(rows[2].cluster_id, 1)
}

// ---- n_clusters is correct -------------------------------------------------

test "bf_cluster: n_clusters matches max cluster_id + 1" {
  let a = [
    @src.BedInterval::new("chr1", 0, 50),
    @src.BedInterval::new("chr1", 100, 150),
    @src.BedInterval::new("chr1", 200, 250),
    @src.BedInterval::new("chr1", 300, 350),
    @src.BedInterval::new("chr1", 400, 450),
  ]
  let (rows, n) = @src.bf_cluster(a)
  assert_eq(n, 5)
  assert_eq(rows[0].cluster_id, 0)
  assert_eq(rows[1].cluster_id, 1)
  assert_eq(rows[2].cluster_id, 2)
  assert_eq(rows[3].cluster_id, 3)
  assert_eq(rows[4].cluster_id, 4)
}
