///|
/// Tests for the Bio.Cluster numerical clustering port (cluster.c).
/// Expected values were generated with Biopython 1.88 (Bio.Cluster) and
/// hand-verified. Note: like the C Clustering Library, the "Euclidean"
/// distance is the weighted mean of SQUARED differences (no square root).

///|
fn bc_approx(a : Double, b : Double, tol : Double) -> Bool {
  (a - b).abs() < tol
}

///|
fn bc_approx_vec(a : Array[Double], b : Array[Double], tol : Double) -> Bool {
  if a.length() != b.length() {
    return false
  }
  for i = 0; i < a.length(); i = i + 1 {
    if (a[i] - b[i]).abs() >= tol {
      return false
    }
  }
  true
}

///|
fn bc_approx_mat(
  a : Array[Array[Double]],
  b : Array[Array[Double]],
  tol : Double,
) -> Bool {
  if a.length() != b.length() {
    return false
  }
  for i = 0; i < a.length(); i = i + 1 {
    if !bc_approx_vec(a[i], b[i], tol) {
      return false
    }
  }
  true
}

///|
/// Distance between two single-element clusters via clusterdistance equals
/// the pairwise metric.

///|
test "bc_euclid_is_weighted_mean_squared_diff" {
  // [0,0] vs [3,4]: (w*9 + w*16)/(2w) = 12.5 (no square root!)
  let data : Array[Array[Double]] = [[0.0, 0.0], [3.0, 4.0]]
  let d = @src.clusterdistance(data, index1=[0], index2=[1], dist="e")
  assert_true(bc_approx(d, 12.5, 1.0e-12))
}

///|
test "bc_cityblock" {
  // (3 + 4)/2 = 3.5
  let data : Array[Array[Double]] = [[0.0, 0.0], [3.0, 4.0]]
  let d = @src.clusterdistance(data, index1=[0], index2=[1], dist="b")
  assert_true(bc_approx(d, 3.5, 1.0e-12))
}

///|
test "bc_correlation" {
  // perfectly correlated -> 0
  let pos : Array[Array[Double]] = [[1.0, 2.0, 3.0], [2.0, 4.0, 6.0]]
  assert_true(
    bc_approx(
      @src.clusterdistance(pos, index1=[0], index2=[1], dist="c"),
      0.0,
      1.0e-12,
    ),
  )
  // perfectly anti-correlated -> 2
  let neg : Array[Array[Double]] = [[1.0, 2.0, 3.0], [3.0, 2.0, 1.0]]
  assert_true(
    bc_approx(
      @src.clusterdistance(neg, index1=[0], index2=[1], dist="c"),
      2.0,
      1.0e-12,
    ),
  )
  // absolute correlation -> 0
  assert_true(
    bc_approx(
      @src.clusterdistance(neg, index1=[0], index2=[1], dist="a"),
      0.0,
      1.0e-12,
    ),
  )
}

///|
test "bc_ucorrelation" {
  // uncentered: r = 10/14 -> 1 - 5/7 = 2/7
  let data : Array[Array[Double]] = [[1.0, 2.0, 3.0], [3.0, 2.0, 1.0]]
  assert_true(
    bc_approx(
      @src.clusterdistance(data, index1=[0], index2=[1], dist="u"),
      0.2857142857142857,
      1.0e-12,
    ),
  )
  // absolute uncentered with y = -x: |r| = 1 -> 0
  let neg : Array[Array[Double]] = [[1.0, 2.0, 3.0], [-1.0, -2.0, -3.0]]
  assert_true(
    bc_approx(
      @src.clusterdistance(neg, index1=[0], index2=[1], dist="x"),
      0.0,
      1.0e-12,
    ),
  )
}

///|
test "bc_spearman" {
  let data : Array[Array[Double]] = [[1.0, 2.0, 3.0], [2.0, 1.0, 3.0]]
  assert_true(
    bc_approx(
      @src.clusterdistance(data, index1=[0], index2=[1], dist="s"),
      0.5,
      1.0e-12,
    ),
  )
  let rev : Array[Array[Double]] = [[1.0, 2.0, 3.0], [3.0, 2.0, 1.0]]
  assert_true(
    bc_approx(
      @src.clusterdistance(rev, index1=[0], index2=[1], dist="s"),
      2.0,
      1.0e-12,
    ),
  )
}

///|
test "bc_kendall" {
  let data : Array[Array[Double]] = [[1.0, 2.0, 3.0], [2.0, 1.0, 3.0]]
  // tau = (2 concordant - 1 discordant)/3 -> 1 - 1/3
  assert_true(
    bc_approx(
      @src.clusterdistance(data, index1=[0], index2=[1], dist="k"),
      0.6666666666666667,
      1.0e-12,
    ),
  )
  let rev : Array[Array[Double]] = [[1.0, 2.0, 3.0], [3.0, 2.0, 1.0]]
  assert_true(
    bc_approx(
      @src.clusterdistance(rev, index1=[0], index2=[1], dist="k"),
      2.0,
      1.0e-12,
    ),
  )
}

///|
test "bc_euclid_missing_values" {
  // y[1] is missing: only x[0], x[2] contribute and they are equal -> 0
  let data : Array[Array[Double]] = [[1.0, 2.0, 3.0], [1.0, 0.0, 3.0]]
  let mask : Array[Array[Int]] = [[1, 1, 1], [1, 0, 1]]
  let d = @src.clusterdistance(data, mask~, index1=[0], index2=[1], dist="e")
  assert_true(bc_approx(d, 0.0, 1.0e-12))
}

///|
test "bc_euclid_weighted" {
  // weight [1,1,4]: (1*0 + 1*0 + 4*4)/6 = 16/6
  let data : Array[Array[Double]] = [[1.0, 2.0, 3.0], [1.0, 2.0, 1.0]]
  let w : Array[Double] = [1.0, 1.0, 4.0]
  let d = @src.clusterdistance(data, weight=w, index1=[0], index2=[1], dist="e")
  assert_true(bc_approx(d, 2.6666666666666665, 1.0e-12))
}

///|
test "bc_distancematrix_lower_triangular" {
  // squared-euclid with weight normalization:
  // d(0,1) = (9+16)/2 = 12.5, d(0,2) = (36+64)/2 = 50, d(1,2) = 12.5
  let data : Array[Array[Double]] = [[0.0, 0.0], [3.0, 4.0], [6.0, 8.0]]
  let dm = @src.distancematrix(data)
  assert_eq(dm.length(), 3)
  assert_eq(dm[0].length(), 0)
  assert_eq(dm[1].length(), 1)
  assert_eq(dm[2].length(), 2)
  assert_true(bc_approx(dm[1][0], 12.5, 1.0e-12))
  assert_true(bc_approx(dm[2][0], 50.0, 1.0e-12))
  assert_true(bc_approx(dm[2][1], 12.5, 1.0e-12))
}

///|
/// Hierarchical clustering on data [[0],[1],[10],[11]] (squared Euclidean).

///|
test "bc_treecluster_four_linkages" {
  let data : Array[Array[Double]] = [[0.0], [1.0], [10.0], [11.0]]
  // single linkage: top = min pairwise squared = 81
  let s = @src.treecluster(data, method="s")
  assert_eq(s.length(), 3)
  assert_true(s.node(0).left == 0 && s.node(0).right == 1)
  assert_true(s.node(1).left == 2 && s.node(1).right == 3)
  assert_true(s.node(2).left == -1 && s.node(2).right == -2)
  assert_true(bc_approx(s.node(2).distance, 81.0, 1.0e-12))
  assert_eq(s.cut(nclusters=2), [0, 0, 1, 1])
  // complete linkage: top = max pairwise squared = 121
  let m = @src.treecluster(data, method="m")
  assert_true(m.node(0).left == 1 && m.node(0).right == 0)
  assert_true(m.node(1).left == 2 && m.node(1).right == 3)
  assert_true(m.node(2).left == -2 && m.node(2).right == -1)
  assert_true(bc_approx(m.node(2).distance, 121.0, 1.0e-12))
  assert_eq(m.cut(nclusters=2), [1, 1, 0, 0])
  // average linkage: top = (100+121+81+100)/4 = 100.5
  let a = @src.treecluster(data, method="a")
  assert_true(bc_approx(a.node(2).distance, 100.5, 1.0e-12))
  assert_eq(a.cut(nclusters=2), [1, 1, 0, 0])
  // centroid linkage: top = (0.5 - 10.5)^2 = 100
  let c = @src.treecluster(data, method="c")
  assert_true(c.node(0).left == 0 && c.node(0).right == 1)
  assert_true(c.node(1).left == 3 && c.node(1).right == 2)
  assert_true(bc_approx(c.node(2).distance, 100.0, 1.0e-12))
  assert_eq(c.cut(nclusters=2), [0, 0, 1, 1])
}

///|
test "bc_treecluster_three_elements" {
  let data : Array[Array[Double]] = [[0.0], [1.0], [5.0]]
  // d(1,0)=1, d(2,0)=25, d(2,1)=16 in squared space
  let s = @src.treecluster(data, method="s")
  assert_true(s.node(0).left == 0 && s.node(0).right == 1)
  assert_true(s.node(1).left == -1 && s.node(1).right == 2)
  assert_true(bc_approx(s.node(1).distance, 16.0, 1.0e-12))
  assert_eq(s.cut(nclusters=2), [0, 0, 1])
  let m = @src.treecluster(data, method="m")
  assert_true(m.node(1).left == 2 && m.node(1).right == -1)
  assert_true(bc_approx(m.node(1).distance, 25.0, 1.0e-12))
  assert_eq(m.cut(nclusters=2), [1, 1, 0])
  let a = @src.treecluster(data, method="a")
  assert_true(bc_approx(a.node(1).distance, 20.5, 1.0e-12))
  assert_eq(a.cut(nclusters=2), [1, 1, 0])
  let c = @src.treecluster(data, method="c")
  assert_true(bc_approx(c.node(1).distance, 20.25, 1.0e-12))
  assert_eq(c.cut(nclusters=2), [0, 0, 1])
}

///|
test "bc_treecluster_from_distance" {
  let dist : Array[Array[Double]] = [[], [1.0], [25.0, 16.0]]
  let s = @src.treecluster_from_distance(dist, method="s")
  assert_true(s.node(0).left == 0 && s.node(0).right == 1)
  assert_true(s.node(1).left == -1 && s.node(1).right == 2)
  assert_true(bc_approx(s.node(1).distance, 16.0, 1.0e-12))
  let m = @src.treecluster_from_distance(dist, method="m")
  assert_true(m.node(1).left == 2 && m.node(1).right == -1)
  assert_true(bc_approx(m.node(1).distance, 25.0, 1.0e-12))
  let a = @src.treecluster_from_distance(dist, method="a")
  assert_true(bc_approx(a.node(1).distance, 20.5, 1.0e-12))
}

///|
test "bc_tree_cut_variants" {
  let data : Array[Array[Double]] = [[0.0], [1.0], [10.0], [11.0]]
  let tree = @src.treecluster(data, method="s")
  assert_eq(tree.cut(nclusters=1), [0, 0, 0, 0])
  assert_eq(tree.cut(nclusters=2), [0, 0, 1, 1])
  assert_eq(tree.cut(nclusters=4), [0, 1, 2, 3])
  // Asking for more clusters than elements clips to singletons
  assert_eq(tree.cut(nclusters=99), [0, 1, 2, 3])
}

///|
test "bc_tree_sort" {
  let data : Array[Array[Double]] = [[0.0], [1.0], [10.0], [11.0]]
  let c = @src.treecluster(data, method="c")
  let idx = c.sort(order=[3.0, 1.0, 2.0, 4.0])
  assert_eq(idx, [1, 0, 2, 3])
  // sorting swaps children of node 0 (order 3 > 1)
  assert_true(c.node(0).left == 1 && c.node(0).right == 0)
  let s = @src.treecluster(data, method="s")
  let idx2 = s.sort(order=[1.0, 4.0, 3.0, 2.0])
  assert_eq(idx2, [0, 1, 3, 2])
  // sorting swaps children of node 1 (order 3 > 2)
  assert_true(s.node(1).left == 3 && s.node(1).right == 2)
}

///|
/// k-means / k-median clustering.

///|
test "bc_kcluster_initialid" {
  let data : Array[Array[Double]] = [
    [0.0, 0.0],
    [1.0, 1.0],
    [10.0, 10.0],
    [11.0, 11.0],
  ]
  let r = @src.kcluster(data, nclusters=2, initialid=[0, 0, 1, 1])
  assert_eq(r.clusterid(), [0, 0, 1, 1])
  assert_true(bc_approx(r.error(), 1.0, 1.0e-12))
  assert_eq(r.ifound(), 1)
  // different initial clustering converges to the same solution
  let r2 = @src.kcluster(data, nclusters=2, initialid=[0, 1, 0, 1])
  assert_eq(r2.clusterid(), [0, 0, 1, 1])
  assert_true(bc_approx(r2.error(), 1.0, 1.0e-12))
}

///|
test "bc_kcluster_median" {
  let data : Array[Array[Double]] = [
    [0.0, 0.0],
    [1.0, 1.0],
    [10.0, 10.0],
    [11.0, 11.0],
  ]
  let r = @src.kcluster(data, nclusters=2, initialid=[0, 0, 1, 1], method="m")
  assert_eq(r.clusterid(), [0, 0, 1, 1])
  assert_true(bc_approx(r.error(), 1.0, 1.0e-12))
}

///|
test "bc_kcluster_transpose" {
  // Cluster the two identical columns; they stay in separate clusters.
  let data : Array[Array[Double]] = [
    [0.0, 0.0],
    [1.0, 1.0],
    [10.0, 10.0],
    [11.0, 11.0],
  ]
  let r = @src.kcluster(data, nclusters=2, transpose=true, initialid=[0, 1])
  assert_eq(r.clusterid(), [0, 1])
  assert_true(bc_approx(r.error(), 0.0, 1.0e-12))
}

///|
test "bc_kcluster_random_passes" {
  @src.cluster_seed(42)
  let data : Array[Array[Double]] = [
    [0.0, 0.0],
    [1.0, 1.0],
    [10.0, 10.0],
    [11.0, 11.0],
  ]
  let r = @src.kcluster(data, nclusters=2, npass=5)
  // the optimal partition groups the two close pairs together
  assert_true(r.clusterid()[0] == r.clusterid()[1])
  assert_true(r.clusterid()[2] == r.clusterid()[3])
  assert_true(r.clusterid()[0] != r.clusterid()[2])
  assert_true(bc_approx(r.error(), 1.0, 1.0e-12))
  assert_true(r.ifound() >= 1)
}

///|
test "bc_kcluster_too_many_clusters" {
  let data : Array[Array[Double]] = [[1.0, 2.0]]
  let r = @src.kcluster(data, nclusters=2)
  assert_eq(r.ifound(), 0)
}

///|
/// k-medoids clustering.

///|
test "bc_kmedoids_initialid" {
  let dist : Array[Array[Double]] = [
    [0.0, 1.0, 4.0],
    [1.0, 0.0, 2.0],
    [4.0, 2.0, 0.0],
  ]
  let r = @src.kmedoids(dist, nclusters=2, initialid=[0, 0, 1])
  // cluster numbers are the item numbers of the medoids
  assert_eq(r.clusterid(), [0, 0, 2])
  assert_true(bc_approx(r.error(), 1.0, 1.0e-12))
  assert_eq(r.ifound(), 1)
  let r2 = @src.kmedoids(dist, nclusters=2, initialid=[0, 1, 0])
  assert_eq(r2.clusterid(), [0, 1, 1])
  assert_true(bc_approx(r2.error(), 2.0, 1.0e-12))
  assert_eq(r2.ifound(), 1)
}

///|
test "bc_kmedoids_random_passes" {
  @src.cluster_seed(7)
  let dist : Array[Array[Double]] = [
    [0.0, 1.0, 9.0, 10.0],
    [1.0, 0.0, 10.0, 11.0],
    [9.0, 10.0, 0.0, 1.0],
    [10.0, 11.0, 1.0, 0.0],
  ]
  let r = @src.kmedoids(dist, nclusters=2, npass=5)
  assert_true(r.error() <= 2.0 + 1.0e-12)
  let c = r.clusterid()
  // optimal: medoids are one of {0,1} and one of {2,3}
  let m0 = c[0]
  let m1 = c[2]
  assert_true(m0 <= 1 && m1 >= 2)
}

///|
/// Cluster centroids.

///|
test "bc_clustercentroids_mean_median" {
  let data : Array[Array[Double]] = [
    [1.0, 2.0],
    [3.0, 4.0],
    [5.0, 6.0],
    [7.0, 8.0],
  ]
  let mean = @src.clustercentroids(data, clusterid=[0, 0, 1, 1], method="a")
  assert_true(bc_approx_mat(mean.cdata(), [[2.0, 3.0], [6.0, 7.0]], 1.0e-12))
  assert_eq(mean.cmask(), [[1, 1], [1, 1]])
  let median = @src.clustercentroids(data, clusterid=[0, 0, 1, 1], method="m")
  assert_true(bc_approx_mat(median.cdata(), [[2.0, 3.0], [6.0, 7.0]], 1.0e-12))
  assert_eq(median.cmask(), [[1, 1], [1, 1]])
}

///|
test "bc_clustercentroids_with_mask" {
  let data : Array[Array[Double]] = [
    [1.0, 2.0],
    [3.0, 4.0],
    [5.0, 6.0],
    [7.0, 8.0],
  ]
  let mask : Array[Array[Int]] = [[1, 1], [1, 1], [1, 0], [1, 1]]
  let r = @src.clustercentroids(data, mask~, clusterid=[0, 1, 0, 1], method="a")
  // cluster 0: mean of [1,2] and [5,missing] = [3, 2]
  // cluster 1: mean of [3,4] and [7,8] = [5, 6]
  assert_true(bc_approx_mat(r.cdata(), [[3.0, 2.0], [5.0, 6.0]], 1.0e-12))
  assert_eq(r.cmask(), [[1, 1], [1, 1]])
}

///|
/// clusterdistance between clusters.

///|
test "bc_clusterdistance_methods" {
  let data : Array[Array[Double]] = [[0.0], [1.0], [5.0]]
  let idx1 : Array[Int] = [0, 1]
  let idx2 : Array[Int] = [2]
  // means: (0.5 - 5)^2 = 20.25
  assert_true(
    bc_approx(
      @src.clusterdistance(data, index1=idx1, index2=idx2, method="a"),
      20.25,
      1.0e-12,
    ),
  )
  // medians: median{0,1} = 0.5 -> 20.25
  assert_true(
    bc_approx(
      @src.clusterdistance(data, index1=idx1, index2=idx2, method="m"),
      20.25,
      1.0e-12,
    ),
  )
  // shortest pairwise: min(25,16) = 16
  assert_true(
    bc_approx(
      @src.clusterdistance(data, index1=idx1, index2=idx2, method="s"),
      16.0,
      1.0e-12,
    ),
  )
  // longest pairwise: max(25,16) = 25
  assert_true(
    bc_approx(
      @src.clusterdistance(data, index1=idx1, index2=idx2, method="x"),
      25.0,
      1.0e-12,
    ),
  )
  // average pairwise: (25+16)/2 = 20.5
  assert_true(
    bc_approx(
      @src.clusterdistance(data, index1=idx1, index2=idx2, method="v"),
      20.5,
      1.0e-12,
    ),
  )
  // invalid index -> -1
  let bad : Array[Int] = [99]
  assert_true(
    @src.clusterdistance(data, index1=bad, index2=idx2, method="a") == -1.0,
  )
}

///|
/// SOM clustering (structural tests; the training uses a deterministic RNG).

///|
test "bc_somcluster_structure" {
  @src.cluster_seed(12345)
  let data : Array[Array[Double]] = [
    [0.0, 0.0],
    [0.5, 0.2],
    [10.0, 10.0],
    [10.5, 10.2],
  ]
  let r = @src.somcluster(data, nxgrid=2, nygrid=1, inittau=0.02, niter=100)
  assert_eq(r.clusterid().length(), 4)
  assert_eq(r.celldata().length(), 2)
  assert_eq(r.celldata()[0].length(), 1)
  assert_eq(r.celldata()[0][0].length(), 2)
  // assignments are valid grid cells
  for i = 0; i < 4; i = i + 1 {
    let (ix, iy) = r.clusterid()[i]
    assert_true(ix >= 0 && ix < 2)
    assert_true(iy >= 0 && iy < 1)
  }
  // the trained cells converge to different directions of the data
  assert_true((r.celldata()[0][0][0] - r.celldata()[1][0][0]).abs() > 1.0e-6)
  assert_true((r.celldata()[0][0][1] - r.celldata()[1][0][1]).abs() > 1.0e-6)
  // each cell vector is normalized to Euclidean length sqrt(ndata)
  // (the C library divides by sqrt(sum/ndata))
  for ix = 0; ix < 2; ix = ix + 1 {
    let norm2 = r.celldata()[ix][0][0] * r.celldata()[ix][0][0] +
      r.celldata()[ix][0][1] * r.celldata()[ix][0][1]
    assert_true((norm2 - 2.0).abs() < 1.0e-9)
  }
  // reproducible for a fixed seed (same result on a second run)
  @src.cluster_seed(12345)
  let r2 = @src.somcluster(data, nxgrid=2, nygrid=1, inittau=0.02, niter=100)
  for i = 0; i < 4; i = i + 1 {
    let (ix, iy) = r.clusterid()[i]
    let (jx, jy) = r2.clusterid()[i]
    assert_true(ix == jx && iy == jy)
  }
  for ix = 0; ix < 2; ix = ix + 1 {
    assert_true(
      bc_approx_vec(r.celldata()[ix][0], r2.celldata()[ix][0], 1.0e-30),
    )
  }
}

///|
/// Principal component analysis.

///|
/// Column-wise comparison of two matrices that are equal up to a sign
/// flip of each column (SVD singular vectors are only defined up to the
/// sign of each component).
fn bc_approx_mat_signed(
  a : Array[Array[Double]],
  b : Array[Array[Double]],
  tol : Double,
) -> Bool {
  if a.length() != b.length() {
    return false
  }
  let ncols = if a.length() > 0 { a[0].length() } else { 0 }
  for j = 0; j < ncols; j = j + 1 {
    let mut direct = true
    let mut flipped = true
    for i = 0; i < a.length(); i = i + 1 {
      if (a[i][j] - b[i][j]).abs() >= tol {
        direct = false
      }
      if (a[i][j] + b[i][j]).abs() >= tol {
        flipped = false
      }
    }
    if !(direct || flipped) {
      return false
    }
  }
  true
}

///|
/// Row-wise comparison of two matrices that are equal up to a sign flip
/// of each row (the principal components carry one singular vector per
/// row, whose sign is arbitrary).
fn bc_approx_mat_signed_rows(
  a : Array[Array[Double]],
  b : Array[Array[Double]],
  tol : Double,
) -> Bool {
  if a.length() != b.length() {
    return false
  }
  for i = 0; i < a.length(); i = i + 1 {
    if !bc_approx_vec(a[i], b[i], tol) && !bc_approx_vec_neg(a[i], b[i], tol) {
      return false
    }
  }
  true
}

///|
fn bc_approx_vec_neg(
  a : Array[Double],
  b : Array[Double],
  tol : Double,
) -> Bool {
  if a.length() != b.length() {
    return false
  }
  for i = 0; i < a.length(); i = i + 1 {
    if (a[i] + b[i]).abs() >= tol {
      return false
    }
  }
  true
}

///|
test "bc_pca_doc_example" {
  let data : Array[Array[Double]] = [
    [0.0, 0.0, 0.0],
    [1.0, 0.0, 0.0],
    [7.0, 3.0, 0.0],
    [4.0, 2.0, 6.0],
  ]
  let r = @src.cluster_pca(data)
  // column means
  assert_true(bc_approx_vec(r.columnmean(), [3.0, 1.25, 1.5], 1.0e-12))
  // eigenvalues (singular values of the centered matrix), largest first
  assert_true(
    bc_approx_vec(
      r.eigenvalues(),
      [6.369049291152344, 4.806117424762438, 0.2940177312113694],
      1.0e-9,
    ),
  )
  // coordinates: one singular vector per column; the sign of each
  // component is arbitrary (values from Biopython 1.88)
  assert_true(
    bc_approx_mat_signed(
      r.coordinates(),
      [
        [-3.565823812757457, -0.24471741725381072, 0.19368511059061072],
        [-2.7783303999811664, 0.21251913785010768, -0.21958097501363233],
        [3.1154058590920206, 3.4067624477985814, 0.026757380374613095],
        [3.2287483536466013, -3.3745641683948797, -0.000861515951591482],
      ],
      1.0e-9,
    ),
  )
  // principal components: one singular vector per row (sign arbitrary)
  assert_true(
    bc_approx_mat_signed_rows(
      r.pc(),
      [
        [0.7874934127762903, 0.38959192747181554, 0.47756911005921093],
        [0.4572365551039173, 0.15027465977499002, -0.8765570485177864],
        [-0.4132660856042432, 0.9086449563379009, -0.05979536613350543],
      ],
      1.0e-9,
    ),
  )
  // reconstruction: data ~= columnmean + coordinates x pc
  for i = 0; i < 4; i = i + 1 {
    for j = 0; j < 3; j = j + 1 {
      let mut v = r.columnmean()[j]
      for k = 0; k < 3; k = k + 1 {
        v = v + r.coordinates()[i][k] * r.pc()[k][j]
      }
      assert_true((v - data[i][j]).abs() < 1.0e-9)
    }
  }
}

///|
test "bc_pca_wide_matrix" {
  let data : Array[Array[Double]] = [[1.0, 2.0, 3.0, 4.0], [5.0, 6.0, 7.0, 8.0]]
  let r = @src.cluster_pca(data)
  assert_true(bc_approx_vec(r.columnmean(), [3.0, 4.0, 5.0, 6.0], 1.0e-12))
  assert_true(bc_approx_vec(r.eigenvalues(), [5.65685424949238, 0.0], 1.0e-9))
  assert_true(bc_approx_mat(r.coordinates(), [[4.0, 0.0], [-4.0, 0.0]], 1.0e-9))
  assert_true(
    bc_approx_mat(
      r.pc(),
      [
        [-0.5, -0.5, -0.5, -0.5],
        [-0.5, 0.8333333333333334, -0.16666666666666666, -0.16666666666666666],
      ],
      1.0e-9,
    ),
  )
}

///|
test "bc_pca_empty_raises" {
  let empty : Array[Array[Double]] = []
  let ok = try {
    ignore(@src.cluster_pca(empty))
    false
  } catch {
    ClusterError(_) => true
  }
  assert_true(ok)
}

///|
/// Record (Cluster/TreeView file format).

///|
test "bc_record_read_full" {
  let content = "UNIQID\tNAME\tGWEIGHT\tGORDER\texp1\texp2\texp3\ngene1\tdesc A\t1.5\t2.0\t2.5\t-1.0\t0.5\ngene2\tdesc B\t1.5\t2.0\t0.0\t3.5\t-2.0\nEWEIGHT\t\t\t\t1.0\t2.0\t1.0\nEORDER\t\t\t\t2.0\t1.0\t3.0\n"
  let rec = @src.bio_cluster_read(content)
  assert_eq(rec.uniqid(), "UNIQID")
  assert_eq(rec.expid(), ["exp1", "exp2", "exp3"])
  assert_eq(rec.geneid(), ["gene1", "gene2"])
  assert_eq(rec.genename(), Some(["desc A", "desc B"]))
  assert_true(
    match rec.gweight() {
      Some(w) => bc_approx_vec(w, [1.5, 1.5], 1.0e-12)
      None => false
    },
  )
  assert_true(
    match rec.gorder() {
      Some(w) => bc_approx_vec(w, [2.0, 2.0], 1.0e-12)
      None => false
    },
  )
  assert_true(
    match rec.eweight() {
      Some(w) => bc_approx_vec(w, [1.0, 2.0, 1.0], 1.0e-12)
      None => false
    },
  )
  assert_true(
    match rec.eorder() {
      Some(w) => bc_approx_vec(w, [2.0, 1.0, 3.0], 1.0e-12)
      None => false
    },
  )
  assert_true(
    bc_approx_mat(rec.data(), [[2.5, -1.0, 0.5], [0.0, 3.5, -2.0]], 1.0e-12),
  )
  // no missing values -> no mask
  assert_eq(rec.mask(), None)
}

///|
test "bc_record_read_missing_values" {
  let content = "UID\texp1\texp2\ng1\t1.0\t\ng2\t2.0\t3.0\n"
  let rec = @src.bio_cluster_read(content)
  assert_eq(rec.uniqid(), "UID")
  assert_eq(rec.expid(), ["exp1", "exp2"])
  assert_eq(rec.geneid(), ["g1", "g2"])
  assert_true(bc_approx_mat(rec.data(), [[1.0, 0.0], [2.0, 3.0]], 1.0e-12))
  assert_eq(rec.mask(), Some([[1, 0], [1, 1]]))
  assert_eq(rec.genename(), None)
  assert_eq(rec.gweight(), None)
  assert_eq(rec.eweight(), None)
}

///|
test "bc_record_read_inconsistent_columns" {
  let content = "UID\texp1\texp2\ng1\t1.0\n"
  let ok = try {
    ignore(@src.bio_cluster_read(content))
    false
  } catch {
    ClusterError(_) => true
  }
  assert_true(ok)
}

///|
test "bc_record_read_empty" {
  let ok = try {
    ignore(@src.bio_cluster_read(""))
    false
  } catch {
    ClusterError(_) => true
  }
  assert_true(ok)
}

///|
test "bc_record_methods" {
  let content = "UID\texp1\texp2\ng1\t0.0\t0.0\ng2\t1.0\t1.0\ng3\t10.0\t10.0\ng4\t11.0\t11.0\n"
  let rec = @src.bio_cluster_read(content)
  let tree = rec.treecluster(method="s")
  assert_eq(tree.length(), 3)
  let kc = rec.kcluster(nclusters=2, initialid=[0, 0, 1, 1])
  assert_eq(kc.clusterid(), [0, 0, 1, 1])
  // distance between the first and last gene: (121 + 121)/2 = 121
  let d = rec.clusterdistance(index1=[0], index2=[3], dist="e")
  assert_true(bc_approx(d, 121.0, 1.0e-9))
  // distance matrix over the rows (default transpose=false): 4 x 4 lower
  // triangular; d(g1,g0) = (1+1)/2 = 1, d(g3,g0) = (121+121)/2 = 121
  let dm = rec.distancematrix()
  assert_eq(dm.length(), 4)
  assert_true(bc_approx(dm[1][0], 1.0, 1.0e-9))
  assert_true(bc_approx(dm[3][0], 121.0, 1.0e-9))
}
