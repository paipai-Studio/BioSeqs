///|
/// Tests for mda_traj_cluster.mbt (RMSD-based trajectory clustering).

///|
test "mda_trajectory_kmeans_two_groups" {
  // Six single-atom frames: three at x=0, three at x=10. With k=2 the
  // initial centroids are frames 0 and 5; the split is 3/3.
  let frames = [
    [[0.0, 0.0, 0.0]],
    [[0.0, 0.0, 0.0]],
    [[0.0, 0.0, 0.0]],
    [[10.0, 0.0, 0.0]],
    [[10.0, 0.0, 0.0]],
    [[10.0, 0.0, 0.0]],
  ]
  let r = @src.mda_trajectory_kmeans(frames, 2)
  assert_eq(r.n_clusters, 2)
  assert_eq(r.assignments[0], 0)
  assert_eq(r.assignments[2], 0)
  assert_eq(r.assignments[3], 1)
  assert_eq(r.assignments[5], 1)
  assert_eq(r.sizes[0], 3)
  assert_eq(r.sizes[1], 3)
  assert_eq(r.centroids.length(), 2)
}

///|
test "mda_trajectory_kmeans_degenerate_inputs" {
  // Empty trajectory -> 0 clusters.
  assert_eq(@src.mda_trajectory_kmeans([], 3).n_clusters, 0)
  // k <= 0 -> 0 clusters.
  assert_eq(@src.mda_trajectory_kmeans([[[0.0, 0.0, 0.0]]], 0).n_clusters, 0)
}

///|
test "mda_trajectory_kmeans_k_exceeds_n" {
  // Three frames at x = 0, 1, 2 with k = 5: actual_k clamps to 3. The
  // integer-spaced initial centroids are frames 0, 0, 1 (i*(n-1)/(k-1)),
  // giving assignments [0, 2, 2] with sizes [1, 0, 2].
  let frames = [[[0.0, 0.0, 0.0]], [[1.0, 0.0, 0.0]], [[2.0, 0.0, 0.0]]]
  let r = @src.mda_trajectory_kmeans(frames, 5)
  assert_eq(r.n_clusters, 3)
  assert_eq(r.assignments[0], 0)
  assert_eq(r.assignments[1], 2)
  assert_eq(r.assignments[2], 2)
  assert_eq(r.sizes[0], 1)
  assert_eq(r.sizes[1], 0)
  assert_eq(r.sizes[2], 2)
}

///|
test "mda_trajectory_cluster_by_cutoff_two_groups" {
  // Frames at x = 0, 0.5, 10, 10.5 with cutoff 1: two tight groups.
  let frames = [
    [[0.0, 0.0, 0.0]],
    [[0.5, 0.0, 0.0]],
    [[10.0, 0.0, 0.0]],
    [[10.5, 0.0, 0.0]],
  ]
  let r = @src.mda_trajectory_cluster_by_cutoff(frames, 1.0)
  assert_eq(r.n_clusters, 2)
  assert_eq(r.assignments[0], 0)
  assert_eq(r.assignments[1], 0)
  assert_eq(r.assignments[2], 1)
  assert_eq(r.assignments[3], 1)
  assert_eq(r.centroids[0], 0)
  assert_eq(r.centroids[1], 2)
  assert_eq(r.sizes[0], 2)
  assert_eq(r.sizes[1], 2)
}

///|
test "mda_trajectory_cluster_by_cutoff_chain_connectivity" {
  // Chained frames at x = 0, 0.8, 1.6 with cutoff 1: everything joins
  // one cluster through single-link connectivity.
  let frames = [[[0.0, 0.0, 0.0]], [[0.8, 0.0, 0.0]], [[1.6, 0.0, 0.0]]]
  let r = @src.mda_trajectory_cluster_by_cutoff(frames, 1.0)
  assert_eq(r.n_clusters, 1)
  assert_eq(r.assignments[0], 0)
  assert_eq(r.assignments[1], 0)
  assert_eq(r.assignments[2], 0)
  assert_eq(r.sizes[0], 3)
}

///|
test "mda_trajectory_cluster_by_cutoff_empty" {
  let r = @src.mda_trajectory_cluster_by_cutoff([], 1.0)
  assert_eq(r.n_clusters, 0)
  assert_eq(r.assignments.length(), 0)
}

///|
test "mda_cluster_populations_fractions" {
  let frames = [
    [[0.0, 0.0, 0.0]],
    [[0.5, 0.0, 0.0]],
    [[10.0, 0.0, 0.0]],
    [[10.5, 0.0, 0.0]],
  ]
  let r = @src.mda_trajectory_cluster_by_cutoff(frames, 1.0)
  let pops = @src.mda_cluster_populations(r)
  assert_eq(pops.length(), 2)
  assert_true(@moonbit.assert_num_close(pops[0], 0.5, 1.0e-12))
  assert_true(@moonbit.assert_num_close(pops[1], 0.5, 1.0e-12))
  // Empty cluster result -> empty populations.
  assert_eq(
    @src.mda_cluster_populations(@src.mda_trajectory_kmeans([], 3)).length(),
    0,
  )
}
