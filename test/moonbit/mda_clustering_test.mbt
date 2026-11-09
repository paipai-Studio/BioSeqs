///|
/// Tests for mda_clustering.mbt (trajectory clustering and averaging).

///|
test "mda_cluster_kmeans_two_groups" {
  // Six single-atom frames: three at x=0, three at x=10. K-means with
  // 2 clusters must split them 3/3 with centroids at frames 0 and 3.
  let frames = [
    [(0.0, 0.0, 0.0)],
    [(0.0, 0.0, 0.0)],
    [(0.0, 0.0, 0.0)],
    [(10.0, 0.0, 0.0)],
    [(10.0, 0.0, 0.0)],
    [(10.0, 0.0, 0.0)],
  ]
  let r = @src.mda_cluster_kmeans(frames, n_clusters=2)
  assert_eq(r.n_clusters, 2)
  assert_eq(r.assignments.length(), 6)
  assert_eq(r.assignments[0], 0)
  assert_eq(r.assignments[2], 0)
  assert_eq(r.assignments[3], 1)
  assert_eq(r.assignments[5], 1)
  assert_eq(r.centroids.length(), 2)
  assert_eq(r.cluster_sizes[0], 3)
  assert_eq(r.cluster_sizes[1], 3)
}

///|
test "mda_cluster_kmeans_empty" {
  let r = @src.mda_cluster_kmeans([])
  assert_eq(r.n_clusters, 0)
  assert_eq(r.assignments.length(), 0)
}

///|
test "mda_cluster_kmeans_single_frame_clamps_clusters" {
  // One frame with n_clusters=3: nc is clamped to 1.
  let r = @src.mda_cluster_kmeans([[(1.0, 2.0, 3.0)]], n_clusters=3)
  assert_eq(r.n_clusters, 1)
  assert_eq(r.assignments.length(), 1)
  assert_eq(r.assignments[0], 0)
  assert_eq(r.cluster_sizes[0], 1)
}

///|
test "mda_average_structure_all_and_subset" {
  let f0 = [(0.0, 0.0, 0.0), (2.0, 0.0, 0.0)]
  let f1 = [(2.0, 0.0, 0.0), (0.0, 0.0, 0.0)]
  // Average over both frames: every coordinate becomes (1, 0, 0).
  let avg = @src.mda_average_structure([f0, f1])
  assert_eq(avg.length(), 2)
  assert_true(@moonbit.assert_num_close(avg[0].0, 1.0, 1.0e-9))
  assert_true(@moonbit.assert_num_close(avg[0].2, 0.0, 1.0e-9))
  assert_true(@moonbit.assert_num_close(avg[1].0, 1.0, 1.0e-9))
  // Restricted to frame 1 only.
  let avg1 = @src.mda_average_structure([f0, f1], frame_indices=[1])
  assert_true(@moonbit.assert_num_close(avg1[0].0, 2.0, 1.0e-9))
  assert_true(@moonbit.assert_num_close(avg1[1].0, 0.0, 1.0e-9))
}

///|
test "mda_average_structure_empty_atoms" {
  let avg = @src.mda_average_structure([[]])
  assert_eq(avg.length(), 0)
}

///|
test "mda_density_peaks_assigns_self_cluster" {
  // Three well-separated single-atom frames: every frame becomes its own
  // center (n=3 -> 3 centers), so each centroid frame maps to its own
  // cluster index regardless of the density ordering.
  let frames = [[(0.0, 0.0, 0.0)], [(0.25, 0.0, 0.0)], [(0.5, 0.0, 0.0)]]
  let r = @src.mda_density_peaks(frames, dc=0.3)
  assert_eq(r.n_clusters, 3)
  assert_eq(r.centroids.length(), 3)
  assert_eq(r.cluster_sizes[0], 1)
  assert_eq(r.cluster_sizes[1], 1)
  assert_eq(r.cluster_sizes[2], 1)
  for c = 0; c < 3; c = c + 1 {
    assert_eq(r.assignments[r.centroids[c]], c)
  }
}

///|
test "mda_density_peaks_empty" {
  let r = @src.mda_density_peaks([])
  assert_eq(r.n_clusters, 0)
}

///|
test "mda_cluster_centroids_from_kmeans" {
  // Reusing the two-group kmeans result: cluster 0 averages to x=0,
  // cluster 1 to x=10.
  let frames = [
    [(0.0, 0.0, 0.0)],
    [(0.0, 0.0, 0.0)],
    [(10.0, 0.0, 0.0)],
    [(10.0, 0.0, 0.0)],
  ]
  let r = @src.mda_cluster_kmeans(frames, n_clusters=2)
  let cents = @src.mda_cluster_centroids(frames, r)
  assert_eq(cents.length(), 2)
  assert_true(@moonbit.assert_num_close(cents[0][0].0, 0.0, 1.0e-9))
  assert_true(@moonbit.assert_num_close(cents[1][0].0, 10.0, 1.0e-9))
}

///|
test "mda_cluster_hierarchical_two_groups" {
  // Frames at x = 0, 0.5, 10, 10.5: single-linkage merges (0,1) and
  // (2,3) first -> assignments [0,0,1,1], sizes [2,2], centroids [0,2].
  let frames = [
    [(0.0, 0.0, 0.0)],
    [(0.5, 0.0, 0.0)],
    [(10.0, 0.0, 0.0)],
    [(10.5, 0.0, 0.0)],
  ]
  let r = @src.mda_cluster_hierarchical(frames, n_clusters=2)
  assert_eq(r.n_clusters, 2)
  assert_eq(r.assignments[0], 0)
  assert_eq(r.assignments[1], 0)
  assert_eq(r.assignments[2], 1)
  assert_eq(r.assignments[3], 1)
  assert_eq(r.cluster_sizes[0], 2)
  assert_eq(r.cluster_sizes[1], 2)
  assert_eq(r.centroids[0], 0)
  assert_eq(r.centroids[1], 2)
}

///|
test "mda_cluster_hierarchical_empty" {
  let r = @src.mda_cluster_hierarchical([])
  assert_eq(r.n_clusters, 0)
}
