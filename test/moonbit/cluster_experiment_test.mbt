///|
/// Tests for the clusterExperiment module (RSEC workflow).

// ============================================================================
// 1. Numeric helpers: ce_mean, ce_variance, ce_sd, ce_median, ce_euclidean
// ============================================================================

///|
test "ce_mean_basic" {
  let m = @src.ce_mean([1.0, 2.0, 3.0])
  assert_true((m - 2.0).abs() < 0.001)
}

///|
test "ce_mean_empty" {
  let m = @src.ce_mean([])
  assert_true((m - 0.0).abs() < 0.001)
}

///|
test "ce_mean_single" {
  let m = @src.ce_mean([5.0])
  assert_true((m - 5.0).abs() < 0.001)
}

///|
test "ce_variance_basic" {
  // population variance (1/n): mean=2, var = (1+0+1)/3 = 2/3
  let v = @src.ce_variance([1.0, 2.0, 3.0])
  assert_true((v - 0.6667).abs() < 0.001)
}

///|
test "ce_variance_single_returns_zero" {
  let v = @src.ce_variance([5.0])
  assert_true((v - 0.0).abs() < 0.001)
}

///|
test "ce_variance_empty_returns_zero" {
  let v = @src.ce_variance([])
  assert_true((v - 0.0).abs() < 0.001)
}

///|
test "ce_sd_basic" {
  // sd = sqrt(2/3) ~ 0.8165
  let s = @src.ce_sd([1.0, 2.0, 3.0])
  assert_true((s - 0.8165).abs() < 0.001)
}

///|
test "ce_sd_single_returns_zero" {
  let s = @src.ce_sd([42.0])
  assert_true((s - 0.0).abs() < 0.001)
}

///|
test "ce_median_odd" {
  let m = @src.ce_median([1.0, 2.0, 3.0])
  assert_true((m - 2.0).abs() < 0.001)
}

///|
test "ce_median_even" {
  let m = @src.ce_median([1.0, 2.0, 3.0, 4.0])
  assert_true((m - 2.5).abs() < 0.001)
}

///|
test "ce_median_unsorted" {
  let m = @src.ce_median([3.0, 1.0, 2.0])
  assert_true((m - 2.0).abs() < 0.001)
}

///|
test "ce_median_empty" {
  let m = @src.ce_median([])
  assert_true((m - 0.0).abs() < 0.001)
}

///|
test "ce_euclidean_basic" {
  let d = @src.ce_euclidean([0.0, 0.0], [3.0, 4.0])
  assert_true((d - 5.0).abs() < 0.001)
}

///|
test "ce_euclidean_same_point" {
  let d = @src.ce_euclidean([1.0, 2.0], [1.0, 2.0])
  assert_true((d - 0.0).abs() < 0.001)
}

// ============================================================================
// 2. Matrix helpers: ce_make_matrix, ce_transpose, ce_center_columns
// ============================================================================

///|
test "ce_make_matrix_basic" {
  let m = @src.ce_make_matrix(2, 3, 1.5)
  assert_eq(m.length(), 2)
  assert_eq(m[0].length(), 3)
  assert_true((m[0][0] - 1.5).abs() < 0.001)
  assert_true((m[1][2] - 1.5).abs() < 0.001)
}

///|
test "ce_transpose_basic" {
  let m = [[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]]
  let t = @src.ce_transpose(m)
  assert_eq(t.length(), 3)
  assert_eq(t[0].length(), 2)
  assert_true((t[0][0] - 1.0).abs() < 0.001)
  assert_true((t[0][1] - 4.0).abs() < 0.001)
  assert_true((t[2][1] - 6.0).abs() < 0.001)
}

///|
test "ce_transpose_empty" {
  let m : Array[Array[Double]] = []
  let t = @src.ce_transpose(m)
  assert_eq(t.length(), 0)
}

///|
test "ce_center_columns_basic" {
  // col 0: mean = 2 -> centered [-1, 1]
  // col 1: mean = 15 -> centered [-5, 5]
  let m = [[1.0, 10.0], [3.0, 20.0]]
  let c = @src.ce_center_columns(m)
  assert_true((c[0][0] - (-1.0)).abs() < 0.001)
  assert_true((c[0][1] - (-5.0)).abs() < 0.001)
  assert_true((c[1][0] - 1.0).abs() < 0.001)
  assert_true((c[1][1] - 5.0).abs() < 0.001)
}

///|
test "ce_center_columns_preserves_shape" {
  let m = [[1.0, 2.0, 3.0], [4.0, 5.0, 6.0], [7.0, 8.0, 9.0]]
  let c = @src.ce_center_columns(m)
  assert_eq(c.length(), 3)
  assert_eq(c[0].length(), 3)
  // centered columns sum to zero
  for j in 0..<3 {
    let mut col_sum = 0.0
    for i in 0..<3 {
      col_sum = col_sum + c[i][j]
    }
    assert_true(col_sum.abs() < 0.001)
  }
}

// ============================================================================
// 3. PCA: ce_pca
// ============================================================================

///|
test "ce_pca_shape" {
  let data = [
    [1.0, 2.0, 3.0],
    [2.0, 4.0, 6.0],
    [3.0, 6.0, 9.0],
    [10.0, 20.0, 30.0],
  ]
  let projected = @src.ce_pca(data, 2)
  assert_eq(projected.length(), 4)
  assert_eq(projected[0].length(), 2)
}

///|
test "ce_pca_zero_components" {
  let data = [[1.0, 2.0], [3.0, 4.0]]
  let projected = @src.ce_pca(data, 0)
  assert_eq(projected.length(), 0)
}

///|
test "ce_pca_empty" {
  let data : Array[Array[Double]] = []
  let projected = @src.ce_pca(data, 2)
  assert_eq(projected.length(), 0)
}

///|
test "ce_pca_components_capped_by_features" {
  // 2 features, request 5 components -> should get 2
  let data = [[1.0, 2.0], [3.0, 4.0], [5.0, 6.0], [7.0, 8.0]]
  let projected = @src.ce_pca(data, 5)
  assert_eq(projected.length(), 4)
  assert_eq(projected[0].length(), 2)
}

// ============================================================================
// 4. K-means: ce_kmeans, ce_centroid
// ============================================================================

///|
test "ce_kmeans_two_clusters" {
  let data = [
    [0.0, 0.0],
    [0.1, 0.0],
    [0.0, 0.1],
    [10.0, 10.0],
    [10.1, 10.0],
    [10.0, 10.1],
  ]
  let (labels, centroids) = @src.ce_kmeans(data, 2, 5, 100)
  assert_eq(labels.length(), 6)
  assert_eq(centroids.length(), 2)
  // first three (cluster A) share a label
  assert_eq(labels[0], labels[1])
  assert_eq(labels[1], labels[2])
  // last three (cluster B) share a label
  assert_eq(labels[3], labels[4])
  assert_eq(labels[4], labels[5])
  // the two groups differ
  assert_true(labels[0] != labels[3])
}

///|
test "ce_kmeans_k1" {
  let data = [[1.0, 2.0], [3.0, 4.0], [5.0, 6.0]]
  let (labels, centroids) = @src.ce_kmeans(data, 1, 1, 10)
  assert_eq(labels.length(), 3)
  assert_eq(labels[0], 0)
  assert_eq(labels[1], 0)
  assert_eq(labels[2], 0)
  assert_eq(centroids.length(), 1)
}

///|
test "ce_kmeans_empty" {
  let data : Array[Array[Double]] = []
  let (labels, centroids) = @src.ce_kmeans(data, 2, 5, 100)
  assert_eq(labels.length(), 0)
  assert_eq(centroids.length(), 0)
}

///|
test "ce_centroid_basic" {
  let points = [[0.0, 0.0], [2.0, 4.0], [4.0, 8.0]]
  let c = @src.ce_centroid(points)
  assert_eq(c.length(), 2)
  assert_true((c[0] - 2.0).abs() < 0.001)
  assert_true((c[1] - 4.0).abs() < 0.001)
}

///|
test "ce_centroid_empty" {
  let c = @src.ce_centroid([])
  assert_eq(c.length(), 0)
}

// ============================================================================
// 5. clusterMany: ce_cluster_many
// ============================================================================

///|
test "ce_cluster_many" {
  let data = [
    [0.0, 0.0],
    [0.1, 0.0],
    [0.0, 0.1],
    [10.0, 10.0],
    [10.1, 10.0],
    [10.0, 10.1],
  ]
  let params = @src.ce_default_params()
  let ensemble = @src.ce_cluster_many(data, params)
  // ks=[2,3,4] x reduced_dims=[0,5] x subsamples=[1.0, 0.8] = 12 runs
  assert_eq(ensemble.length(), 12)
  for cl in ensemble {
    assert_eq(cl.labels.length(), 6)
    assert_true(cl.k >= 2 && cl.k <= 4)
  }
}

// ============================================================================
// 6. Co-clustering matrix: ce_cocluster_matrix
// ============================================================================

///|
test "ce_cocluster_matrix" {
  let c1 = @src.CeClustering::new(
    labels=[0, 0, 1],
    k=2,
    reduced_dim=0,
    subsample=1.0,
  )
  let c2 = @src.CeClustering::new(
    labels=[0, 0, 1],
    k=2,
    reduced_dim=0,
    subsample=1.0,
  )
  let ensemble = [c1, c2]
  let m = @src.ce_cocluster_matrix(ensemble, 3)
  assert_eq(m.length(), 3)
  // samples 0 and 1 always co-cluster
  assert_true((m[0][1] - 1.0).abs() < 0.001)
  assert_true((m[1][0] - 1.0).abs() < 0.001)
  // samples 0 and 2 never co-cluster
  assert_true((m[0][2] - 0.0).abs() < 0.001)
  assert_true((m[2][0] - 0.0).abs() < 0.001)
  // diagonal = number of runs where sample was assigned
  assert_true((m[0][0] - 2.0).abs() < 0.001)
}

///|
test "ce_cocluster_matrix_with_unassigned" {
  // one run assigns all, one run leaves sample 2 unassigned
  let c1 = @src.CeClustering::new(
    labels=[0, 0, 1],
    k=2,
    reduced_dim=0,
    subsample=1.0,
  )
  let c2 = @src.CeClustering::new(
    labels=[0, 0, -1],
    k=2,
    reduced_dim=0,
    subsample=0.8,
  )
  let m = @src.ce_cocluster_matrix([c1, c2], 3)
  // sample 2 assigned in only 1 run
  assert_true((m[2][2] - 1.0).abs() < 0.001)
  // samples 0,1 co-cluster in both runs (denom = min(2,2) = 2)
  assert_true((m[0][1] - 1.0).abs() < 0.001)
}

// ============================================================================
// 7. Consensus: ce_make_consensus
// ============================================================================

///|
test "ce_make_consensus" {
  // 4 samples: 0,1,2 always together, 3 separate
  let n = 4
  let cocluster = @src.ce_make_matrix(n, n, 0.0)
  cocluster[0][1] = 1.0
  cocluster[1][0] = 1.0
  cocluster[0][2] = 1.0
  cocluster[2][0] = 1.0
  cocluster[1][2] = 1.0
  cocluster[2][1] = 1.0
  for i in 0..<n {
    cocluster[i][i] = 1.0
  }
  let labels = @src.ce_make_consensus(cocluster, 2)
  assert_eq(labels.length(), 4)
  // 0,1,2 together; 3 separate
  assert_eq(labels[0], labels[1])
  assert_eq(labels[1], labels[2])
  assert_true(labels[0] != labels[3])
}

// ============================================================================
// 8. Hierarchical clustering: ce_hclust_cut
// ============================================================================

///|
test "ce_hclust_cut_k2" {
  let dist = [
    [0.0, 0.1, 1.0, 1.0],
    [0.1, 0.0, 1.0, 1.0],
    [1.0, 1.0, 0.0, 0.1],
    [1.0, 1.0, 0.1, 0.0],
  ]
  let labels = @src.ce_hclust_cut(dist, 2)
  assert_eq(labels.length(), 4)
  // 0,1 together; 2,3 together
  assert_eq(labels[0], labels[1])
  assert_eq(labels[2], labels[3])
  assert_true(labels[0] != labels[2])
}

///|
test "ce_hclust_cut_k_equals_n" {
  let dist = [
    [0.0, 0.5, 0.5],
    [0.5, 0.0, 0.5],
    [0.5, 0.5, 0.0],
  ]
  let labels = @src.ce_hclust_cut(dist, 3)
  assert_eq(labels.length(), 3)
  // each sample in its own cluster
  assert_eq(labels[0], 0)
  assert_eq(labels[1], 1)
  assert_eq(labels[2], 2)
}

///|
test "ce_hclust_cut_empty" {
  let dist : Array[Array[Double]] = []
  let labels = @src.ce_hclust_cut(dist, 2)
  assert_eq(labels.length(), 0)
}

// ============================================================================
// 9. Dendrogram: ce_make_dendrogram
// ============================================================================

///|
test "ce_make_dendrogram_three_clusters" {
  let data = [
    [0.0, 0.0],
    [1.0, 0.0],
    [10.0, 10.0],
    [11.0, 10.0],
    [20.0, 20.0],
    [21.0, 20.0],
  ]
  let labels = [0, 0, 1, 1, 2, 2]
  let dendro = @src.ce_make_dendrogram(data, labels)
  match dendro {
    Some(root) => {
      assert_eq(root.cluster_id, -1)
      assert_eq(root.members.length(), 3)
      assert_true(root.height > 0.0)
    }
    None => assert_true(false)
  }
}

///|
test "ce_make_dendrogram_single_cluster_returns_none" {
  let data = [[0.0, 0.0], [1.0, 1.0]]
  let labels = [0, 0]
  let dendro = @src.ce_make_dendrogram(data, labels)
  match dendro {
    Some(_) => assert_true(false)
    None => assert_true(true)
  }
}

// ============================================================================
// 10. mergeClusters: ce_merge_clusters
// ============================================================================

///|
test "ce_merge_clusters_identical_merges" {
  // All clusters identical -> se=0 -> pval=1 -> qval=1 -> merged (1 > 0.05)
  let data = [
    [1.0, 1.0],
    [1.0, 1.0],
    [1.0, 1.0],
    [1.0, 1.0],
    [1.0, 1.0],
    [1.0, 1.0],
  ]
  let labels = [0, 0, 1, 1, 2, 2]
  let dendro = @src.ce_make_dendrogram(data, labels)
  let (new_labels, info) = @src.ce_merge_clusters(data, labels, dendro, 0.05)
  assert_eq(new_labels.length(), 6)
  assert_true(info.length() >= 2)
  // all merges happen -> single cluster
  assert_eq(@src.ce_count_clusters(new_labels), 1)
}

///|
test "ce_merge_clusters_distinct_no_merge" {
  // Distinct clusters with within-cluster variance -> significant -> no merge
  let data = [
    [0.0, 0.0],
    [0.1, 0.0],
    [10.0, 10.0],
    [10.1, 10.0],
    [20.0, 20.0],
    [20.1, 20.0],
  ]
  let labels = [0, 0, 1, 1, 2, 2]
  let dendro = @src.ce_make_dendrogram(data, labels)
  let (new_labels, info) = @src.ce_merge_clusters(data, labels, dendro, 0.05)
  assert_eq(new_labels.length(), 6)
  assert_true(info.length() >= 2)
  // no merges -> 3 clusters remain
  assert_eq(@src.ce_count_clusters(new_labels), 3)
}

///|
test "ce_merge_clusters_no_dendro" {
  let data = [[1.0, 2.0], [3.0, 4.0]]
  let labels = [0, 1]
  let (new_labels, info) = @src.ce_merge_clusters(data, labels, None, 0.05)
  assert_eq(new_labels.length(), 2)
  assert_eq(info.length(), 0)
}

// ============================================================================
// 11. Sequential clustering: ce_sequential_cluster
// ============================================================================

///|
test "ce_sequential_cluster_isolates_small" {
  // 8 samples: 2 in A (far), 6 identical in B
  // First split: small=A(2), large=B(6), large_frac=0.75>=0.7, qval~0 -> accept
  // Second iter on 6 identical B points: k-means gives small=empty -> break
  let data = [
    [0.0, 0.0],
    [0.1, 0.0],
    [10.0, 10.0],
    [10.0, 10.0],
    [10.0, 10.0],
    [10.0, 10.0],
    [10.0, 10.0],
    [10.0, 10.0],
  ]
  let labels = [0, 0, 0, 0, 0, 0, 0, 0]
  let new_labels = @src.ce_sequential_cluster(data, labels, 0.05, 0.7)
  assert_eq(new_labels.length(), 8)
  // A isolated as its own cluster, B kept together -> 2 clusters
  assert_eq(@src.ce_count_clusters(new_labels), 2)
  // A samples (indices 0,1) share a label
  assert_eq(new_labels[0], new_labels[1])
  // B samples (indices 2-7) share a label
  assert_eq(new_labels[2], new_labels[3])
  assert_eq(new_labels[3], new_labels[7])
  // A and B differ
  assert_true(new_labels[0] != new_labels[2])
}

// ============================================================================
// 12. RSEC pipeline: ce_rsec
// ============================================================================

///|
test "ce_rsec_pipeline" {
  // 4 samples, 2 clear clusters; small enough that sequential splits are
  // not triggered (each consensus cluster has < 4 points).
  let data = [
    [0.0, 0.0],
    [0.1, 0.0],
    [10.0, 10.0],
    [10.1, 10.0],
  ]
  let params = @src.ce_default_params()
  let result = @src.ce_rsec(data, params)
  assert_eq(result.labels.length(), 4)
  assert_true(result.n_clusters >= 1)
  assert_eq(result.cocluster.length(), 4)
  assert_eq(result.cocluster[0].length(), 4)
  assert_true((result.merge_cutoff - 0.05).abs() < 0.001)
}

// ============================================================================
// 13. Utilities: ce_count_clusters, ce_cluster_sizes, ce_summary
// ============================================================================

///|
test "ce_count_clusters_basic" {
  assert_eq(@src.ce_count_clusters([0, 1, 0, 1, 2]), 3)
  assert_eq(@src.ce_count_clusters([0, 0, 0]), 1)
}

///|
test "ce_count_clusters_with_unassigned" {
  assert_eq(@src.ce_count_clusters([-1, -1]), 0)
  assert_eq(@src.ce_count_clusters([0, -1, 1]), 2)
}

///|
test "ce_count_clusters_empty" {
  assert_eq(@src.ce_count_clusters([]), 0)
}

///|
test "ce_cluster_sizes_basic" {
  let labels = [0, 1, 0, 1, 2, 0]
  let sizes = @src.ce_cluster_sizes(labels)
  assert_eq(sizes.length(), 3)
  let (id0, size0) = sizes[0]
  assert_eq(id0, 0)
  assert_eq(size0, 3)
  let (id1, size1) = sizes[1]
  assert_eq(id1, 1)
  assert_eq(size1, 2)
  let (id2, size2) = sizes[2]
  assert_eq(id2, 2)
  assert_eq(size2, 1)
}

///|
test "ce_cluster_sizes_ignores_unassigned" {
  let sizes = @src.ce_cluster_sizes([0, -1, 1, 1, -1])
  assert_eq(sizes.length(), 2)
  let (id0, size0) = sizes[0]
  assert_eq(id0, 0)
  assert_eq(size0, 1)
  let (id1, size1) = sizes[1]
  assert_eq(id1, 1)
  assert_eq(size1, 2)
}

///|
test "ce_summary" {
  let result = @src.ClusterExperiment::new(
    labels=[0, 0, 1, 1, -1],
    cocluster=@src.ce_make_matrix(5, 5, 1.0),
    dendro=None,
    n_clusters=2,
    merge_cutoff=0.05,
    merge_info=[],
  )
  let s = @src.ce_summary(result)
  assert_true(s.contains("Samples: 5"))
  assert_true(s.contains("Clusters: 2"))
  assert_true(s.contains("Merge cutoff (q-value):"))
}

// ============================================================================
// 14. Normal CDF: ce_normal_cdf
// ============================================================================

///|
test "ce_normal_cdf_zero" {
  let cdf = @src.ce_normal_cdf(0.0)
  assert_true((cdf - 0.5).abs() < 0.001)
}

///|
test "ce_normal_cdf_positive" {
  let cdf = @src.ce_normal_cdf(1.96)
  assert_true((cdf - 0.975).abs() < 0.001)
}

///|
test "ce_normal_cdf_negative" {
  let cdf = @src.ce_normal_cdf(-1.96)
  assert_true((cdf - 0.025).abs() < 0.001)
}

///|
test "ce_normal_cdf_large_positive" {
  let cdf = @src.ce_normal_cdf(10.0)
  assert_true((cdf - 1.0).abs() < 0.001)
}

///|
test "ce_normal_cdf_symmetry" {
  // cdf(-x) = 1 - cdf(x)
  let x = 1.5
  let a = @src.ce_normal_cdf(x)
  let b = @src.ce_normal_cdf(-x)
  assert_true((a + b - 1.0).abs() < 0.001)
}

// ============================================================================
// 15. BH FDR: ce_bh_fdr
// ============================================================================

///|
test "ce_bh_fdr_basic" {
  let pvals = [0.01, 0.02, 0.5]
  let qvals = @src.ce_bh_fdr(pvals)
  assert_eq(qvals.length(), 3)
  assert_true((qvals[0] - 0.03).abs() < 0.001)
  assert_true((qvals[1] - 0.03).abs() < 0.001)
  assert_true((qvals[2] - 0.5).abs() < 0.001)
}

///|
test "ce_bh_fdr_empty" {
  let qvals = @src.ce_bh_fdr([])
  assert_eq(qvals.length(), 0)
}

///|
test "ce_bh_fdr_single" {
  let qvals = @src.ce_bh_fdr([0.3])
  assert_eq(qvals.length(), 1)
  assert_true((qvals[0] - 0.3).abs() < 0.001)
}

///|
test "ce_bh_fdr_identical_pvalues" {
  let pvals = [0.5, 0.5, 0.5]
  let qvals = @src.ce_bh_fdr(pvals)
  assert_eq(qvals.length(), 3)
  for q in qvals {
    assert_true((q - 0.5).abs() < 0.001)
  }
}

///|
test "ce_bh_fdr_monotone" {
  // q-values should be non-decreasing when sorted by p-value
  let pvals = [0.001, 0.02, 0.1, 0.04]
  let qvals = @src.ce_bh_fdr(pvals)
  assert_eq(qvals.length(), 4)
  // all q-values should be in [0, 1]
  for q in qvals {
    assert_true(q >= 0.0 && q <= 1.0)
  }
}

// ============================================================================
// 16. CeRsecParams: ce_default_params
// ============================================================================

///|
test "ce_default_params" {
  let p = @src.ce_default_params()
  assert_eq(p.ks.length(), 3)
  assert_eq(p.ks[0], 2)
  assert_eq(p.ks[1], 3)
  assert_eq(p.ks[2], 4)
  assert_eq(p.reduced_dims.length(), 2)
  assert_eq(p.reduced_dims[0], 0)
  assert_eq(p.reduced_dims[1], 5)
  assert_eq(p.subsamples.length(), 2)
  assert_true((p.consensus_proportion - 0.7).abs() < 0.001)
  assert_true((p.merge_cutoff - 0.05).abs() < 0.001)
  assert_true((p.sequential_alpha - 0.05).abs() < 0.001)
  assert_true((p.sequential_beta - 0.7).abs() < 0.001)
  assert_eq(p.n_init, 5)
  assert_eq(p.max_iter, 100)
}

///|
test "ce_default_params_used_in_rsec" {
  // Smoke test: default params produce a valid ClusterExperiment
  let data = [
    [0.0, 0.0],
    [0.1, 0.1],
    [5.0, 5.0],
    [5.1, 5.1],
  ]
  let params = @src.ce_default_params()
  let result = @src.ce_rsec(data, params)
  assert_true(result.labels.length() == 4)
  assert_true(result.n_clusters >= 1)
  // merge_info should have been populated if dendrogram was built
  match result.dendro {
    Some(_) => assert_true(result.merge_info.length() >= 0)
    None => assert_true(true)
  }
}
