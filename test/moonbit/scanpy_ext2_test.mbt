///|
/// Tests for the scanpy extension batch #102:
/// neighbors metrics (#1), neighbors->obsp pipeline (#2), multilevel Louvain
/// (#3), rank_genes_groups method dispatch + logreg (#4), tl.kmeans (#5),
/// sctransform (#6), seurat_v3 VST HVG (#7), recipe_weinreb17 (#8), PHATE
/// (#9), TriMAP (#10), draw_graph layouts (#11), scIB suite (#12) and the
/// scirpy-style IR object model (#13).

///|
/// Euclidean distance between two vectors.
fn ext2_dist(a : Array[Double], b : Array[Double]) -> Double {
  let mut s = 0.0
  let d = if a.length() < b.length() { a.length() } else { b.length() }
  for t = 0; t < d; t = t + 1 {
    let diff = a[t] - b[t]
    s = s + diff * diff
  }
  s.sqrt()
}

///|
/// Edge weight of a SparseGraph (0.0 when the edge is absent).
fn ext2_edge(g : @src.SparseGraph, i : Int32, j : Int32) -> Double {
  if i >= g.n || j >= g.n {
    return 0.0
  }
  let mut w = 0.0
  for nb in g.adj[i] {
    if nb.0 == j {
      w = nb.1
    }
  }
  w
}

///|
test "scanpy_neighbors_metric_cosine" {
  // (1,0) is nearly parallel to (0.9,0.1); orthogonal and antiparallel pairs
  // are clamped to distance 1.0
  let x = [[1.0, 0.0], [0.9, 0.1], [0.0, 1.0], [-1.0, 0.0]]
  let res = @src.scanpy_neighbors_metric(x, 3, metric="cosine")
  // nearest neighbour of point 0 must be point 1
  assert_true(ext2_edge(res.distances, 0, 1) < 0.1)
  let w02 = ext2_edge(res.distances, 0, 2)
  let w03 = ext2_edge(res.distances, 0, 3)
  assert_true((w02 - 1.0).abs() < 1.0e-9)
  assert_true((w03 - 1.0).abs() < 1.0e-9)
}

///|
test "scanpy_neighbors_metric_manhattan" {
  let x = [[0.0, 0.0], [3.0, 0.0], [1.0, 1.0]]
  let res = @src.scanpy_neighbors_metric(x, 2, metric="manhattan")
  // manhattan distance (0,0)-(1,1) = 2 < (0,0)-(3,0) = 3
  assert_true((ext2_edge(res.distances, 0, 2) - 2.0).abs() < 1.0e-9)
  assert_true((ext2_edge(res.distances, 0, 1) - 3.0).abs() < 1.0e-9)
}

///|
test "scanpy_neighbors_metric_correlation" {
  // perfectly correlated vectors have correlation distance 0
  let x = [[1.0, 2.0, 3.0], [2.0, 4.0, 6.0], [1.0, 2.0, 2.0]]
  let res = @src.scanpy_neighbors_metric(x, 2, metric="correlation")
  assert_true(ext2_edge(res.distances, 0, 1) < 1.0e-9)
}

///|
test "scanpy_neighbors_metric_euclidean_matches_reference" {
  // the euclidean/umap path must be bit-identical to the reference pipeline
  let x = [
    [0.0, 0.0],
    [0.0, 1.0],
    [1.0, 0.0],
    [5.0, 5.0],
    [5.0, 6.0],
    [6.0, 5.0],
  ]
  let res = @src.scanpy_neighbors_metric(x, 3)
  let reference = @src.scanpy_neighbors_with_distances(x, 3)
  for i = 0; i < 6; i = i + 1 {
    for j = 0; j < 6; j = j + 1 {
      let a = ext2_edge(res.connectivities, i, j)
      let b = ext2_edge(reference.connectivities, i, j)
      assert_true((a - b).abs() < 1.0e-12)
    }
  }
}

///|
test "scanpy_neighbors_metric_gauss_matches_reference" {
  let x = [[0.0, 0.0], [0.0, 1.0], [1.0, 0.0], [5.0, 5.0]]
  let res = @src.scanpy_neighbors_metric(x, 2, method_="gauss")
  let gauss_ref = @src.scanpy_neighbors_gauss(x, 2)
  for i = 0; i < 4; i = i + 1 {
    for j = 0; j < 4; j = j + 1 {
      let a = ext2_edge(res.connectivities, i, j)
      let b = ext2_edge(gauss_ref, i, j)
      assert_true((a - b).abs() < 1.0e-12)
    }
  }
}

///|
test "scanpy_neighbors_anndata_obsp_pipeline" {
  let x = [
    [0.0, 0.0],
    [0.0, 1.0],
    [1.0, 0.0],
    [5.0, 5.0],
    [5.0, 6.0],
    [6.0, 5.0],
  ]
  let obs_names = ["c0", "c1", "c2", "c3", "c4", "c5"]
  let adata = @src.AnnData::new(x, obs_names, ["g0", "g1"])
  let res = @src.scanpy_neighbors_anndata(adata, 3)
  assert_eq(res.n_obs, 6)
  assert_true(adata.obsp.contains("distances"))
  assert_true(adata.obsp.contains("connectivities"))
  let conn = adata.obsp["connectivities"]
  assert_eq(conn.length(), 6)
  // diagonal self-weight 1 like scanpy
  assert_true((conn[0][0] - 1.0).abs() < 1.0e-12)
  // graph round-trip from obsp
  let graph = @src.scanpy_graph_from_obsp(adata, key="connectivities")
  assert_eq(graph.num_vertices(), 6)
  // two separated blobs -> two leiden communities written to obs
  let n_comm = @src.scanpy_leiden_anndata(adata)
  assert_eq(n_comm, 2)
  let labels = adata.obs_column("leiden").unwrap()
  assert_eq(labels.length(), 6)
  assert_eq(labels[0], labels[1])
  assert_eq(labels[0], labels[2])
  assert_eq(labels[3], labels[4])
  assert_true(labels[0] != labels[3])
  // louvain wrapper as well
  let n_comm2 = @src.scanpy_louvain_anndata(adata, obs_key="louvain")
  assert_eq(n_comm2, 2)
}

///|
test "scanpy_louvain_multilevel_two_blobs" {
  // two 4-cliques joined by a single weak bridge
  let n = 8
  let adj : Array[Array[(Int32, Double)]] = Array(capacity=n)
  let diag = Array::make(n, 0.0)
  let pairs = [
    (0, 1),
    (0, 2),
    (0, 3),
    (1, 2),
    (1, 3),
    (2, 3),
    (4, 5),
    (4, 6),
    (4, 7),
    (5, 6),
    (5, 7),
    (6, 7),
    (3, 4),
  ]
  for i = 0; i < n; i = i + 1 {
    adj.push([])
  }
  for p in pairs {
    adj[p.0].push((p.1, 1.0))
    adj[p.1].push((p.0, 1.0))
  }
  let graph = @src.SparseGraph::new(n, diag, adj)
  let res = @src.scanpy_louvain_multilevel(graph)
  assert_eq(res.n_communities, 2)
  assert_eq(res.partition.length(), 8)
  assert_eq(res.partition[0], res.partition[1])
  assert_eq(res.partition[0], res.partition[3])
  assert_eq(res.partition[4], res.partition[7])
  assert_true(res.partition[0] != res.partition[4])
  assert_true(res.n_levels >= 1)
  assert_true(res.modularity > 0.3)
  // modularity of the found partition is consistent with the standalone fn
  let q = @src.scanpy_modularity(graph, res.partition)
  assert_true((q - res.modularity).abs() < 1.0e-9)
}

///|
test "scanpy_rank_genes_groups_method_dispatch" {
  // gene0 high in A, gene1 high in B, gene2 flat
  let x = [
    [10.0, 1.0, 5.0],
    [10.0, 1.0, 5.0],
    [10.0, 1.0, 5.0],
    [1.0, 10.0, 5.0],
    [1.0, 10.0, 5.0],
    [1.0, 10.0, 5.0],
  ]
  let groups = ["A", "A", "A", "B", "B", "B"]
  let t_res = @src.scanpy_rank_genes_groups_method(x, groups, "B", "A")
  assert_eq(t_res.method_, "t-test")
  assert_eq(t_res.gene_indices[0], 0)
  assert_eq(t_res.gene_indices.length(), 3)
  assert_true(t_res.logfoldchanges[0] > 0.0)
  let w_res = @src.scanpy_rank_genes_groups_method(
    x,
    groups,
    "B",
    "A",
    method_="wilcoxon",
  )
  assert_eq(w_res.method_, "wilcoxon")
  assert_eq(w_res.gene_indices[0], 0)
  assert_true(w_res.pvals_adj[0] >= w_res.pvals[0] - 1.0e-12)
  let l_res = @src.scanpy_rank_genes_groups_method(
    x,
    groups,
    "B",
    "A",
    method_="logreg",
  )
  assert_eq(l_res.method_, "logreg")
  assert_eq(l_res.gene_indices[0], 0)
  // logreg reports no p-values
  assert_true((l_res.pvals[0] - 1.0).abs() < 1.0e-12)
  // n_top truncation
  let top1 = @src.scanpy_rank_genes_groups_method(
    x,
    groups,
    "B",
    "A",
    method_="t-test",
    n_top=1,
  )
  assert_eq(top1.gene_indices.length(), 1)
}

///|
test "scanpy_kmeans_two_blobs" {
  let x = [
    [0.0, 0.0],
    [0.0, 1.0],
    [1.0, 0.0],
    [1.0, 1.0],
    [10.0, 10.0],
    [10.0, 11.0],
    [11.0, 10.0],
    [11.0, 11.0],
  ]
  let res = @src.scanpy_kmeans(x, 2, n_init=3)
  assert_eq(res.n_clusters, 2)
  assert_eq(res.labels.length(), 8)
  assert_eq(res.labels[0], res.labels[1])
  assert_eq(res.labels[0], res.labels[3])
  assert_eq(res.labels[4], res.labels[7])
  assert_true(res.labels[0] != res.labels[4])
  assert_true(res.inertia < 4.1)
  assert_eq(res.centers.length(), 2)
  // AnnData wrapper writes obs + uns
  let adata = @src.AnnData::new(x, ["a", "b", "c", "d", "e", "f", "g", "h"], [
    "g0", "g1",
  ])
  let km = @src.scanpy_kmeans_anndata(adata, 2, n_init=3)
  assert_eq(km.labels.length(), 8)
  let labels = adata.obs_column("kmeans").unwrap()
  assert_eq(labels.length(), 8)
  match adata.get_uns("kmeans") {
    Some(@src.JsonObj(m)) => assert_true(m.contains("method"))
    Some(_) => fail("uns[kmeans] should be a JsonObj")
    None => fail("uns[kmeans] missing")
  }
}

///|
test "scanpy_sctransform_residuals" {
  // gene0 high, gene1 depth-dependent, gene2 sparse (filtered), gene3 medium
  let x = [
    [10.0, 5.0, 0.0, 3.0],
    [12.0, 10.0, 0.0, 4.0],
    [11.0, 15.0, 1.0, 3.0],
    [13.0, 20.0, 0.0, 5.0],
    [10.0, 25.0, 0.0, 4.0],
    [12.0, 30.0, 1.0, 3.0],
    [11.0, 35.0, 0.0, 5.0],
    [9.0, 40.0, 0.0, 4.0],
  ]
  let res = @src.scanpy_sctransform(x, min_cells=3)
  assert_eq(res.n_cells, 8)
  // gene2 expressed in only 2 cells -> filtered out
  assert_eq(res.kept_genes.length(), 3)
  let mut has_gene2 = false
  for g in res.kept_genes {
    if g == 2 {
      has_gene2 = true
    }
  }
  assert_true(!has_gene2)
  assert_eq(res.residuals.length(), 8)
  assert_eq(res.residuals_clipped[0].length(), 3)
  // all parameters positive; clipped residuals within the clip range
  for t = 0; t < res.theta.length(); t = t + 1 {
    assert_true(res.theta[t] > 0.0)
  }
  for row in res.residuals_clipped {
    for r in row {
      assert_true(r <= 10.0 && r >= -10.0)
    }
  }
  // regularized parameters shrink toward the pooled estimates
  assert_true(res.pooled_beta != 0.0)
}

///|
test "scanpy_hvg_seurat_v3_vst" {
  // gene0: bimodal highly variable; gene1: constant; gene2: mild;
  // gene3: filler balancing every cell to the same count depth
  let x = [
    [0.0, 5.0, 3.0, 52.0],
    [0.0, 5.0, 4.0, 51.0],
    [0.0, 5.0, 3.0, 52.0],
    [30.0, 5.0, 4.0, 21.0],
    [30.0, 5.0, 3.0, 22.0],
    [30.0, 5.0, 4.0, 21.0],
  ]
  let res = @src.scanpy_hvg_seurat_v3(x, n_top_genes=2, min_cells=2)
  assert_eq(res.kept_genes.length(), 4)
  // the two non-constant genes form the top-2 (index-order tie)
  assert_true(res.highly_variable[0])
  assert_true(!res.highly_variable[1])
  assert_true(res.highly_variable[2])
  // constant genes are excluded from the loess fit (stabilized variance 0)
  assert_true((res.residual_variance[1] - 0.0).abs() < 1.0e-12)
}

///|
test "scanpy_recipe_weinreb17_pipeline" {
  let x = [
    [10.0, 0.0, 5.0, 1.0],
    [12.0, 1.0, 5.0, 0.0],
    [11.0, 0.0, 6.0, 1.0],
    [30.0, 2.0, 5.0, 0.0],
    [28.0, 1.0, 4.0, 1.0],
    [31.0, 2.0, 5.0, 0.0],
  ]
  let res = @src.scanpy_recipe_weinreb17(
    x,
    min_cells=2,
    n_top_genes=2,
    n_comps=2,
    n_neighbors=3,
  )
  assert_eq(res.scaled.length(), 6)
  assert_eq(res.scaled[0].length(), 2)
  assert_eq(res.hvg_indices.length(), 2)
  assert_eq(res.pca.coords.length(), 6)
  assert_eq(res.pca.coords[0].length(), 2)
  assert_eq(res.diffusion.eigen_basis.length(), 6)
  assert_true(res.counts_median > 0.0)
}

///|
test "scanpy_phate_two_blobs" {
  let x = [
    [0.0, 0.0],
    [0.0, 1.0],
    [1.0, 0.0],
    [1.0, 1.0],
    [10.0, 10.0],
    [10.0, 11.0],
    [11.0, 10.0],
    [11.0, 11.0],
  ]
  let res = @src.scanpy_phate(x, n_neighbors=3, mds_iter=100)
  assert_eq(res.embedding.length(), 8)
  assert_eq(res.embedding[0].length(), 2)
  assert_true(res.t >= 1)
  // finite coordinates
  for row in res.embedding {
    for v in row {
      assert_true(v == v && v < 1.0e6 && v > -1.0e6)
    }
  }
  // the two blobs are separated in the embedding: cross-cluster distance
  // exceeds within-cluster distance
  let mut cross_sum = 0.0
  let mut within_sum = 0.0
  for i = 0; i < 4; i = i + 1 {
    for j = 4; j < 8; j = j + 1 {
      cross_sum = cross_sum + ext2_dist(res.embedding[i], res.embedding[j])
    }
    for j = 0; j < 4; j = j + 1 {
      if i != j {
        within_sum = within_sum + ext2_dist(res.embedding[i], res.embedding[j])
      }
    }
  }
  let cross_mean = cross_sum / 16.0
  let within_mean = within_sum / 12.0
  assert_true(cross_mean > within_mean)
}

///|
test "scanpy_trimap_two_blobs" {
  let x = [
    [0.0, 0.0],
    [0.0, 1.0],
    [1.0, 0.0],
    [10.0, 10.0],
    [10.0, 11.0],
    [11.0, 10.0],
  ]
  let res = @src.scanpy_trimap(x, n_iter=400, seed=42)
  assert_eq(res.embedding.length(), 6)
  assert_eq(res.n_inliers, 10)
  assert_eq(res.n_iter, 400)
  for row in res.embedding {
    for v in row {
      assert_true(v == v && v < 1.0e6 && v > -1.0e6)
    }
  }
  // blobs separated in the embedding
  let mut cross_sum = 0.0
  let mut within_sum = 0.0
  for i = 0; i < 3; i = i + 1 {
    for j = 3; j < 6; j = j + 1 {
      cross_sum = cross_sum + ext2_dist(res.embedding[i], res.embedding[j])
    }
    for j = 0; j < 3; j = j + 1 {
      if i != j {
        within_sum = within_sum + ext2_dist(res.embedding[i], res.embedding[j])
      }
    }
  }
  assert_true(cross_sum / 9.0 > within_sum / 6.0)
}

///|
test "scanpy_draw_graph_sparse_layouts" {
  // two triangles joined by a bridge (2-3)
  let n = 6
  let adj : Array[Array[(Int32, Double)]] = Array(capacity=n)
  let diag = Array::make(n, 0.0)
  let pairs = [(0, 1), (0, 2), (1, 2), (3, 4), (3, 5), (4, 5), (2, 3)]
  for i = 0; i < n; i = i + 1 {
    adj.push([])
  }
  for p in pairs {
    adj[p.0].push((p.1, 1.0))
    adj[p.1].push((p.0, 1.0))
  }
  let graph = @src.SparseGraph::new(n, diag, adj)
  for layout in ["fr", "kk", "fa2"] {
    let pos = @src.scanpy_draw_graph_sparse(graph, layout~, n_iter=300)
    assert_eq(pos.length(), 6)
    for row in pos {
      assert_eq(row.length(), 2)
      for v in row {
        assert_true(v == v)
      }
    }
    // same-cluster points are closer than cross-cluster points
    let mut within = 0.0
    let mut cross = 0.0
    for i = 0; i < 3; i = i + 1 {
      for j = 3; j < 6; j = j + 1 {
        cross = cross + ext2_dist(pos[i], pos[j])
      }
      for j = 0; j < 3; j = j + 1 {
        if i != j {
          within = within + ext2_dist(pos[i], pos[j])
        }
      }
    }
    assert_true(within / 6.0 < cross / 9.0)
  }
}

///|
test "scanpy_scib_metrics_and_correction" {
  // mixed batches along a line: alternating batch ids keep every kNN
  // neighbourhood balanced, so kBET sits near its ideal value
  let mixed = [
    [0.0, 0.0],
    [0.1, 0.0],
    [0.2, 0.0],
    [0.3, 0.0],
    [0.4, 0.0],
    [0.5, 0.0],
  ]
  let mixed_batch = [0, 1, 0, 1, 0, 1]
  let m = @src.scanpy_scib_metrics(mixed, mixed_batch, [1, 1, 0, 0, 1, 1], k=5)
  assert_true((m.graph_connectivity - 1.0).abs() < 1.0e-9)
  assert_true(m.kbet > 0.8)
  // separated batches: neighbourhoods are batch-pure, so kBET collapses
  let split = [
    [0.0, 0.0],
    [0.1, 0.0],
    [0.2, 0.0],
    [5.0, 0.0],
    [5.1, 0.0],
    [5.2, 0.0],
  ]
  let split_batch = [0, 0, 0, 1, 1, 1]
  let s = @src.scanpy_scib_metrics(split, split_batch, [0, 0, 0, 1, 1, 1], k=3)
  assert_true(s.kbet < m.kbet)
  // label silhouette high for the label-aligned split embedding
  assert_true(s.silhouette_label > m.silhouette_label)
  // correction dispatch keeps the shape
  let none = @src.scanpy_scib_correct(split, split_batch, method_="none")
  assert_eq(none.method_, "none")
  assert_eq(none.corrected.length(), 6)
  for i = 0; i < 6; i = i + 1 {
    for d = 0; d < 2; d = d + 1 {
      assert_true((none.corrected[i][d] - split[i][d]).abs() < 1.0e-12)
    }
  }
  let har = @src.scanpy_scib_correct(split, split_batch, method_="harmony")
  assert_eq(har.method_, "harmony")
  assert_eq(har.n_batches, 2)
  assert_eq(har.corrected.length(), 6)
  assert_eq(har.corrected[0].length(), 2)
  let sc = @src.scanpy_scib_correct(split, split_batch, method_="scanorama")
  assert_eq(sc.method_, "scanorama")
  assert_eq(sc.corrected.length(), 6)
  // scanorama centers each batch: batch means collapse toward zero
  let mut mean0 = 0.0
  let mut mean1 = 0.0
  for d = 0; d < 2; d = d + 1 {
    mean0 = mean0 +
      (sc.corrected[0][d] + sc.corrected[1][d] + sc.corrected[2][d]) / 3.0
    mean1 = mean1 +
      (sc.corrected[3][d] + sc.corrected[4][d] + sc.corrected[5][d]) / 3.0
  }
  assert_true(mean0.abs() < 0.5 && mean1.abs() < 0.5)
  // silhouette helper sanity: perfectly separated groups score 1
  let sil = @src.scanpy_scib_silhouette(split, split_batch)
  assert_true(sil > 0.9)
}

///|
test "scanpy_ir_pair_chains_and_network" {
  let chains : Array[@src.ScanpyIrChain] = [
    {
      cell_id: "c1",
      locus: "TRA",
      v_call: "TRAV1",
      j_call: "TRAJ1",
      cdr3_aa: "CAAS",
      productive: true,
    },
    {
      cell_id: "c1",
      locus: "TRB",
      v_call: "TRBV1",
      j_call: "TRBJ1",
      cdr3_aa: "CASS",
      productive: true,
    },
    {
      cell_id: "c2",
      locus: "TRA",
      v_call: "TRAV1",
      j_call: "TRAJ1",
      cdr3_aa: "CAAS",
      productive: true,
    },
    {
      cell_id: "c2",
      locus: "TRB",
      v_call: "TRBV1",
      j_call: "TRBJ1",
      cdr3_aa: "CASS",
      productive: true,
    },
    {
      cell_id: "c3",
      locus: "TRA",
      v_call: "TRAV2",
      j_call: "TRAJ1",
      cdr3_aa: "CAAT",
      productive: true,
    },
  ]
  let data = @src.scanpy_ir_pair_chains(chains)
  assert_eq(data.cells.length(), 3)
  assert_eq(data.clonotypes.length(), 2)
  // c1 and c2 share both primary keys -> same clonotype
  assert_eq(data.cells[0].clonotype_id, data.cells[1].clonotype_id)
  assert_true(data.cells[2].clonotype_id != data.cells[0].clonotype_id)
  assert_eq(data.cells[0].primary_vj.unwrap().cdr3_aa, "CAAS")
  assert_eq(data.cells[0].primary_vdj.unwrap().locus, "TRB")
  match data.cells[0].secondary_vj {
    Some(_) => fail("c1 has a single TRA chain")
    None => ()
  }
  // exact network: the two clonotypes share no CDR3 -> no edge
  let net = @src.scanpy_ir_clonotype_network(data)
  assert_eq(net.num_vertices(), 2)
  assert_true(ext2_edge(net, 0, 1) == 0.0)
  // hamming network: CASS vs CAAT differ at 2 positions
  let net1 = @src.scanpy_ir_clonotype_network(
    data,
    mode="hamming",
    max_distance=1,
  )
  assert_true(ext2_edge(net1, 0, 1) == 0.0)
  let net2 = @src.scanpy_ir_clonotype_network(
    data,
    mode="hamming",
    max_distance=2,
  )
  assert_true(ext2_edge(net2, 0, 1) == 1.0)
  // abundance table
  let sizes = @src.scanpy_ir_clonotype_sizes(data)
  assert_eq(sizes.length(), 2)
  assert_true(
    (sizes[0].1 == 2 && sizes[1].1 == 1) || (sizes[0].1 == 1 && sizes[1].1 == 2),
  )
}

///|
test "scanpy_ir_multichain_and_locus_routing" {
  let chains : Array[@src.ScanpyIrChain] = [
    {
      cell_id: "c1",
      locus: "TRA",
      v_call: "TRAV1",
      j_call: "TRAJ1",
      cdr3_aa: "CAAA",
      productive: true,
    },
    {
      cell_id: "c1",
      locus: "TRG",
      v_call: "TRGV1",
      j_call: "TRGJ1",
      cdr3_aa: "CAGG",
      productive: true,
    },
    {
      cell_id: "c1",
      locus: "TRA",
      v_call: "TRAV9",
      j_call: "TRAJ9",
      cdr3_aa: "CATT",
      productive: false,
    },
    {
      cell_id: "c1",
      locus: "TRA",
      v_call: "TRAV8",
      j_call: "TRAJ8",
      cdr3_aa: "CACG",
      productive: true,
    },
    {
      cell_id: "c2",
      locus: "IGH",
      v_call: "IGHV1",
      j_call: "IGHJ1",
      cdr3_aa: "CAR",
      productive: true,
    },
  ]
  assert_true(@src.scanpy_ir_is_vj("TRG"))
  assert_true(!@src.scanpy_ir_is_vj("IGH"))
  let data = @src.scanpy_ir_pair_chains(chains)
  assert_eq(data.cells.length(), 2)
  // c1 has 3 TRA chains -> multichain, excluded from clonotypes
  assert_true(data.cells[0].multichain)
  assert_eq(data.cells[0].clonotype_id, -1)
  // primary VJ chain is the first productive one
  assert_eq(data.cells[0].primary_vj.unwrap().v_call, "TRAV1")
  assert_eq(data.cells[0].secondary_vj.unwrap().v_call, "TRAV8")
  // c2 has a single IGH (VDJ) chain and gets its own clonotype
  assert_true(!data.cells[1].multichain)
  assert_eq(data.cells[1].clonotype_id, 0)
  assert_eq(data.clonotypes.length(), 1)
  assert_eq(data.clonotypes[0].vdj_key, "IGHV1:IGHJ1:CAR")
  assert_eq(data.clonotypes[0].vj_key, "")
}
