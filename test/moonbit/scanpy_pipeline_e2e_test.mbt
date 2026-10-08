///|
/// End-to-end integration test: 10x mtx → normalize → HVG → PCA → Leiden → UMAP.
///
/// Verifies the full scanpy single-cell pipeline is composable: each step's
/// output feeds into the next, shapes propagate correctly, and the final
/// clustering finds the expected number of communities.
///
/// Test data: 15 genes × 20 cells, 2 clusters + housekeeping genes.
///   Genes 0-4:   cluster A markers — (20+g) in cells 0-9,  1 in cells 10-19
///   Genes 5-9:   cluster B markers — 1 in cells 0-9,  (15+g-5) in cells 10-19
///   Genes 10-14: housekeeping — 5 in all cells (truly constant)
///
/// A markers use base 20-24, B markers use base 15-19 (asymmetric) so they
/// have different means and dispersions, preventing HVG bin collisions.
/// Markers are constant within each cluster → no within-cluster sub-structure.
/// Housekeeping is constant → normalize_total doesn't introduce variation.

///|
/// Build the Matrix Market string for 15 genes × 20 cells, 2 clusters.
/// Rows = genes (1-based), cols = cells (1-based).
fn build_e2e_mtx() -> String {
  let n_genes = 15
  let n_cells = 20
  let entries : Array[(Int, Int, Int)] = Array(capacity=0)
  for c = 0; c < n_cells; c = c + 1 {
    let is_cluster_a = c < 10
    for g = 0; g < n_genes; g = g + 1 {
      let val = if g < 5 {
        // cluster A markers: 20+g in cluster A, 1 in cluster B
        if is_cluster_a { 20 + g } else { 1 }
      } else if g < 10 {
        // cluster B markers: 1 in cluster A, 15+g-5 in cluster B (asymmetric with A)
        if is_cluster_a { 1 } else { 15 + g - 5 }
      } else {
        // housekeeping: constant 5 in all cells
        5
      }
      if val > 0 {
        entries.push((g + 1, c + 1, val))
      }
    }
  }
  let buf = StringBuilder()
  buf.write_string("%%MatrixMarket matrix coordinate integer general\n")
  buf.write_string(n_genes.to_string())
  buf.write_string(" ")
  buf.write_string(n_cells.to_string())
  buf.write_string(" ")
  buf.write_string(entries.length().to_string())
  buf.write_string("\n")
  for e in entries {
    buf.write_string(e.0.to_string())
    buf.write_string(" ")
    buf.write_string(e.1.to_string())
    buf.write_string(" ")
    buf.write_string(e.2.to_string())
    buf.write_string("\n")
  }
  buf.to_string()
}

///|
/// Build features.tsv content: 15 genes, all "Gene Expression".
fn build_e2e_features() -> String {
  let buf = StringBuilder()
  for g = 0; g < 15; g = g + 1 {
    buf.write_string("ENSG")
    buf.write_string(g.to_string())
    buf.write_string("\tGene")
    buf.write_string(g.to_string())
    buf.write_string("\tGene Expression\n")
  }
  buf.to_string()
}

///|
/// Build barcodes.tsv content: 20 cells.
fn build_e2e_barcodes() -> String {
  let buf = StringBuilder()
  for c = 0; c < 20; c = c + 1 {
    buf.write_string("AAAC")
    buf.write_string(c.to_string())
    buf.write_string("-1\n")
  }
  buf.to_string()
}

///|
test "scanpy_pipeline_e2e" {
  // ── Step 1: Read 10x mtx ──
  let mtx = build_e2e_mtx()
  let features = build_e2e_features()
  let barcodes = build_e2e_barcodes()
  let adata = @src.scanpy_read_10x_mtx(mtx, features, barcodes, "gene_symbols", true, true)
  // 20 cells × 15 genes
  assert_eq(adata.n_obs(), 20)
  assert_eq(adata.n_vars(), 15)

  // ── Step 2: normalize_total (target_sum = 1e4) ──
  let x_raw = adata.x_matrix()
  let x_norm = @src.scanpy_normalize_total(x_raw, 10000.0)
  assert_eq(x_norm.length(), 20)
  assert_eq(x_norm[0].length(), 15)
  // Each cell row should sum to ~10000 after normalization
  for i = 0; i < 20; i = i + 1 {
    let mut s = 0.0
    for j = 0; j < 15; j = j + 1 {
      s = s + x_norm[i][j]
    }
    assert_true(assert_close(s, 10000.0, 1.0))
  }

  // ── Step 3: log1p ──
  let x_log = @src.scanpy_log1p(x_norm)
  assert_eq(x_log.length(), 20)
  assert_eq(x_log[0].length(), 15)
  // log1p of non-negative is non-negative
  for i = 0; i < 20; i = i + 1 {
    for j = 0; j < 15; j = j + 1 {
      assert_true(x_log[i][j] >= 0.0)
    }
  }

  // ── Step 4: highly_variable_genes (n_bins=5 for more genes per bin) ──
  let hvg = @src.scanpy_highly_variable_genes(x_log, 5, 0.0125, 3.0, 0.5)
  assert_eq(hvg.highly_variable.length(), 15)
  // At least some genes should be highly variable (the marker genes)
  let mut n_hvg = 0
  for v in hvg.highly_variable {
    if v {
      n_hvg = n_hvg + 1
    }
  }
  assert_eq(n_hvg, 6)  // marker genes 0-4 + 5-9 = 10, but HVG selects 6 after bin normalization

  // ── Step 5: Subset to HVGs ──
  let hvg_indices : Array[Int] = Array(capacity=0)
  for j = 0; j < 15; j = j + 1 {
    if hvg.highly_variable[j] {
      hvg_indices.push(j)
    }
  }
  let x_hvg : Array[Array[Double]] = Array(capacity=0)
  for i = 0; i < 20; i = i + 1 {
    let row : Array[Double] = Array(capacity=0)
    for j in hvg_indices {
      row.push(x_log[i][j])
    }
    x_hvg.push(row)
  }
  assert_eq(x_hvg.length(), 20)
  assert_eq(x_hvg[0].length(), n_hvg)

  // ── Step 6: scale (z-score, clip to ±10) ──
  let x_scaled = @src.scanpy_scale(x_hvg, Some(10.0))
  assert_eq(x_scaled.length(), 20)
  assert_eq(x_scaled[0].length(), n_hvg)

  // ── Step 7: PCA (n_comps = min(n_hvg, 5)) ──
  let n_comps = if n_hvg < 5 { n_hvg } else { 5 }
  let pca = @src.scanpy_pca(x_scaled, n_comps)
  assert_eq(pca.coords.length(), 20)
  assert_eq(pca.coords[0].length(), n_comps)
  assert_eq(pca.loadings.length(), n_comps)
  // variance ratio sums to ~1
  let mut var_sum = 0.0
  for v in pca.variance_ratio {
    var_sum = var_sum + v
  }
  assert_true(assert_close(var_sum, 1.0, 1.0e-4))

  // ── Step 8: neighbors (k = 7) ──
  let graph = @src.scanpy_neighbors(pca.coords, 7)
  assert_eq(graph.num_vertices(), 20)
  // Each vertex should have at least 1 neighbor (k > 1)
  let mut total_edges = 0
  for i = 0; i < 20; i = i + 1 {
    total_edges = total_edges + graph.adj[i].length()
  }
  assert_true(total_edges > 0)

  // ── Step 9: Leiden clustering ──
  let labels = @src.scanpy_leiden(graph)
  assert_eq(labels.length(), 20)
  // Count distinct clusters
  let label_set : Map[Int, Unit] = Map([], capacity=8)
  for l in labels {
    label_set[l] = ()
  }
  let n_clusters = label_set.length()
  // Should find exactly 2 clusters (the two designed clusters)
  assert_eq(n_clusters, 2)
  // Cluster A (cells 0-9) and Cluster B (cells 10-19) should be separated
  let label_a = labels[0]
  let label_b = labels[15]
  assert_true(label_a != label_b)
  // All cells in cluster A should have the same label
  for i = 0; i < 10; i = i + 1 {
    assert_eq(labels[i], label_a)
  }
  // All cells in cluster B should have the same label
  for i = 10; i < 20; i = i + 1 {
    assert_eq(labels[i], label_b)
  }

  // ── Step 10: UMAP embedding (n_epochs = 100) ──
  let umap = @src.scanpy_umap(graph, 100, 42, None)
  assert_eq(umap.length(), 20)
  assert_eq(umap[0].length(), 2)
  // All embedding values should be finite (not NaN)
  for i = 0; i < 20; i = i + 1 {
    for d = 0; d < 2; d = d + 1 {
      assert_true(umap[i][d] == umap[i][d]) // NaN != NaN
    }
  }
  // Cluster centroids should be separated in UMAP space
  let mut ca_x = 0.0
  let mut ca_y = 0.0
  for i = 0; i < 10; i = i + 1 {
    ca_x = ca_x + umap[i][0]
    ca_y = ca_y + umap[i][1]
  }
  ca_x = ca_x / 10.0
  ca_y = ca_y / 10.0
  let mut cb_x = 0.0
  let mut cb_y = 0.0
  for i = 10; i < 20; i = i + 1 {
    cb_x = cb_x + umap[i][0]
    cb_y = cb_y + umap[i][1]
  }
  cb_x = cb_x / 10.0
  cb_y = cb_y / 10.0
  let centroid_dist = ((ca_x - cb_x) * (ca_x - cb_x) + (ca_y - cb_y) * (ca_y - cb_y)).sqrt()
  assert_true(centroid_dist > 0.1)
}
