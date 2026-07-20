///|
/// Seurat Demo - Comprehensive single-cell RNA-seq analysis.
/// Inspired by Bioconductor Seurat package.

///|
fn main {
  println("=== Seurat Single-cell RNA-seq Analysis Demo ===")
  println("Inspired by Bioconductor Seurat package")
  println("")
  
  println("1. Creating simulated single-cell data...")
  let obj = @bio.seurat_create_example_data()
  println("   Cells: \{obj.col_names.length().to_string()}")
  println("   Genes: \{obj.row_names.length().to_string()}")
  println("")
  
  println("2. Normalizing data (LogNormalize)...")
  let normalized = @bio.normalize_total(obj)
  println("   Normalization completed")
  println("")
  
  println("3. Finding highly variable genes...")
  let with_hvg = @bio.find_variable_features(normalized)
  println("   Found \{with_hvg.var_features.length().to_string()} highly variable genes")
  println("   First 5 HVGs:")
  let mut i = 0
  while i < 5 && i < with_hvg.var_features.length() {
    println("     \{with_hvg.var_features[i]}")
    i = i + 1
  }
  println("")
  
  println("4. Scaling data...")
  let scaled = @bio.scale_data(with_hvg)
  println("   Scaling completed")
  println("")
  
  println("5. Running PCA...")
  let with_pca = @bio.run_pca(scaled)
  println("   Number of PCs: \{with_pca.pca[0].length().to_string()}")
  println("   Variance explained by top 5 PCs:")
  let mut j = 0
  while j < 5 && j < with_pca.var_explained.length() {
    println("     PC\{(j + 1).to_string()}: \{(with_pca.var_explained[j] * 100.0).to_string()}%")
    j = j + 1
  }
  println("")
  
  println("6. Finding nearest neighbors...")
  let with_neighbors = @bio.find_neighbors(with_pca)
  println("   Neighbor graph built with \{with_neighbors.neighbors[0].length().to_string()} neighbors per cell")
  println("")
  
  println("7. Finding clusters...")
  let with_clusters = @bio.find_clusters(with_neighbors)
  
  let unique_clusters : Array[Int] = Array::new()
  for c in with_clusters.clusters {
    let mut exists = false
    for u in unique_clusters {
      if u == c {
        exists = true
        break
      }
    }
    if !exists {
      unique_clusters.push(c)
    }
  }
  
  println("   Number of clusters found: \{unique_clusters.length().to_string()}")
  
  let cluster_sizes : Array[Int] = Array::make(unique_clusters.length(), 0)
  for c in with_clusters.clusters {
    let mut idx = 0
    while idx < unique_clusters.length() {
      if unique_clusters[idx] == c {
        cluster_sizes[idx] = cluster_sizes[idx] + 1
        break
      }
      idx = idx + 1
    }
  }
  
  let mut k = 0
  while k < unique_clusters.length() {
    println("     Cluster \{unique_clusters[k].to_string()}: \{cluster_sizes[k].to_string()} cells")
    k = k + 1
  }
  println("")
  
  println("8. Running UMAP...")
  let with_umap = @bio.run_umap(with_pca)
  println("   UMAP completed with 2 dimensions")
  println("   First 3 UMAP coordinates:")
  let mut l = 0
  while l < 3 {
    println("     Cell \{l.to_string()}: (\{with_umap.umap[l][0].to_string()}, \{with_umap.umap[l][1].to_string()})")
    l = l + 1
  }
  println("")
  
  println("9. Finding marker genes...")
  let markers = @bio.find_all_markers(with_clusters)
  println("   Found \{markers.length().to_string()} marker genes")
  
  if markers.length() > 0 {
    println("   Top 5 markers by adjusted p-value:")
    let top_markers = @bio.get_top_markers(markers, 5)
    let mut o = 0
    while o < top_markers.length() {
      let marker = top_markers[o]
      println("     \{marker.gene}: log2FC=\{marker.avg_log2_fc.to_string()}, padj=\{marker.adjusted_p_value.to_string()}")
      o = o + 1
    }
  }
  println("")
  
  println("=== Demo completed successfully! ===")
  println("")
  println("Key takeaways:")
  println("- Seurat provides a complete single-cell analysis workflow")
  println("- Normalization -> HVG detection -> Scaling -> PCA -> Clustering -> UMAP -> Marker identification")
  println("- PCA reduces dimensionality while preserving biological variation")
  println("- Clustering identifies cell populations")
  println("- Marker genes define cell type identities")
}
