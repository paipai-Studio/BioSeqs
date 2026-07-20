///|
/// Tests for Seurat module

///|
test "seurat_create_example_data" {
  let obj = @src.seurat_create_example_data()
  assert_eq(obj.row_names.length(), 500)
  assert_eq(obj.col_names.length(), 200)
  assert_eq(obj.counts.length(), 500)
  assert_eq(obj.counts[0].length(), 200)
}

///|
test "seurat_simple_pipeline" {
  let obj = @src.seurat_create_example_data()
  let normalized = @src.normalize_total(obj)
  let with_hvg = @src.find_variable_features(normalized)
  let scaled = @src.scale_data(with_hvg)
  let with_pca = @src.run_pca(scaled)
  let with_neighbors = @src.find_neighbors(with_pca)
  let with_clusters = @src.find_clusters(with_neighbors)
  let with_umap = @src.run_umap(with_pca)
  
  assert_eq(with_clusters.clusters.length(), 200)
  assert_eq(with_umap.umap.length(), 200)
}

///|
test "normalize_total" {
  let obj = @src.seurat_create_example_data()
  let normalized = @src.normalize_total(obj)
  assert_eq(normalized.data.length(), 500)
  assert_eq(normalized.data[0].length(), 200)
  
  let mut has_positive = false
  let mut i = 0
  while i < normalized.data.length() && !has_positive {
    let gene = normalized.data[i]
    let mut j = 0
    while j < gene.length() && !has_positive {
      if gene[j] > 0.0 {
        has_positive = true
      }
      j = j + 1
    }
    i = i + 1
  }
  assert_true(has_positive)
}

///|
test "find_variable_features" {
  let obj = @src.seurat_create_example_data()
  let normalized = @src.normalize_total(obj)
  let with_hvg = @src.find_variable_features(normalized)
  assert_true(with_hvg.var_features.length() > 0)
  assert_true(with_hvg.var_features.length() <= 500)
}

///|
test "scale_data" {
  let obj = @src.seurat_create_example_data()
  let normalized = @src.normalize_total(obj)
  let scaled = @src.scale_data(normalized)
  assert_eq(scaled.scaled_data.length(), 500)
  assert_eq(scaled.scaled_data[0].length(), 200)
}

///|
test "run_pca" {
  let obj = @src.seurat_create_example_data()
  let normalized = @src.normalize_total(obj)
  let with_hvg = @src.find_variable_features(normalized)
  let scaled = @src.scale_data(with_hvg)
  let with_pca = @src.run_pca(scaled)
  
  assert_true(with_pca.pca.length() > 0)
  assert_true(with_pca.pca[0].length() > 0)
  assert_eq(with_pca.pca.length(), 200)
  assert_true(with_pca.var_explained.length() > 0)
  
  let mut sum_var = 0.0
  let mut i = 0
  while i < with_pca.var_explained.length() {
    sum_var = sum_var + with_pca.var_explained[i]
    i = i + 1
  }
  assert_true(sum_var > 0.0)
}

///|
test "find_neighbors" {
  let obj = @src.seurat_create_example_data()
  let normalized = @src.normalize_total(obj)
  let with_hvg = @src.find_variable_features(normalized)
  let scaled = @src.scale_data(with_hvg)
  let with_pca = @src.run_pca(scaled)
  let with_neighbors = @src.find_neighbors(with_pca)
  
  assert_eq(with_neighbors.neighbors.length(), 200)
  assert_true(with_neighbors.neighbors[0].length() > 0)
}

///|
test "find_clusters" {
  let obj = @src.seurat_create_example_data()
  let normalized = @src.normalize_total(obj)
  let with_hvg = @src.find_variable_features(normalized)
  let scaled = @src.scale_data(with_hvg)
  let with_pca = @src.run_pca(scaled)
  let with_neighbors = @src.find_neighbors(with_pca)
  let with_clusters = @src.find_clusters(with_neighbors)
  
  assert_eq(with_clusters.clusters.length(), 200)
  
  let unique_clusters : Array[Int] = Array::new()
  let mut i = 0
  while i < with_clusters.clusters.length() {
    let c = with_clusters.clusters[i]
    let mut exists = false
    let mut j = 0
    while j < unique_clusters.length() {
      if unique_clusters[j] == c {
        exists = true
        break
      }
      j = j + 1
    }
    if !exists {
      unique_clusters.push(c)
    }
    i = i + 1
  }
  assert_true(unique_clusters.length() > 0)
}

///|
test "run_umap" {
  let obj = @src.seurat_create_example_data()
  let normalized = @src.normalize_total(obj)
  let with_hvg = @src.find_variable_features(normalized)
  let scaled = @src.scale_data(with_hvg)
  let with_pca = @src.run_pca(scaled)
  let with_umap = @src.run_umap(with_pca)
  
  assert_eq(with_umap.umap.length(), 200)
  assert_eq(with_umap.umap[0].length(), 2)
}

///|
test "find_all_markers" {
  let obj = @src.seurat_create_example_data()
  let normalized = @src.normalize_total(obj)
  let with_hvg = @src.find_variable_features(normalized)
  let scaled = @src.scale_data(with_hvg)
  let with_pca = @src.run_pca(scaled)
  let with_neighbors = @src.find_neighbors(with_pca)
  let with_clusters = @src.find_clusters(with_neighbors)
  
  let markers = @src.find_all_markers(with_clusters)
  
  assert_true(markers.length() >= 0)
  
  if markers.length() > 0 {
    let m = markers[0]
    assert_true(m.gene.length() > 0)
  }
}
