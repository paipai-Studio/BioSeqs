///|
fn main {
  println("=== scran Demo ===")

  let (counts, gene_names) = @src.create_example_scran_data()

  println("\n1. Compute sum factors (scran normalization):")
  let sizes = [10]
  let sf = @src.compute_sum_factors(counts, sizes)
  println("  Number of cells: " + sf.factors.length().to_string())
  println("  First 5 size factors:")
  let mut i = 0
  while i < 5 && i < sf.factors.length() {
    println("    Cell " + i.to_string() + ": " + sf.factors[i].to_string())
    i = i + 1
  }

  println("\n2. Log normalize counts:")
  let normalized = @src.scran_log_normalize(counts, sf.factors)
  println(
    "  Normalized matrix shape: " +
    normalized.length().to_string() +
    " cells x " +
    normalized[0].length().to_string() +
    " genes",
  )

  println("\n3. Build SNN graph:")
  let (adjacency, weights) = @src.build_snn_graph(normalized, 10)
  println("  Graph nodes: " + adjacency.length().to_string())
  println(
    "  Average neighbors per cell: " +
    ((
      adjacency[0].length().to_double() +
      adjacency[1].length().to_double() +
      adjacency[2].length().to_double()
    ) /
    3.0).to_string(),
  )

  println("\n4. Leiden clustering:")
  let clusters = @src.leiden_clustering(adjacency, weights, 1.0)
  println("  Number of clusters: " + clusters.n_clusters.to_string())

  let cluster_counts : Map[Int, Int] = Map([], capacity=clusters.n_clusters)
  i = 0
  while i < clusters.clusters.length() {
    let c = clusters.clusters[i]
    let current = cluster_counts.get(c)
    if current.is_some() {
      cluster_counts[c] = current.unwrap() + 1
    } else {
      cluster_counts[c] = 1
    }
    i = i + 1
  }

  println("  Cluster sizes:")
  let keys = Array::make(clusters.n_clusters, 0)
  let mut k = 0
  while k < clusters.n_clusters {
    keys[k] = k
    k = k + 1
  }
  k = 0
  while k < keys.length() {
    let v = cluster_counts.get(keys[k])
    if v.is_some() {
      println(
        "    Cluster " +
        keys[k].to_string() +
        ": " +
        v.unwrap().to_string() +
        " cells",
      )
    }
    k = k + 1
  }

  println("\n5. Score marker genes:")
  let markers = @src.score_markers(counts, clusters.clusters, gene_names)
  println(
    "  Marker results for " +
    markers.gene.length().to_string() +
    " gene-cluster pairs",
  )

  println("\n6. Top markers for cluster 0 (by log2FC):")
  let top_markers : Array[(Double, String)] = Array::new()
  i = 0
  while i < markers.gene.length() {
    if markers.cluster[i] == 0 && markers.fdr[i] < 0.05 {
      top_markers.push((markers.avg_log2_fc[i], markers.gene[i]))
    }
    i = i + 1
  }

  let mut j = 0
  while j < top_markers.length() - 1 {
    let mut k = j + 1
    while k < top_markers.length() {
      if top_markers[j].0 < top_markers[k].0 {
        let tmp = top_markers[j]
        top_markers[j] = top_markers[k]
        top_markers[k] = tmp
      }
      k = k + 1
    }
    j = j + 1
  }

  j = 0
  while j < 5 && j < top_markers.length() {
    println(
      "    " + top_markers[j].1 + ": log2FC=" + top_markers[j].0.to_string(),
    )
    j = j + 1
  }

  println("\nscran demo completed!")
}
