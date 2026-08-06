///|
fn bluster_demo_ints(values : Array[Int]) -> String {
  values.map(fn(value) { value.to_string() }).join(", ")
}

///|
fn bluster_demo_doubles(values : Array[Double]) -> String {
  values.map(fn(value) { value.to_string() }).join(", ")
}

///|
fn bluster_demo_optional(value : Double?) -> String {
  match value {
    Some(actual) => actual.to_string()
    None => "NA"
  }
}

///|
fn main {
  println("=== Bioconductor bluster 1.23.0 Demo ===")

  let coordinates = [
    [0.0, 0.0],
    [0.1, 0.0],
    [0.2, 0.1],
    [0.3, 0.1],
    [10.0, 10.0],
    [10.1, 10.0],
    [10.2, 10.1],
    [10.3, 10.1],
  ]

  let kmeans_config = @bio.BlusterKmeansConfig::create(2, starts=8, seed=101) catch {
    _ => abort("failed to create the k-means configuration")
  }
  let kmeans = @bio.bluster_cluster_rows_kmeans(coordinates, kmeans_config) catch {
    _ => abort("k-means clustering failed")
  }
  println("\n1. Deterministic multi-start K-means")
  println("  clusters: " + bluster_demo_ints(kmeans.clusters))
  println(
    "  total within sum of squares: " +
    kmeans.total_within_sum_squares.to_string(),
  )

  let graph = @bio.bluster_make_snn_graph(
    coordinates,
    k=2,
    weighting=@bio.bluster_rank_weight(),
  ) catch {
    _ => abort("SNN graph construction failed")
  }
  let graph_config = @bio.BlusterGraphConfig::create(
    k=2,
    resolution=0.5,
    seed=103,
  ) catch {
    _ => abort("failed to create the graph configuration")
  }
  let graph_clusters = @bio.bluster_cluster_rows_graph(
    coordinates, graph_config,
  ) catch {
    _ => abort("graph clustering failed")
  }
  println("\n2. Exact SNN graph and Louvain-style optimization")
  println("  graph edges: " + graph.edge_count().to_string())
  println("  clusters: " + bluster_demo_ints(graph_clusters.clusters))
  println("  modularity: " + graph_clusters.modularity.to_string())

  let two_step_config = @bio.BlusterTwoStepConfig::create(
    first_centers=4,
    first_starts=5,
    second_k=1,
    resolution=0.5,
    seed=107,
  ) catch {
    _ => abort("failed to create the two-step configuration")
  }
  let two_step = @bio.bluster_cluster_rows_two_step(
    coordinates, two_step_config,
  ) catch {
    _ => abort("two-step clustering failed")
  }
  println("\n3. Two-step vector quantization")
  println("  centroids: " + two_step.centroids.length().to_string())
  println("  clusters: " + bluster_demo_ints(two_step.clusters))

  let rand = @bio.bluster_pairwise_rand(
    kmeans.clusters,
    graph_clusters.clusters,
  ) catch {
    _ => abort("Rand comparison failed")
  }
  let silhouette = @bio.bluster_approx_silhouette(
    coordinates,
    graph_clusters.clusters,
  ) catch {
    _ => abort("silhouette approximation failed")
  }
  let purity = @bio.bluster_neighbor_purity(
    coordinates,
    graph_clusters.clusters,
    k=2,
  ) catch {
    _ => abort("neighbor purity failed")
  }
  let rmsd = @bio.bluster_cluster_rmsd(coordinates, graph_clusters.clusters) catch {
    _ => abort("cluster RMSD failed")
  }
  println("\n4. Cluster diagnostics")
  println("  adjusted Rand index: " + bluster_demo_optional(rand.index))
  println("  silhouette widths: " + bluster_demo_doubles(silhouette.widths))
  println("  neighbor purities: " + bluster_demo_doubles(purity.purity))
  println("  first cluster RMSD: " + bluster_demo_optional(rmsd.values[0]))

  let stability = @bio.bluster_bootstrap_kmeans_stability(
    coordinates,
    kmeans_config,
    iterations=8,
    seed=109,
  ) catch {
    _ => abort("bootstrap stability failed")
  }
  println("\n5. Bootstrap stability")
  println(
    "  first cluster coherence: " +
    bluster_demo_optional(stability.ratios[0][0]),
  )
  println(
    "  cluster separation: " + bluster_demo_optional(stability.ratios[0][1]),
  )

  let experiment = @bio.SingleCellExperiment::new(
    [
      [8.0, 7.0, 9.0, 8.0, 1.0, 2.0, 1.0, 2.0],
      [1.0, 2.0, 1.0, 2.0, 8.0, 7.0, 9.0, 8.0],
    ],
    ["MarkerA", "MarkerB"],
    ["C1", "C2", "C3", "C4", "C5", "C6", "C7", "C8"],
  )
  experiment.reduced_dims["PCA"] = coordinates
  experiment.metadata["source"] = "bluster_demo"
  let clustered = @bio.bluster_cluster_sce(
    experiment,
    graph_config,
    output_column="community",
  ) catch {
    _ => abort("SCE clustering failed")
  }
  println("\n6. Immutable SingleCellExperiment integration")
  println(
    "  written labels: " + clustered.experiment.col_data["community"].join(", "),
  )
  println(
    "  original unchanged: " +
    (!experiment.col_data.contains("community")).to_string(),
  )
  println("\nbluster demo completed.")
}
