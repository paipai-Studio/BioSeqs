// Bioconductor Banksy-inspired spatial clustering workflow.

///|
fn banksy_demo_config() -> @src.BanksyComputeConfig {
  @src.BanksyComputeConfig::create(
    max_harmonic=1,
    k_geom=[4, 4],
    spatial_mode=@src.banksy_knn_median(),
    seed=7,
  ) catch {
    BanksyError(message) => abort("invalid BANKSY configuration: " + message)
  }
}

///|
fn banksy_demo_cluster_config() -> @src.BanksyClusterConfig {
  @src.BanksyClusterConfig::create(
    n_clusters=2,
    n_starts=8,
    max_iterations=100,
    tolerance=1.0e-8,
    seed=3,
  ) catch {
    BanksyError(message) =>
      abort("invalid BANKSY cluster configuration: " + message)
  }
}

///|
fn main {
  println("=== Bioconductor Banksy Demo ===")
  let (expression, coordinates, gene_names, spot_names) = @src.banksy_example_data()
  let config = banksy_demo_config()
  let model = @src.banksy_compute_with_config(
    expression,
    coordinates,
    config,
    gene_names~,
    spot_names~,
  ) catch {
    BanksyError(message) =>
      abort("BANKSY feature computation failed: " + message)
  }

  println("\n1. Compute spatial neighborhood harmonics")
  println("  " + model.summary())
  println(
    "  first spot left_marker: expression=" +
    model.expression[0][0].to_string() +
    ", H0=" +
    model.harmonics[0][0][0].to_string() +
    ", H1=" +
    model.harmonics[1][0][0].to_string(),
  )
  let mut weight_sum = 0.0
  for neighbor in model.neighbors[0][0] {
    weight_sum = weight_sum + neighbor.weight
  }
  println(
    "  H0 neighbors=" +
    model.neighbors[0][0].length().to_string() +
    ", normalized weight sum=" +
    weight_sum.to_string(),
  )

  println("\n2. Compare cell-typing and spatial-domain lambda values")
  let typing = @src.banksy_run(
    model,
    lambda=0.2,
    n_components=3,
    cluster_config=banksy_demo_cluster_config(),
    max_harmonic=1,
    scale=true,
    smooth=true,
    smoothing_k=3,
    smoothing_threshold=0.5,
    smoothing_iterations=10,
  ) catch {
    BanksyError(message) => abort("BANKSY cell-typing run failed: " + message)
  }
  let domains = @src.banksy_run(
    model,
    lambda=0.8,
    n_components=3,
    cluster_config=banksy_demo_cluster_config(),
    max_harmonic=1,
    scale=true,
    smooth=true,
    smoothing_k=3,
    smoothing_threshold=0.5,
    smoothing_iterations=10,
  ) catch {
    BanksyError(message) => abort("BANKSY domain run failed: " + message)
  }
  println("  lambda=0.2: " + typing.summary())
  println("  lambda=0.8: " + domains.summary())
  println(
    "  PCA dimensions=" +
    typing.pca.n_spots().to_string() +
    "x" +
    typing.pca.n_components.to_string() +
    ", cluster sizes=" +
    typing.clustering.cluster_sizes.to_string(),
  )
  let ari = @src.banksy_adjusted_rand_index(
    typing.final_labels(),
    domains.final_labels(),
  ) catch {
    BanksyError(message) => abort("BANKSY ARI failed: " + message)
  }
  println("  adjusted Rand index between lambda runs=" + ari.to_string())

  println("\n3. Sweep spatial weights")
  let sweep = @src.banksy_parameter_sweep(
    model,
    [0.0, 0.2, 0.5, 0.8],
    n_components=3,
    cluster_config=banksy_demo_cluster_config(),
    max_harmonic=1,
    scale=true,
  ) catch {
    BanksyError(message) => abort("BANKSY parameter sweep failed: " + message)
  }
  for run in sweep {
    println(
      "  lambda=" +
      run.lambda.to_string() +
      ", silhouette=" +
      run.clustering.silhouette.to_string() +
      ", spatial agreement=" +
      run.spatial_agreement.to_string(),
    )
  }

  println("\n4. Enrich an immutable SpatialExperiment")
  let experiment = @src.SpatialExperiment::new()
  ignore(@src.se_add_assay(experiment, "logcounts", expression))
  for gene_name in gene_names {
    ignore(@src.se_add_row(experiment, Map([("gene_name", gene_name)])))
  }
  for spot in 0..<coordinates.length() {
    let section = if coordinates[spot][0] < 3.0 { "left" } else { "right" }
    ignore(
      @src.se_add_col(
        experiment,
        Map([("spot_id", spot_names[spot]), ("section", section)]),
      ),
    )
    ignore(
      @src.se_add_spatial_coord(
        experiment,
        @src.SpatialCoord::new_2d(coordinates[spot][0], coordinates[spot][1]),
      ),
    )
  }
  let integrated = @src.banksy_spatial_experiment(
    experiment,
    config~,
    lambda=0.8,
    n_components=3,
    cluster_config=banksy_demo_cluster_config(),
    max_harmonic=1,
    group_key="section",
    smooth=true,
    smoothing_k=3,
  ) catch {
    BanksyError(message) =>
      abort("BANKSY SpatialExperiment integration failed: " + message)
  }
  println(
    "  assays include H0=" +
    integrated.experiment.assay.contains("H0").to_string() +
    ", H1=" +
    integrated.experiment.assay.contains("H1").to_string(),
  )
  println(
    "  first cluster=" +
    integrated.experiment.col_data[0]["banksy_cluster"] +
    ", input remains unchanged=" +
    (!experiment.assay.contains("H0")).to_string(),
  )
  println("=== Demo Complete ===")
}
