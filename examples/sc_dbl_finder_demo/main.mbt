///|
fn main {
  let (data, clusters, samples, known_doublets) = @bio.scdf_create_advanced_example() catch {
    _ => abort("failed to create the scDblFinder example")
  }
  let config = @bio.ScDblFinderConfig::create(
    expected_doublet_rate=0.15,
    artificial_doublets=48,
    selected_features=18,
    dimensions=5,
    neighbors=8,
    iterations=2,
    classifier_steps=120,
    seed=17,
  ) catch {
    _ => abort("failed to create the scDblFinder configuration")
  }
  let result = @bio.sc_dbl_finder(
    data,
    clusters~,
    samples~,
    known_doublets~,
    config~,
  ) catch {
    _ => abort("scDblFinder fitting failed")
  }

  println(result.summary())
  println("\nCapture thresholds")
  for capture in 0..<result.sample_ids.length() {
    let threshold = result.thresholds[capture]
    println(
      result.sample_ids[capture] +
      "\tthreshold=" +
      threshold.threshold.to_string() +
      "\tcalled_rate=" +
      threshold.called_rate.to_string() +
      "\tFNR=" +
      threshold.false_negative_rate.to_string(),
    )
  }

  println("\nHighest-scoring cells")
  let top = result.top_doublets(limit=8) catch {
    _ => abort("failed to rank doublets")
  }
  for cell in top {
    println(
      cell.cell_name +
      "\tscore=" +
      cell.score.to_string() +
      "\tclass=" +
      (if cell.is_doublet { "doublet" } else { "singlet" }) +
      "\torigin=" +
      cell.origin,
    )
  }

  println("\nPairwise origin enrichment")
  let enrichment = @bio.scdf_pairwise_enrichment(result) catch {
    _ => abort("failed to compute origin enrichment")
  }
  for row in enrichment {
    println(
      row.combination +
      "\tobserved=" +
      row.observed.to_string() +
      "\texpected=" +
      row.expected.to_string() +
      "\tadjusted_p=" +
      row.adjusted_p_value.to_string(),
    )
  }

  let singlets = result.filter_singlets(data) catch {
    _ => abort("failed to filter singlets")
  }
  println(
    "\nRetained " +
    singlets.cell_names.length().to_string() +
    " singlets from " +
    data.cell_names.length().to_string() +
    " cells",
  )

  let automatic_clusters = @bio.scdf_fast_cluster(
    data,
    3,
    dimensions=5,
    selected_features=18,
  ) catch {
    _ => abort("automatic clustering failed")
  }
  println("Automatic clusters for first five cells")
  for cell in 0..<5 {
    println(data.cell_names[cell] + "\t" + automatic_clusters[cell])
  }

  let assay : Array[Array[Double]] = []
  for gene in 0..<data.gene_names.length() {
    let row : Array[Double] = []
    for cell in 0..<data.cell_names.length() {
      row.push(data.counts[cell][gene])
    }
    assay.push(row)
  }
  let experiment = @bio.SingleCellExperiment::new(
    assay,
    data.gene_names,
    data.cell_names,
  )
  experiment.col_data["cluster"] = clusters
  experiment.col_data["capture"] = samples
  let truth_labels : Array[String] = []
  for value in known_doublets {
    truth_labels.push(if value { "doublet" } else { "singlet" })
  }
  experiment.col_data["known_doublet"] = truth_labels
  let output = @bio.sc_dbl_finder_single_cell_experiment(
    experiment,
    cluster_column="cluster",
    sample_column="capture",
    known_doublet_column="known_doublet",
    config~,
  ) catch {
    _ => abort("SingleCellExperiment integration failed")
  }
  println(
    "\nSingleCellExperiment write-back: score=" +
    output.experiment.col_data["scDblFinder.score"][0] +
    ", original_unchanged=" +
    (!experiment.col_data.contains("scDblFinder.score")).to_string(),
  )
}
