///|
fn single_r_demo_optional_label(value : String?) -> String {
  match value {
    Some(label) => label
    None => "NA"
  }
}

///|
fn main {
  println("=== Bioconductor SingleR 2.15.2 Advanced Demo ===")
  let (reference, data) = @bio.single_r_advanced_example() catch {
    _ => abort("failed to construct the SingleR example")
  }
  println(
    "\n1. Gene x sample/cell inputs: reference=" +
    reference.n_genes.to_string() +
    "x" +
    reference.n_samples.to_string() +
    ", test=" +
    data.n_genes.to_string() +
    "x" +
    data.n_cells.to_string(),
  )

  let training = @bio.single_r_advanced_train(reference, data.gene_names) catch {
    _ => abort("SingleR marker training failed")
  }
  let t_vs_b = training.markers_between("T cell", "B cell") catch {
    _ => abort("failed to query pairwise markers")
  }
  println("\n2. Directed classic marker training")
  println(
    "  shared genes=" +
    training.common_gene_names.length().to_string() +
    ", selected markers=" +
    training.marker_gene_names().length().to_string(),
  )
  println("  T cell > B cell markers: " + t_vs_b.join(", "))

  let config = @bio.SingleRAdvancedConfig::create(
    quantile=0.8,
    fine_tune=true,
    prune=true,
    nmads=3.0,
  ) catch {
    _ => abort("failed to construct the SingleR configuration")
  }
  let result = @bio.single_r_advanced_classify(data, training, config~) catch {
    _ => abort("SingleR cell annotation failed")
  }
  println("\n3. Quantile scores, iterative fine-tuning and pruning")
  for cell in 0..<result.n_cells() {
    let label = result.labels[cell]
    let score = result.scores[cell][result.label_index(label)]
    println(
      "  " +
      result.cell_names[cell] +
      ": label=" +
      label +
      ", pruned=" +
      single_r_demo_optional_label(result.pruned_labels[cell]) +
      ", score=" +
      score.to_string() +
      ", delta.next=" +
      result.delta_next[cell].to_string() +
      ", fine.tuned=" +
      result.fine_tuned[cell].to_string(),
    )
  }

  let clusters = ["lymphoid", "lymphoid", "myeloid", "cytotoxic", "lymphoid"]
  let clustered = @bio.single_r_advanced_clusters(
    data,
    clusters,
    reference,
    config~,
  ) catch {
    _ => abort("SingleR cluster annotation failed")
  }
  println("\n4. Cluster-level annotation from summed profiles")
  for cluster in 0..<clustered.n_cells() {
    println(
      "  " + clustered.cell_names[cluster] + ": " + clustered.labels[cluster],
    )
  }

  let alternate = @bio.SingleRAdvancedReference::create(
    reference.expression,
    [
      "T lymphocyte", "T lymphocyte", "B lymphocyte", "B lymphocyte", "Myeloid",
      "Myeloid", "Cytotoxic", "Cytotoxic",
    ],
    gene_names=reference.gene_names,
    sample_names=reference.sample_names,
    reference_name="immune_reference_two",
  ) catch {
    _ => abort("failed to construct the second reference")
  }
  let combined = @bio.single_r_advanced_combine(
    data,
    [reference, alternate],
    config~,
  ) catch {
    _ => abort("SingleR multi-reference integration failed")
  }
  println("\n5. Multi-reference marker-space recomputation")
  for cell in 0..<combined.n_cells() {
    println(
      "  " +
      combined.cell_names[cell] +
      ": " +
      combined.labels[cell] +
      " from " +
      combined.selected_references[cell],
    )
  }

  let experiment = @bio.SingleCellExperiment::new(
    data.expression,
    data.gene_names,
    data.cell_names,
  )
  experiment.assays["logcounts"] = data.expression.map(fn(row) { row.copy() })
  experiment.col_data["cluster"] = clusters
  experiment.metadata["source"] = "single_r_demo"
  let output = @bio.single_r_advanced_sce(
    experiment,
    reference,
    cluster_column="cluster",
    output_prefix="SingleR215",
    config~,
  ) catch {
    _ => abort("SingleR SingleCellExperiment integration failed")
  }
  println("\n6. Immutable SingleCellExperiment write-back")
  println(
    "  mode=" +
    output.experiment.metadata["SingleR215.mode"] +
    ", labels=" +
    output.experiment.col_data["SingleR215.labels"].join(", "),
  )
  println(
    "  original unchanged=" +
    (!experiment.col_data.contains("SingleR215.labels")).to_string(),
  )
  println("\nSingleR advanced demo completed.")
}
