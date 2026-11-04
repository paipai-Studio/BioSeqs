// Bioconductor nnSVG-inspired spatially variable gene analysis.

///|
fn main {
  let (expression, coordinates, gene_names) = @src.nnsvg_example_data()
  let config = @src.NnsvgConfig::create(
    n_neighbors=4,
    length_scale_grid=8,
    proportion_grid=8,
    refinement_steps=2,
    fdr_threshold=0.1,
  ) catch {
    NnsvgError(message) => abort("invalid nnSVG configuration: " + message)
  }

  let result = @src.nnsvg_with_config(
    expression,
    coordinates,
    config,
    gene_names~,
  ) catch {
    NnsvgError(message) => abort("nnSVG fit failed: " + message)
  }

  println("=== Direct nnSVG analysis ===")
  println(result.summary())
  println("rank\tgene\tLR\tp-value\tFDR\tprop_sv\tlength_scale")
  for gene in result.top(result.n_genes()) {
    println(
      gene.rank.to_string() +
      "\t" +
      gene.gene_name +
      "\t" +
      gene.likelihood_ratio.to_string() +
      "\t" +
      gene.p_value.to_string() +
      "\t" +
      gene.adjusted_p_value.to_string() +
      "\t" +
      gene.proportion_spatial_variance.to_string() +
      "\t" +
      gene.length_scale.to_string(),
    )
  }

  let counts = [
    [0.0, 0.0, 4.0, 5.0, 0.0, 0.0],
    [8.0, 9.0, 7.0, 8.0, 9.0, 8.0],
    [5.0, 5.0, 5.0, 5.0, 5.0, 5.0],
  ]
  let filter = @src.nnsvg_filter_genes(
    counts,
    ["low_expression", "spatial_candidate", "MT-ND1"],
    minimum_count=3.0,
    minimum_spot_percentage=50.0,
  ) catch {
    NnsvgError(message) => abort("nnSVG filtering failed: " + message)
  }
  println("\n=== Upstream-compatible gene filtering ===")
  println(
    "kept=" +
    filter.kept_indices.length().to_string() +
    ", low-expression=" +
    filter.removed_low_expression.length().to_string() +
    ", mitochondrial=" +
    filter.removed_mitochondrial.length().to_string(),
  )

  let experiment = @src.SpatialExperiment::new()
  ignore(@src.se_add_assay(experiment, "logcounts", expression))
  for gene_name in gene_names {
    ignore(@src.se_add_row(experiment, Map([("gene_name", gene_name)])))
  }
  for spot in 0..<coordinates.length() {
    ignore(
      @src.se_add_col(
        experiment,
        Map([("spot_id", "spot_" + (spot + 1).to_string())]),
      ),
    )
    ignore(
      @src.se_add_spatial_coord(
        experiment,
        @src.SpatialCoord::new_2d(coordinates[spot][0], coordinates[spot][1]),
      ),
    )
  }
  let integrated = @src.nnsvg_spatial_experiment(experiment, config~) catch {
    NnsvgError(message) =>
      abort("nnSVG SpatialExperiment integration failed: " + message)
  }
  println("\n=== SpatialExperiment integration ===")
  println(
    "rowData fields written for " +
    integrated.experiment.row_data.length().to_string() +
    " genes; first-gene padj=" +
    integrated.experiment.row_data[0]["nnsvg_padj"],
  )
}
