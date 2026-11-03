///|
fn muscat_demo_config() -> @bio.MuscatAdvancedConfig {
  @bio.MuscatAdvancedConfig::create(
    min_cells=3,
    min_count=1.0,
    min_samples=2,
    max_iterations=100,
    tolerance=1.0e-8,
    dispersion_prior_df=10.0,
    ridge=1.0e-6,
    fdr_threshold=0.1,
    detection_filter=0.9,
  ) catch {
    _ => abort("failed to create muscat advanced configuration")
  }
}

///|
fn muscat_demo_result(
  result : @bio.MuscatAdvancedResults,
  genes : Array[String],
  clusters : Array[String],
  contrast : String,
) -> Unit {
  println("\n" + result.summary())
  for cluster in clusters {
    println("  cluster " + cluster)
    for gene in genes {
      match result.get(gene, cluster, contrast) {
        Some(value) =>
          println(
            "    " +
            gene +
            ": log2FC=" +
            value.log_fc.to_string() +
            ", Wald=" +
            value.statistic.to_string() +
            ", p=" +
            value.p_value.to_string() +
            ", local FDR=" +
            value.local_fdr.to_string() +
            ", global FDR=" +
            value.global_fdr.to_string(),
          )
        None => ()
      }
    }
  }
}

///|
fn muscat_demo_sce(data : @bio.MuscatAdvancedData) -> @bio.SingleCellExperiment {
  let experiment = @bio.SingleCellExperiment::new(
    data.counts,
    data.gene_names,
    data.cell_names,
  )
  experiment.col_data["sample"] = data.sample_ids.copy()
  experiment.col_data["cluster"] = data.cluster_ids.copy()
  experiment.col_data["group"] = data.group_ids.copy()
  experiment.row_data["symbol"] = data.gene_names.copy()
  experiment.metadata["source"] = "muscat advanced demo"
  experiment
}

///|
fn main {
  println("=== Bioconductor muscat 1.27.4 Advanced Demo ===")
  let data = @bio.muscat_advanced_example() catch {
    _ => abort("failed to construct muscat example data")
  }
  let config = muscat_demo_config()

  println("\n1. Strict gene x cell data contract")
  println(
    "  " +
    data.gene_names.length().to_string() +
    " genes, " +
    data.cell_names.length().to_string() +
    " cells, " +
    data.sample_names.length().to_string() +
    " samples, " +
    data.cluster_names.length().to_string() +
    " clusters",
  )

  let sums = @bio.muscat_aggregate_advanced(
    data,
    aggregation=@bio.muscat_sum_counts(),
  ) catch {
    _ => abort("sum-count aggregation failed")
  }
  let detections = @bio.muscat_aggregate_advanced(
    data,
    aggregation=@bio.muscat_number_detected(),
  ) catch {
    _ => abort("number-detected aggregation failed")
  }
  println("\n2. Cluster-sample pseudobulk aggregation")
  for cluster in 0..<sums.cluster_names.length() {
    println(
      "  " +
      sums.cluster_names[cluster] +
      ": cells/sample=" +
      sums.cell_counts[cluster][0].to_string() +
      ", first library=" +
      sums.library_sizes[cluster][0].to_string(),
    )
  }

  let design = @bio.muscat_group_design_advanced(sums, reference="ctrl") catch {
    _ => abort("group design construction failed")
  }
  let contrasts = @bio.muscat_default_contrasts_advanced(design) catch {
    _ => abort("default contrast construction failed")
  }
  let contrast = contrasts[0].name
  println("\n3. Replicated design and contrast")
  println("  coefficients: " + design.coefficient_names.join(", "))
  println("  contrast: " + contrast)

  let ds = @bio.muscat_pbds_advanced(sums, design, contrasts, config~) catch {
    _ => abort("differential-state model failed")
  }
  let dd = @bio.muscat_pbdd_advanced(detections, design, contrasts, config~) catch {
    _ => abort("differential-detection model failed")
  }
  muscat_demo_result(ds, data.gene_names, data.cluster_names, contrast)
  muscat_demo_result(dd, data.gene_names, data.cluster_names, contrast)

  let stagewise = @bio.muscat_stagewise_ds_dd(ds, dd, alpha=0.1) catch {
    _ => abort("stagewise DS/DD testing failed")
  }
  println("\n4. Harmonic-mean screening and two-stage confirmation")
  for result in stagewise.results {
    if result.classification != "none" {
      println(
        "  " +
        result.cluster_name +
        "/" +
        result.gene_name +
        ": " +
        result.classification +
        ", screen FDR=" +
        result.screen_fdr.to_string(),
      )
    }
  }

  let sce_output = @bio.muscat_advanced_sce(
    muscat_demo_sce(data),
    "sample",
    "cluster",
    "group",
    reference="ctrl",
    output_prefix="muscat",
    config~,
  ) catch {
    _ => abort("SingleCellExperiment integration failed")
  }
  println("\n5. Immutable SingleCellExperiment write-back")
  println(
    "  metadata version: " + sce_output.experiment.metadata["muscat.version"],
  )
  println(
    "  row annotations: muscat.A.dsLogFC, muscat.A.dsFdr, " +
    "muscat.A.ddLogFC, muscat.A.ddFdr, muscat.A.stageClass",
  )
}
