///|
fn scuttle_demo_doubles(values : Array[Double]) -> String {
  values.map(fn(value) { value.to_string() }).join(", ")
}

///|
fn main {
  println("=== Bioconductor scuttle 1.23.1 Demo ===")

  let outlier_config = @bio.ScuttleOutlierConfig::create(
    nmads=2.0,
    direction=@bio.scuttle_outlier_higher(),
    batches=["A", "A", "A", "B", "B", "B"],
  ) catch {
    _ => abort("failed to create the outlier configuration")
  }
  let outliers = @bio.scuttle_is_outlier(
    [9.0, 10.0, 11.0, 98.0, 100.0, 150.0],
    config=outlier_config,
  ) catch {
    _ => abort("batch-aware outlier detection failed")
  }
  println("\n1. Batch-aware MAD outlier detection")
  println("  batches: " + outliers.batch_names.join(", "))
  println("  upper thresholds: " + scuttle_demo_doubles(outliers.higher))
  println("  discarded observations: " + outliers.discarded_count().to_string())

  let counts = [
    [8.0, 4.0, 16.0, 8.0],
    [2.0, 6.0, 4.0, 12.0],
    [0.0, 2.0, 0.0, 4.0],
    [10.0, 8.0, 20.0, 16.0],
  ]
  let controls = @bio.ScuttleNamedCellSubset::create("controls", [0, 1]) catch {
    _ => abort("failed to create the control-cell subset")
  }
  let qc = @bio.scuttle_per_feature_qc(counts, subsets=[controls]) catch {
    _ => abort("per-feature QC failed")
  }
  println("\n2. Per-feature QC with a named cell subset")
  println("  means: " + scuttle_demo_doubles(qc.means))
  println(
    "  detected percentages: " + scuttle_demo_doubles(qc.detected_percent),
  )
  println(
    "  control/global ratios: " + scuttle_demo_doubles(qc.subset_ratios[0]),
  )

  let pathway_a = @bio.ScuttleFeatureSet::create("PathwayA", [0, 1]) catch {
    _ => abort("failed to create PathwayA")
  }
  let pathway_b = @bio.ScuttleFeatureSet::create("PathwayB", [1, 2, 3]) catch {
    _ => abort("failed to create PathwayB")
  }
  let aggregated = @bio.scuttle_aggregate_feature_sets(counts, [
    pathway_a, pathway_b,
  ]) catch {
    _ => abort("feature-set aggregation failed")
  }
  println("\n3. Overlapping feature-set aggregation")
  for index in 0..<aggregated.names.length() {
    println(
      "  " +
      aggregated.names[index] +
      " (n=" +
      aggregated.sizes[index].to_string() +
      "): " +
      scuttle_demo_doubles(aggregated.values[index]),
    )
  }

  let balanced = @bio.scuttle_downsample_batches(
    counts,
    ["A", "A", "B", "B"],
    seed=42,
  ) catch {
    _ => abort("batch coverage equalization failed")
  }
  println("\n4. Exact batch coverage equalization")
  println("  batch summaries: " + scuttle_demo_doubles(balanced.summaries))
  println("  proportions: " + scuttle_demo_doubles(balanced.proportions))

  let experiment = @bio.SingleCellExperiment::new(
    counts,
    ["G1", "G2", "G3", "G4"],
    ["C1", "C2", "C3", "C4"],
  )
  experiment.col_data["batch"] = ["A", "A", "B", "B"]
  experiment.metadata["source"] = "scuttle_demo"
  let qc_output = @bio.scuttle_per_feature_qc_sce(experiment, subsets=[controls]) catch {
    _ => abort("SCE per-feature QC failed")
  }
  let downsampled = @bio.scuttle_downsample_sce(
    qc_output.experiment,
    0.5,
    output_assay="half",
    seed=42,
  ) catch {
    _ => abort("SCE downsampling failed")
  }
  let pathway_sce = @bio.scuttle_aggregate_feature_sets_sce(
    downsampled,
    [pathway_a, pathway_b],
    assay_names=["counts", "half"],
  ) catch {
    _ => abort("SCE feature aggregation failed")
  }
  println("\n5. Immutable SingleCellExperiment integration")
  println(
    "  QC rows written: " +
    qc_output.experiment.row_data["scuttle.mean"].length().to_string(),
  )
  println(
    "  downsampled assay rows: " +
    downsampled.assays["half"].length().to_string(),
  )
  println("  aggregated rows: " + pathway_sce.row_names.join(", "))
  println(
    "  original unchanged: " +
    (!experiment.row_data.contains("scuttle.mean") &&
    !experiment.assays.contains("half")).to_string(),
  )
  println("\nscuttle demo completed.")
}
