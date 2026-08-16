///|
fn flowsom_demo_ints(values : Array[Int]) -> String {
  values.map(fn(value) { value.to_string() }).join(", ")
}

///|
fn flowsom_demo_optional(value : Double?) -> String {
  match value {
    Some(actual) => actual.to_string()
    None => "NA"
  }
}

///|
fn flowsom_demo_bool_count(values : Array[Bool]) -> Int {
  let mut total = 0
  for value in values {
    if value {
      total = total + 1
    }
  }
  total
}

///|
fn main {
  println("=== Bioconductor FlowSOM 2.21.0 Demo ===")

  // FlowSOM uses cell x marker matrices.
  let cells = [
    [8.0, 1.0, 5.0],
    [8.2, 1.1, 5.1],
    [7.8, 0.9, 4.9],
    [8.1, 1.2, 5.0],
    [1.0, 8.0, 5.0],
    [1.1, 8.2, 5.2],
    [0.9, 7.8, 4.8],
    [1.2, 8.1, 5.1],
    [4.0, 4.0, 9.0],
    [4.2, 4.1, 9.2],
    [3.8, 3.9, 8.8],
    [4.1, 4.2, 9.1],
  ]
  let config = @bio.FlowSomConfig::create(
    xdim=3,
    ydim=2,
    rlen=5,
    mst_runs=2,
    alpha_start=0.08,
    alpha_end=0.01,
    radius_start=2.0,
    radius_end=0.0,
    importance=[1.0, 1.0, 0.75],
    metaclusters=3,
    meta_starts=6,
    outlier_mad=3.0,
    seed=101,
  ) catch {
    _ => abort("failed to create the FlowSOM configuration")
  }
  let model = @bio.flowsom_train(cells, config, marker_names=[
    "CD3", "CD19", "CD45",
  ]) catch {
    _ => abort("FlowSOM training failed")
  }

  println("\n1. Topology-aware SOM training")
  println(@bio.flowsom_summary(model))
  println("  node counts: " + flowsom_demo_ints(model.node_counts))
  println(
    "  cell meta-clusters: " +
    flowsom_demo_ints(model.cell_meta_clusters.map(fn(value) { value + 1 })),
  )

  println("\n2. Minimum spanning tree")
  for edge in model.mst_edges {
    println(
      "  node " +
      (edge.from + 1).to_string() +
      " -> " +
      (edge.to + 1).to_string() +
      ", weight=" +
      edge.weight.to_string(),
    )
  }

  let mut occupied = 0
  while occupied < model.node_counts.length() &&
        model.node_counts[occupied] == 0 {
    occupied = occupied + 1
  }
  println("\n3. Node statistics on unweighted marker values")
  println("  first occupied node: " + (occupied + 1).to_string())
  println(
    "  CD3 median fluorescence: " +
    flowsom_demo_optional(model.node_medians[occupied][0]),
  )
  println(
    "  CD3 coefficient of variation: " +
    flowsom_demo_optional(model.node_cvs[occupied][0]),
  )
  let positive = @bio.flowsom_node_positive_percentages(model, [5.0, 5.0, 7.0]) catch {
    _ => abort("node positivity calculation failed")
  }
  println(
    "  CD3-positive fraction: " + flowsom_demo_optional(positive[occupied][0]),
  )
  println(
    "  training outliers: " +
    flowsom_demo_bool_count(model.outliers.per_cell).to_string(),
  )

  let projected = @bio.flowsom_map_new(model, [
    [model.codes[0][0], model.codes[0][1], model.codes[0][2] / 0.75],
    [20.0, 20.0, 20.0],
  ]) catch {
    _ => abort("new data projection failed")
  }
  println("\n4. New-data projection and MAD outliers")
  println(
    "  mapped nodes: " +
    flowsom_demo_ints(projected.mapping.clusters.map(fn(value) { value + 1 })),
  )
  println(
    "  mapped meta-clusters: " +
    flowsom_demo_ints(projected.meta_clusters.map(fn(value) { value + 1 })),
  )
  println(
    "  outlier flags: " +
    projected.outliers.map(fn(value) { value.to_string() }).join(", "),
  )

  let frame = @bio.FlowFrame::new(cells, [
    @bio.ParameterDescription::new("CD3", 25.0, 0.0, 25.0),
    @bio.ParameterDescription::new("CD19", 25.0, 0.0, 25.0),
    @bio.ParameterDescription::new("CD45", 25.0, 0.0, 25.0),
  ])
  let frame_config = @bio.FlowSomConfig::create(
    xdim=3,
    ydim=2,
    rlen=5,
    mst_runs=2,
    importance=[1.0, 0.75],
    metaclusters=3,
    meta_starts=6,
    seed=103,
  ) catch {
    _ => abort("failed to create the FlowFrame configuration")
  }
  let frame_model = @bio.flowsom_train_flow_frame(frame, frame_config, marker_indices=[
    0, 2,
  ]) catch {
    _ => abort("FlowFrame integration failed")
  }
  println("\n5. flowCore FlowFrame integration")
  println("  selected markers: " + frame_model.marker_names.join(", "))
  println(
    "  mapped events: " + frame_model.mapping.clusters.length().to_string(),
  )

  // SingleCellExperiment assays use marker x cell orientation.
  let experiment = @bio.SingleCellExperiment::new(
    [
      [8.0, 8.2, 7.8, 8.1, 1.0, 1.1, 0.9, 1.2, 4.0, 4.2, 3.8, 4.1],
      [1.0, 1.1, 0.9, 1.2, 8.0, 8.2, 7.8, 8.1, 4.0, 4.1, 3.9, 4.2],
      [5.0, 5.1, 4.9, 5.0, 5.0, 5.2, 4.8, 5.1, 9.0, 9.2, 8.8, 9.1],
    ],
    ["CD3", "CD19", "CD45"],
    ["C1", "C2", "C3", "C4", "C5", "C6", "C7", "C8", "C9", "C10", "C11", "C12"],
  )
  experiment.metadata["source"] = "flowsom_demo"
  let clustered = @bio.flowsom_cluster_sce(experiment, config) catch {
    _ => abort("SingleCellExperiment integration failed")
  }
  println("\n6. Immutable SingleCellExperiment integration")
  println(
    "  SOM labels: " +
    clustered.experiment.col_data["FlowSOM.cluster"].join(", "),
  )
  println(
    "  meta-cluster labels: " +
    clustered.experiment.col_data["FlowSOM.metacluster"].join(", "),
  )
  println(
    "  original unchanged: " +
    (!experiment.col_data.contains("FlowSOM.cluster")).to_string(),
  )
  println("\nFlowSOM demo completed.")
}
