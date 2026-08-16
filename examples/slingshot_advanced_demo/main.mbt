///|
fn slingshot_demo_coordinates() -> Array[Array[Double]] {
  [
    [-0.1, 0.0],
    [0.0, 0.1],
    [0.1, -0.1],
    [0.9, 0.0],
    [1.0, 0.1],
    [1.1, -0.1],
    [1.9, 0.0],
    [2.0, 0.1],
    [2.1, -0.1],
    [2.9, 0.9],
    [3.0, 1.0],
    [3.1, 1.1],
    [2.9, -0.9],
    [3.0, -1.0],
    [3.1, -1.1],
  ]
}

///|
fn slingshot_demo_labels() -> Array[String] {
  ["A", "A", "A", "B", "B", "B", "C", "C", "C", "D", "D", "D", "E", "E", "E"]
}

///|
fn slingshot_demo_lineage(
  lineage : Array[Int],
  cluster_names : Array[String],
) -> String {
  lineage.map(fn(index) { cluster_names[index] }).join(" -> ")
}

///|
fn slingshot_demo_doubles(values : Array[Double]) -> String {
  values.map(fn(value) { value.to_string() }).join(", ")
}

///|
fn main {
  println("=== Bioconductor slingshot 2.21.0 Advanced Demo ===")
  let coordinates = slingshot_demo_coordinates()
  let labels = slingshot_demo_labels()
  let config = @bio.SlingshotAdvancedConfig::create(
    start_clusters=["A"],
    end_clusters=["D", "E"],
    distance=@bio.slingshot_scaled_full(),
    extension=@bio.slingshot_extend_line(),
    shrink=1.0,
    reweight=true,
    reassign=true,
    max_iterations=8,
    tolerance=1.0e-4,
    smoother_span=0.4,
    curve_points=30,
  ) catch {
    _ => abort("failed to create the advanced slingshot configuration")
  }
  let result = @bio.slingshot_advanced(coordinates, labels, config) catch {
    _ => abort("advanced slingshot fitting failed")
  }

  println("\n1. Constrained cluster minimum spanning tree")
  for edge in result.lineage_model.edges {
    println(
      "  " +
      result.lineage_model.cluster_names[edge.from] +
      " -> " +
      result.lineage_model.cluster_names[edge.to] +
      ", distance=" +
      edge.distance.to_string(),
    )
  }

  println("\n2. Root-to-leaf lineages and simultaneous curves")
  for index in 0..<result.lineage_model.lineages.length() {
    println(
      "  " +
      result.curves[index].name +
      ": " +
      slingshot_demo_lineage(
        result.lineage_model.lineages[index],
        result.lineage_model.cluster_names,
      ) +
      ", arc length=" +
      result.curves[index].length.to_string(),
    )
  }
  println("  iterations: " + result.iterations.to_string())
  println("  converged: " + result.converged.to_string())

  println("\n3. Cell x lineage pseudotime weights")
  for cell in [0, 7, 10, 13] {
    println(
      "  cell " +
      (cell + 1).to_string() +
      ": branch=" +
      result.branch_ids[cell] +
      ", average pseudotime=" +
      result.average_pseudotime[cell].to_string() +
      ", weights=[" +
      slingshot_demo_doubles(result.weights[cell]) +
      "]",
    )
  }

  let projected = @bio.slingshot_predict(result, [[2.8, 0.8], [2.8, -0.8]]) catch {
    _ => abort("new-cell trajectory projection failed")
  }
  println("\n4. New-cell projection")
  for cell in 0..<projected.weights.length() {
    println(
      "  query " +
      (cell + 1).to_string() +
      ": branch=" +
      projected.branch_ids[cell] +
      ", average pseudotime=" +
      projected.average_pseudotime[cell].to_string() +
      ", weights=[" +
      slingshot_demo_doubles(projected.weights[cell]) +
      "]",
    )
  }

  // SingleCellExperiment assays use feature x cell orientation.
  let experiment = @bio.SingleCellExperiment::new(
    [
      [
        1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0, 13.0, 14.0,
        15.0,
      ],
    ],
    ["gene"],
    [
      "c1", "c2", "c3", "c4", "c5", "c6", "c7", "c8", "c9", "c10", "c11", "c12",
      "c13", "c14", "c15",
    ],
  )
  experiment.reduced_dims["PCA"] = coordinates
  experiment.col_data["cluster"] = labels
  let sce_output = @bio.slingshot_advanced_sce(experiment, config) catch {
    _ => abort("SingleCellExperiment trajectory integration failed")
  }
  println("\n5. Immutable SingleCellExperiment integration")
  println("  output columns: slingshot.branch, slingshot.pseudotime")
  println(
    "  pseudotime matrix: " +
    sce_output.experiment.reduced_dims["slingshot.pseudotime"]
    .length()
    .to_string() +
    " cells x " +
    sce_output.experiment.reduced_dims["slingshot.pseudotime"][0]
    .length()
    .to_string() +
    " lineages",
  )
  println(
    "  original unchanged: " +
    (!experiment.col_data.contains("slingshot.branch")).to_string(),
  )
  println("\n" + @bio.slingshot_advanced_summary(result))
}
