///|
fn format_ints(values : Array[Int]) -> String {
  let mut output = "["
  for index in 0..<values.length() {
    if index > 0 {
      output = output + ", "
    }
    output = output + values[index].to_string()
  }
  output + "]"
}

///|
fn format_doubles(values : Array[Double]) -> String {
  let mut output = "["
  for index in 0..<values.length() {
    if index > 0 {
      output = output + ", "
    }
    output = output + values[index].to_string()
  }
  output + "]"
}

///|
fn main {
  println("=== Bioconductor miloR Demo ===")
  let data = @src.milo_sample_data()

  println("\n1. Exact KNN graph from reduced dimensions")
  let graph = @src.milo_build_graph(
    data.coordinates,
    cell_names=data.cell_names,
    k=8,
    dimensions=2,
  ) catch {
    _ => abort("failed to build the Milo KNN graph")
  }
  println("   " + graph.summary())
  println("   first directed KNN: " + format_ints(graph.knn_indices[0]))
  println(
    "   first undirected neighbors: " + format_ints(graph.graph_neighbors[0]),
  )

  println("\n2. Refined, overlapping neighborhoods")
  let neighborhoods = graph.make_neighborhoods(
    proportion=0.5,
    refined=true,
    seed=19,
  ) catch {
    _ => abort("failed to sample Milo neighborhoods")
  }
  println(
    "   refined representatives: " +
    neighborhoods.neighborhood_indices.length().to_string(),
  )
  println(
    "   first neighborhood: " + format_ints(neighborhoods.neighborhoods[0]),
  )

  println("\n3. Neighborhood-by-sample counts and mean expression")
  let counted = neighborhoods.count_cells(data.sample_ids) catch {
    _ => abort("failed to count cells by sample")
  }
  let feature_x : Array[Double] = []
  let feature_y : Array[Double] = []
  for coordinate in data.coordinates {
    feature_x.push(coordinate[0])
    feature_y.push(coordinate[1])
  }
  let means = counted.neighborhood_expression([feature_x, feature_y]) catch {
    _ => abort("failed to aggregate neighborhood expression")
  }
  println("   sample order: " + counted.sample_names.to_string())
  println("   first count row: " + format_ints(counted.neighborhood_counts[0]))
  println("   mean x by neighborhood: " + format_doubles(means[0]))

  println("\n4. Fixed-effect negative-binomial differential abundance")
  let design = @src.milo_binary_design(data.sample_conditions, reference="A") catch {
    _ => abort("failed to construct the Milo design matrix")
  }
  let results = counted.test_neighborhoods(design.matrix, 1, cell_sizes=[
    4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0,
  ]) catch {
    _ => abort("failed to test Milo neighborhoods")
  }
  println("   coefficient: " + design.coefficient_name)
  let mut left_depleted = 0
  let mut right_enriched = 0
  for result in results {
    let x = counted.coordinates[result.index_cell][0]
    if x < 2.5 && result.log_fc < 0.0 {
      left_depleted = left_depleted + 1
    } else if x >= 2.5 && result.log_fc > 0.0 {
      right_enriched = right_enriched + 1
    }
    println(
      "   nhood " +
      (result.neighborhood + 1).to_string() +
      ": log2FC=" +
      result.log_fc.to_string() +
      ", FDR=" +
      result.fdr.to_string() +
      ", spatial FDR=" +
      result.spatial_fdr.to_string(),
    )
  }
  println(
    "   state directions: left depleted=" +
    left_depleted.to_string() +
    ", right enriched=" +
    right_enriched.to_string(),
  )

  println("\n5. Alternative graph-aware FDR and SCE integration")
  let p_values : Array[Double] = []
  for result in results {
    p_values.push(result.p_value)
  }
  let overlap_fdr = @src.milo_graph_spatial_fdr(
    counted,
    p_values,
    weighting=@src.milo_graph_overlap_weighting(),
  ) catch {
    _ => abort("failed to compute graph-overlap FDR")
  }
  println("   graph-overlap FDR: " + format_doubles(overlap_fdr))

  let sce = @src.SingleCellExperiment::new(
    [Array::make(data.cell_names.length(), 1.0)],
    ["marker"],
    data.cell_names,
  )
  let with_pca = @src.sce_set_reduced_dim(sce, "PCA", data.coordinates)
  let sce_graph = @src.milo_from_single_cell_experiment(
    with_pca,
    reduced_dim="PCA",
    k=8,
  ) catch {
    _ => abort("failed to construct Milo from SingleCellExperiment")
  }
  println("   SCE-backed graph: " + sce_graph.summary())

  println(
    "\nScope: portable fixed-effect NB-GLM; GLMM, edgeR TMM/RLE/QL backends and plotting are not included.",
  )
  println("\n=== Demo Complete ===")
}
