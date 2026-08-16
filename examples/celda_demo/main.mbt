// Bioconductor celda_CG-inspired joint cell and feature clustering workflow.

///|
fn celda_demo_config() -> @src.CeldaCGConfig {
  @src.CeldaCGConfig::create(
    3,
    4,
    max_iterations=20,
    stop_iterations=4,
    chains=3,
    seed=2026,
  ) catch {
    _ => abort("invalid celda_CG configuration")
  }
}

///|
fn main {
  println("=== Bioconductor celda_CG Demo ===")
  let (counts, features, cells, samples) = @src.celda_cg_example_data()

  println("\n1. Jointly cluster cells and feature modules")
  let result = @src.celda_cg_fit(
    counts,
    samples,
    celda_demo_config(),
    feature_names=features,
    cell_names=cells,
  ) catch {
    _ => abort("celda_CG model fitting failed")
  }
  println("   " + result.summary())
  println(
    "   chains=" +
    result.diagnostics.chain_scores.length().to_string() +
    ", selected chain=" +
    (result.diagnostics.best_chain + 1).to_string(),
  )
  for population in 0..<result.config.cell_populations {
    let members = result.cells_in_population(population) catch {
      _ => abort("population query failed")
    }
    println(
      "   population " +
      (population + 1).to_string() +
      ": " +
      members.length().to_string() +
      " cells",
    )
  }

  println("\n2. Inspect feature modules and their strongest markers")
  for module_ in 0..<result.config.feature_modules {
    let top = result.top_features(module_, limit=2) catch {
      _ => abort("feature module query failed")
    }
    let mut marker_text = ""
    for index in 0..<top.length() {
      if index > 0 {
        marker_text = marker_text + ", "
      }
      marker_text = marker_text +
        top[index].feature_name +
        " (" +
        top[index].probability.to_string() +
        ")"
    }
    println("   module " + (module_ + 1).to_string() + ": " + marker_text)
  }

  println("\n3. Predict populations for new cells")
  let new_counts : Array[Array[Double]] = []
  for feature in 0..<counts.length() {
    new_counts.push([counts[feature][0], counts[feature][4], counts[feature][8]])
  }
  let prediction = @src.celda_cg_predict_cells(result, new_counts) catch {
    _ => abort("celda_CG prediction failed")
  }
  for cell in 0..<prediction.assignments.length() {
    println(
      "   query " +
      (cell + 1).to_string() +
      " -> population " +
      (prediction.assignments[cell] + 1).to_string(),
    )
  }

  println("\n4. Compare candidate K values by BIC")
  let grid_config = @src.CeldaCGConfig::create(
    1,
    1,
    max_iterations=8,
    stop_iterations=2,
    chains=1,
    seed=31,
  ) catch {
    _ => abort("invalid grid-search configuration")
  }
  let grid = @src.celda_cg_grid_search(
    counts,
    samples,
    [2, 3],
    [4],
    grid_config,
    feature_names=features,
    cell_names=cells,
  ) catch {
    _ => abort("celda_CG grid search failed")
  }
  for entry in grid.entries {
    println(
      "   K=" +
      entry.cell_populations.to_string() +
      ", L=" +
      entry.feature_modules.to_string() +
      ", BIC=" +
      entry.score.to_string(),
    )
  }
  println(
    "   selected K=" + grid.best_result().config.cell_populations.to_string(),
  )

  println("\n5. Add fitted counts and labels to an SCE copy")
  let sce = @src.SingleCellExperiment::new(counts, features, cells)
  sce.col_data["donor"] = samples
  let output = @src.celda_cg_sce(
    sce,
    celda_demo_config(),
    sample_column="donor",
  ) catch {
    _ => abort("celda_CG SCE integration failed")
  }
  println(
    "   fitted assay rows=" +
    @src.sce_get_assay(output.experiment, "celda_fitted").length().to_string(),
  )
  println(
    "   cell labels=" +
    @src.sce_get_col_data(output.experiment, "celda_cell_population")
    .length()
    .to_string() +
    ", source unchanged=" +
    (@src.sce_get_assay(sce, "celda_fitted").length() == 0).to_string(),
  )
  println("\n=== Demo Complete ===")
}
