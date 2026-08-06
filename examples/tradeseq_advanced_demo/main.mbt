///|
fn tradeseq_demo_counts() -> Array[Array[Double]] {
  let flat : Array[Double] = []
  let association : Array[Double] = []
  let endpoint : Array[Double] = []
  let early_branch : Array[Double] = []
  let zero : Array[Double] = []
  for cell in 0..<24 {
    flat.push((10 + cell % 3).to_double())
    if cell < 8 {
      association.push((3 + cell).to_double())
      endpoint.push((6 + cell % 2).to_double())
      early_branch.push((8 + cell % 2).to_double())
    } else if cell < 16 {
      let step = cell - 8
      association.push((12 + 3 * step).to_double())
      endpoint.push((10 + 4 * step).to_double())
      early_branch.push(
        if step < 4 {
          (26 + step).to_double()
        } else {
          (14 + step % 2).to_double()
        },
      )
    } else {
      let step = cell - 16
      association.push((12 + 3 * step).to_double())
      endpoint.push((10 + step / 3).to_double())
      early_branch.push(
        if step < 4 {
          (4 + step).to_double()
        } else {
          (14 + step % 2).to_double()
        },
      )
    }
    zero.push(0.0)
  }
  [flat, association, endpoint, early_branch, zero]
}

///|
fn tradeseq_demo_trajectory() -> (Array[Array[Double]], Array[Array[Double]]) {
  let pseudotime : Array[Array[Double]] = []
  let weights : Array[Array[Double]] = []
  for cell in 0..<24 {
    if cell < 8 {
      let time = cell.to_double() * 0.05
      pseudotime.push([time, time])
      weights.push([0.5, 0.5])
    } else if cell < 16 {
      let time = 0.4 + (cell - 8).to_double() * 0.08
      pseudotime.push([time, 0.0])
      weights.push([1.0, 0.0])
    } else {
      let time = 0.4 + (cell - 16).to_double() * 0.08
      pseudotime.push([0.0, time])
      weights.push([0.0, 1.0])
    }
  }
  (pseudotime, weights)
}

///|
fn tradeseq_demo_config() -> @bio.TradeSeqAdvancedConfig {
  @bio.TradeSeqAdvancedConfig::create(
    n_knots=4,
    smoothing_penalty=1.0,
    max_iterations=80,
    tolerance=1.0e-5,
    ridge=1.0e-6,
    fdr_threshold=0.1,
    test_points=6,
  ) catch {
    _ => abort("failed to create tradeSeq configuration")
  }
}

///|
fn tradeseq_demo_print_test(table : @bio.TradeSeqAdvancedTestTable) -> Unit {
  println("\n" + table.test_name)
  for result in table.results {
    println(
      "  " +
      result.gene_id +
      ": Wald=" +
      result.wald_statistic.to_string() +
      ", df=" +
      result.degrees_freedom.to_string() +
      ", p=" +
      result.p_value.to_string() +
      ", FDR=" +
      result.adjusted_p_value.to_string() +
      ", max log2FC=" +
      result.log2_fold_change.to_string(),
    )
  }
}

///|
fn tradeseq_demo_slingshot() -> @bio.SlingshotAdvancedResult {
  let coordinates = [
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
  let labels = [
    "A", "A", "A", "B", "B", "B", "C", "C", "C", "D", "D", "D", "E", "E", "E",
  ]
  let config = @bio.SlingshotAdvancedConfig::create(
    start_clusters=["A"],
    end_clusters=["D", "E"],
    distance=@bio.slingshot_center_euclidean(),
    extension=@bio.slingshot_extend_none(),
    max_iterations=6,
    curve_points=24,
  ) catch {
    _ => abort("failed to create Slingshot configuration")
  }
  @bio.slingshot_advanced(coordinates, labels, config) catch {
    _ => abort("Slingshot fitting failed")
  }
}

///|
fn main {
  println("=== Bioconductor tradeSeq 1.27.0 Advanced Demo ===")
  let counts = tradeseq_demo_counts()
  let (pseudotime, weights) = tradeseq_demo_trajectory()
  let gene_names = ["flat", "association", "endpoint", "early", "zero"]
  let config = tradeseq_demo_config()
  let fit = @bio.tradeseq_fit_advanced(
    counts,
    pseudotime,
    weights,
    gene_names~,
    lineage_names=["left", "right"],
    offsets=Array::make(24, 0.0),
    config~,
  ) catch {
    _ => abort("tradeSeq NB-GAM fitting failed")
  }

  println("\n1. Multi-lineage negative-binomial GAM")
  println(@bio.tradeseq_advanced_summary(fit))
  for model in fit.models {
    println(
      "  " +
      model.gene_id +
      ": dispersion=" +
      model.dispersion.to_string() +
      ", AIC=" +
      model.aic.to_string() +
      ", converged=" +
      model.converged.to_string(),
    )
  }

  println("\n2. Trajectory Wald tests with BH-FDR")
  let association = @bio.tradeseq_association_test_advanced(fit) catch {
    _ => abort("associationTest failed")
  }
  let start_end = @bio.tradeseq_start_vs_end_test_advanced(fit) catch {
    _ => abort("startVsEndTest failed")
  }
  let diff_end = @bio.tradeseq_diff_end_test_advanced(fit) catch {
    _ => abort("diffEndTest failed")
  }
  let pattern = @bio.tradeseq_pattern_test_advanced(fit, n_points=5) catch {
    _ => abort("patternTest failed")
  }
  let early = @bio.tradeseq_early_de_test_advanced(fit, 0.35, 0.75, n_points=5) catch {
    _ => abort("earlyDETest failed")
  }
  tradeseq_demo_print_test(association)
  tradeseq_demo_print_test(start_end)
  tradeseq_demo_print_test(diff_end)
  tradeseq_demo_print_test(pattern)
  tradeseq_demo_print_test(early)

  let prediction = @bio.tradeseq_predict_smooth_advanced(
    fit,
    "endpoint",
    n_points=8,
  ) catch {
    _ => abort("smooth prediction failed")
  }
  println("\n3. Smooth prediction with delta-method uncertainty")
  for lineage in 0..<fit.n_lineages {
    println(
      "  " +
      fit.lineage_names[lineage] +
      ": start=" +
      prediction.fitted[lineage][0].to_string() +
      " +/- " +
      prediction.standard_errors[lineage][0].to_string() +
      ", end=" +
      prediction.fitted[lineage][7].to_string() +
      " +/- " +
      prediction.standard_errors[lineage][7].to_string(),
    )
  }

  let evaluation = @bio.tradeseq_evaluate_k_advanced(
    [counts[0], counts[1], counts[2]],
    pseudotime,
    weights,
    [3, 4],
    gene_names=["flat", "association", "endpoint"],
    lineage_names=["left", "right"],
    offsets=Array::make(24, 0.0),
    config~,
  ) catch {
    _ => abort("knot evaluation failed")
  }
  println("\n4. AIC knot evaluation")
  for index in 0..<evaluation.candidates.length() {
    println(
      "  knots=" +
      evaluation.candidates[index].to_string() +
      ", mean AIC=" +
      evaluation.mean_aic[index].to_string(),
    )
  }
  println("  selected knots=" + evaluation.selected_n_knots.to_string())

  let slingshot = tradeseq_demo_slingshot()
  let trend : Array[Double] = []
  for cell in 0..<15 {
    trend.push((cell + 2).to_double())
  }
  let slingshot_fit = @bio.tradeseq_fit_from_slingshot(
    [trend],
    ["slingshotTrend"],
    slingshot,
    offsets=Array::make(15, 0.0),
    config=@bio.TradeSeqAdvancedConfig::create(
      n_knots=3,
      max_iterations=50,
      ridge=1.0e-5,
      test_points=4,
    ),
  ) catch {
    _ => abort("Slingshot to tradeSeq integration failed")
  }
  println("\n5. Direct Slingshot integration")
  println(
    "  " +
    slingshot_fit.n_cells.to_string() +
    " cells, lineages=" +
    slingshot_fit.lineage_names.join(","),
  )

  let cell_names : Array[String] = []
  for cell in 0..<24 {
    cell_names.push("cell" + (cell + 1).to_string())
  }
  let experiment = @bio.SingleCellExperiment::new(
    counts, gene_names, cell_names,
  )
  experiment.reduced_dims["slingshot.pseudotime"] = pseudotime
  experiment.reduced_dims["slingshot.weights"] = weights
  let sce_output = @bio.tradeseq_advanced_sce(experiment, config~) catch {
    _ => abort("SingleCellExperiment tradeSeq integration failed")
  }
  println("\n6. Immutable SingleCellExperiment write-back")
  println(
    "  fitted assay: " +
    sce_output.experiment.assays["tradeSeq.fitted"].length().to_string() +
    " genes x " +
    sce_output.experiment.assays["tradeSeq.fitted"][0].length().to_string() +
    " cells",
  )
  println(
    "  original unchanged: " +
    (!experiment.assays.contains("tradeSeq.fitted")).to_string(),
  )
}
