///|
fn example_counts() -> Array[Array[Int]] {
  [
    [40, 3, 2],
    [35, 4, 1],
    [42, 2, 3],
    [38, 5, 2],
    [45, 3, 1],
    [36, 2, 4],
    [2, 40, 3],
    [4, 35, 2],
    [3, 43, 1],
    [5, 38, 2],
    [1, 44, 3],
    [3, 36, 4],
  ]
}

///|
fn example_groups() -> Array[String] {
  ["A", "A", "A", "A", "A", "A", "B", "B", "B", "B", "B", "B"]
}

///|
fn main {
  println("=== Bioconductor DirichletMultinomial Demo ===")
  let counts = example_counts()
  let groups = example_groups()
  let taxa = ["taxon_a", "taxon_b", "taxon_c"]
  let config = @src.DmnConfig::create(
    components=2,
    max_iterations=50,
    optimizer_iterations=120,
    soft_kmeans_iterations=60,
    seed=2026,
  ) catch {
    DirichletMultinomialError(message) =>
      abort("failed to create DirichletMultinomial configuration: " + message)
  }

  println("\n1. Fit a finite Dirichlet-multinomial mixture")
  let fit = @src.dirichlet_multinomial_fit(counts, taxon_names=taxa, config~) catch {
    DirichletMultinomialError(message) =>
      abort("failed to fit Dirichlet-multinomial mixture: " + message)
  }
  println("   " + fit.summary())
  println("   mixture weights=" + fit.mixture_weights().to_string())
  println("   component proportions=" + fit.component_proportions().to_string())
  println("   assignments=" + fit.assignments().to_string())

  println("\n2. Compare K=1 and K=2 with the upstream Laplace criterion")
  let selection = @src.dirichlet_multinomial_select(
    counts,
    2,
    criterion=@src.DmnLaplace,
    taxon_names=taxa,
    config=config.with_components(1),
  ) catch {
    DirichletMultinomialError(message) =>
      abort("failed to select Dirichlet-multinomial model: " + message)
  }
  println("   Laplace scores=" + selection.scores.to_string())
  println(
    "   selected components=" + selection.best_component_count().to_string(),
  )

  println("\n3. Train a dmngroup-style generative classifier")
  let classifier = @src.dirichlet_multinomial_group_fit(
    counts,
    groups,
    components_by_group=[1, 1],
    taxon_names=taxa,
    config=config.with_components(1),
  ) catch {
    DirichletMultinomialError(message) =>
      abort("failed to fit Dirichlet-multinomial classifier: " + message)
  }
  let novel = [[50, 1, 1], [1, 50, 1]]
  let probabilities = classifier.predict(novel) catch {
    DirichletMultinomialError(message) =>
      abort("failed to classify novel samples: " + message)
  }
  let predictions = classifier.predict_assignments(novel) catch {
    DirichletMultinomialError(message) =>
      abort("failed to assign novel samples: " + message)
  }
  println("   class probabilities=" + probabilities.to_string())
  println("   predictions=" + predictions.to_string())

  println("\n4. Run deterministic stratified cross-validation and ROC")
  let cross_validation = @src.dirichlet_multinomial_cross_validate(
    counts,
    groups,
    3,
    components_by_group=[1, 1],
    taxon_names=taxa,
    config=config.with_components(1),
  ) catch {
    DirichletMultinomialError(message) =>
      abort("failed to cross-validate classifier: " + message)
  }
  let truth : Array[Bool] = []
  let scores : Array[Double] = []
  for sample in 0..<groups.length() {
    truth.push(groups[sample] == "B")
    scores.push(cross_validation.probabilities[sample][1])
  }
  let roc = @src.dirichlet_multinomial_roc(truth, scores) catch {
    DirichletMultinomialError(message) =>
      abort("failed to compute ROC: " + message)
  }
  println(
    "   accuracy=" +
    cross_validation.accuracy.to_string() +
    ", AUC=" +
    roc.auc.to_string() +
    ", best threshold=" +
    roc.best_threshold.to_string(),
  )

  println("\n5. Fit a feature x sample SummarizedExperiment assay")
  let assays : Map[String, Array[Array[Double]]] = Map([
    (
      "counts",
      [
        [40.0, 35.0, 42.0, 2.0, 4.0, 3.0],
        [3.0, 4.0, 2.0, 40.0, 35.0, 43.0],
        [2.0, 1.0, 3.0, 3.0, 2.0, 1.0],
      ],
    ),
  ])
  let experiment = @src.summarized_experiment(assays, [], [], Map([]))
  let se_fit = @src.dirichlet_multinomial_fit_se(
    experiment,
    "counts",
    taxon_names=taxa,
    config~,
  ) catch {
    DirichletMultinomialError(message) =>
      abort("failed to fit SummarizedExperiment assay: " + message)
  }
  println(
    "   transposed dimensions=" +
    se_fit.sample_count().to_string() +
    " samples x " +
    se_fit.taxon_count().to_string() +
    " taxa",
  )

  println("\n6. Invalid count matrices produce explicit diagnostics")
  let rejected = try {
    ignore(@src.dirichlet_multinomial_fit([[3, -1], [2, 4]]))
    false
  } catch {
    DirichletMultinomialError(message) => {
      println("   " + message)
      true
    }
  }
  println("   malformed counts rejected=" + rejected.to_string())
  println("\n=== Demo Complete ===")
}
