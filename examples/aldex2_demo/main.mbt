///|
fn print_feature(result : @src.Aldex2Result, name : String) -> Unit {
  let effect = match result.effects.find_feature(name) {
    Some(value) => value
    None => abort("ALDEx2 example feature not found: " + name)
  }
  let tests = match result.tests.find_feature(name) {
    Some(value) => value
    None => abort("ALDEx2 example test result not found: " + name)
  }
  println(
    "   " +
    name +
    ": effect=" +
    effect.effect.to_string() +
    ", overlap=" +
    effect.overlap.to_string() +
    ", Welch eBH=" +
    tests.welch_adjusted.to_string() +
    ", Wilcoxon eBH=" +
    tests.wilcoxon_adjusted.to_string(),
  )
}

///|
fn main {
  println("=== Bioconductor ALDEx2 Demo ===")
  let (counts, conditions, features, samples) = @src.aldex2_example_data()
  let config = @src.Aldex2Config::create(
    mc_samples=64,
    denominator=@src.Aldex2Iqlr,
    seed=2026,
  ) catch {
    _ => abort("failed to create ALDEx2 configuration")
  }

  println("\n1. Draw deterministic Dirichlet instances and apply IQLR")
  let clr = @src.aldex2_clr(
    counts,
    conditions,
    config~,
    feature_names=features,
    sample_names=samples,
  ) catch {
    Aldex2Error(message) => abort("failed to construct ALDEx2 CLR: " + message)
  }
  println(
    "   retained features=" +
    clr.feature_count().to_string() +
    ", removed all-zero features=" +
    clr.removed_feature_count().to_string() +
    ", samples=" +
    clr.sample_count().to_string() +
    ", Monte Carlo instances=" +
    clr.mc_sample_count().to_string(),
  )

  println("\n2. Estimate posterior tests and standardized effects")
  let result = @src.aldex2(
    counts,
    conditions,
    config~,
    feature_names=features,
    sample_names=samples,
  ) catch {
    Aldex2Error(message) => abort("failed to run ALDEx2: " + message)
  }
  print_feature(result, "increased")
  print_feature(result, "decreased")
  print_feature(result, "stable_high")
  print_feature(result, "rare_increased")

  println("\n3. Rank, select, and summarize compositional effects")
  let ranked = result.ranked()
  let summary = result.summary(maximum_adjusted=0.1, minimum_effect=1.0)
  println(
    "   largest absolute effect=" +
    ranked[0].feature_name +
    " (" +
    ranked[0].effect.to_string() +
    ")",
  )
  println(
    "   Welch eBH <= 0.10=" +
    summary.significant_welch_count.to_string() +
    ", |effect| >= 1=" +
    summary.large_effect_count.to_string() +
    ", jointly selected=" +
    result.select(maximum_adjusted=0.1, minimum_effect=1.0).length().to_string(),
  )

  println("\n4. Compute posterior expected Aitchison distances")
  let distances = result.clr.expected_distance()
  println(
    "   distance(C1, T1)=" +
    distances[0][4].to_string() +
    ", symmetric=" +
    (distances[0][4] == distances[4][0]).to_string(),
  )

  println("\n5. Add posterior summaries to an immutable container copy")
  let assays : Map[String, Array[Array[Double]]] = Map([("counts", counts)])
  let metadata : Map[String, String] = Map([("source", "aldex2_demo")])
  let experiment = @src.summarized_experiment(assays, [], [], metadata)
  let output = @src.aldex2_summarized_experiment(
    experiment,
    conditions,
    config~,
    feature_names=features,
    sample_names=samples,
  ) catch {
    Aldex2Error(message) =>
      abort("failed to integrate ALDEx2 with SummarizedExperiment: " + message)
  }
  let effect_assay = match @src.se_assay(output.experiment, "aldex2_effect") {
    Some(value) => value
    None => abort("ALDEx2 effect assay missing")
  }
  println(
    "   effect rows=" +
    effect_assay.length().to_string() +
    ", source object unchanged=" +
    (@src.se_assay(experiment, "aldex2_effect") is None).to_string(),
  )

  println("\n6. Export an ALDEx2-compatible tabular summary")
  println(result.to_tsv())

  println("7. Invalid count matrices produce explicit diagnostics")
  let rejected = try {
    ignore(
      @src.aldex2([[2.0, -1.0, 3.0, 4.0], [5.0, 6.0, 7.0, 8.0]], [
        "A", "A", "B", "B",
      ]),
    )
    false
  } catch {
    Aldex2Error(message) => {
      println("   " + message)
      true
    }
  }
  println("   malformed counts rejected=" + rejected.to_string())
  println("\n=== Demo Complete ===")
}
