///|
fn main {
  let data = @bio.gsva_create_example_data()
  let gene_sets = @bio.gsva_create_example_gene_sets()

  let params = @bio.GSVAParams::new()
  let scores = @bio.gsva_run(data, gene_sets, params)

  println(
    "GSVA Scores computed for " +
    scores.scores.length().to_string() +
    " gene sets",
  )
  println("Number of samples: " + scores.sample_names.length().to_string())

  let summary = @bio.gsva_summarize(scores)
  println(summary)

  println("\n=== 1. Enrichment Map ===")
  let map = @bio.gsva_enrichment_map(scores, gene_sets, 0.1)
  println(map)

  println("\n=== 2. Phenotype Correlation ===")
  let phenotypes = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0]
  let corr_report = @bio.gsva_phenotype_report(scores, phenotypes, "Age")
  println(corr_report)

  println("\n=== 3. Survival Analysis ===")
  let survival_time = [10.0, 20.0, 30.0, 40.0, 50.0, 60.0]
  let event = [1, 1, 1, 0, 0, 1]
  let surv_report = @bio.gsva_survival_report(scores, survival_time, event)
  println(surv_report)

  println("\n=== 4. Score Distribution ===")
  let dist = @bio.gsva_score_distribution(scores)
  println(dist)
}
