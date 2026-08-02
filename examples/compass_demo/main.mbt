///|
/// COMPASS module demo - profile-profile alignment output parser.
///
/// COMPASS (Comparison of Multiple Protein Alignments with Assessment of
/// Statistical Significance) is a tool for profile-profile alignment.
///
/// Run: moon run examples/compass_demo
fn main {
  println("=== COMPASS Module Demo ===")
  println("Profile-profile alignment output parser")
  println(
    "Reference: Sadreyev RI et al. (2003) Protein Sci 12(10):2262-72",
  )
  println("")

  // ------------------------------------------------------------------
  // 1. Generate and parse sample COMPASS output
  // ------------------------------------------------------------------
  println("=== 1. Generate and parse sample COMPASS output ===")
  let text = @bio.compass_sample_data()
  let records = @bio.parse_compass(text)
  println("Parsed " + records.length().to_string() + " COMPASS records")
  println("")

  // ------------------------------------------------------------------
  // 2. Display version information
  // ------------------------------------------------------------------
  println("=== 2. Version information ===")
  let version = @bio.compass_version(text)
  println("COMPASS version: " + version)
  println("")

  // ------------------------------------------------------------------
  // 3. Display all parsed records
  // ------------------------------------------------------------------
  println("=== 3. Parsed records ===")
  for r in records {
    println("Query: " + r.query_name())
    println("  Template: " + r.template_name())
    println(
      "  Sequences: query=" +
      r.query_n_seqs().to_string() +
      " template=" +
      r.template_n_seqs().to_string(),
    )
    println(
      "  Columns: query=" +
      r.query_n_cols().to_string() +
      " template=" +
      r.template_n_cols().to_string(),
    )
    println("  Smith-Waterman score: " + r.sw_score().to_string())
    println("  E-value: " + r.e_value().to_string())
    println(
      "  Percentage identity: " + r.percentage_identity().to_string() + "%",
    )
    println(
      "  Query range: " +
      r.query_start().to_string() +
      "-" +
      r.query_end().to_string(),
    )
    println(
      "  Template range: " +
      r.template_start().to_string() +
      "-" +
      r.template_end().to_string(),
    )
    println("  Aligned query: " + r.aligned_query())
    println("  Consensus:     " + r.consensus_line())
    println("  Aligned template: " + r.aligned_template())
    println(
      "  Alignment length: " +
      @bio.compass_alignment_length(r).to_string(),
    )
    println("")
  }

  // ------------------------------------------------------------------
  // 4. Filter records by E-value
  // ------------------------------------------------------------------
  println("=== 4. Filter by E-value <= 1e-06 ===")
  let significant = @bio.compass_filter_by_evalue(records, 1.0e-06)
  println(
    "Found " +
    significant.length().to_string() +
    " significant alignment(s):",
  )
  for r in significant {
    println(
      "  " +
      r.query_name() +
      " vs " +
      r.template_name() +
      " (E=" +
      r.e_value().to_string() +
      ")",
    )
  }
  println("")

  // ------------------------------------------------------------------
  // 5. Filter records by percentage identity
  // ------------------------------------------------------------------
  println("=== 5. Filter by identity >= 20.0% ===")
  let high_identity = @bio.compass_filter_by_identity(records, 20.0)
  println(
    "Found " +
    high_identity.length().to_string() +
    " alignment(s) with identity >= 20%:",
  )
  for r in high_identity {
    println(
      "  " +
      r.query_name() +
      " vs " +
      r.template_name() +
      " (identity=" +
      r.percentage_identity().to_string() +
      "%)",
    )
  }
  println("")

  // ------------------------------------------------------------------
  // 6. Generate summary
  // ------------------------------------------------------------------
  println("=== 6. Full summary ===")
  println(@bio.compass_summary(records))

  println("=== Demo complete ===")
}
