///|
fn codeml_demo_control() -> @bio.CodemlControl {
  @bio.codeml_parse_control(@bio.codeml_example_control_text()) catch {
    CodemlError(message) => abort("failed to parse CODEML control: " + message)
  }
}

///|
fn codeml_demo_results() -> @bio.CodemlResults {
  @bio.codeml_parse_results(@bio.codeml_example_results_text()) catch {
    CodemlError(message) => abort("failed to parse CODEML results: " + message)
  }
}

///|
fn codeml_demo_ints(values : Array[Int]) -> String {
  "[" + values.map(fn(value) { value.to_string() }).join(", ") + "]"
}

///|
fn main {
  println("=== Biopython Bio.Phylo.PAML.codeml Offline Demo ===")
  let control = codeml_demo_control()
  println("\n1. Strict CODEML control parsing")
  println("  sequence file: " + control.sequence_file())
  println("  tree file: " + control.tree_file())
  println("  NSsites: " + codeml_demo_ints(control.ns_sites()))
  let written = @bio.codeml_write_control(control)
  let reparsed = @bio.codeml_parse_control(written) catch {
    CodemlError(message) => abort("control round-trip failed: " + message)
  }
  println(
    "  canonical round-trip: " +
    (reparsed.ns_sites() == control.ns_sites()).to_string(),
  )

  let results = codeml_demo_results()
  let models = results.models()
  let null_model = models[0]
  let selection_model = models[1]
  println("\n2. CODONML metadata and multiple NSsites models")
  println(
    "  program=" +
    results.program() +
    ", PAML=" +
    results.version() +
    ", sequences=" +
    results.sequence_count().to_string() +
    ", sites=" +
    results.site_count().to_string(),
  )
  println(
    "  M0 lnL=" +
    null_model.ln_likelihood().to_string() +
    ", parameters=" +
    null_model.parameter_count().to_string(),
  )
  println(
    "  M2 lnL=" +
    selection_model.ln_likelihood().to_string() +
    ", parameters=" +
    selection_model.parameter_count().to_string(),
  )

  let branches = null_model.branches()
  let site_classes = selection_model.site_classes()
  let positive_sites = selection_model.positive_sites()
  println("\n3. Branch, site-class, and BEB estimates")
  println(
    "  branches=" +
    branches.length().to_string() +
    ", first branch=" +
    branches[0].branch(),
  )
  println(
    "  site classes=" +
    site_classes.length().to_string() +
    ", selected class proportion=" +
    site_classes[2].proportion().to_string(),
  )
  println(
    "  BEB site " +
    positive_sites[0].position().to_string() +
    positive_sites[0].amino_acid() +
    ", probability=" +
    positive_sites[0].probability().to_string() +
    positive_sites[0].significance(),
  )

  let bic = selection_model.bic(results.site_count()) catch {
    CodemlError(message) => abort("BIC calculation failed: " + message)
  }
  let comparison = @bio.codeml_likelihood_ratio(null_model, selection_model) catch {
    CodemlError(message) => abort("likelihood-ratio test failed: " + message)
  }
  println("\n4. Information criteria and nested-model LRT")
  println(
    "  M2 AIC=" + selection_model.aic().to_string() + ", BIC=" + bic.to_string(),
  )
  println(
    "  statistic=" +
    comparison.statistic().to_string() +
    ", df=" +
    comparison.degrees_of_freedom().to_string() +
    ", p=" +
    comparison.p_value().to_string() +
    ", significant=" +
    comparison.significant().to_string(),
  )

  println(
    "\nThe demo parses existing CODEML artifacts; it does not launch PAML.",
  )
}
