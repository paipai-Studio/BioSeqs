///|
fn baseml_demo_control() -> @bio.BasemlControl {
  @bio.baseml_parse_control(@bio.baseml_example_control_text()) catch {
    BasemlError(message) => abort("failed to parse BASEML control: " + message)
  }
}

///|
fn baseml_demo_results() -> @bio.BasemlResults {
  @bio.baseml_parse_results(@bio.baseml_example_results_text()) catch {
    BasemlError(message) => abort("failed to parse BASEML results: " + message)
  }
}

///|
fn baseml_demo_null_model() -> @bio.BasemlResults {
  @bio.baseml_parse_results(
    "BASEML (in paml version 4.7, January 2013) alignment.phylip JC69 dGamma\n" +
    "lnL(ntime: 7 np: 8): -330.000000 +0.000000\n" +
    "0.01 0.01 0.01 0.01 0.01 0.01 0.01 1.0\n" +
    "tree length = 0.018\n" +
    "((A:0.004,B:0.004):0.002,C:0.008);\n",
  ) catch {
    BasemlError(message) => abort("failed to parse null model: " + message)
  }
}

///|
fn main {
  println("=== Biopython Bio.Phylo.PAML.baseml Offline Demo ===")
  let control = baseml_demo_control()
  let model = match control.model_number() {
    Some(value) => value
    None => abort("control has no model")
  }
  println("\n1. Strict BASEML control parsing")
  println("  sequence file: " + control.sequence_file())
  println("  tree file: " + control.tree_file())
  println(
    "  model: " +
    model.to_string() +
    " (" +
    (@bio.baseml_model_name(model) catch {
      BasemlError(message) => abort(message)
    }) +
    ")",
  )
  let canonical = @bio.baseml_write_control(control)
  let reparsed = @bio.baseml_parse_control(canonical) catch {
    BasemlError(message) => abort("control round-trip failed: " + message)
  }
  println("  canonical round-trip: " + (reparsed == control).to_string())

  let result = baseml_demo_results()
  println("\n2. PAML result, parameters, and standard errors")
  println(
    "  PAML=" +
    result.version() +
    ", model=" +
    result.model_description() +
    ", lnL=" +
    result.ln_likelihood().to_string(),
  )
  println(
    "  parameters=" +
    result.parameter_count().to_string() +
    ", SE values=" +
    result.standard_errors().length().to_string() +
    ", tree length=" +
    result.tree_length().to_string(),
  )

  println("\n3. REV rate matrix and discrete gamma")
  match result.q_matrix() {
    Some(matrix) => {
      let average = match matrix.average_ts_tv() {
        Some(value) => value.to_string()
        None => "not reported"
      }
      println(
        "  Q dimensions=" +
        matrix.rows().length().to_string() +
        "x" +
        matrix.rows()[0].length().to_string() +
        ", average Ts/Tv=" +
        average,
      )
    }
    None => abort("result has no Q matrix")
  }
  let alpha_text = match result.alpha() {
    Some(value) => value.to_string()
    None => "not reported"
  }
  println(
    "  rate categories=" +
    result.rates().length().to_string() +
    ", alpha=" +
    alpha_text,
  )

  let comparison = @bio.baseml_likelihood_ratio(
    baseml_demo_null_model(),
    result,
  ) catch {
    BasemlError(message) => abort("likelihood-ratio test failed: " + message)
  }
  let bic = result.bic(222) catch {
    BasemlError(message) => abort("BIC calculation failed: " + message)
  }
  println("\n4. Information criteria and nested-model LRT")
  println("  REV AIC=" + result.aic().to_string() + ", BIC=" + bic.to_string())
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
    "\nThe demo parses existing BASEML artifacts; it does not launch PAML.",
  )
}
