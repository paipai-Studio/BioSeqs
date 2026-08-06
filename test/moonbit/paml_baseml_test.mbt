///|
fn baseml_test_close(
  actual : Double,
  expected : Double,
  tolerance : Double,
) -> Unit {
  if (actual - expected).abs() > tolerance {
    abort(
      "BASEML value " +
      actual.to_string() +
      " differs from " +
      expected.to_string(),
    )
  }
}

///|
fn baseml_test_expect_control_error(text : String) -> Unit {
  let raised = try {
    ignore(@src.baseml_parse_control(text))
    false
  } catch {
    @src.BasemlError(_) => true
  }
  if !raised {
    abort("expected BasemlError")
  }
}

///|
fn baseml_test_expect_result_error(text : String) -> Unit {
  let raised = try {
    ignore(@src.baseml_parse_results(text))
    false
  } catch {
    @src.BasemlError(_) => true
  }
  if !raised {
    abort("expected BasemlError")
  }
}

///|
fn baseml_test_minimal_result() -> String {
  "BASEML (in paml version 4.1, August 2008) alignment.phy K80\n" +
  "lnL(ntime: 1 np: 2): -12.5 +0.0\n" +
  "0.1 2.5\n" +
  "tree length = 0.25\n" +
  "(A:0.10,B:0.15);\n"
}

///|
fn baseml_test_kappa_result() -> String {
  baseml_test_minimal_result() +
  "Parameters (kappa) in the rate matrix (TN93):\n" +
  "3.25 0.75\n"
}

///|
fn baseml_test_branch_result() -> String {
  baseml_test_minimal_result() +
  "Parameters (kappa) in the rate matrix (F84):\n" +
  "1..2 0.10 2.50 3.50 1.25\n" +
  "2..3 0.20 4.50 5.50 2.25\n"
}

///|
fn baseml_test_auto_gamma_result() -> String {
  baseml_test_minimal_result() +
  "Parameters (kappa) in the rate matrix (K80):\n" +
  "2.5\n" +
  "alpha (gamma, K=2) = 0.25\n" +
  "rate: 0.5 1.5\n" +
  "freq: 0.4 0.6\n" +
  "rho for the auto-discrete-gamma model: 0.9\n" +
  "transition probabilities between rate categories:\n" +
  "0.8 0.2\n" +
  "0.1 0.9\n"
}

///|
fn baseml_test_node_result() -> String {
  baseml_test_minimal_result() +
  "(frequency parameters for branches) [frequencies at nodes]\n" +
  "Node #1 ( 0.10 0.20 0.30 0.40 )\n" +
  "Node #2 ( 0.11 0.21 0.31 0.37 0.25 0.25 0.30 0.20 )\n" +
  "Note: node 2 is root.\n"
}

///|
fn baseml_test_old_frequency_result() -> String {
  baseml_test_minimal_result() +
  "base frequency parameters\n" +
  "0.20 0.30 0.10 0.40\n"
}

///|
fn baseml_test_null_result() -> String {
  "BASEML (in paml version 4.7, January 2013) alignment.phy JC69\n" +
  "lnL(ntime: 1 np: 2): -105.0 +0.0\n" +
  "0.1 0.2\n" +
  "tree length = 0.2\n" +
  "(A:0.1,B:0.1);\n"
}

///|
fn baseml_test_alternative_result() -> String {
  "BASEML (in paml version 4.7, January 2013) alignment.phy HKY85\n" +
  "lnL(ntime: 1 np: 4): -100.0 +0.0\n" +
  "0.1 0.2 0.3 0.4\n" +
  "tree length = 0.2\n" +
  "(A:0.1,B:0.1);\n"
}

///|
test "baseml control parses required paths" {
  let control = @src.baseml_parse_control(@src.baseml_example_control_text())
  assert_eq(control.sequence_file(), "alignment.phylip")
  assert_eq(control.output_file(), "baseml.out")
  assert_eq(control.tree_file(), "species.tree")
}

///|
test "baseml control parses all official options" {
  let control = @src.baseml_parse_control(@src.baseml_example_control_text())
  assert_eq(control.options().length(), 21)
  assert_eq(control.option("runmode"), Some("0"))
  assert_eq(control.option("Small_Diff"), Some("0.000007"))
}

///|
test "baseml control exposes model number" {
  let control = @src.baseml_parse_control(@src.baseml_example_control_text())
  assert_eq(control.model_number(), Some(7))
  assert_eq(control.model_options(), "")
}

///|
test "baseml control parses model 9 options" {
  let control = @src.baseml_parse_control(
    "seqfile=a\noutfile=b\ntreefile=c\nmodel=9 [1 (TC CT AG GA)]\n",
  )
  assert_eq(control.model_number(), Some(9))
  assert_eq(control.model_options(), "[1 (TC CT AG GA)]")
}

///|
test "baseml control parses model 10 options" {
  let control = @src.baseml_parse_control(
    "seqfile=a\noutfile=b\ntreefile=c\nmodel=10 [5 (AC CA) (AG GA)]\n",
  )
  assert_eq(control.model_number(), Some(10))
  assert_eq(control.model_options(), "[5 (AC CA) (AG GA)]")
}

///|
test "baseml control strips comments and CRLF" {
  let control = @src.baseml_parse_control(
    "seqfile = a.phy * alignment\r\noutfile = out\r\ntreefile = t.nwk\r\nmodel = 6 * TN93\r\n",
  )
  assert_eq(control.sequence_file(), "a.phy")
  assert_eq(control.option("model"), Some("6"))
}

///|
test "baseml control canonical round trip" {
  let first = @src.baseml_parse_control(@src.baseml_example_control_text())
  let written = @src.baseml_write_control(first)
  let second = @src.baseml_parse_control(written)
  assert_eq(first, second)
  assert_true(written.has_prefix("seqfile = alignment.phylip\n"))
}

///|
test "baseml control preserves option order" {
  let control = @src.baseml_parse_control(
    "seqfile=a\noutfile=b\ntreefile=c\nalpha=0.5\nmodel=4\n",
  )
  let options = control.options()
  assert_eq(options[0].name(), "alpha")
  assert_eq(options[1].name(), "model")
}

///|
test "baseml control options are defensive copies" {
  let control = @src.baseml_parse_control(@src.baseml_example_control_text())
  let options = control.options()
  ignore(options.pop())
  assert_eq(control.options().length(), 21)
}

///|
test "baseml maps all model names" {
  let names : Array[String] = []
  for model in 0..<=10 {
    names.push(@src.baseml_model_name(model))
  }
  assert_eq(names, [
    "JC69", "K80", "F81", "F84", "HKY85", "T92", "TN93", "REV", "UNREST", "REVu",
    "UNRESTu",
  ])
}

///|
test "baseml control rejects missing seqfile" {
  baseml_test_expect_control_error("outfile=b\ntreefile=c\n")
}

///|
test "baseml control rejects missing outfile" {
  baseml_test_expect_control_error("seqfile=a\ntreefile=c\n")
}

///|
test "baseml control rejects missing treefile" {
  baseml_test_expect_control_error("seqfile=a\noutfile=b\n")
}

///|
test "baseml control rejects malformed lines" {
  baseml_test_expect_control_error(
    "seqfile=a\noutfile=b\ntreefile=c\nmodel 7\n",
  )
}

///|
test "baseml control rejects multiple equals signs" {
  baseml_test_expect_control_error(
    "seqfile=a\noutfile=b\ntreefile=c\nmodel=7=8\n",
  )
}

///|
test "baseml control rejects duplicate keys" {
  baseml_test_expect_control_error(
    "seqfile=a\noutfile=b\ntreefile=c\nmodel=7\nmodel=8\n",
  )
}

///|
test "baseml control rejects unknown options" {
  baseml_test_expect_control_error(
    "seqfile=a\noutfile=b\ntreefile=c\nunknown=1\n",
  )
}

///|
test "baseml control rejects empty values" {
  baseml_test_expect_control_error("seqfile=a\noutfile=b\ntreefile=c\nmodel=\n")
}

///|
test "baseml control rejects invalid integers" {
  baseml_test_expect_control_error(
    "seqfile=a\noutfile=b\ntreefile=c\nncatG=five\n",
  )
}

///|
test "baseml control rejects invalid doubles" {
  baseml_test_expect_control_error(
    "seqfile=a\noutfile=b\ntreefile=c\nalpha=fast\n",
  )
}

///|
test "baseml control rejects unsupported model" {
  baseml_test_expect_control_error(
    "seqfile=a\noutfile=b\ntreefile=c\nmodel=11\n",
  )
}

///|
test "baseml control rejects options on standard model" {
  baseml_test_expect_control_error(
    "seqfile=a\noutfile=b\ntreefile=c\nmodel=8 [custom]\n",
  )
}

///|
test "baseml control rejects non-positive ncatG" {
  baseml_test_expect_control_error(
    "seqfile=a\noutfile=b\ntreefile=c\nncatG=0\n",
  )
}

///|
test "baseml control rejects invalid nparK" {
  baseml_test_expect_control_error(
    "seqfile=a\noutfile=b\ntreefile=c\nnparK=5\n",
  )
}

///|
test "baseml control rejects invalid nhomo" {
  baseml_test_expect_control_error(
    "seqfile=a\noutfile=b\ntreefile=c\nnhomo=-1\n",
  )
}

///|
test "baseml control rejects invalid binary option" {
  baseml_test_expect_control_error(
    "seqfile=a\noutfile=b\ntreefile=c\ngetSE=2\n",
  )
}

///|
test "baseml model name rejects out of range value" {
  let raised = try {
    ignore(@src.baseml_model_name(-1))
    false
  } catch {
    @src.BasemlError(_) => true
  }
  assert_true(raised)
}

///|
test "baseml result parses version and model description" {
  let result = @src.baseml_parse_results(@src.baseml_example_results_text())
  assert_eq(result.version(), "4.7")
  assert_eq(result.model_description(), "REV dGamma (ncatG=5)")
}

///|
test "baseml result parses likelihoods" {
  let result = @src.baseml_parse_results(@src.baseml_example_results_text())
  baseml_test_close(result.ln_likelihood(), -319.774741, 1.0e-9)
  match result.ln_l_max() {
    Some(value) => baseml_test_close(value, -316.049385, 1.0e-9)
    None => abort("missing ln Lmax")
  }
}

///|
test "baseml result parses parameter vector" {
  let result = @src.baseml_parse_results(@src.baseml_example_results_text())
  assert_eq(result.parameter_count(), 13)
  assert_eq(result.parameter_list().length(), 13)
  baseml_test_close(result.parameter_list()[7], 998.99998, 1.0e-8)
}

///|
test "baseml result parses standard errors" {
  let result = @src.baseml_parse_results(@src.baseml_example_results_text())
  assert_eq(result.standard_errors().length(), 13)
  baseml_test_close(result.standard_errors()[12], 4.0, 1.0e-12)
}

///|
test "baseml result parses tree" {
  let result = @src.baseml_parse_results(@src.baseml_example_results_text())
  baseml_test_close(result.tree_length(), 0.01826, 1.0e-12)
  assert_true(result.tree().has_prefix("(((Homo_sapie:"))
}

///|
test "baseml result parses REV rate parameters" {
  let result = @src.baseml_parse_results(@src.baseml_example_results_text())
  assert_eq(result.rate_parameters().length(), 5)
  baseml_test_close(result.rate_parameters()[1], 130.94908, 1.0e-8)
}

///|
test "baseml result parses base frequencies" {
  let result = @src.baseml_parse_results(@src.baseml_example_results_text())
  match result.base_frequencies() {
    Some(frequencies) => {
      baseml_test_close(frequencies.thymine(), 0.20090, 1.0e-12)
      baseml_test_close(frequencies.cytosine(), 0.16306, 1.0e-12)
      baseml_test_close(frequencies.adenine(), 0.37027, 1.0e-12)
      baseml_test_close(frequencies.guanine(), 0.26577, 1.0e-12)
    }
    None => abort("missing base frequencies")
  }
}

///|
test "baseml result parses Q matrix" {
  let result = @src.baseml_parse_results(@src.baseml_example_results_text())
  match result.q_matrix() {
    Some(matrix) => {
      assert_eq(matrix.rows().length(), 4)
      assert_eq(matrix.rows()[0].length(), 4)
      baseml_test_close(matrix.rows()[2][3], 0.003044, 1.0e-12)
      match matrix.average_ts_tv() {
        Some(value) => baseml_test_close(value, 3.3698, 1.0e-12)
        None => abort("missing average Ts/Tv")
      }
    }
    None => abort("missing Q matrix")
  }
}

///|
test "baseml result parses gamma rates" {
  let result = @src.baseml_parse_results(@src.baseml_example_results_text())
  match result.alpha() {
    Some(value) => baseml_test_close(value, 200.95753, 1.0e-8)
    None => abort("missing alpha")
  }
  assert_eq(result.rates().length(), 5)
  assert_eq(result.rate_frequencies(), [0.2, 0.2, 0.2, 0.2, 0.2])
}

///|
test "baseml result parses scalar and multiple kappas" {
  let result = @src.baseml_parse_results(baseml_test_kappa_result())
  assert_eq(result.kappas(), [3.25, 0.75])
}

///|
test "baseml result parses branch-specific kappas" {
  let result = @src.baseml_parse_results(baseml_test_branch_result())
  let branches = result.branch_parameters()
  assert_eq(branches.length(), 2)
  assert_eq(branches[0].branch(), "1..2")
  baseml_test_close(branches[0].time(), 0.10, 1.0e-12)
  baseml_test_close(branches[0].kappa(), 2.50, 1.0e-12)
  baseml_test_close(branches[0].transitions(), 3.50, 1.0e-12)
  baseml_test_close(branches[0].transversions(), 1.25, 1.0e-12)
}

///|
test "baseml result parses auto discrete gamma" {
  let result = @src.baseml_parse_results(baseml_test_auto_gamma_result())
  match result.rho() {
    Some(value) => baseml_test_close(value, 0.9, 1.0e-12)
    None => abort("missing rho")
  }
  assert_eq(result.rates(), [0.5, 1.5])
  assert_eq(result.rate_frequencies(), [0.4, 0.6])
  assert_eq(result.transition_probabilities().length(), 2)
  baseml_test_close(result.transition_probabilities()[1][1], 0.9, 1.0e-12)
}

///|
test "baseml result parses nonhomogeneous nodes" {
  let result = @src.baseml_parse_results(baseml_test_node_result())
  let nodes = result.nodes()
  assert_eq(nodes.length(), 2)
  assert_eq(nodes[0].node(), 1)
  assert_false(nodes[0].is_root())
  assert_true(nodes[1].is_root())
  assert_eq(nodes[0].frequency_parameters(), [0.1, 0.2, 0.3, 0.4])
}

///|
test "baseml result parses realized node base frequencies" {
  let result = @src.baseml_parse_results(baseml_test_node_result())
  match result.nodes()[1].base_frequencies() {
    Some(frequencies) => {
      baseml_test_close(frequencies.thymine(), 0.25, 1.0e-12)
      baseml_test_close(frequencies.guanine(), 0.20, 1.0e-12)
    }
    None => abort("missing node base frequencies")
  }
}

///|
test "baseml result parses PAML 4.1 base frequency heading" {
  let result = @src.baseml_parse_results(baseml_test_old_frequency_result())
  match result.base_frequencies() {
    Some(frequencies) =>
      baseml_test_close(frequencies.cytosine(), 0.30, 1.0e-12)
    None => abort("missing old-style base frequencies")
  }
}

///|
test "baseml result arrays are defensive copies" {
  let result = @src.baseml_parse_results(@src.baseml_example_results_text())
  let parameters = result.parameter_list()
  parameters[0] = 999.0
  let rates = result.rates()
  ignore(rates.pop())
  baseml_test_close(result.parameter_list()[0], 0.000004, 1.0e-12)
  assert_eq(result.rates().length(), 5)
}

///|
test "baseml Q matrix rows are defensive copies" {
  let result = @src.baseml_parse_results(@src.baseml_example_results_text())
  match result.q_matrix() {
    Some(matrix) => {
      let rows = matrix.rows()
      rows[0][0] = 100.0
      match result.q_matrix() {
        Some(again) => baseml_test_close(again.rows()[0][0], -2.483179, 1.0e-12)
        None => abort("missing copied Q matrix")
      }
    }
    None => abort("missing Q matrix")
  }
}

///|
test "baseml node arrays are defensive copies" {
  let result = @src.baseml_parse_results(baseml_test_node_result())
  let nodes = result.nodes()
  let parameters = nodes[0].frequency_parameters()
  parameters[0] = 9.0
  baseml_test_close(result.nodes()[0].frequency_parameters()[0], 0.1, 1.0e-12)
}

///|
test "baseml result rejects empty input" {
  baseml_test_expect_result_error("")
}

///|
test "baseml result rejects missing header" {
  baseml_test_expect_result_error(
    "lnL(ntime: 1 np: 1): -1\n0.1\ntree length = 0.1\n(A:0.1);\n",
  )
}

///|
test "baseml result rejects missing likelihood" {
  baseml_test_expect_result_error(
    "BASEML (in paml version 4.7, January 2013) a.phy JC69\n",
  )
}

///|
test "baseml result rejects malformed likelihood" {
  baseml_test_expect_result_error(
    "BASEML (in paml version 4.7, January 2013) a.phy JC69\n" +
    "lnL(ntime: 1 np: x): bad\n",
  )
}

///|
test "baseml result rejects parameter vector mismatch" {
  baseml_test_expect_result_error(
    "BASEML (in paml version 4.7, January 2013) a.phy JC69\n" +
    "lnL(ntime: 1 np: 2): -1\n" +
    "0.1\n" +
    "tree length = 0.1\n" +
    "(A:0.1);\n",
  )
}

///|
test "baseml result rejects standard error mismatch" {
  baseml_test_expect_result_error(
    baseml_test_minimal_result() + "SEs for parameters:\n0.1\n",
  )
}

///|
test "baseml result rejects missing tree length" {
  baseml_test_expect_result_error(
    "BASEML (in paml version 4.7, January 2013) a.phy JC69\n" +
    "lnL(ntime: 1 np: 1): -1\n" +
    "0.1\n" +
    "(A:0.1);\n",
  )
}

///|
test "baseml result rejects missing branch-length tree" {
  baseml_test_expect_result_error(
    "BASEML (in paml version 4.7, January 2013) a.phy JC69\n" +
    "lnL(ntime: 1 np: 1): -1\n" +
    "0.1\n" +
    "tree length = 0.1\n" +
    "(A,B);\n",
  )
}

///|
test "baseml result rejects incomplete Q matrix" {
  baseml_test_expect_result_error(
    baseml_test_minimal_result() +
    "Rate matrix Q, Average Ts/Tv = 2.0\n" +
    "1.0 0.0 0.0 0.0\n",
  )
}

///|
test "baseml result rejects transition matrix row count mismatch" {
  baseml_test_expect_result_error(
    baseml_test_minimal_result() +
    "rate: 0.5 1.5\n" +
    "transition probabilities between rate categories:\n" +
    "0.8 0.2\n",
  )
}

///|
test "baseml result rejects transition matrix column mismatch" {
  baseml_test_expect_result_error(
    baseml_test_minimal_result() +
    "rate: 0.5 1.5\n" +
    "transition probabilities between rate categories:\n" +
    "0.8 0.2\n" +
    "1.0\n",
  )
}

///|
test "baseml result rejects duplicate nodes" {
  baseml_test_expect_result_error(
    baseml_test_minimal_result() +
    "Node #1 ( 0.1 0.2 0.3 0.4 )\n" +
    "Node #1 ( 0.2 0.2 0.2 0.4 )\n",
  )
}

///|
test "baseml computes AIC" {
  let result = @src.baseml_parse_results(baseml_test_null_result())
  baseml_test_close(result.aic(), 214.0, 1.0e-12)
}

///|
test "baseml computes BIC" {
  let result = @src.baseml_parse_results(baseml_test_null_result())
  baseml_test_close(result.bic(100), 2.0 * @math.ln(100.0) + 210.0, 1.0e-12)
}

///|
test "baseml BIC rejects non-positive observations" {
  let result = @src.baseml_parse_results(baseml_test_null_result())
  let raised = try {
    ignore(result.bic(0))
    false
  } catch {
    @src.BasemlError(_) => true
  }
  assert_true(raised)
}

///|
test "baseml computes nested likelihood ratio" {
  let null_model = @src.baseml_parse_results(baseml_test_null_result())
  let alternative = @src.baseml_parse_results(baseml_test_alternative_result())
  let comparison = @src.baseml_likelihood_ratio(null_model, alternative)
  baseml_test_close(comparison.statistic(), 10.0, 1.0e-12)
  assert_eq(comparison.degrees_of_freedom(), 2)
  baseml_test_close(comparison.p_value(), 0.006737946999, 1.0e-9)
  assert_eq(comparison.alpha(), 0.05)
  assert_true(comparison.significant())
}

///|
test "baseml LRT rejects reversed nesting" {
  let null_model = @src.baseml_parse_results(baseml_test_null_result())
  let alternative = @src.baseml_parse_results(baseml_test_alternative_result())
  let raised = try {
    ignore(@src.baseml_likelihood_ratio(alternative, null_model))
    false
  } catch {
    @src.BasemlError(_) => true
  }
  assert_true(raised)
}

///|
test "baseml LRT rejects lower alternative likelihood" {
  let null_model = @src.baseml_parse_results(baseml_test_minimal_result())
  let alternative = @src.baseml_parse_results(
    "BASEML (in paml version 4.7, January 2013) a.phy HKY85\n" +
    "lnL(ntime: 1 np: 3): -20.0\n" +
    "0.1 0.2 0.3\n" +
    "tree length = 0.2\n" +
    "(A:0.1,B:0.1);\n",
  )
  let raised = try {
    ignore(@src.baseml_likelihood_ratio(null_model, alternative))
    false
  } catch {
    @src.BasemlError(_) => true
  }
  assert_true(raised)
}

///|
test "baseml LRT rejects invalid alpha" {
  let null_model = @src.baseml_parse_results(baseml_test_null_result())
  let alternative = @src.baseml_parse_results(baseml_test_alternative_result())
  let raised = try {
    ignore(@src.baseml_likelihood_ratio(null_model, alternative, alpha=1.0))
    false
  } catch {
    @src.BasemlError(_) => true
  }
  assert_true(raised)
}
