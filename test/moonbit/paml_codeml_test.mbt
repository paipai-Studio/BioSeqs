///|
fn codeml_test_close(
  actual : Double,
  expected : Double,
  tolerance : Double,
) -> Unit {
  if (actual - expected).abs() > tolerance {
    abort(
      "CODEML value " +
      actual.to_string() +
      " differs from " +
      expected.to_string(),
    )
  }
}

///|
fn codeml_test_expect_control_error(text : String) -> Unit {
  let raised = try {
    ignore(@src.codeml_parse_control(text))
    false
  } catch {
    @src.CodemlError(_) => true
  }
  if !raised {
    abort("expected CodemlError")
  }
}

///|
fn codeml_test_expect_result_error(text : String) -> Unit {
  let raised = try {
    ignore(@src.codeml_parse_results(text))
    false
  } catch {
    @src.CodemlError(_) => true
  }
  if !raised {
    abort("expected CodemlError")
  }
}

///|
fn codeml_test_branch_site_text() -> String {
  "CODONML (in paml version 4.7, January 2013) alignment.phylip\n" +
  "Model: branch-site model A\n" +
  "Codon frequency model: F3x4\n" +
  "Site-class models: PositiveSelection\n" +
  "ns = 5 ls = 74\n" +
  "lnL(ntime: 7 np: 12): -308.031579 +0.000000\n" +
  "0.1 0.2 0.3 0.4 0.5 0.6 0.7 1.5 0.9 0.1 0.2 2.0\n" +
  "tree length = 0.05472\n" +
  "((A:0.01, B:0.02):0.03, C:0.04);\n" +
  "proportion       0.70000  0.10000  0.15000  0.05000\n" +
  "background w     0.10000  1.00000  0.10000  1.00000\n" +
  "foreground w     0.10000  1.00000  4.50000  4.50000\n" +
  "Bayes Empirical Bayes (BEB) analysis\n" +
  " 21 R 0.981* 3.750 +- 0.820\n" +
  " 45 K 0.999** 8.250 +- 1.400\n"
}

///|
fn codeml_test_clade_text() -> String {
  "CODONML (in paml version 4.9j, February 2020) alignment.phylip\n" +
  "Model: Clade model C\n" +
  "Codon frequency model: F3x4\n" +
  "Site-class models: PositiveSelection\n" +
  "ns = 5 ls = 74\n" +
  "lnL(ntime: 7 np: 13): -307.500000 +0.000000\n" +
  "0.1 0.2 0.3 0.4 0.5 0.6 0.7 1.5 0.9 0.1 0.2 2.0 3.0\n" +
  "proportion        0.60000   0.30000   0.10000\n" +
  "branch type 0:    0.10000   1.00000   2.50000\n" +
  "branch type 1:    0.20000   1.00000   4.00000\n"
}

///|
fn codeml_test_free_ratio_text() -> String {
  "CODONML (in paml version 4.7, January 2013) alignment.phylip\n" +
  "Model: free ratios for branches\n" +
  "Codon frequency model: F3x4\n" +
  "ns = 3 ls = 90\n" +
  "lnL(ntime: 4 np: 8): -120.500000 +0.000000\n" +
  "0.1 0.2 0.3 0.4 1.1 0.8 0.4 0.2\n" +
  "tree length = 0.12000\n" +
  "((A:0.01, B:0.02):0.03, C:0.04);\n" +
  "kappa (ts/tv) = 2.10000\n" +
  "w (dN/dS) for branches: 1.10000 0.80000 0.40000 0.20000\n" +
  " 4..5 0.010 100.0 50.0 1.1000 0.0050 0.0045 0.5 0.2\n" +
  " 5..1 0.020 100.0 50.0 0.8000 0.0040 0.0050 0.4 0.3\n" +
  "dS tree:\n" +
  "((A:0.004, B:0.005):0.006, C:0.007);\n" +
  "dN tree:\n" +
  "((A:0.005, B:0.004):0.003, C:0.002);\n" +
  "w ratios as labels for TreeView:\n" +
  "((A #1.1, B #0.8) #0.4, C #0.2);\n"
}

///|
fn codeml_test_pairwise_text() -> String {
  "CODONML (in paml version 4.7, January 2013) alignment.phylip\n" +
  "Model: One dN/dS ratio for branches,\n" +
  "Codon frequency model: F3x4\n" +
  "ns = 3 ls = 74\n" +
  "2 (Pan_troglo) ... 1 (Homo_sapie)\n" +
  "lnL = -291.465693\n" +
  "0.01262 999.00000 0.00100\n" +
  "t= 0.0126 S=81.4 N=140.6 dN/dS=0.0010 dN=0.0000 dS=0.0115\n" +
  "3 (Gorilla_go) ... 1 (Homo_sapie)\n" +
  "lnL = -290.441129\n" +
  "0.01265 999.00000 0.00100\n" +
  "t=0.0127 S=81.7 N=140.3 dN/dS=0.0010 dN=0.0000 dS=0.0114\n" +
  "3 (Gorilla_go) ... 2 (Pan_troglo)\n" +
  "lnL = -296.416525\n" +
  "0.02582 999.00000 0.00100\n" +
  "t=0.0258 S=81.5 N=140.5 dN/dS=0.0010 dN=0.0000 dS=0.0234\n"
}

///|
fn codeml_test_aaml_text() -> String {
  "AAML (in paml version 4.7, January 2013) aa_alignment.phylip\n" +
  "Model: Poisson for branches, ns = 3 ls = 74\n" +
  "ln Lmax (unconstrained) = -192.057157\n" +
  "AA distances (raw proportions of different sites)\n" +
  "Human\n" +
  "Chimp 0.0100\n" +
  "Gorilla 0.0200 0.0300\n" +
  "\n" +
  "ML distances of aa seqs.\n" +
  "Human\n" +
  "Chimp 0.0110\n" +
  "Gorilla 0.0220 0.0330\n" +
  "\n"
}

///|
fn codeml_test_se_text() -> String {
  "CODONML (in paml version 4.6, August 2010) alignment.phylip\n" +
  "Model: one-ratio\n" +
  "Codon frequency model: F3x4\n" +
  "ns = 3 ls = 60\n" +
  "lnL(ntime: 4 np: 6): -100.000000 +0.000000\n" +
  "0.1 0.2 0.3 0.4 2.0 0.5\n" +
  "SEs for parameters:\n" +
  "0.01 0.02 0.03 0.04 0.20 0.05\n"
}

///|
fn codeml_test_multigene_text() -> String {
  "CODONML (in paml version 4.7, January 2013) genes.phylip\n" +
  "Model: One dN/dS ratio for branches, (2 genes: separate data)\n" +
  "Codon frequency model: F3x4\n" +
  "Site-class models: one-ratio\n" +
  "ns = 4 ls = 120\n" +
  "Gene 1\n" +
  "lnL(ntime: 3 np: 5): -80.000000 +0.000000\n" +
  "0.1 0.2 0.3 2.0 0.4\n" +
  "tree length = 0.10\n" +
  "Gene 2\n" +
  "lnL(ntime: 3 np: 5): -90.000000 +0.000000\n" +
  "0.1 0.2 0.3 2.1 0.5\n" +
  "tree length = 0.20\n"
}

///|
fn codeml_test_joint_gene_text() -> String {
  "CODONML (in paml version 4.7, January 2013) genes.phylip\n" +
  "Model: One dN/dS ratio for branches, (2 genes: joint data)\n" +
  "Codon frequency model: F3x4\n" +
  "ns = 4 ls = 120\n" +
  "lnL(ntime: 7 np: 10): -170.000000 +0.000000\n" +
  "0.1 0.2 0.3 0.4 0.5 0.6 2.0 0.3 1.0 2.5\n" +
  "rates for 2 genes: 1 2.50000\n" +
  "gene # 1: kappa = 1.70000 omega = 0.30000\n" +
  "gene # 2: kappa = 1.90000 omega = 1.20000\n"
}

///|
test "codeml control parses required paths" {
  let control = @src.codeml_parse_control(@src.codeml_example_control_text())
  assert_eq(control.sequence_file(), "alignment.phylip")
  assert_eq(control.output_file(), "codeml.out")
  assert_eq(control.tree_file(), "species.tree")
}

///|
test "codeml control parses options" {
  let control = @src.codeml_parse_control(@src.codeml_example_control_text())
  assert_eq(control.option("CodonFreq"), Some("2"))
  assert_eq(control.option("missing"), None)
}

///|
test "codeml control parses NSsites integer list" {
  let control = @src.codeml_parse_control(@src.codeml_example_control_text())
  assert_eq(control.ns_sites(), [0, 1, 2])
}

///|
test "codeml control strips star comments" {
  let control = @src.codeml_parse_control(
    "seqfile=a.phy * alignment\noutfile=x.out\ntreefile=t.tree\nomega=1 * start\n",
  )
  assert_eq(control.sequence_file(), "a.phy")
  assert_eq(control.option("omega"), Some("1"))
}

///|
test "codeml control accepts CRLF" {
  let control = @src.codeml_parse_control(
    "seqfile=a.phy\r\noutfile=x.out\r\ntreefile=t.tree\r\nrunmode=-2\r\n",
  )
  assert_eq(control.option("runmode"), Some("-2"))
}

///|
test "codeml control canonical round trip" {
  let control = @src.codeml_parse_control(@src.codeml_example_control_text())
  let written = @src.codeml_write_control(control)
  let reparsed = @src.codeml_parse_control(written)
  assert_eq(reparsed, control)
}

///|
test "codeml control preserves option order" {
  let control = @src.codeml_parse_control(
    "seqfile=a\noutfile=b\ntreefile=c\nomega=1\nkappa=2\nmodel=0\n",
  )
  let options = control.options()
  assert_eq(options[0].name(), "omega")
  assert_eq(options[1].name(), "kappa")
  assert_eq(options[2].name(), "model")
}

///|
test "codeml control option accessor is defensive" {
  let control = @src.codeml_parse_control(@src.codeml_example_control_text())
  let options = control.options()
  options[0] = options[1]
  assert_eq(control.options()[0].name(), "verbose")
}

///|
test "codeml control rejects missing seqfile" {
  codeml_test_expect_control_error("outfile=x\ntreefile=t\n")
}

///|
test "codeml control rejects missing outfile" {
  codeml_test_expect_control_error("seqfile=a\ntreefile=t\n")
}

///|
test "codeml control rejects missing treefile" {
  codeml_test_expect_control_error("seqfile=a\noutfile=x\n")
}

///|
test "codeml control rejects malformed line" {
  codeml_test_expect_control_error(
    "seqfile=a\noutfile=x\ntreefile=t\nmodel 0\n",
  )
}

///|
test "codeml control rejects unknown option" {
  codeml_test_expect_control_error(
    "seqfile=a\noutfile=x\ntreefile=t\nunknown=1\n",
  )
}

///|
test "codeml control rejects duplicate path" {
  codeml_test_expect_control_error(
    "seqfile=a\nseqfile=b\noutfile=x\ntreefile=t\n",
  )
}

///|
test "codeml control rejects duplicate option" {
  codeml_test_expect_control_error(
    "seqfile=a\noutfile=x\ntreefile=t\nomega=1\nomega=2\n",
  )
}

///|
test "codeml control rejects invalid integer" {
  codeml_test_expect_control_error(
    "seqfile=a\noutfile=x\ntreefile=t\nmodel=zero\n",
  )
}

///|
test "codeml control rejects invalid double" {
  codeml_test_expect_control_error(
    "seqfile=a\noutfile=x\ntreefile=t\nomega=one\n",
  )
}

///|
test "codeml control rejects negative NSsites" {
  codeml_test_expect_control_error(
    "seqfile=a\noutfile=x\ntreefile=t\nNSsites=0 -1\n",
  )
}

///|
test "codeml results parse header metadata" {
  let results = @src.codeml_parse_results(@src.codeml_example_results_text())
  assert_eq(results.program(), "CODONML")
  assert_eq(results.version(), "4.7")
  assert_eq(results.codon_frequency_model(), "F3x4")
  assert_eq(results.sequence_count(), 5)
  assert_eq(results.site_count(), 74)
}

///|
test "codeml results parse multiple NSsites models" {
  let results = @src.codeml_parse_results(@src.codeml_example_results_text())
  assert_eq(results.models().length(), 2)
  assert_eq(results.models()[0].number(), 0)
  assert_eq(results.models()[1].number(), 2)
}

///|
test "codeml model parses likelihood and parameter count" {
  let model = @src.codeml_parse_results(@src.codeml_example_results_text()).models()[0]
  codeml_test_close(model.ln_likelihood(), -308.032865, 1.0e-9)
  assert_eq(model.parameter_count(), 9)
  assert_true(model.parameter_list().length() > 0)
}

///|
test "codeml model parses tree and length" {
  let model = @src.codeml_parse_results(@src.codeml_example_results_text()).models()[0]
  codeml_test_close(model.tree_length().unwrap(), 0.05472, 1.0e-10)
  assert_true(model.tree().length() > 20)
}

///|
test "codeml model parses kappa and omega" {
  let model = @src.codeml_parse_results(@src.codeml_example_results_text()).models()[0]
  codeml_test_close(model.kappa().unwrap(), 1.51919, 1.0e-10)
  codeml_test_close(model.omegas()[0], 0.00010, 1.0e-10)
}

///|
test "codeml model parses branch table" {
  let model = @src.codeml_parse_results(@src.codeml_example_results_text()).models()[0]
  assert_eq(model.branches().length(), 2)
  assert_eq(model.branches()[1].branch(), "8..2")
  codeml_test_close(model.branches()[1].ds().unwrap(), 0.0183, 1.0e-10)
}

///|
test "codeml model parses dN and dS tree lengths" {
  let model = @src.codeml_parse_results(@src.codeml_example_results_text()).models()[0]
  codeml_test_close(model.dn_tree_length().unwrap(), 0.0, 1.0e-12)
  codeml_test_close(model.ds_tree_length().unwrap(), 0.0746, 1.0e-10)
}

///|
test "codeml model parses site classes" {
  let model = @src.codeml_parse_results(@src.codeml_example_results_text()).models()[1]
  assert_eq(model.site_classes().length(), 3)
  codeml_test_close(model.site_classes()[2].proportion(), 0.02, 1.0e-10)
  codeml_test_close(model.site_classes()[2].omega().unwrap(), 5.0, 1.0e-10)
}

///|
test "codeml model parses BEB positive site" {
  let model = @src.codeml_parse_results(@src.codeml_example_results_text()).models()[1]
  let sites = model.positive_sites()
  assert_eq(sites.length(), 1)
  assert_eq(sites[0].analysis_method(), "BEB")
  assert_eq(sites[0].position(), 42)
  assert_eq(sites[0].amino_acid(), "K")
  assert_eq(sites[0].significance(), "**")
  codeml_test_close(sites[0].probability(), 0.997, 1.0e-10)
  codeml_test_close(sites[0].posterior_mean_omega().unwrap(), 5.2, 1.0e-10)
}

///|
test "codeml model lookup finds requested NSsites class" {
  let results = @src.codeml_parse_results(@src.codeml_example_results_text())
  assert_eq(results.model_by_number(2).unwrap().number(), 2)
  assert_eq(results.model_by_number(8), None)
}

///|
test "codeml result model accessor is deeply defensive" {
  let results = @src.codeml_parse_results(@src.codeml_example_results_text())
  let models = results.models()
  let classes = models[1].site_classes()
  classes[0] = classes[1]
  assert_eq(results.models()[1].site_classes()[0].index(), 0)
}

///|
test "codeml branch-site A parses four classes" {
  let model = @src.codeml_parse_results(codeml_test_branch_site_text()).models()[0]
  assert_eq(model.number(), 2)
  assert_eq(model.site_classes().length(), 4)
}

///|
test "codeml branch-site A parses foreground and background omega" {
  let classes = @src.codeml_parse_results(codeml_test_branch_site_text()).models()[0].site_classes()
  assert_eq(classes[2].branch_types().length(), 2)
  assert_eq(classes[2].branch_types()[0].name(), "background")
  codeml_test_close(classes[2].branch_types()[1].value(), 4.5, 1.0e-10)
}

///|
test "codeml branch-site A parses multiple selected sites" {
  let sites = @src.codeml_parse_results(codeml_test_branch_site_text()).models()[0].positive_sites()
  assert_eq(sites.length(), 2)
  assert_eq(sites[0].significance(), "*")
  assert_eq(sites[1].significance(), "**")
}

///|
test "codeml clade model C parses branch types" {
  let classes = @src.codeml_parse_results(codeml_test_clade_text()).models()[0].site_classes()
  assert_eq(classes.length(), 3)
  assert_eq(classes[0].branch_types()[0].name(), "0")
  assert_eq(classes[0].branch_types()[1].name(), "1")
  codeml_test_close(classes[2].branch_types()[1].value(), 4.0, 1.0e-10)
}

///|
test "codeml free-ratio parses branch omegas" {
  let model = @src.codeml_parse_results(codeml_test_free_ratio_text()).models()[0]
  assert_eq(model.omegas().length(), 4)
  codeml_test_close(model.omegas()[0], 1.1, 1.0e-10)
}

///|
test "codeml free-ratio parses specialized trees" {
  let model = @src.codeml_parse_results(codeml_test_free_ratio_text()).models()[0]
  assert_true(model.dn_tree().length() > 0)
  assert_true(model.ds_tree().length() > 0)
  assert_true(model.omega_tree().length() > 0)
}

///|
test "codeml free-ratio parses branch rows" {
  let branches = @src.codeml_parse_results(codeml_test_free_ratio_text()).models()[0].branches()
  assert_eq(branches.length(), 2)
  codeml_test_close(branches[0].omega().unwrap(), 1.1, 1.0e-10)
}

///|
test "codeml pairwise parses all unordered comparisons" {
  let pairwise = @src.codeml_parse_results(codeml_test_pairwise_text()).pairwise()
  assert_eq(pairwise.length(), 3)
  assert_eq(pairwise[0].first(), "Pan_troglo")
  assert_eq(pairwise[0].second(), "Homo_sapie")
}

///|
test "codeml pairwise parses likelihood" {
  let pair = @src.codeml_parse_results(codeml_test_pairwise_text()).pairwise()[0]
  codeml_test_close(pair.ln_likelihood(), -291.465693, 1.0e-9)
}

///|
test "codeml pairwise parses site counts and rates" {
  let pair = @src.codeml_parse_results(codeml_test_pairwise_text()).pairwise()[0]
  codeml_test_close(pair.time(), 0.0126, 1.0e-10)
  codeml_test_close(pair.synonymous_sites(), 81.4, 1.0e-10)
  codeml_test_close(pair.nonsynonymous_sites(), 140.6, 1.0e-10)
  codeml_test_close(pair.omega(), 0.001, 1.0e-10)
  codeml_test_close(pair.dn(), 0.0, 1.0e-10)
  codeml_test_close(pair.ds(), 0.0115, 1.0e-10)
}

///|
test "codeml AAML parses program and unconstrained likelihood" {
  let results = @src.codeml_parse_results(codeml_test_aaml_text())
  assert_eq(results.program(), "AAML")
  codeml_test_close(results.ln_l_max().unwrap(), -192.057157, 1.0e-9)
}

///|
test "codeml AAML parses raw distance triangle" {
  let distances = @src.codeml_parse_results(codeml_test_aaml_text()).distances()
  assert_eq(distances.length(), 6)
  assert_eq(distances[0].kind(), "raw")
  assert_eq(distances[0].first(), "Chimp")
  assert_eq(distances[0].second(), "Human")
  codeml_test_close(distances[0].value(), 0.01, 1.0e-10)
}

///|
test "codeml AAML parses ML distance triangle" {
  let distances = @src.codeml_parse_results(codeml_test_aaml_text()).distances()
  assert_eq(distances[3].kind(), "ml")
  codeml_test_close(distances[5].value(), 0.033, 1.0e-10)
}

///|
test "codeml parses parameter standard errors" {
  let model = @src.codeml_parse_results(codeml_test_se_text()).models()[0]
  assert_true(model.standard_errors().length() > 0)
  assert_eq(model.parameter_count(), 6)
}

///|
test "codeml parses separate multi-gene sections" {
  let models = @src.codeml_parse_results(codeml_test_multigene_text()).models()
  assert_eq(models.length(), 2)
  assert_eq(models[0].gene(), 1)
  assert_eq(models[1].gene(), 2)
  codeml_test_close(models[1].tree_length().unwrap(), 0.2, 1.0e-10)
}

///|
test "codeml parses joint multi-gene relative rates" {
  let model = @src.codeml_parse_results(codeml_test_joint_gene_text()).models()[0]
  assert_eq(model.rates(), [1.0, 2.5])
}

///|
test "codeml parses joint multi-gene kappa and omega" {
  let genes = @src.codeml_parse_results(codeml_test_joint_gene_text()).models()[0].gene_parameters()
  assert_eq(genes.length(), 2)
  assert_eq(genes[1].gene(), 2)
  codeml_test_close(genes[1].kappa(), 1.9, 1.0e-10)
  codeml_test_close(genes[1].omega(), 1.2, 1.0e-10)
}

///|
test "codeml AIC matches definition" {
  let model = @src.codeml_parse_results(@src.codeml_example_results_text()).models()[0]
  codeml_test_close(model.aic(), 634.06573, 1.0e-6)
}

///|
test "codeml BIC matches definition" {
  let model = @src.codeml_parse_results(@src.codeml_example_results_text()).models()[0]
  let expected = 9.0 * @math.ln(74.0) + 616.06573
  codeml_test_close(model.bic(74), expected, 1.0e-6)
}

///|
test "codeml likelihood-ratio statistic and degrees" {
  let models = @src.codeml_parse_results(@src.codeml_example_results_text()).models()
  let comparison = @src.codeml_likelihood_ratio(models[0], models[1])
  codeml_test_close(comparison.statistic(), 12.002572, 1.0e-6)
  assert_eq(comparison.degrees_of_freedom(), 3)
}

///|
test "codeml likelihood-ratio chi-square tail" {
  let models = @src.codeml_parse_results(@src.codeml_example_results_text()).models()
  let comparison = @src.codeml_likelihood_ratio(models[0], models[1])
  assert_true(comparison.p_value() > 0.007)
  assert_true(comparison.p_value() < 0.008)
  assert_true(comparison.significant())
  codeml_test_close(comparison.alpha(), 0.05, 1.0e-12)
}

///|
test "codeml best AIC selects improved model" {
  let models = @src.codeml_parse_results(@src.codeml_example_results_text()).models()
  assert_eq(@src.codeml_best_aic(models).number(), 2)
}

///|
test "codeml BIC rejects non-positive observations" {
  let model = @src.codeml_parse_results(@src.codeml_example_results_text()).models()[0]
  let raised = try {
    ignore(model.bic(0))
    false
  } catch {
    @src.CodemlError(_) => true
  }
  assert_true(raised)
}

///|
test "codeml LRT rejects reversed nesting" {
  let models = @src.codeml_parse_results(@src.codeml_example_results_text()).models()
  let raised = try {
    ignore(@src.codeml_likelihood_ratio(models[1], models[0]))
    false
  } catch {
    @src.CodemlError(_) => true
  }
  assert_true(raised)
}

///|
test "codeml LRT rejects invalid alpha" {
  let models = @src.codeml_parse_results(@src.codeml_example_results_text()).models()
  let raised = try {
    ignore(@src.codeml_likelihood_ratio(models[0], models[1], alpha=1.0))
    false
  } catch {
    @src.CodemlError(_) => true
  }
  assert_true(raised)
}

///|
test "codeml best AIC rejects empty input" {
  let raised = try {
    ignore(@src.codeml_best_aic([]))
    false
  } catch {
    @src.CodemlError(_) => true
  }
  assert_true(raised)
}

///|
test "codeml results reject empty text" {
  codeml_test_expect_result_error("")
}

///|
test "codeml results reject missing PAML header" {
  codeml_test_expect_result_error(
    "Model: one-ratio\nlnL(ntime: 1 np: 1): -1.0\n",
  )
}

///|
test "codeml results reject header without estimates" {
  codeml_test_expect_result_error(
    "CODONML (in paml version 4.7, January 2013) a.phy\nModel: x\n",
  )
}

///|
test "codeml results reject malformed model likelihood" {
  codeml_test_expect_result_error(
    "CODONML (in paml version 4.7, January 2013) a.phy\n" +
    "Model 0: one-ratio\n" +
    "lnL(ntime: 1 np: 1): not-a-number\n",
  )
}

///|
test "codeml results reject incomplete pairwise estimate" {
  codeml_test_expect_result_error(
    "CODONML (in paml version 4.7, January 2013) a.phy\n" +
    "2 (B) ... 1 (A)\n" +
    "lnL = -1.0\n" +
    "t=0.1 S=10 N=20 dN/dS=0.5 dN=0.01\n",
  )
}
