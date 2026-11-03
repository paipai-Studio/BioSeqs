///|
fn yn00_test_close(
  actual : Double,
  expected : Double,
  tolerance : Double,
) -> Unit {
  if (actual - expected).abs() > tolerance {
    abort(
      "YN00 value " +
      actual.to_string() +
      " differs from " +
      expected.to_string(),
    )
  }
}

///|
fn yn00_test_expect_control_error(text : String) -> Unit {
  let raised = try {
    ignore(@src.yn00_parse_control(text))
    false
  } catch {
    @src.Yn00Error(_) => true
  }
  if !raised {
    abort("expected Yn00Error")
  }
}

///|
fn yn00_test_expect_result_error(text : String) -> Unit {
  let raised = try {
    ignore(@src.yn00_parse_results(text))
    false
  } catch {
    @src.Yn00Error(_) => true
  }
  if !raised {
    abort("expected Yn00Error")
  }
}

///|
fn yn00_test_two_sequence_result(
  first : String,
  second : String,
  ng_row : String,
  lwl85_omega : String,
  modified_ds : String,
) -> String {
  "YN00 pairwise.phy\n" +
  "ns = 2 ls = 100\n" +
  "(A) Nei-Gojobori (1986) method\n" +
  first +
  "\n" +
  ng_row +
  "\n" +
  "(B) Yang & Nielsen (2000) method\n" +
  "2 1 30.0 70.0 0.1 2.0 0.5 0.01 +- 0.001 0.02 +- 0.002\n" +
  "(C) LWL85, LPB93 & LWLm methods\n" +
  "2 (" +
  second +
  ") vs. 1 (" +
  first +
  ")\n" +
  "LWL85: dS = 0.0200 dN = 0.0100 w = " +
  lwl85_omega +
  " S = 30.0 N = 70.0\n" +
  "LWL85m: dS = " +
  modified_ds +
  " dN = -nan w = -nan S = -nan N = -nan (rho = -nan)\n" +
  "LPB93: dS = 0.0200 dN = 0.0100 w = 0.5000\n"
}

///|
test "yn00 control parses required paths" {
  let control = @src.yn00_parse_control(@src.yn00_example_control_text())
  assert_eq(control.sequence_file(), "alignment.phylip")
  assert_eq(control.output_file(), "yn00.out")
}

///|
test "yn00 control parses all official options" {
  let control = @src.yn00_parse_control(@src.yn00_example_control_text())
  assert_eq(control.option("verbose"), Some("1"))
  assert_eq(control.option("icode"), Some("0"))
  assert_eq(control.option("weighting"), Some("0"))
  assert_eq(control.option("commonf3x4"), Some("0"))
  assert_eq(control.option("ndata"), Some("1"))
}

///|
test "yn00 control ignores comments and blank lines" {
  let control = @src.yn00_parse_control(
    "* header\n\nseqfile = a.phy * input\noutfile = out.txt * output\n",
  )
  assert_eq(control.sequence_file(), "a.phy")
  assert_eq(control.output_file(), "out.txt")
}

///|
test "yn00 control accepts CRLF" {
  let text = @src.yn00_example_control_text().replace_all(old="\n", new="\r\n")
  let control = @src.yn00_parse_control(text)
  assert_eq(control.options().length(), 5)
}

///|
test "yn00 control canonical round trip" {
  let control = @src.yn00_parse_control(@src.yn00_example_control_text())
  let canonical = @src.yn00_write_control(control)
  assert_eq(
    @src.yn00_write_control(@src.yn00_parse_control(canonical)),
    canonical,
  )
}

///|
test "yn00 control preserves option order" {
  let control = @src.yn00_parse_control(
    "seqfile=a\noutfile=b\nndata=2\nverbose=0\nicode=3\n",
  )
  let options = control.options()
  assert_eq(options[0].name(), "ndata")
  assert_eq(options[1].name(), "verbose")
  assert_eq(options[2].name(), "icode")
}

///|
test "yn00 control options are defensive copies" {
  let control = @src.yn00_parse_control(@src.yn00_example_control_text())
  let options = control.options()
  ignore(options.pop())
  assert_eq(control.options().length(), 5)
}

///|
test "yn00 control reports absent option" {
  let control = @src.yn00_parse_control("seqfile=a\noutfile=b\n")
  assert_eq(control.option("icode"), None)
}

///|
test "yn00 control normalizes integers" {
  let control = @src.yn00_parse_control(
    "seqfile=a\noutfile=b\nverbose=+1\nicode=03\nndata=002\n",
  )
  assert_eq(control.option("verbose"), Some("1"))
  assert_eq(control.option("icode"), Some("3"))
  assert_eq(control.option("ndata"), Some("2"))
}

///|
test "yn00 control rejects missing seqfile" {
  yn00_test_expect_control_error("outfile=out\n")
}

///|
test "yn00 control rejects missing outfile" {
  yn00_test_expect_control_error("seqfile=in\n")
}

///|
test "yn00 control rejects line without equals" {
  yn00_test_expect_control_error("seqfile=in\noutfile out\n")
}

///|
test "yn00 control rejects multiple equals" {
  yn00_test_expect_control_error("seqfile=in=x\noutfile=out\n")
}

///|
test "yn00 control rejects empty key" {
  yn00_test_expect_control_error("=in\noutfile=out\n")
}

///|
test "yn00 control rejects empty value" {
  yn00_test_expect_control_error("seqfile=\noutfile=out\n")
}

///|
test "yn00 control rejects duplicate key" {
  yn00_test_expect_control_error("seqfile=in\noutfile=out\nseqfile=again\n")
}

///|
test "yn00 control rejects unknown option" {
  yn00_test_expect_control_error("seqfile=in\noutfile=out\nunknown=1\n")
}

///|
test "yn00 control rejects non-integer option" {
  yn00_test_expect_control_error("seqfile=in\noutfile=out\nicode=1.5\n")
}

///|
test "yn00 control rejects verbose outside binary range" {
  yn00_test_expect_control_error("seqfile=in\noutfile=out\nverbose=2\n")
}

///|
test "yn00 control rejects weighting outside binary range" {
  yn00_test_expect_control_error("seqfile=in\noutfile=out\nweighting=-1\n")
}

///|
test "yn00 control rejects commonf3x4 outside binary range" {
  yn00_test_expect_control_error("seqfile=in\noutfile=out\ncommonf3x4=3\n")
}

///|
test "yn00 control rejects negative genetic code" {
  yn00_test_expect_control_error("seqfile=in\noutfile=out\nicode=-1\n")
}

///|
test "yn00 control rejects unsupported genetic code" {
  yn00_test_expect_control_error("seqfile=in\noutfile=out\nicode=11\n")
}

///|
test "yn00 control rejects non-positive data count" {
  yn00_test_expect_control_error("seqfile=in\noutfile=out\nndata=0\n")
}

///|
test "yn00 results parse alignment metadata" {
  let result = @src.yn00_parse_results(@src.yn00_example_results_text())
  assert_eq(result.alignment_file(), "alignment.phylip")
  assert_eq(result.sequence_count(), 3)
  assert_eq(result.codon_count(), 74)
}

///|
test "yn00 results preserve sequence order" {
  let result = @src.yn00_parse_results(@src.yn00_example_results_text())
  assert_eq(result.sequence_names(), ["Homo_sapie", "Pan_troglo", "Gorilla_go"])
}

///|
test "yn00 results contain complete triangular pairs" {
  let result = @src.yn00_parse_results(@src.yn00_example_results_text())
  assert_eq(result.pairs().length(), 3)
}

///|
test "yn00 pair lookup is symmetric" {
  let result = @src.yn00_parse_results(@src.yn00_example_results_text())
  let forward = result.pair("Homo_sapie", "Pan_troglo").unwrap()
  let reverse = result.pair("Pan_troglo", "Homo_sapie").unwrap()
  assert_eq(forward, reverse)
}

///|
test "yn00 pair lookup rejects diagonal" {
  let result = @src.yn00_parse_results(@src.yn00_example_results_text())
  assert_eq(result.pair("Homo_sapie", "Homo_sapie"), None)
}

///|
test "yn00 pair lookup reports unknown names" {
  let result = @src.yn00_parse_results(@src.yn00_example_results_text())
  assert_eq(result.pair("Homo_sapie", "Unknown"), None)
}

///|
test "yn00 parses NG86 estimate" {
  let pair = @src.yn00_parse_results(@src.yn00_example_results_text())
    .pair("Homo_sapie", "Pan_troglo")
    .unwrap()
  yn00_test_close(pair.ng86().omega(), -1.0, 1.0e-12)
  yn00_test_close(pair.ng86().dn(), 0.0, 1.0e-12)
  yn00_test_close(pair.ng86().ds(), 0.0207, 1.0e-12)
}

///|
test "yn00 parses Yang-Nielsen site counts" {
  let estimate = @src.yn00_parse_results(@src.yn00_example_results_text())
    .pair("Homo_sapie", "Pan_troglo")
    .unwrap()
    .yang_nielsen()
  yn00_test_close(estimate.synonymous_sites(), 67.3, 1.0e-12)
  yn00_test_close(estimate.nonsynonymous_sites(), 154.7, 1.0e-12)
}

///|
test "yn00 parses Yang-Nielsen model parameters" {
  let estimate = @src.yn00_parse_results(@src.yn00_example_results_text())
    .pair("Homo_sapie", "Pan_troglo")
    .unwrap()
    .yang_nielsen()
  yn00_test_close(estimate.time(), 0.0136, 1.0e-12)
  yn00_test_close(estimate.kappa(), 3.6564, 1.0e-12)
  yn00_test_close(estimate.omega(), 0.0, 1.0e-12)
}

///|
test "yn00 parses Yang-Nielsen rates and errors" {
  let estimate = @src.yn00_parse_results(@src.yn00_example_results_text())
    .pair("Homo_sapie", "Pan_troglo")
    .unwrap()
    .yang_nielsen()
  yn00_test_close(estimate.dn(), 0.0, 1.0e-12)
  yn00_test_close(estimate.dn_standard_error(), 0.0, 1.0e-12)
  yn00_test_close(estimate.ds(), 0.015, 1.0e-12)
  yn00_test_close(estimate.ds_standard_error(), 0.0151, 1.0e-12)
}

///|
test "yn00 parses LWL85 estimate" {
  let estimate = @src.yn00_parse_results(@src.yn00_example_results_text())
    .pair("Homo_sapie", "Pan_troglo")
    .unwrap()
    .lwl85()
  assert_eq(estimate.method_name(), "LWL85")
  yn00_test_close(estimate.ds().unwrap(), 0.0227, 1.0e-12)
  assert_eq(estimate.dn(), Some(0.0))
  assert_eq(estimate.omega(), Some(0.0))
  assert_eq(estimate.synonymous_sites(), Some(45.0))
  assert_eq(estimate.nonsynonymous_sites(), Some(177.0))
}

///|
test "yn00 maps undefined modified LWL85 values to none" {
  let estimate = @src.yn00_parse_results(@src.yn00_example_results_text())
    .pair("Homo_sapie", "Pan_troglo")
    .unwrap()
    .lwl85_modified()
  assert_eq(estimate.ds(), None)
  assert_eq(estimate.dn(), None)
  assert_eq(estimate.omega(), None)
  assert_eq(estimate.synonymous_sites(), None)
  assert_eq(estimate.nonsynonymous_sites(), None)
  assert_eq(estimate.rho(), None)
}

///|
test "yn00 parses defined modified LWL85 values" {
  let estimate = @src.yn00_parse_results(@src.yn00_example_results_text())
    .pair("Homo_sapie", "Gorilla_go")
    .unwrap()
    .lwl85_modified()
  yn00_test_close(estimate.ds().unwrap(), 0.018, 1.0e-12)
  assert_eq(estimate.dn(), Some(0.011))
  assert_eq(estimate.omega(), Some(0.6111))
  assert_eq(estimate.rho(), Some(0.507))
}

///|
test "yn00 parses LPB93 estimate without site counts" {
  let estimate = @src.yn00_parse_results(@src.yn00_example_results_text())
    .pair("Homo_sapie", "Pan_troglo")
    .unwrap()
    .lpb93()
  assert_eq(estimate.method_name(), "LPB93")
  assert_eq(estimate.ds(), Some(0.0129))
  assert_eq(estimate.dn(), Some(0.0))
  assert_eq(estimate.omega(), Some(0.0))
  assert_eq(estimate.synonymous_sites(), None)
  assert_eq(estimate.nonsynonymous_sites(), None)
}

///|
test "yn00 parses PAML 4.1 compact NG86 rows" {
  let text = yn00_test_two_sequence_result(
    "Alpha", "Beta", "Beta0.5000(0.0100 0.0200)", "0.5000", "-nan",
  )
  let pair = @src.yn00_parse_results(text).pair("Alpha", "Beta").unwrap()
  assert_eq(pair.ng86().omega(), 0.5)
  assert_eq(pair.ng86().dn(), 0.01)
}

///|
test "yn00 parses PAML 4.8 long names glued to negative omega" {
  let first = "patient1_1000_326|1M|XXX|XXX|2"
  let second = "patient1_1000_16|1M|XXX|XXX|20"
  let text = yn00_test_two_sequence_result(
    first,
    second,
    second + "-1.0000 (0.0028 0.0000)",
    "inf",
    "-nan",
  )
  let result = @src.yn00_parse_results(text)
  assert_eq(result.sequence_names(), [first, second])
  assert_eq(result.pair(first, second).unwrap().ng86().dn(), 0.0028)
}

///|
test "yn00 parses PAML 4.9i dotted numeric names" {
  let first = "patient1.1000.326"
  let second = "patient1.1000.16"
  let text = yn00_test_two_sequence_result(
    first,
    second,
    second + " -1.0000 (0.0028 0.0000)",
    "inf",
    "-nan",
  )
  let result = @src.yn00_parse_results(text)
  assert_eq(result.sequence_names(), [first, second])
}

///|
test "yn00 maps infinity to none" {
  let text = yn00_test_two_sequence_result(
    "Alpha", "Beta", "Beta -1.0000 (0.0100 0.0000)", "inf", "-nan",
  )
  let estimate = @src.yn00_parse_results(text)
    .pair("Alpha", "Beta")
    .unwrap()
    .lwl85()
  assert_eq(estimate.omega(), None)
}

///|
test "yn00 maps Windows indefinite values to none" {
  let text = yn00_test_two_sequence_result(
    "Alpha", "Beta", "Beta -1.0000 (0.0100 0.0000)", "0.5000", "-1.#IND",
  )
  let estimate = @src.yn00_parse_results(text)
    .pair("Alpha", "Beta")
    .unwrap()
    .lwl85_modified()
  assert_eq(estimate.ds(), None)
}

///|
test "yn00 results accept CRLF" {
  let text = @src.yn00_example_results_text().replace_all(old="\n", new="\r\n")
  assert_eq(@src.yn00_parse_results(text).pairs().length(), 3)
}

///|
test "yn00 sequence names are defensive copies" {
  let result = @src.yn00_parse_results(@src.yn00_example_results_text())
  let names = result.sequence_names()
  names[0] = "Changed"
  assert_eq(result.sequence_names()[0], "Homo_sapie")
}

///|
test "yn00 pairs are defensive copies" {
  let result = @src.yn00_parse_results(@src.yn00_example_results_text())
  let pairs = result.pairs()
  ignore(pairs.pop())
  assert_eq(result.pairs().length(), 3)
}

///|
test "yn00 builds symmetric NG86 matrix" {
  let result = @src.yn00_parse_results(@src.yn00_example_results_text())
  let matrix = result.matrix("NG86", "dS")
  assert_eq(matrix.names(), result.sequence_names())
  assert_eq(matrix.value("Homo_sapie", "Pan_troglo"), Some(0.0207))
  assert_eq(matrix.value("Pan_troglo", "Homo_sapie"), Some(0.0207))
  assert_eq(matrix.value("Homo_sapie", "Homo_sapie"), None)
}

///|
test "yn00 builds Yang-Nielsen standard error matrix" {
  let result = @src.yn00_parse_results(@src.yn00_example_results_text())
  let matrix = result.matrix("YN00", "dS SE")
  yn00_test_close(
    matrix.value("Homo_sapie", "Pan_troglo").unwrap(),
    0.0151,
    1.0e-12,
  )
}

///|
test "yn00 matrix preserves undefined values" {
  let result = @src.yn00_parse_results(@src.yn00_example_results_text())
  let matrix = result.matrix("LWL85m", "dS")
  assert_eq(matrix.value("Homo_sapie", "Pan_troglo"), None)
  yn00_test_close(
    matrix.value("Homo_sapie", "Gorilla_go").unwrap(),
    0.018,
    1.0e-12,
  )
}

///|
test "yn00 matrix values are defensive copies" {
  let result = @src.yn00_parse_results(@src.yn00_example_results_text())
  let matrix = result.matrix("NG86", "dS")
  let values = matrix.values()
  values[0][1] = None
  assert_eq(matrix.value("Homo_sapie", "Pan_troglo"), Some(0.0207))
}

///|
test "yn00 matrix reports unknown name" {
  let result = @src.yn00_parse_results(@src.yn00_example_results_text())
  let matrix = result.matrix("NG86", "dS")
  assert_eq(matrix.value("Unknown", "Homo_sapie"), None)
}

///|
test "yn00 computes mean over all pairs" {
  let result = @src.yn00_parse_results(@src.yn00_example_results_text())
  yn00_test_close(
    result.mean("YN00", "dS"),
    (0.015 + 0.02 + 0.0303) / 3.0,
    1.0e-12,
  )
}

///|
test "yn00 mean skips undefined estimates" {
  let result = @src.yn00_parse_results(@src.yn00_example_results_text())
  yn00_test_close(result.mean("LWL85m", "dS"), 0.018, 1.0e-12)
}

///|
test "yn00 matrix rejects unknown method" {
  let result = @src.yn00_parse_results(@src.yn00_example_results_text())
  let raised = try {
    ignore(result.matrix("UNKNOWN", "dS"))
    false
  } catch {
    @src.Yn00Error(_) => true
  }
  assert_true(raised)
}

///|
test "yn00 matrix rejects unknown statistic" {
  let result = @src.yn00_parse_results(@src.yn00_example_results_text())
  let raised = try {
    ignore(result.matrix("YN00", "bad"))
    false
  } catch {
    @src.Yn00Error(_) => true
  }
  assert_true(raised)
}

///|
test "yn00 mean rejects wholly undefined statistic" {
  let text = yn00_test_two_sequence_result(
    "Alpha", "Beta", "Beta -1.0000 (0.0100 0.0000)", "0.5000", "-nan",
  )
  let result = @src.yn00_parse_results(text)
  let raised = try {
    ignore(result.mean("LWL85m", "rho"))
    false
  } catch {
    @src.Yn00Error(_) => true
  }
  assert_true(raised)
}

///|
test "yn00 rejects empty results" {
  yn00_test_expect_result_error("")
}

///|
test "yn00 rejects missing NG86 section" {
  yn00_test_expect_result_error(
    @src.yn00_example_results_text().replace_all(
      old="(A) Nei-Gojobori (1986) method",
      new="",
    ),
  )
}

///|
test "yn00 rejects missing Yang-Nielsen section" {
  yn00_test_expect_result_error(
    @src.yn00_example_results_text().replace_all(
      old="(B) Yang & Nielsen (2000) method",
      new="",
    ),
  )
}

///|
test "yn00 rejects missing counting section" {
  yn00_test_expect_result_error(
    @src.yn00_example_results_text().replace_all(
      old="(C) LWL85, LPB93 & LWLm methods",
      new="",
    ),
  )
}

///|
test "yn00 rejects duplicate method section" {
  yn00_test_expect_result_error(
    @src.yn00_example_results_text() + "(A) Nei-Gojobori (1986) method\n",
  )
}

///|
test "yn00 rejects missing alignment header" {
  yn00_test_expect_result_error(
    @src.yn00_example_results_text().replace_all(
      old="YN00 alignment.phylip\n",
      new="",
    ),
  )
}

///|
test "yn00 rejects duplicate alignment header" {
  yn00_test_expect_result_error(
    "YN00 duplicate.phy\n" + @src.yn00_example_results_text(),
  )
}

///|
test "yn00 rejects malformed dimensions" {
  yn00_test_expect_result_error(
    @src.yn00_example_results_text().replace_all(
      old="ns = 3 ls = 74",
      new="ns = bad ls = 74",
    ),
  )
}

///|
test "yn00 rejects fewer than two sequences" {
  yn00_test_expect_result_error(
    @src.yn00_example_results_text().replace_all(
      old="ns = 3 ls = 74",
      new="ns = 1 ls = 74",
    ),
  )
}

///|
test "yn00 rejects non-positive codon count" {
  yn00_test_expect_result_error(
    @src.yn00_example_results_text().replace_all(
      old="ns = 3 ls = 74",
      new="ns = 3 ls = 0",
    ),
  )
}

///|
test "yn00 rejects invalid comparison index" {
  yn00_test_expect_result_error(
    @src.yn00_example_results_text().replace_all(
      old="2 (Pan_troglo) vs. 1 (Homo_sapie)",
      new="4 (Pan_troglo) vs. 1 (Homo_sapie)",
    ),
  )
}

///|
test "yn00 rejects inconsistent indexed name" {
  yn00_test_expect_result_error(
    @src.yn00_example_results_text().replace_all(
      old="3 (Gorilla_go) vs. 2 (Pan_troglo)",
      new="3 (Different) vs. 2 (Pan_troglo)",
    ),
  )
}

///|
test "yn00 rejects duplicate sequence names" {
  yn00_test_expect_result_error(
    @src.yn00_example_results_text().replace_all(
      old="2 (Pan_troglo) vs. 1 (Homo_sapie)",
      new="2 (Homo_sapie) vs. 1 (Homo_sapie)",
    ),
  )
}

///|
test "yn00 rejects duplicate pair header" {
  let duplicate = "2 (Pan_troglo) vs. 1 (Homo_sapie)\n"
  yn00_test_expect_result_error(@src.yn00_example_results_text() + duplicate)
}

///|
test "yn00 rejects counting estimate before pair header" {
  yn00_test_expect_result_error(
    @src.yn00_example_results_text().replace_all(
      old="(C) LWL85, LPB93 & LWLm methods\n",
      new="(C) LWL85, LPB93 & LWLm methods\n" +
        "LWL85: dS = 0.1 dN = 0.1 w = 1.0 S = 1.0 N = 1.0\n",
    ),
  )
}

///|
test "yn00 rejects missing counting method" {
  yn00_test_expect_result_error(
    @src.yn00_example_results_text().replace_all(
      old="LPB93: dS = 0.0129 dN = 0.0000 w = 0.0000\n",
      new="",
    ),
  )
}

///|
test "yn00 rejects missing counting statistic" {
  yn00_test_expect_result_error(
    @src.yn00_example_results_text().replace_all(
      old="LWL85: dS = 0.0227 dN = 0.0000 w = 0.0000 S = 45.0 N = 177.0",
      new="LWL85: dS = 0.0227 dN = 0.0000 S = 45.0 N = 177.0",
    ),
  )
}

///|
test "yn00 rejects duplicate counting estimate" {
  yn00_test_expect_result_error(
    @src.yn00_example_results_text().replace_all(
      old="LWL85: dS = 0.0227 dN = 0.0000 w = 0.0000 S = 45.0 N = 177.0\n",
      new="LWL85: dS = 0.0227 dN = 0.0000 w = 0.0000 S = 45.0 N = 177.0\n" +
        "LWL85: dS = 0.0227 dN = 0.0000 w = 0.0000 S = 45.0 N = 177.0\n",
    ),
  )
}

///|
test "yn00 rejects malformed Yang-Nielsen row" {
  yn00_test_expect_result_error(
    @src.yn00_example_results_text().replace_all(
      old="2 1 67.3 154.7 0.0136 3.6564 0.0000 0.0000 +- 0.0000 0.0150 +- 0.0151",
      new="2 1 67.3 154.7 0.0136 3.6564",
    ),
  )
}

///|
test "yn00 rejects invalid Yang-Nielsen index" {
  yn00_test_expect_result_error(
    @src.yn00_example_results_text().replace_all(
      old="2 1 67.3 154.7",
      new="4 1 67.3 154.7",
    ),
  )
}

///|
test "yn00 rejects duplicate Yang-Nielsen pair" {
  let row = "2 1 67.3 154.7 0.0136 3.6564 0.0000 0.0000 +- 0.0000 0.0150 +- 0.0151\n"
  yn00_test_expect_result_error(
    @src.yn00_example_results_text().replace_all(
      old="(C) LWL85, LPB93 & LWLm methods",
      new=row + "(C) LWL85, LPB93 & LWLm methods",
    ),
  )
}

///|
test "yn00 rejects missing NG86 row" {
  yn00_test_expect_result_error(
    @src.yn00_example_results_text().replace_all(
      old="Gorilla_go 0.5000 (0.0100 0.0200)-1.0000 (0.0000 0.0421)\n",
      new="",
    ),
  )
}

///|
test "yn00 rejects incomplete NG86 row" {
  yn00_test_expect_result_error(
    @src.yn00_example_results_text().replace_all(
      old="Gorilla_go 0.5000 (0.0100 0.0200)-1.0000 (0.0000 0.0421)",
      new="Gorilla_go 0.5000 (0.0100 0.0200)",
    ),
  )
}

///|
test "yn00 rejects malformed NG86 estimate" {
  yn00_test_expect_result_error(
    @src.yn00_example_results_text().replace_all(
      old="Pan_troglo -1.0000 (0.0000 0.0207)",
      new="Pan_troglo -1.0000 (0.0000)",
    ),
  )
}
