// ===========================================================================
// cobra3_io_test.mbt — JSON / YAML / TSV / SBML(FBC v2) roundtrip
// ===========================================================================

///|
test "cobra3 io json roundtrip" {
  let m = cobra3_toy()
  let text = @src.cobra3_model_to_json(m)
  let back = @src.cobra3_model_from_json(text).unwrap()
  assert_eq(back.id, "toy")
  assert_eq(back.metabolites.length(), 2)
  assert_eq(back.reactions.length(), 3)
  let r1 = @src.cobra3_model_get_reaction(back, "R1").unwrap()
  assert_eq(r1.lower_bound, 0.0)
  assert_eq(r1.upper_bound, 1000.0)
  let gpr = @src.cobra3_gpr_to_string(r1.gpr.unwrap())
  assert_true(gpr.contains("g1"))
  assert_true(gpr.contains("g2"))
  let r2 = @src.cobra3_model_get_reaction(back, "R2").unwrap()
  assert_true((r2.objective - 1.0).abs() < 1.0e-9)
}

///|
test "cobra3 io json parser nested" {
  let doc = "{\"a\":[1,2,{\"b\":true}],\"c\":\"x\\ny\",\"d\":null}"
  let v = @src.cobra3_json_parse(doc).unwrap()
  match v {
    JObj(_) => assert_true(true)
    _ => fail("expected object")
  }
}

///|
test "cobra3 io json unicode escape" {
  match @src.cobra3_json_parse("\"\\u4e2d\"") {
    Some(@src.JStr(s)) => assert_eq(s, "中")
    _ => fail("expect string")
  }
}

///|
test "cobra3 io legacy sbml roundtrip keeps notes bounds" {
  let m = legacy_branch_model()
  let xml = @src.cobra_write_sbml(m)
  let back = @src.cobra_read_sbml(xml).unwrap()
  let mut lb = -999.0
  let mut ub = -1.0
  let mut obj = 0.0
  for r in back.reactions {
    if r.id == "R1" {
      lb = r.lower_bound
      ub = r.upper_bound
    }
    if r.id == "BIO" {
      obj = r.objective
    }
  }
  assert_true((lb - 0.0).abs() < 1.0e-9)
  assert_true((ub - 10.0).abs() < 1.0e-9)
  assert_true((obj - 1.0).abs() < 1.0e-9)
}

///|
test "cobra3 io json rejects trailing content" {
  match @src.cobra3_json_parse("{}x") {
    None => assert_true(true)
    Some(_) => fail("trailing content should be rejected")
  }
}

///|
test "cobra3 io json escapes control characters" {
  let mut m = @src.cobra3_model_new("ctrl")
  let ch = (1).to_char().unwrap_or('?')
  let name = "a" + Char::to_string(ch) + "b"
  let met : @src.Cobra3Metabolite = {
    id: "m1",
    name,
    formula: "",
    charge: 0.0,
    compartment: "c",
    boundary: false,
    annotation: Map([], capacity=1),
  }
  m = @src.cobra3_model_add_metabolite(m, met)
  let text = @src.cobra3_model_to_json(m)
  assert_true(text.contains("\\u0001"))
}

///|
test "cobra3 io yaml roundtrip" {
  let m = cobra3_toy()
  let text = @src.cobra3_model_to_yaml(m)
  let back = @src.cobra3_model_from_yaml(text).unwrap()
  assert_eq(back.reactions.length(), 3)
}

///|
test "cobra3 io tsv roundtrip" {
  let m = cobra3_toy()
  let text = @src.cobra3_model_to_tsv(m)
  let back = @src.cobra3_model_from_tsv(text).unwrap()
  assert_eq(back.metabolites.length(), 2)
  assert_eq(back.reactions.length(), 3)
  let r1 = @src.cobra3_model_get_reaction(back, "R1").unwrap()
  assert_true((r1.metabolites["A"] + 1.0).abs() < 1.0e-9)
  assert_true((r1.metabolites["B"] - 1.0).abs() < 1.0e-9)
}

///|
test "cobra3 io sbml fbc roundtrip" {
  let m = cobra3_toy()
  let xml = @src.cobra3_model_to_sbml(m)
  assert_true(xml.contains("<fbc:listOfGeneProducts>"))
  assert_true(xml.contains("fbc:lowerFluxBound"))
  assert_true(xml.contains("<fbc:geneProductAssociation>"))
  let back = @src.cobra3_read_sbml(xml).unwrap()
  assert_eq(back.metabolites.length(), 2)
  assert_eq(back.reactions.length(), 3)
  let r1 = @src.cobra3_model_get_reaction(back, "R1").unwrap()
  assert_true((r1.lower_bound - 0.0).abs() < 1.0e-9)
  assert_true((r1.upper_bound - 1000.0).abs() < 1.0e-9)
  let r2 = @src.cobra3_model_get_reaction(back, "R2").unwrap()
  assert_true((r2.objective - 1.0).abs() < 1.0e-9)
  let gpr = @src.cobra3_gpr_to_string(r1.gpr.unwrap())
  assert_true(gpr.contains("g1"))
}

///|
test "cobra3 io sbml preserves formula and charge" {
  let mut m = @src.cobra3_model_new("fbc")
  let met : @src.Cobra3Metabolite = {
    id: "atp_c",
    name: "ATP",
    formula: "C10H16N5O13P3",
    charge: -4.0,
    compartment: "c",
    boundary: false,
    annotation: Map([], capacity=1),
  }
  m = @src.cobra3_model_add_metabolite(m, met)
  m = @src.cobra3_model_add_reaction_simple(
    m,
    "R",
    [("atp_c", -1.0)],
    0.0,
    1000.0,
  )
  let xml = @src.cobra3_model_to_sbml(m)
  let back = @src.cobra3_read_sbml(xml).unwrap()
  let bm = @src.cobra3_model_get_metabolite(back, "atp_c").unwrap()
  assert_eq(bm.formula, "C10H16N5O13P3")
  assert_true((bm.charge + 4.0).abs() < 1.0e-9)
}
