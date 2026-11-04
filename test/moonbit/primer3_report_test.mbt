///|
test "primer3 report json parseable" {
  let template = "ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAGATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAGATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG"
  let report = @src.primer3_design_report_json(template)
  let parsed = @src.parse_json(report)
  assert_true(parsed is Some(_))
  match parsed {
    Some(JsonObj(m)) => {
      match m.get("module") {
        Some(JsonStr(s)) => assert_eq(s, "primer3_report")
        _ => assert_true(false)
      }
      assert_true(m.contains("template"))
      assert_true(m.contains("forward"))
      assert_true(m.contains("reverse"))
      assert_true(m.contains("pairs"))
    }
    _ => assert_true(false)
  }
}

///|
test "primer3 report template stats" {
  let template = "ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAGATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAGATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG"
  let v = @src.primer3_design_report_value(
    template,
    @src.primer3_default_config(),
    50,
  )
  match v {
    JsonObj(m) =>
      match m.get("template") {
        Some(JsonObj(t)) =>
          match (t.get("length"), t.get("gc")) {
            (Some(JsonNum(l)), Some(JsonNum(g))) => {
              assert_eq(l, template.length().to_double())
              assert_true(g > 0.4 && g < 0.8)
            }
            _ => assert_true(false)
          }
        _ => assert_true(false)
      }
    _ => assert_true(false)
  }
}

///|
test "primer3 report candidates present" {
  let template = "ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAGATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAGATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG"
  let design = @src.primer3_design_primers(
    template,
    @src.primer3_default_config(),
  )
  let v = @src.primer3_design_report_value(
    template,
    @src.primer3_default_config(),
    50,
  )
  match v {
    JsonObj(m) => {
      match m.get("forward") {
        Some(JsonArr(f)) => assert_eq(f.length(), design.n_forward())
        _ => assert_true(false)
      }
      match m.get("reverse") {
        Some(JsonArr(r)) => assert_eq(r.length(), design.n_reverse())
        _ => assert_true(false)
      }
    }
    _ => assert_true(false)
  }
}

///|
test "primer3 report pairs product length" {
  let template = "ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAGATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAGATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG"
  let v = @src.primer3_design_report_value(
    template,
    @src.primer3_default_config(),
    50,
  )
  match v {
    JsonObj(m) =>
      match m.get("pairs") {
        Some(JsonArr(pairs)) =>
          if pairs.length() > 0 {
            match pairs[0] {
              JsonObj(p) => {
                match p.get("product_len") {
                  Some(JsonNum(pl)) => assert_true(pl > 0.0 && pl <= 5000.0)
                  _ => assert_true(false)
                }
                assert_true(p.contains("forward_index"))
                assert_true(p.contains("reverse_index"))
                assert_true(p.contains("pair_score"))
              }
              _ => assert_true(false)
            }
          }
        _ => assert_true(false)
      }
    _ => assert_true(false)
  }
}

///|
test "primer3 report max pairs limit" {
  let template = "ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAGATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAGATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG"
  let v = @src.primer3_design_report_value(
    template,
    @src.primer3_default_config(),
    3,
  )
  match v {
    JsonObj(m) =>
      match m.get("pairs") {
        Some(JsonArr(pairs)) => assert_true(pairs.length() <= 3)
        _ => assert_true(false)
      }
    _ => assert_true(false)
  }
}

///|
test "primer3 report short template empty lists" {
  let v = @src.primer3_design_report_value(
    "ACGTACGTAC",
    @src.primer3_default_config(),
    50,
  )
  match v {
    JsonObj(m) => {
      match m.get("forward") {
        Some(JsonArr(f)) => assert_eq(f.length(), 0)
        _ => assert_true(false)
      }
      match m.get("reverse") {
        Some(JsonArr(r)) => assert_eq(r.length(), 0)
        _ => assert_true(false)
      }
      match m.get("pairs") {
        Some(JsonArr(p)) => assert_eq(p.length(), 0)
        _ => assert_true(false)
      }
    }
    _ => assert_true(false)
  }
}

///|
test "primer3 result to json fields" {
  let template = "ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAGATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAGATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG"
  let design = @src.primer3_design_primers(
    template,
    @src.primer3_default_config(),
  )
  if design.n_forward() > 0 {
    let j = @src.primer3_result_to_json(design.forward()[0])
    match j {
      JsonObj(m) => {
        assert_true(m.contains("sequence"))
        assert_true(m.contains("start"))
        assert_true(m.contains("length"))
        assert_true(m.contains("tm"))
        assert_true(m.contains("gc"))
        assert_true(m.contains("self_complement"))
        assert_true(m.contains("is_forward"))
        assert_true(m.contains("score"))
        match m.get("is_forward") {
          Some(JsonBool(b)) => assert_true(b)
          _ => assert_true(false)
        }
      }
      _ => assert_true(false)
    }
  }
}
