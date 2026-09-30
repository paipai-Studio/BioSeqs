///|
test "expr basic arithmetic" {
  let f = @src.hts_filter_init("1 + 2 * 3")
  let symbols : Map[String, @src.ExprVal] = {}
  let r = @src.hts_filter_eval(f, symbols)
  assert_true((r.d - 7.0).abs() < 0.000001)
}

///|
test "expr comparison" {
  let f = @src.hts_filter_init("5 > 3")
  let symbols : Map[String, @src.ExprVal] = {}
  let r = @src.hts_filter_eval(f, symbols)
  assert_true((r.d - 1.0).abs() < 0.000001)
}

///|
test "expr boolean and" {
  let f1 = @src.hts_filter_init("1 && 1")
  let f2 = @src.hts_filter_init("1 && 0")
  let f3 = @src.hts_filter_init("0 && 1")
  let symbols : Map[String, @src.ExprVal] = {}
  assert_true((@src.hts_filter_eval(f1, symbols).d - 1.0).abs() < 0.000001)
  assert_true((@src.hts_filter_eval(f2, symbols).d - 0.0).abs() < 0.000001)
  assert_true((@src.hts_filter_eval(f3, symbols).d - 0.0).abs() < 0.000001)
}

///|
test "expr boolean or" {
  let f1 = @src.hts_filter_init("1 || 0")
  let f2 = @src.hts_filter_init("0 || 0")
  let symbols : Map[String, @src.ExprVal] = {}
  assert_true((@src.hts_filter_eval(f1, symbols).d - 1.0).abs() < 0.000001)
  assert_true((@src.hts_filter_eval(f2, symbols).d - 0.0).abs() < 0.000001)
}

///|
test "expr variable lookup" {
  let f = @src.hts_filter_init("QUAL > 30")
  let symbols : Map[String, @src.ExprVal] = {
    "QUAL": @src.ExprVal::num(40.0)
  }
  let r = @src.hts_filter_eval(f, symbols)
  assert_true((r.d - 1.0).abs() < 0.000001)
}

///|
test "expr variable filter" {
  let f = @src.hts_filter_init("QUAL > 30 && DP > 10")
  let symbols_pass : Map[String, @src.ExprVal] = {
    "QUAL": @src.ExprVal::num(40.0),
    "DP": @src.ExprVal::num(15.0)
  }
  let symbols_fail : Map[String, @src.ExprVal] = {
    "QUAL": @src.ExprVal::num(20.0),
    "DP": @src.ExprVal::num(15.0)
  }
  assert_true(@src.hts_filter_test(f, symbols_pass))
  assert_false(@src.hts_filter_test(f, symbols_fail))
}

///|
test "expr bitwise and" {
  let f = @src.hts_filter_init("FLAG & 4")
  let symbols_unmapped : Map[String, @src.ExprVal] = {
    "FLAG": @src.ExprVal::num(4.0)
  }
  let symbols_mapped : Map[String, @src.ExprVal] = {
    "FLAG": @src.ExprVal::num(0.0)
  }
  assert_true(@src.hts_filter_test(f, symbols_unmapped))
  assert_false(@src.hts_filter_test(f, symbols_mapped))
}

///|
test "expr equality" {
  let f1 = @src.hts_filter_init("5 == 5")
  let f2 = @src.hts_filter_init("5 == 3")
  let f3 = @src.hts_filter_init("5 != 3")
  let symbols : Map[String, @src.ExprVal] = {}
  assert_true((@src.hts_filter_eval(f1, symbols).d - 1.0).abs() < 0.000001)
  assert_true((@src.hts_filter_eval(f2, symbols).d - 0.0).abs() < 0.000001)
  assert_true((@src.hts_filter_eval(f3, symbols).d - 1.0).abs() < 0.000001)
}

///|
test "expr unary minus" {
  let f = @src.hts_filter_init("-5 + 10")
  let symbols : Map[String, @src.ExprVal] = {}
  let r = @src.hts_filter_eval(f, symbols)
  assert_true((r.d - 5.0).abs() < 0.000001)
}

///|
test "expr unary not" {
  let f1 = @src.hts_filter_init("!0")
  let f2 = @src.hts_filter_init("!1")
  let symbols : Map[String, @src.ExprVal] = {}
  assert_true((@src.hts_filter_eval(f1, symbols).d - 1.0).abs() < 0.000001)
  assert_true((@src.hts_filter_eval(f2, symbols).d - 0.0).abs() < 0.000001)
}

///|
test "expr parentheses" {
  let f = @src.hts_filter_init("(1 + 2) * 3")
  let symbols : Map[String, @src.ExprVal] = {}
  let r = @src.hts_filter_eval(f, symbols)
  assert_true((r.d - 9.0).abs() < 0.000001)
}

///|
test "expr division and modulo" {
  let f1 = @src.hts_filter_init("10 / 3")
  let f2 = @src.hts_filter_init("10 % 3")
  let symbols : Map[String, @src.ExprVal] = {}
  let r1 = @src.hts_filter_eval(f1, symbols)
  let r2 = @src.hts_filter_eval(f2, symbols)
  assert_true((r1.d - 3.333333).abs() < 0.001)
  assert_true((r2.d - 1.0).abs() < 0.000001)
}

///|
test "expr float arithmetic" {
  let f = @src.hts_filter_init("3.14 * 2.0")
  let symbols : Map[String, @src.ExprVal] = {}
  let r = @src.hts_filter_eval(f, symbols)
  assert_true((r.d - 6.28).abs() < 0.001)
}

///|
test "expr string literal" {
  let f = @src.hts_filter_init("\"hello\"")
  let symbols : Map[String, @src.ExprVal] = {}
  let r = @src.hts_filter_eval(f, symbols)
  assert_eq(r.s, "hello")
}

///|
test "expr func sqrt" {
  let f = @src.hts_filter_init("sqrt(16)")
  let symbols : Map[String, @src.ExprVal] = {}
  let r = @src.hts_filter_eval(f, symbols)
  assert_true((r.d - 4.0).abs() < 0.001)
}

///|
test "expr func pow" {
  let f = @src.hts_filter_init("pow(2, 10)")
  let symbols : Map[String, @src.ExprVal] = {}
  let r = @src.hts_filter_eval(f, symbols)
  assert_true((r.d - 1024.0).abs() < 0.001)
}

///|
test "expr func length" {
  let f = @src.hts_filter_init("length(SEQ)")
  let symbols : Map[String, @src.ExprVal] = {
    "SEQ": @src.ExprVal::str("ACGTACGT")
  }
  let r = @src.hts_filter_eval(f, symbols)
  assert_true((r.d - 8.0).abs() < 0.000001)
}

///|
test "expr complex filter" {
  // Real-world samtools filter: MAPQ>=20 && FLAG&4==0
  let f = @src.hts_filter_init("MAPQ >= 20 && FLAG & 4 == 0")
  let symbols_pass : Map[String, @src.ExprVal] = {
    "MAPQ": @src.ExprVal::num(30.0),
    "FLAG": @src.ExprVal::num(0.0)
  }
  let symbols_fail : Map[String, @src.ExprVal] = {
    "MAPQ": @src.ExprVal::num(30.0),
    "FLAG": @src.ExprVal::num(4.0)
  }
  assert_true(@src.hts_filter_test(f, symbols_pass))
  assert_false(@src.hts_filter_test(f, symbols_fail))
}

///|
test "ExprVal constructors" {
  let n = @src.ExprVal::num(42.0)
  assert_false(n.is_str)
  assert_true((n.d - 42.0).abs() < 0.000001)

  let s = @src.ExprVal::str("hello")
  assert_true(s.is_str)
  assert_eq(s.s, "hello")

  let u = @src.ExprVal::undef()
  assert_eq(u.is_true, -1)
}

