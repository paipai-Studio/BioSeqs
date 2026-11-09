///|
/// Tests for pandas-style query expressions over obs/var.

///|
fn q_build() -> @src.AnnData {
  @src.AnnData::new(
    [[1.0, 10.0], [2.0, 20.0], [3.0, 30.0], [4.0, 40.0]],
    ["c1", "c2", "c3", "c4"],
    ["g1", "g2"],
  )
  .set_obs("batch", ["a", "b", "a", "b"])
  .set_obs("n_genes", ["100", "500", "900", "500"])
  .set_var("highly_variable", ["true", "false"])
}

///|
test "query_tokenize_basic" {
  match @src.anndata_query_tokenize("batch == 'a'") {
    Some(toks) => {
      match toks[0] {
        @src.QIdent(s) => assert_eq(s, "batch")
        _ => fail("tok0 should be ident")
      }
      match toks[1] {
        @src.QOp(s) => assert_eq(s, "==")
        _ => fail("tok1 should be op")
      }
      match toks[2] {
        @src.QStr(s) => assert_eq(s, "a")
        _ => fail("tok2 should be string")
      }
      match toks[3] {
        @src.QEnd => ()
        _ => fail("tok3 should be end")
      }
    }
    None => fail("tokenize failed")
  }
  // stray characters rejected
  match @src.anndata_query_tokenize("batch # a") {
    None => ()
    Some(_) => fail("stray chars should not tokenize")
  }
}

///|
test "query_parse_precedence" {
  // & binds tighter than |
  match
    @src.anndata_query_parse("batch == 'a' | batch == 'b' & n_genes > 100") {
    Some(@src.QOr(_, @src.QAnd(_, _))) => ()
    _ => fail("precedence wrong")
  }
  // parentheses override
  match
    @src.anndata_query_parse("(batch == 'a' | batch == 'b') & n_genes > 100") {
    Some(@src.QAnd(@src.QOr(_, _), _)) => ()
    _ => fail("parenthesized shape wrong")
  }
  // not
  match @src.anndata_query_parse("! (batch == 'a')") {
    Some(@src.QNot(_)) => ()
    _ => fail("not shape wrong")
  }
  // parse errors
  match @src.anndata_query_parse("batch ==") {
    None => ()
    Some(_) => fail("dangling op should not parse")
  }
  match @src.anndata_query_parse("(batch == 'a'") {
    None => ()
    Some(_) => fail("unclosed paren should not parse")
  }
}

///|
test "query_obs_string_equality" {
  let ad = q_build()
  let hits = @src.anndata_obs_query(ad, "batch == 'a'")
  assert_eq(hits.length(), 2)
  assert_eq(hits[0], 0)
  assert_eq(hits[1], 2)
}

///|
test "query_obs_numeric_comparison" {
  let ad = q_build()
  let hits = @src.anndata_obs_query(ad, "n_genes > 100")
  assert_eq(hits.length(), 3)
  assert_eq(hits[0], 1)
  let between = @src.anndata_obs_query(ad, "n_genes >= 500 & n_genes <= 500")
  assert_eq(between.length(), 2)
}

///|
test "query_obs_combined" {
  let ad = q_build()
  let hits = @src.anndata_obs_query(ad, "batch == 'b' & n_genes > 100")
  assert_eq(hits.length(), 2)
  assert_eq(hits[0], 1)
  assert_eq(hits[1], 3)
  let either = @src.anndata_obs_query(ad, "batch == 'a' | n_genes > 800")
  assert_eq(either.length(), 2)
  assert_eq(either[0], 0)
  assert_eq(either[1], 2)
  let neg = @src.anndata_obs_query(ad, "! (batch == 'a')")
  assert_eq(neg.length(), 2)
}

///|
test "query_var_axis" {
  let ad = q_build()
  let hits = @src.anndata_var_query(ad, "highly_variable == 'true'")
  assert_eq(hits.length(), 1)
  assert_eq(hits[0], 0)
}

///|
test "query_subset" {
  let ad = q_build()
  let sub = @src.anndata_query_subset(ad, "batch == 'a'")
  assert_eq(sub.n_obs(), 2)
  assert_eq(sub.obs_names[0], "c1")
  assert_eq(sub.obs_names[1], "c3")
  assert_eq(sub.get_x(0, 0), 1.0)
  let sub_var = @src.anndata_query_subset_var(ad, "highly_variable == 'false'")
  assert_eq(sub_var.n_vars(), 1)
  assert_eq(sub_var.var_names[0], "g2")
}

///|
test "query_missing_column_and_parse_error" {
  let ad = q_build()
  // missing column yields no matches (empty cell "")
  let hits = @src.anndata_obs_query(ad, "missing == 'a'")
  assert_eq(hits.length(), 0)
  // parse error yields empty index list
  assert_eq(@src.anndata_obs_query(ad, "batch === a").length(), 0)
}
