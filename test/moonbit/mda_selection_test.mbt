///|
/// Tests for mda_selection.mbt (legacy selection tokenizer/evaluator).
///
/// The evaluator (mda_sel_evaluate) requires MdaAtomInfo values, but that
/// struct is a read-only pub struct with no public constructor exposed, so
/// from an external test package only the tokenizer can be exercised. These
/// tests therefore cover mda_sel_tokenize token sequences.

///|
test "mda_sel_tokenize_simple_words" {
  let toks = @src.mda_sel_tokenize("name CA")
  assert_eq(toks.length(), 2)
  match toks[0] {
    MdaSelWord(w) => assert_eq(w, "name")
    _ => fail("expected a word token")
  }
  match toks[1] {
    MdaSelWord(w) => assert_eq(w, "CA")
    _ => fail("expected a word token")
  }
}

///|
test "mda_sel_tokenize_operators_and_parens" {
  let toks = @src.mda_sel_tokenize("protein and ( not resname ALA )")
  assert_eq(toks.length(), 7)
  match toks[0] {
    MdaSelWord(w) => assert_eq(w, "protein")
    _ => fail("expected word 'protein'")
  }
  match toks[1] {
    MdaSelOp(op) => assert_eq(op, "and")
    _ => fail("expected op 'and'")
  }
  assert_true(toks[2] is MdaSelLParen)
  match toks[3] {
    MdaSelOp(op) => assert_eq(op, "not")
    _ => fail("expected op 'not'")
  }
  match toks[4] {
    MdaSelWord(w) => assert_eq(w, "resname")
    _ => fail("expected word 'resname'")
  }
  match toks[5] {
    MdaSelWord(w) => assert_eq(w, "ALA")
    _ => fail("expected word 'ALA'")
  }
  assert_true(toks[6] is MdaSelRParen)
}

///|
test "mda_sel_tokenize_numbers_and_colon" {
  // The tokenizer does not split "1:50"; a bare colon token only appears
  // when it is written as a separate word.
  // parse_int parses the leading digits of "1:50" as 1, so the compound
  // word becomes a number token.
  let toks = @src.mda_sel_tokenize("resid 1:50")
  assert_eq(toks.length(), 2)
  match toks[0] {
    MdaSelWord(w) => assert_eq(w, "resid")
    _ => fail("expected word 'resid'")
  }
  match toks[1] {
    MdaSelNumber(v) => assert_true(@moonbit.assert_num_close(v, 1.0, 1.0e-12))
    _ => fail("expected number 1 from '1:50'")
  }
  // Spaced form produces numbers and a colon token.
  let toks2 = @src.mda_sel_tokenize("resid 1 : 50")
  assert_eq(toks2.length(), 4)
  match toks2[1] {
    MdaSelNumber(v) => assert_true(@moonbit.assert_num_close(v, 1.0, 1.0e-12))
    _ => fail("expected number 1")
  }
  assert_true(toks2[2] is MdaSelColon)
  match toks2[3] {
    MdaSelNumber(v) => assert_true(@moonbit.assert_num_close(v, 50.0, 1.0e-12))
    _ => fail("expected number 50")
  }
}

///|
test "mda_sel_tokenize_empty_and_whitespace" {
  assert_eq(@src.mda_sel_tokenize("").length(), 0)
  // Extra spaces collapse away; surrounding whitespace is trimmed.
  let toks = @src.mda_sel_tokenize("  name   CA  ")
  assert_eq(toks.length(), 2)
  match toks[0] {
    MdaSelWord(w) => assert_eq(w, "name")
    _ => fail("expected word 'name'")
  }
  match toks[1] {
    MdaSelWord(w) => assert_eq(w, "CA")
    _ => fail("expected word 'CA'")
  }
}

///|
test "mda_sel_tokenize_or_operator" {
  let toks = @src.mda_sel_tokenize("water or backbone")
  assert_eq(toks.length(), 3)
  match toks[0] {
    MdaSelWord(w) => assert_eq(w, "water")
    _ => fail("expected word 'water'")
  }
  match toks[1] {
    MdaSelOp(op) => assert_eq(op, "or")
    _ => fail("expected op 'or'")
  }
  match toks[2] {
    MdaSelWord(w) => assert_eq(w, "backbone")
    _ => fail("expected word 'backbone'")
  }
}
