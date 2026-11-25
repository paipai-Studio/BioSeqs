///|
/// Tests for biotite_dotbracket.mbt — RNA dot-bracket parsing.
test "biotite_dotbracket_balanced" {
  let db = @src.BiotiteDotBracket::from_string("((..))")
  assert_eq(db.num_base_pairs(), 2)
  assert_true(!db.has_pseudoknot())
  assert_true(db.is_paired(0))
  assert_true(!db.is_paired(2))
  match db.get_partner(0) {
    Some(p) => assert_eq(p, 5)
    None => fail("position 0 should be paired")
  }
  assert_eq(db.unpaired_positions().length(), 2)
  assert_eq(db.to_string(), "((..))")
}

///|
test "biotite_dotbracket_pseudoknot" {
  let db = @src.BiotiteDotBracket::from_string("([)]")
  assert_eq(db.num_base_pairs(), 2)
  assert_true(db.has_pseudoknot())
}

///|
test "biotite_dotbracket_unpaired" {
  let db = @src.BiotiteDotBracket::from_string("....")
  assert_eq(db.num_base_pairs(), 0)
  assert_eq(db.unpaired_positions().length(), 4)
}

///|
test "biotite_dotbracket_unmatched_brackets" {
  // Unmatched closing and opening brackets do not form pairs.
  let db = @src.BiotiteDotBracket::from_string(").(")
  assert_eq(db.num_base_pairs(), 0)
}
