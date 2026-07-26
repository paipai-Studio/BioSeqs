///|
test "compound_location_join_creation" {
  let loc1 = @src.FeatureLocation::new_exact(1, 100)
  let loc2 = @src.FeatureLocation::new_exact(200, 300)
  let compound = @src.CompoundLocation::join([loc1, loc2])
  assert_eq(compound.operator(), "join")
  assert_eq(compound.parts().length(), 2)
}

///|
test "compound_location_order_creation" {
  let loc1 = @src.FeatureLocation::new_exact(1, 100)
  let loc2 = @src.FeatureLocation::new_exact(200, 300)
  let compound = @src.CompoundLocation::order([loc1, loc2])
  assert_eq(compound.operator(), "order")
  assert_eq(compound.parts().length(), 2)
}

///|
test "compound_location_start_end" {
  let loc1 = @src.FeatureLocation::new_exact(10, 50)
  let loc2 = @src.FeatureLocation::new_exact(100, 200)
  let compound = @src.CompoundLocation::join([loc1, loc2])
  let start = compound.compound_start()
  let end = compound.compound_end()
  assert_eq(start.to_int(), 10)
  assert_eq(end.to_int(), 200)
}

///|
test "compound_location_start_end_empty" {
  let compound = @src.CompoundLocation::join([])
  let start = compound.compound_start()
  let end = compound.compound_end()
  assert_true(start is @src.Position::UnknownPosition)
  assert_true(end is @src.Position::UnknownPosition)
}

///|
test "compound_location_to_string_join" {
  let loc1 = @src.FeatureLocation::new_exact(1, 100)
  let loc2 = @src.FeatureLocation::new_exact(200, 300)
  let compound = @src.CompoundLocation::join([loc1, loc2])
  let s = compound.compound_to_string()
  assert_eq(s, "join(1..100,200..300)")
}

///|
test "compound_location_to_string_order" {
  let loc1 = @src.FeatureLocation::new_exact(1, 100)
  let loc2 = @src.FeatureLocation::new_exact(200, 300)
  let compound = @src.CompoundLocation::order([loc1, loc2])
  let s = compound.compound_to_string()
  assert_eq(s, "order(1..100,200..300)")
}

///|
test "compound_location_to_string_fuzzy" {
  let loc1 = @src.FeatureLocation::new_safe(
    @src.Position::new_before(1),
    @src.Position::new_exact(100),
  )
  let loc2 = @src.FeatureLocation::new_safe(
    @src.Position::new_exact(200),
    @src.Position::new_after(300),
  )
  let compound = @src.CompoundLocation::join([loc1, loc2])
  let s = compound.compound_to_string()
  assert_eq(s, "join(<1..100,200..>300)")
}

///|
test "compound_location_to_string_with_strand" {
  let loc1 = @src.FeatureLocation::new_exact(1, 100, strand=Some(-1))
  let loc2 = @src.FeatureLocation::new_exact(200, 300, strand=Some(-1))
  let compound = @src.CompoundLocation::join([loc1, loc2])
  let s = compound.compound_to_string()
  assert_true(s.contains("complement(join("))
  assert_true(s.contains("1..100"))
  assert_true(s.contains("200..300"))
}

///|
test "parse_position_string_exact" {
  let pos = @src.parse_position_string("100")
  assert_eq(pos.to_int(), 100)
}

///|
test "parse_position_string_before" {
  let pos = @src.parse_position_string("<50")
  assert_true(pos is @src.Position::BeforePosition(_))
  assert_eq(pos.to_int(), 50)
}

///|
test "parse_position_string_after" {
  let pos = @src.parse_position_string(">200")
  assert_true(pos is @src.Position::AfterPosition(_))
  assert_eq(pos.to_int(), 200)
}

///|
test "parse_position_string_uncertain" {
  let pos = @src.parse_position_string("?75?")
  assert_true(pos is @src.Position::UncertainPosition(_))
  assert_eq(pos.to_int(), 75)
}

///|
test "parse_position_string_unknown" {
  let pos = @src.parse_position_string("?")
  assert_true(pos is @src.Position::UnknownPosition)
}

///|
test "parse_position_string_empty_raises" {
  try {
    ignore(@src.parse_position_string(""))
    fail("expected LocationParserError for empty position string")
  } catch {
    _ => ()
  }
}

///|
test "parse_simple_location_basic" {
  let loc = @src.parse_simple_location("1..100")
  assert_eq(loc.start_int(), 1)
  assert_eq(loc.end_int(), 100)
  assert_true(loc.strand() is None)
}

///|
test "parse_simple_location_fuzzy_start" {
  let loc = @src.parse_simple_location("<1..100")
  assert_true(loc.start() is @src.Position::BeforePosition(_))
  assert_eq(loc.start_int(), 1)
  assert_eq(loc.end_int(), 100)
}

///|
test "parse_simple_location_fuzzy_end" {
  let loc = @src.parse_simple_location("1..>100")
  assert_eq(loc.start_int(), 1)
  assert_true(loc.end() is @src.Position::AfterPosition(_))
  assert_eq(loc.end_int(), 100)
}

///|
test "parse_simple_location_both_fuzzy" {
  let loc = @src.parse_simple_location("<1..>100")
  assert_true(loc.start() is @src.Position::BeforePosition(_))
  assert_true(loc.end() is @src.Position::AfterPosition(_))
  assert_eq(loc.start_int(), 1)
  assert_eq(loc.end_int(), 100)
}

///|
test "parse_location_string_simple" {
  let parsed = @src.parse_location_string("1..100")
  assert_true(parsed.is_feature_location())
  let loc = parsed.to_feature_location()
  assert_true(loc is Some(_))
  assert_eq(loc.unwrap().start_int(), 1)
  assert_eq(loc.unwrap().end_int(), 100)
}

///|
test "parse_location_string_fuzzy" {
  let parsed = @src.parse_location_string("<1..100")
  assert_true(parsed.is_feature_location())
  let loc = parsed.to_feature_location()
  assert_true(loc is Some(_))
  assert_true(loc.unwrap().start() is @src.Position::BeforePosition(_))
}

///|
test "parse_location_string_complement" {
  let parsed = @src.parse_location_string("complement(1..100)")
  assert_true(parsed.is_feature_location())
  let loc = parsed.to_feature_location()
  assert_true(loc is Some(_))
  assert_eq(loc.unwrap().strand(), Some(-1))
  assert_eq(loc.unwrap().start_int(), 1)
  assert_eq(loc.unwrap().end_int(), 100)
}

///|
test "parse_location_string_join" {
  let parsed = @src.parse_location_string("join(1..100,200..300)")
  assert_true(parsed.is_compound_location())
  let comp = parsed.to_compound_location()
  assert_true(comp is Some(_))
  assert_eq(comp.unwrap().operator(), "join")
  assert_eq(comp.unwrap().parts().length(), 2)
}

///|
test "parse_location_string_order" {
  let parsed = @src.parse_location_string("order(1..100,200..300)")
  assert_true(parsed.is_compound_location())
  let comp = parsed.to_compound_location()
  assert_true(comp is Some(_))
  assert_eq(comp.unwrap().operator(), "order")
  assert_eq(comp.unwrap().parts().length(), 2)
}

///|
test "parse_location_string_nested_complement_join" {
  let parsed = @src.parse_location_string("complement(join(1..100,200..300))")
  assert_true(parsed.is_compound_location())
  let comp = parsed.to_compound_location()
  assert_true(comp is Some(_))
  assert_eq(comp.unwrap().operator(), "join")
  assert_eq(comp.unwrap().parts().length(), 2)
  for part in comp.unwrap().parts() {
    assert_eq(part.strand(), Some(-1))
  }
}

///|
test "parse_location_string_empty_raises" {
  try {
    ignore(@src.parse_location_string(""))
    fail("expected LocationParserError for empty location string")
  } catch {
    _ => ()
  }
}

///|
test "parse_location_string_invalid_format_raises" {
  try {
    ignore(@src.parse_location_string("not_valid"))
    fail("expected LocationParserError for invalid format")
  } catch {
    _ => ()
  }
}

///|
test "parsed_location_helpers" {
  let parsed = @src.parse_location_string("1..100")
  assert_true(parsed.is_feature_location())
  assert_false(parsed.is_compound_location())
  assert_true(parsed.to_feature_location() is Some(_))
  assert_true(parsed.to_compound_location() is None)
}

///|
test "parsed_location_compound_helpers" {
  let parsed = @src.parse_location_string("join(1..100,200..300)")
  assert_true(parsed.is_compound_location())
  assert_false(parsed.is_feature_location())
  assert_true(parsed.to_compound_location() is Some(_))
  assert_true(parsed.to_feature_location() is None)
}

///|
test "seq_feature_extended_creation" {
  let loc = @src.FeatureLocation::new_exact(10, 200)
  let feat = @src.SeqFeatureExtended::new(
    loc,
    "gene",
    id="gene1",
    description="Test gene",
  )
  assert_eq(feat.type_(), "gene")
  assert_eq(feat.id(), "gene1")
  assert_eq(feat.description(), "Test gene")
  assert_eq(feat.location().start_int(), 10)
  assert_eq(feat.location().end_int(), 200)
}

///|
test "seq_feature_extended_defaults" {
  let loc = @src.FeatureLocation::new_exact(0, 100)
  let feat = @src.SeqFeatureExtended::new(loc, "CDS")
  assert_eq(feat.id(), "")
  assert_eq(feat.description(), "")
  assert_eq(feat.qualifiers().length(), 0)
}

///|
test "seq_feature_extended_add_qualifier" {
  let loc = @src.FeatureLocation::new_exact(0, 100)
  let feat = @src.SeqFeatureExtended::new(loc, "gene")
  let feat2 = feat.add_qualifier("gene", "BRCA1")
  assert_eq(feat2.get_qualifier("gene").length(), 1)
  assert_eq(feat2.get_qualifier("gene")[0], "BRCA1")
}

///|
test "seq_feature_extended_add_multiple_qualifiers" {
  let loc = @src.FeatureLocation::new_exact(0, 100)
  let feat = @src.SeqFeatureExtended::new(loc, "gene")
  let feat2 = feat
    .add_qualifier("gene", "BRCA1")
    .add_qualifier("gene", "BRCA2")
    .add_qualifier("chromosome", "17")
  assert_eq(feat2.get_qualifier("gene").length(), 2)
  assert_eq(feat2.get_qualifier("gene")[0], "BRCA1")
  assert_eq(feat2.get_qualifier("gene")[1], "BRCA2")
  assert_eq(feat2.get_qualifier("chromosome").length(), 1)
  assert_eq(feat2.get_qualifier("chromosome")[0], "17")
}

///|
test "seq_feature_extended_get_qualifier_missing" {
  let loc = @src.FeatureLocation::new_exact(0, 100)
  let feat = @src.SeqFeatureExtended::new(loc, "gene")
  let vals = feat.get_qualifier("missing")
  assert_eq(vals.length(), 0)
}

///|
test "seq_feature_extended_to_string" {
  let loc = @src.FeatureLocation::new_exact(1, 100)
  let feat = @src.SeqFeatureExtended::new(loc, "gene")
  let feat2 = feat.add_qualifier("gene", "BRCA1")
  let s = feat2.seq_feature_to_string()
  assert_true(s.contains("gene 1..100"))
  assert_true(s.contains("/gene=\"BRCA1\""))
}

///|
test "location_to_string_simple" {
  let loc = @src.FeatureLocation::new_exact(1, 100)
  let s = @src.location_to_string(loc)
  assert_eq(s, "1..100")
}

///|
test "location_to_string_fuzzy" {
  let loc = @src.FeatureLocation::new_safe(
    @src.Position::new_before(1),
    @src.Position::new_after(100),
  )
  let s = @src.location_to_string(loc)
  assert_eq(s, "<1..>100")
}

///|
test "location_to_string_complement" {
  let loc = @src.FeatureLocation::new_exact(1, 100, strand=Some(-1))
  let s = @src.location_to_string(loc)
  assert_eq(s, "complement(1..100)")
}

///|
test "feature_location_complement_swap" {
  let loc_plus = @src.FeatureLocation::new_exact(1, 100, strand=Some(1))
  let loc_minus = @src.feature_location_complement(loc_plus)
  assert_eq(loc_minus.strand(), Some(-1))
  assert_eq(loc_minus.start_int(), 1)
  assert_eq(loc_minus.end_int(), 100)
}

///|
test "feature_location_complement_minus_to_plus" {
  let loc_minus = @src.FeatureLocation::new_exact(1, 100, strand=Some(-1))
  let loc_plus = @src.feature_location_complement(loc_minus)
  assert_eq(loc_plus.strand(), Some(1))
}

///|
test "feature_location_complement_none_strand" {
  let loc = @src.FeatureLocation::new_safe(
    @src.Position::new_exact(1),
    @src.Position::new_exact(100),
  )
  let complemented = @src.feature_location_complement(loc)
  assert_true(complemented.strand() is None)
}

///|
test "feature_location_shift_positive" {
  let loc = @src.FeatureLocation::new_exact(10, 50)
  let shifted = @src.feature_location_shift(loc, 100)
  assert_eq(shifted.start_int(), 110)
  assert_eq(shifted.end_int(), 150)
}

///|
test "feature_location_shift_negative" {
  let loc = @src.FeatureLocation::new_exact(100, 200)
  let shifted = @src.feature_location_shift(loc, -50)
  assert_eq(shifted.start_int(), 50)
  assert_eq(shifted.end_int(), 150)
}

///|
test "feature_location_shift_preserves_strand" {
  let loc = @src.FeatureLocation::new_exact(10, 50, strand=Some(-1))
  let shifted = @src.feature_location_shift(loc, 10)
  assert_eq(shifted.strand(), Some(-1))
  assert_eq(shifted.start_int(), 20)
  assert_eq(shifted.end_int(), 60)
}

///|
test "feature_location_shift_fuzzy" {
  let loc = @src.FeatureLocation::new_safe(
    @src.Position::new_before(10),
    @src.Position::new_after(50),
  )
  let shifted = @src.feature_location_shift(loc, 5)
  assert_true(shifted.start() is @src.Position::BeforePosition(_))
  assert_true(shifted.end() is @src.Position::AfterPosition(_))
  assert_eq(shifted.start_int(), 15)
  assert_eq(shifted.end_int(), 55)
}

///|
test "merge_locations_overlapping" {
  let loc1 = @src.FeatureLocation::new_exact(1, 50)
  let loc2 = @src.FeatureLocation::new_exact(30, 100)
  let merged = @src.merge_locations(loc1, loc2)
  assert_eq(merged.start_int(), 1)
  assert_eq(merged.end_int(), 100)
}

///|
test "merge_locations_adjacent" {
  let loc1 = @src.FeatureLocation::new_exact(1, 50)
  let loc2 = @src.FeatureLocation::new_exact(51, 100)
  let merged = @src.merge_locations(loc1, loc2)
  assert_eq(merged.start_int(), 1)
  assert_eq(merged.end_int(), 100)
}

///|
test "merge_locations_same_strand" {
  let loc1 = @src.FeatureLocation::new_exact(1, 50, strand=Some(1))
  let loc2 = @src.FeatureLocation::new_exact(30, 100, strand=Some(1))
  let merged = @src.merge_locations(loc1, loc2)
  assert_eq(merged.strand(), Some(1))
}

///|
test "merge_locations_incompatible_strand_raises" {
  let loc1 = @src.FeatureLocation::new_exact(1, 50, strand=Some(1))
  let loc2 = @src.FeatureLocation::new_exact(30, 100, strand=Some(-1))
  try {
    ignore(@src.merge_locations(loc1, loc2))
    fail("expected LocationParserError for incompatible strands")
  } catch {
    _ => ()
  }
}

///|
test "merge_locations_one_strand_none" {
  let loc1 = @src.FeatureLocation::new_exact(1, 50, strand=Some(1))
  let loc2 = @src.FeatureLocation::new_exact(30, 100)
  let merged = @src.merge_locations(loc1, loc2)
  assert_eq(merged.strand(), Some(1))
}

///|
test "compound_location_to_string_empty_parts" {
  let compound = @src.CompoundLocation::join([])
  let s = compound.compound_to_string()
  assert_eq(s, "join()")
}

///|
test "parse_location_string_with_whitespace" {
  let parsed = @src.parse_location_string("  1..100  ")
  assert_true(parsed.is_feature_location())
  let loc = parsed.to_feature_location()
  assert_true(loc is Some(_))
  assert_eq(loc.unwrap().start_int(), 1)
  assert_eq(loc.unwrap().end_int(), 100)
}

///|
test "parse_simple_location_with_whitespace" {
  let loc = @src.parse_simple_location("  <1..100  ")
  assert_true(loc.start() is @src.Position::BeforePosition(_))
  assert_eq(loc.start_int(), 1)
}

///|
test "compound_location_to_string_multiple_parts" {
  let loc1 = @src.FeatureLocation::new_exact(1, 100)
  let loc2 = @src.FeatureLocation::new_exact(200, 300)
  let loc3 = @src.FeatureLocation::new_exact(400, 500)
  let compound = @src.CompoundLocation::join([loc1, loc2, loc3])
  let s = compound.compound_to_string()
  assert_eq(s, "join(1..100,200..300,400..500)")
}

///|
test "location_to_string_with_strand_minus" {
  let loc = @src.FeatureLocation::new_exact(50, 200, strand=Some(-1))
  let s = @src.location_to_string(loc)
  assert_eq(s, "complement(50..200)")
}

///|
test "parse_complement_location_simple" {
  let (parsed, strand) = @src.parse_complement_location("complement(1..100)")
  assert_eq(strand, -1)
  assert_true(parsed.is_feature_location())
  let loc = parsed.to_feature_location()
  assert_true(loc is Some(_))
  assert_eq(loc.unwrap().start_int(), 1)
  assert_eq(loc.unwrap().end_int(), 100)
}

///|
test "parse_complement_location_compound" {
  let (parsed, strand) = @src.parse_complement_location(
    "complement(join(1..100,200..300))",
  )
  assert_eq(strand, -1)
  assert_true(parsed.is_compound_location())
}

///|
test "parse_complement_location_missing_prefix_raises" {
  try {
    ignore(@src.parse_complement_location("1..100"))
    fail("expected LocationParserError for missing complement prefix")
  } catch {
    _ => ()
  }
}

///|
test "parse_complement_location_missing_paren_raises" {
  try {
    ignore(@src.parse_complement_location("complement(1..100"))
    fail("expected LocationParserError for missing closing paren")
  } catch {
    _ => ()
  }
}

///|
test "parse_join_location_basic" {
  let compound = @src.parse_join_location("join(1..100,200..300)")
  assert_eq(compound.operator(), "join")
  assert_eq(compound.parts().length(), 2)
}

///|
test "parse_join_location_order" {
  let compound = @src.parse_join_location("order(1..100,200..300)")
  assert_eq(compound.operator(), "order")
  assert_eq(compound.parts().length(), 2)
}

///|
test "parse_join_location_missing_paren_raises" {
  try {
    ignore(@src.parse_join_location("join(1..100,200..300"))
    fail("expected LocationParserError for missing closing paren")
  } catch {
    _ => ()
  }
}

///|
test "parse_join_location_invalid_prefix_raises" {
  try {
    ignore(@src.parse_join_location("invalid(1..100)"))
    fail("expected LocationParserError for invalid operator")
  } catch {
    _ => ()
  }
}

///|
test "seq_feature_extended_qualifiers_independent" {
  let loc = @src.FeatureLocation::new_exact(0, 100)
  let feat1 = @src.SeqFeatureExtended::new(loc, "gene")
  let feat2 = feat1.add_qualifier("gene", "BRCA1")
  assert_eq(feat1.qualifiers().length(), 0)
  assert_eq(feat2.qualifiers().length(), 1)
}