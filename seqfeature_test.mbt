///|
test "feature_location_basic" {
  let loc = FeatureLocation::new(
    Position::new_exact(5),
    Position::new_exact(10),
  )
  assert_eq(loc.start_int(), 5)
  assert_eq(loc.end_int(), 10)
  assert_true(loc.strand() is None)
}

///|
test "feature_location_strand" {
  let loc_plus = FeatureLocation::new_exact(7, 110, strand=Some(1))
  assert_eq(loc_plus.strand(), Some(1))
  let loc_minus = FeatureLocation::new_exact(9, 108, strand=Some(-1))
  assert_eq(loc_minus.strand(), Some(-1))
}

///|
test "seqfeature_basic" {
  let f1 = SeqFeature::new(
    location=Some(FeatureLocation::new_exact(5, 10)),
    feature_type="domain",
  )
  assert_eq(f1.feature_type(), "domain")
  assert_true(f1.strand() is None)
}

///|
test "seqfeature_extract" {
  let seq = Seq::new("MKQHKAMIVALIVICITAVVAAL")
  let f = SeqFeature::new(
    location=Some(FeatureLocation::new_exact(8, 15)),
    feature_type="domain",
  )
  let extracted = f.extract(seq)
  assert_eq(extracted.to_string(), "VALIVIC")
}

///|
test "seqfeature_extract_minus_strand" {
  let seq = Seq::new("ATGCCG")
  let f = SeqFeature::new(
    location=Some(FeatureLocation::new_exact(0, 6)),
    feature_type="CDS",
    strand=Some(-1),
  )
  let extracted = f.extract(seq)
  assert_eq(extracted.to_string(), "CGGCAT")
}

///|
test "seqfeature_qualifiers" {
  let f = SeqFeature::new(
    location=Some(FeatureLocation::new_exact(0, 10)),
    feature_type="gene",
  )
  let f2 = f.add_qualifier("gene", "abc")
  let qualifiers = f2.qualifiers()
  assert_true(qualifiers.get("gene") != None)
}

///|
test "compound_location" {
  let loc1 = FeatureLocation::new_exact(0, 10, strand=Some(1))
  let loc2 = FeatureLocation::new_exact(20, 30, strand=Some(1))
  let compound = CompoundLocation::new(
    [loc1, loc2],
    operator="join",
  )
  assert_eq(compound.operator(), "join")
  assert_eq(compound.parts().length(), 2)
}

///|
test "position_types" {
  let exact = Position::new_exact(10)
  assert_eq(exact.to_int(), 10)
  let within = Position::new_within(5, 5, 15)
  assert_eq(within.to_int(), 5)
}