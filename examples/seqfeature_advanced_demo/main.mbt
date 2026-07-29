///|
/// SeqFeature Advanced Demo
/// Demonstrates Bio.SeqFeatureAdvanced functionality:
///   - CompoundLocation with join and order operators
///   - GenBank-format string generation
///   - LocationParser for various location string formats
///   - SeqFeatureExtended with qualifiers
///   - Location utilities (complement, shift, merge)
///
/// Run: moon run examples/seqfeature_advanced_demo/

fn main {
  println("=== SeqFeature Advanced Demo ===")

  // =========================================================================
  // 1. CompoundLocation with join and order operators
  // =========================================================================
  println("")
  println("--- 1. CompoundLocation Constructors ---")

  let loc1 = @bio.FeatureLocation::new_safe(
    @bio.Position::new_exact(1),
    @bio.Position::new_exact(100),
  )
  let loc2 = @bio.FeatureLocation::new_safe(
    @bio.Position::new_exact(200),
    @bio.Position::new_exact(300),
  )
  let loc3 = @bio.FeatureLocation::new_safe(
    @bio.Position::new_exact(400),
    @bio.Position::new_exact(500),
  )

  println("Created three FeatureLocations:")
  println("  loc1: \{@bio.location_to_string(loc1)}")
  println("  loc2: \{@bio.location_to_string(loc2)}")
  println("  loc3: \{@bio.location_to_string(loc3)}")

  let join_loc = @bio.CompoundLocation::join([loc1, loc2, loc3])
  println("")
  println("Join compound location:")
  println("  operator: \{join_loc.operator()}")
  println("  parts count: \{join_loc.parts().length()}")
  println("  GenBank string: \{join_loc.compound_to_string()}")
  println("  start: \{join_loc.compound_start().repr()}")
  println("  end: \{join_loc.compound_end().repr()}")

  let order_loc = @bio.CompoundLocation::order([loc1, loc2])
  println("")
  println("Order compound location:")
  println("  operator: \{order_loc.operator()}")
  println("  GenBank string: \{order_loc.compound_to_string()}")

  // =========================================================================
  // 2. Compound locations with strand
  // =========================================================================
  println("")
  println("--- 2. CompoundLocation with Strand ---")

  let fwd1 = @bio.FeatureLocation::new_safe(
    @bio.Position::new_exact(1),
    @bio.Position::new_exact(100),
    strand=Some(1),
  )
  let fwd2 = @bio.FeatureLocation::new_safe(
    @bio.Position::new_exact(200),
    @bio.Position::new_exact(300),
    strand=Some(1),
  )
  let join_fwd = @bio.CompoundLocation::join([fwd1, fwd2])
  println("Forward strand join: \{join_fwd.compound_to_string()}")

  let rev1 = @bio.FeatureLocation::new_safe(
    @bio.Position::new_exact(1),
    @bio.Position::new_exact(100),
    strand=Some(-1),
  )
  let rev2 = @bio.FeatureLocation::new_safe(
    @bio.Position::new_exact(200),
    @bio.Position::new_exact(300),
    strand=Some(-1),
  )
  let join_rev = @bio.CompoundLocation::join([rev1, rev2])
  println("Reverse strand join: \{join_rev.compound_to_string()}")

  // =========================================================================
  // 3. LocationParser - parsing location strings
  // =========================================================================
  println("")
  println("--- 3. LocationParser ---")

  try {
    let parsed1 = @bio.parse_location_string("1..100")
    println("")
    println("Parsed: \"1..100\"")
    println("  is_feature_location: \{parsed1.is_feature_location()}")
    match parsed1.to_feature_location() {
      Some(loc) => {
        println("  start: \{loc.start().repr()}")
        println("  end: \{loc.end().repr()}")
        println("  strand: \{loc.strand()}")
      }
      None => println("  (unexpected compound location)")
    }
  } catch {
    _ => println("  Failed to parse '1..100'")
  }

  try {
    let parsed2 = @bio.parse_location_string("<1..100")
    println("")
    println("Parsed: \"<1..100\"")
    println("  is_feature_location: \{parsed2.is_feature_location()}")
    match parsed2.to_feature_location() {
      Some(loc) => {
        println("  start: \{loc.start().repr()}")
        println("  end: \{loc.end().repr()}")
      }
      None => println("  (unexpected compound location)")
    }
  } catch {
    _ => println("  Failed to parse '<1..100'")
  }

  try {
    let parsed3 = @bio.parse_location_string("complement(1..100)")
    println("")
    println("Parsed: \"complement(1..100)\"")
    println("  is_feature_location: \{parsed3.is_feature_location()}")
    match parsed3.to_feature_location() {
      Some(loc) => {
        println("  start: \{loc.start().repr()}")
        println("  end: \{loc.end().repr()}")
        println("  strand: \{loc.strand()}")
      }
      None => println("  (unexpected compound location)")
    }
  } catch {
    _ => println("  Failed to parse 'complement(1..100)'")
  }

  try {
    let parsed4 = @bio.parse_location_string("join(1..100,200..300)")
    println("")
    println("Parsed: \"join(1..100,200..300)\"")
    println("  is_compound_location: \{parsed4.is_compound_location()}")
    match parsed4.to_compound_location() {
      Some(comp) => {
        println("  operator: \{comp.operator()}")
        println("  parts count: \{comp.parts().length()}")
        println("  GenBank: \{comp.compound_to_string()}")
      }
      None => println("  (unexpected feature location)")
    }
  } catch {
    _ => println("  Failed to parse 'join(1..100,200..300)'")
  }

  try {
    let parsed5 = @bio.parse_location_string(
      "complement(join(1..100,200..300))",
    )
    println("")
    println("Parsed: \"complement(join(1..100,200..300))\"")
    println("  is_compound_location: \{parsed5.is_compound_location()}")
    match parsed5.to_compound_location() {
      Some(comp) => {
        println("  operator: \{comp.operator()}")
        println("  parts count: \{comp.parts().length()}")
        println("  GenBank: \{comp.compound_to_string()}")
      }
      None => println("  (unexpected feature location)")
    }
  } catch {
    _ => println("  Failed to parse 'complement(join(1..100,200..300))'")
  }

  // =========================================================================
  // 4. SeqFeatureExtended with qualifiers
  // =========================================================================
  println("")
  println("--- 4. SeqFeatureExtended ---")

  let gene_loc = @bio.FeatureLocation::new_safe(
    @bio.Position::new_exact(100),
    @bio.Position::new_exact(500),
    strand=Some(1),
  )
  let gene_feature = @bio.SeqFeatureExtended::new(
    gene_loc,
    "gene",
    id="gene_001",
    description="BRCA1 gene",
  )
  println("Created gene feature:")
  println("  type: \{gene_feature.type_()}")
  println("  id: \{gene_feature.id()}")
  println("  description: \{gene_feature.description()}")
  println("  location: \{@bio.location_to_string(gene_feature.location())}")

  let with_qualifiers = gene_feature
    .add_qualifier("gene", "BRCA1")
    .add_qualifier("product", "BRCA1 protein")
    .add_qualifier("note", "Breast cancer susceptibility")
    .add_qualifier("note", "DNA repair gene")
  println("")
  println("After adding qualifiers:")
  println("  qualifiers count: \{with_qualifiers.qualifiers().length()}")

  let notes = with_qualifiers.get_qualifier("note")
  println("  'note' values (\{notes.length()}):")
  for note in notes {
    println("    - \{note}")
  }

  println("")
  println("GenBank format:")
  println(with_qualifiers.seq_feature_to_string())

  // =========================================================================
  // 5. Location utilities
  // =========================================================================
  println("")
  println("--- 5. Location Utilities ---")

  let test_loc = @bio.FeatureLocation::new_safe(
    @bio.Position::new_exact(1),
    @bio.Position::new_exact(100),
  )
  println("location_to_string(test_loc): \{@bio.location_to_string(test_loc)}")

  let fuzzy_loc = @bio.FeatureLocation::new_safe(
    @bio.Position::new_before(1),
    @bio.Position::new_exact(100),
  )
  println("location_to_string(fuzzy): \{@bio.location_to_string(fuzzy_loc)}")

  let comp_loc = @bio.feature_location_complement(test_loc)
  println("")
  println("Complement of test_loc:")
  println("  original strand: \{test_loc.strand()}")
  println("  complement strand: \{comp_loc.strand()}")
  println("  complement string: \{@bio.location_to_string(comp_loc)}")

  let shifted = @bio.feature_location_shift(test_loc, 50)
  println("")
  println("Shift test_loc by +50:")
  println("  original: \{@bio.location_to_string(test_loc)}")
  println("  shifted: \{@bio.location_to_string(shifted)}")

  let shifted_back = @bio.feature_location_shift(test_loc, -10)
  println("Shift test_loc by -10:")
  println("  shifted: \{@bio.location_to_string(shifted_back)}")

  // Merge
  println("")
  println("Merge overlapping locations:")
  let merge_a = @bio.FeatureLocation::new_safe(
    @bio.Position::new_exact(1),
    @bio.Position::new_exact(200),
    strand=Some(1),
  )
  let merge_b = @bio.FeatureLocation::new_safe(
    @bio.Position::new_exact(150),
    @bio.Position::new_exact(300),
    strand=Some(1),
  )
  println("  loc_a: \{@bio.location_to_string(merge_a)}")
  println("  loc_b: \{@bio.location_to_string(merge_b)}")
  try {
    let merged = @bio.merge_locations(merge_a, merge_b)
    println("  merged: \{@bio.location_to_string(merged)}")
  } catch {
    _ => println("  Failed to merge locations")
  }

  println("")
  println("=== Demo complete ===")
}