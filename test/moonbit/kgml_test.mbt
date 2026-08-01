///|
/// Tests for Bio.KEGG.KGML pathway map parser.

// ---------------------------------------------------------------------------
// Data structure creation and accessors
// ---------------------------------------------------------------------------

test "kgml_graphics_creation" {
  let g = @src.KgmlGraphics::new(
    100.0, 200.0, 46.0, 17.0, "rect", "K00844", "#000000", "#BFFFBF",
  )
  assert_eq(g.x(), 100.0)
  assert_eq(g.y(), 200.0)
  assert_eq(g.width(), 46.0)
  assert_eq(g.height(), 17.0)
  assert_eq(g.gtype(), "rect")
  assert_eq(g.name(), "K00844")
  assert_eq(g.fgcolor(), "#000000")
  assert_eq(g.bgcolor(), "#BFFFBF")
}

test "kgml_entry_creation" {
  let e = @src.KgmlEntry::new(1, "ko:K00844", "gene")
  assert_eq(e.id(), 1)
  assert_eq(e.name(), "ko:K00844")
  assert_eq(e.etype(), "gene")
  assert_eq(e.reaction(), "")
  assert_eq(e.link(), "")
  assert_true(e.graphics().is_empty())
  assert_eq(e.components().length(), 0)
}

test "kgml_subtype_creation" {
  let st = @src.KgmlSubType::new("compound", "C00118")
  assert_eq(st.name(), "compound")
  assert_eq(st.value(), "C00118")
}

test "kgml_relation_creation" {
  let r = @src.KgmlRelation::new(1, 2, "ECrel")
  assert_eq(r.entry1(), 1)
  assert_eq(r.entry2(), 2)
  assert_eq(r.rtype(), "ECrel")
  assert_eq(r.subtypes().length(), 0)
}

test "kgml_reaction_creation" {
  let rxn = @src.KgmlReaction::new("rn:R01786", "irreversible")
  assert_eq(rxn.name(), "rn:R01786")
  assert_eq(rxn.rtype(), "irreversible")
  assert_eq(rxn.substrates().length(), 0)
  assert_eq(rxn.products().length(), 0)
}

test "kgml_pathway_creation" {
  let p = @src.KgmlPathway::new("path:ko00010", "ko", "00010", "Glycolysis")
  assert_eq(p.name(), "path:ko00010")
  assert_eq(p.org(), "ko")
  assert_eq(p.number(), "00010")
  assert_eq(p.title(), "Glycolysis")
  assert_eq(p.image(), "")
  assert_eq(p.link(), "")
  assert_eq(p.entries().length(), 0)
  assert_eq(p.relations().length(), 0)
  assert_eq(p.reactions().length(), 0)
}

// ---------------------------------------------------------------------------
// XML parsing
// ---------------------------------------------------------------------------

test "kgml_parse_basic_pathway" {
  let xml = @src.kgml_sample_pathway()
  let pw = @src.parse_kgml(xml)
  assert_true(pw.is_some())
  let p = pw.unwrap()
  assert_eq(p.name(), "path:ko00010")
  assert_eq(p.org(), "ko")
  assert_eq(p.number(), "00010")
  assert_eq(p.title(), "Glycolysis / Gluconeogenesis")
  assert_true(p.image().length() > 0)
  assert_true(p.link().length() > 0)
}

test "kgml_parse_entries" {
  let xml = @src.kgml_sample_pathway()
  let p = @src.parse_kgml(xml).unwrap()
  assert_eq(p.entries().length(), 5)
  // Entry 1: gene K00844
  let e1 = p.entries()[0]
  assert_eq(e1.id(), 1)
  assert_eq(e1.name(), "ko:K00844")
  assert_eq(e1.etype(), "gene")
  assert_eq(e1.reaction(), "rn:R01786")
  // Entry 3: compound
  let e3 = p.entries()[2]
  assert_eq(e3.id(), 3)
  assert_eq(e3.name(), "cpd:C00118")
  assert_eq(e3.etype(), "compound")
  // Entry 5: map
  let e5 = p.entries()[4]
  assert_eq(e5.id(), 5)
  assert_eq(e5.etype(), "map")
}

test "kgml_parse_entry_graphics" {
  let xml = @src.kgml_sample_pathway()
  let p = @src.parse_kgml(xml).unwrap()
  let e1 = p.entries()[0]
  let g = e1.graphics()
  assert_true(g.is_some())
  let gfx = g.unwrap()
  assert_eq(gfx.x(), 100.0)
  assert_eq(gfx.y(), 200.0)
  assert_eq(gfx.width(), 46.0)
  assert_eq(gfx.height(), 17.0)
  assert_eq(gfx.gtype(), "rect")
  assert_eq(gfx.name(), "K00844")
  assert_eq(gfx.fgcolor(), "#000000")
  assert_eq(gfx.bgcolor(), "#BFFFBF")
}

test "kgml_parse_relations" {
  let xml = @src.kgml_sample_pathway()
  let p = @src.parse_kgml(xml).unwrap()
  assert_eq(p.relations().length(), 3)
  let r1 = p.relations()[0]
  assert_eq(r1.entry1(), 1)
  assert_eq(r1.entry2(), 2)
  assert_eq(r1.rtype(), "ECrel")
  assert_eq(r1.subtypes().length(), 1)
  assert_eq(r1.subtypes()[0].name(), "compound")
  assert_eq(r1.subtypes()[0].value(), "C00118")
}

test "kgml_parse_reactions" {
  let xml = @src.kgml_sample_pathway()
  let p = @src.parse_kgml(xml).unwrap()
  assert_eq(p.reactions().length(), 2)
  let r1 = p.reactions()[0]
  assert_eq(r1.name(), "rn:R01786")
  assert_eq(r1.rtype(), "irreversible")
  assert_eq(r1.substrates().length(), 1)
  assert_eq(r1.substrates()[0], "cpd:C00118")
  assert_eq(r1.products().length(), 1)
  assert_eq(r1.products()[0], "cpd:C00197")
  // Second reaction has multiple substrates and products
  let r2 = p.reactions()[1]
  assert_eq(r2.name(), "rn:R01070")
  assert_eq(r2.rtype(), "reversible")
  assert_eq(r2.substrates().length(), 2)
  assert_eq(r2.products().length(), 2)
}

test "kgml_parse_invalid_returns_none" {
  let result = @src.parse_kgml("not a kgml file")
  assert_true(result.is_none())
}

test "kgml_parse_empty_returns_none" {
  let result = @src.parse_kgml("")
  assert_true(result.is_none())
}

// ---------------------------------------------------------------------------
// Query methods
// ---------------------------------------------------------------------------

test "kgml_get_entry_by_id" {
  let xml = @src.kgml_sample_pathway()
  let p = @src.parse_kgml(xml).unwrap()
  let e = p.get_entry(3)
  assert_true(e.is_some())
  assert_eq(e.unwrap().name(), "cpd:C00118")
  let none = p.get_entry(999)
  assert_true(none.is_none())
}

test "kgml_get_entries_by_type" {
  let xml = @src.kgml_sample_pathway()
  let p = @src.parse_kgml(xml).unwrap()
  let genes = p.get_entries_by_type("gene")
  assert_eq(genes.length(), 3)
  let compounds = p.get_entries_by_type("compound")
  assert_eq(compounds.length(), 1)
  let maps = p.get_entries_by_type("map")
  assert_eq(maps.length(), 1)
}

test "kgml_get_relations_for_entry" {
  let xml = @src.kgml_sample_pathway()
  let p = @src.parse_kgml(xml).unwrap()
  // Entry 1 is in relation 1 (1->2)
  let r1 = p.get_relations_for_entry(1)
  assert_eq(r1.length(), 1)
  // Entry 2 is in relation 1 (1->2) and relation 2 (2->4)
  let r2 = p.get_relations_for_entry(2)
  assert_eq(r2.length(), 2)
  // Entry 4 is in relation 2 (2->4) and relation 3 (4->5)
  let r4 = p.get_relations_for_entry(4)
  assert_eq(r4.length(), 2)
}

// ---------------------------------------------------------------------------
// Serialization
// ---------------------------------------------------------------------------

test "kgml_to_string_roundtrip" {
  let xml = @src.kgml_sample_pathway()
  let p = @src.parse_kgml(xml).unwrap()
  let serialized = @src.kgml_to_string(p)
  // Re-parse the serialized output
  let p2 = @src.parse_kgml(serialized)
  assert_true(p2.is_some())
  let p2u = p2.unwrap()
  assert_eq(p2u.name(), p.name())
  assert_eq(p2u.org(), p.org())
  assert_eq(p2u.number(), p.number())
  assert_eq(p2u.title(), p.title())
  assert_eq(p2u.entries().length(), p.entries().length())
  assert_eq(p2u.relations().length(), p.relations().length())
  assert_eq(p2u.reactions().length(), p.reactions().length())
}

test "kgml_to_string_contains_pathway_tag" {
  let xml = @src.kgml_sample_pathway()
  let p = @src.parse_kgml(xml).unwrap()
  let s = @src.kgml_to_string(p)
  assert_true(s.contains("<pathway"))
  assert_true(s.contains("name=\"path:ko00010\""))
  assert_true(s.contains("</pathway>"))
}

test "kgml_to_string_contains_entries" {
  let xml = @src.kgml_sample_pathway()
  let p = @src.parse_kgml(xml).unwrap()
  let s = @src.kgml_to_string(p)
  assert_true(s.contains("<entry"))
  assert_true(s.contains("ko:K00844"))
  assert_true(s.contains("cpd:C00118"))
}

test "kgml_to_string_contains_relations" {
  let xml = @src.kgml_sample_pathway()
  let p = @src.parse_kgml(xml).unwrap()
  let s = @src.kgml_to_string(p)
  assert_true(s.contains("<relation"))
  assert_true(s.contains("ECrel"))
  assert_true(s.contains("<subtype"))
}

test "kgml_to_string_contains_reactions" {
  let xml = @src.kgml_sample_pathway()
  let p = @src.parse_kgml(xml).unwrap()
  let s = @src.kgml_to_string(p)
  assert_true(s.contains("<reaction"))
  assert_true(s.contains("rn:R01786"))
  assert_true(s.contains("<substrate"))
  assert_true(s.contains("<product"))
}
