///|
test "compound_new" {
  let c = @src.Compound::new("C0001", "Glucose")
  assert_eq(c.get_id(), "C0001")
  assert_eq(c.get_name(), "Glucose")
  assert_eq(c.get_formula(), "")
  assert_eq(c.get_charge(), 0)
}

///|
test "compound_with_chemical" {
  let c = @src.Compound::with_chemical("C0001", "Glucose", "C6H12O6", 0, "C(C1C(C(C(C(O1)O)O)O)O)O")
  assert_eq(c.get_id(), "C0001")
  assert_eq(c.get_formula(), "C6H12O6")
  assert_eq(c.get_charge(), 0)
  assert_eq(c.get_smiles(), "C(C1C(C(C(C(O1)O)O)O)O)O")
}

///|
test "compound_setters" {
  let c = @src.Compound::new("C0001", "Test")
  let c2 = c.set_formula("H2O")
  assert_eq(c2.get_formula(), "H2O")
  let c3 = c2.set_charge(1)
  assert_eq(c3.get_charge(), 1)
  let c4 = c3.set_smiles("O")
  assert_eq(c4.get_smiles(), "O")
}

///|
test "compound_aliases" {
  let c = @src.Compound::new("C0001", "Glucose")
  let c2 = c.add_alias("D-Glucose")
  let c3 = c2.add_alias("Dextrose")
  assert_eq(c3.get_aliases().length(), 2)
  assert_true(c3.has_alias("D-Glucose"))
  assert_true(c3.has_alias("dextrose")) // case insensitive
  assert_true(!c3.has_alias("Fructose"))
  assert_true(c3.has_alias("Glucose")) // matches name
}

///|
test "compound_pathways" {
  let c = @src.Compound::new("C0001", "Glucose")
  let c2 = c.add_pathway("Glycolysis")
  let c3 = c2.add_pathway("Pentose Phosphate")
  assert_eq(c3.get_pathways().length(), 2)
}

///|
test "compound_reactions" {
  let c = @src.Compound::new("C0001", "Glucose")
  let c2 = c.add_reaction("R0001")
  assert_eq(c2.get_reactions().length(), 1)
}

///|
test "compound_to_string" {
  let c = @src.Compound::with_chemical("C0001", "Glucose", "C6H12O6", 0, "CCO")
  let s = c.to_string()
  assert_true(s.contains("C0001"))
  assert_true(s.contains("Glucose"))
}

///|
test "parse_formula_simple" {
  let counts = @src.parse_formula("H2O")
  assert_eq(counts.get("H").unwrap(), 2)
  assert_eq(counts.get("O").unwrap(), 1)
}

///|
test "parse_formula_complex" {
  let counts = @src.parse_formula("C6H12O6")
  assert_eq(counts.get("C").unwrap(), 6)
  assert_eq(counts.get("H").unwrap(), 12)
  assert_eq(counts.get("O").unwrap(), 6)
}

///|
test "parse_formula_with_metals" {
  let counts = @src.parse_formula("Fe2O3")
  assert_eq(counts.get("Fe").unwrap(), 2)
  assert_eq(counts.get("O").unwrap(), 3)
}

///|
test "parse_formula_empty" {
  let counts = @src.parse_formula("")
  assert_eq(counts.size(), 0)
}

///|
test "compound_molecular_weight" {
  let mw = @src.compound_molecular_weight("H2O")
  let expected = 2.0 * 1.008 + 15.999
  assert_true((mw - expected).abs() < 0.01)
}

///|
test "compound_molecular_weight_glucose" {
  let mw = @src.compound_molecular_weight("C6H12O6")
  let expected = 6.0 * 12.011 + 12.0 * 1.008 + 6.0 * 15.999
  assert_true((mw - expected).abs() < 0.01)
}

///|
test "compound_molecular_weight_empty" {
  let mw = @src.compound_molecular_weight("")
  assert_true((mw - 0.0).abs() < 0.01)
}

///|
test "reaction_new" {
  let r = @src.ChemicalReaction::new("R0001", "Glucose uptake")
  assert_eq(r.get_id(), "R0001")
  assert_eq(r.get_name(), "Glucose uptake")
  assert_eq(r.get_direction(), "forward")
}

///|
test "reaction_substrates_products" {
  let r = @src.ChemicalReaction::new("R0001", "Test")
  let r2 = r.add_substrate("Glucose")
  let r3 = r2.add_substrate("ATP")
  let r4 = r3.add_product("G6P")
  assert_eq(r4.get_substrates().length(), 2)
  assert_eq(r4.get_products().length(), 1)
}

///|
test "reaction_enzymes" {
  let r = @src.ChemicalReaction::new("R0001", "Test")
  let r2 = r.add_enzyme("HK1")
  assert_eq(r2.get_enzymes().length(), 1)
}

///|
test "reaction_direction" {
  let r = @src.ChemicalReaction::new("R0001", "Test")
  let r2 = r.set_direction("reversible")
  assert_eq(r2.get_direction(), "reversible")
}

///|
test "reaction_equation_forward" {
  let r = @src.ChemicalReaction::new("R0001", "Test")
  let r2 = r.add_substrate("A").add_substrate("B").add_product("C")
  let eq = r2.get_equation()
  assert_true(eq.contains("A + B"))
  assert_true(eq.contains("->"))
  assert_true(eq.contains("C"))
}

///|
test "reaction_equation_reverse" {
  let r = @src.ChemicalReaction::new("R0001", "Test")
  let r2 = r.add_substrate("A").add_product("B").set_direction("reverse")
  let eq = r2.get_equation()
  assert_true(eq.contains("<-"))
}

///|
test "reaction_equation_reversible" {
  let r = @src.ChemicalReaction::new("R0001", "Test")
  let r2 = r.add_substrate("A").add_product("B").set_direction("reversible")
  let eq = r2.get_equation()
  assert_true(eq.contains("<=>"))
}

///|
test "reaction_to_string" {
  let r = @src.ChemicalReaction::new("R0001", "Test")
  let s = r.to_string()
  assert_true(s.contains("R0001"))
  assert_true(s.contains("Test"))
}

///|
test "pathway_map_new" {
  let p = @src.CompoundPathwayMap::new("Glycolysis")
  assert_eq(p.get_name(), "Glycolysis")
  assert_eq(p.get_compounds().length(), 0)
  assert_eq(p.get_reactions().length(), 0)
}

///|
test "pathway_map_add_compound" {
  let p = @src.CompoundPathwayMap::new("Glycolysis")
  let c = @src.Compound::with_chemical("C0003", "Glucose", "C6H12O6", 0, "CCO")
  let p2 = p.add_compound(c)
  assert_eq(p2.get_compounds().length(), 1)
}

///|
test "pathway_map_add_reaction" {
  let p = @src.CompoundPathwayMap::new("Glycolysis")
  let r = @src.ChemicalReaction::new("R0001", "HK1 reaction")
  let p2 = p.add_reaction(r)
  assert_eq(p2.get_reactions().length(), 1)
}

///|
test "pathway_map_find_compound" {
  let p = @src.CompoundPathwayMap::new("Glycolysis")
  let c = @src.Compound::with_chemical("C0003", "Glucose", "C6H12O6", 0, "CCO")
  let p2 = p.add_compound(c)
  let found = p2.find_compound("C0003")
  assert_true(found is Option::Some)
  match found {
    Some(comp) => assert_eq(comp.get_name(), "Glucose")
    None => assert_true(false)
  }
  let not_found = p2.find_compound("C9999")
  assert_true(not_found is Option::None)
}

///|
test "pathway_map_find_reaction" {
  let p = @src.CompoundPathwayMap::new("Glycolysis")
  let r = @src.ChemicalReaction::new("R0001", "HK1 reaction")
  let p2 = p.add_reaction(r)
  let found = p2.find_reaction("R0001")
  assert_true(found is Option::Some)
  let not_found = p2.find_reaction("R9999")
  assert_true(not_found is Option::None)
}

///|
test "pathway_map_summary" {
  let p = @src.CompoundPathwayMap::new("Glycolysis")
  let c = @src.Compound::with_chemical("C0003", "Glucose", "C6H12O6", 0, "CCO")
  let r = @src.ChemicalReaction::new("R0001", "HK1 reaction")
  let p2 = p.add_compound(c).add_reaction(r)
  let s = p2.summary()
  assert_true(s.contains("Glycolysis"))
  assert_true(s.contains("Glucose"))
}

///|
test "compound_templates" {
  let glucose = @src.compound_glucose()
  assert_eq(glucose.get_id(), "C0003")
  assert_eq(glucose.get_formula(), "C6H12O6")

  let fructose = @src.compound_fructose()
  assert_eq(fructose.get_id(), "C0005")

  let pyruvate = @src.compound_pyruvate()
  assert_eq(pyruvate.get_id(), "C0024")
  assert_eq(pyruvate.get_charge(), -1)

  let acetate = @src.compound_acetate()
  assert_eq(acetate.get_id(), "C0009")

  let citrate = @src.compound_citrate()
  assert_eq(citrate.get_id(), "C0015")
  assert_eq(citrate.get_charge(), -3)
}
