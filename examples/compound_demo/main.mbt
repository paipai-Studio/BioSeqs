///|
/// Bio.Compound module example - Chemical compound data structures.
///
/// Run: moon run examples/compound_demo/main.mbt
fn main {
  println("=== Bio.Compound Module Demo ===")
  println("")

  // Compound creation
  println("=== Compound Creation ===")
  let c1 = @bio.Compound::new("C0001", "Glucose")
  println("Compound: \{c1.to_string()}")
  println("ID: \{c1.get_id()}")
  println("Name: \{c1.get_name()}")
  println("Formula: '\{c1.get_formula()}'")
  println("Charge: \{c1.get_charge()}")
  println("")

  // Compound with chemical properties
  println("=== Compound With Chemical Properties ===")
  let c2 = @bio.Compound::with_chemical("C0001", "Glucose", "C6H12O6", 0, "C(C1C(C(C(C(O1)O)O)O)O)O")
  println("Compound: \{c2.to_string()}")
  println("Formula: \{c2.get_formula()}")
  println("SMILES: \{c2.get_smiles()}")
  println("")

  // Compound setters
  println("=== Compound Setters ===")
  let c3 = @bio.Compound::new("C0001", "Test")
  let c3b = c3.set_formula("H2O").set_charge(1).set_smiles("O")
  println("After setters: formula='\{c3b.get_formula()}', charge=\{c3b.get_charge()}, smiles='\{c3b.get_smiles()}'")
  println("")

  // Compound aliases
  println("=== Compound Aliases ===")
  let c4 = @bio.Compound::new("C0001", "Glucose")
  let c4b = c4.add_alias("D-Glucose").add_alias("Dextrose")
  println("Aliases: \{c4b.get_aliases()}")
  println("Has 'D-Glucose': \{c4b.has_alias("D-Glucose")}")
  println("Has 'dextrose' (case insensitive): \{c4b.has_alias("dextrose")}")
  println("Has 'Glucose' (matches name): \{c4b.has_alias("Glucose")}")
  println("Has 'Fructose': \{c4b.has_alias("Fructose")}")
  println("")

  // Compound pathways
  println("=== Compound Pathways ===")
  let c5 = @bio.Compound::new("C0001", "Glucose")
  let c5b = c5.add_pathway("Glycolysis").add_pathway("Pentose Phosphate")
  println("Pathways: \{c5b.get_pathways()}")
  println("")

  // Compound reactions
  println("=== Compound Reactions ===")
  let c6 = @bio.Compound::new("C0001", "Glucose")
  let c6b = c6.add_reaction("R0001")
  println("Reactions: \{c6b.get_reactions()}")
  println("")

  // Formula parsing
  println("=== Formula Parsing ===")
  let f1 = @bio.parse_formula("H2O")
  println("H2O: H=\{f1.get("H").unwrap()}, O=\{f1.get("O").unwrap()}")

  let f2 = @bio.parse_formula("C6H12O6")
  println("C6H12O6: C=\{f2.get("C").unwrap()}, H=\{f2.get("H").unwrap()}, O=\{f2.get("O").unwrap()}")

  let f3 = @bio.parse_formula("Fe2O3")
  println("Fe2O3: Fe=\{f3.get("Fe").unwrap()}, O=\{f3.get("O").unwrap()}")

  let f4 = @bio.parse_formula("")
  println("Empty formula: size=\{f4.size()}")
  println("")

  // Molecular weight
  println("=== Molecular Weight ===")
  let mw1 = @bio.compound_molecular_weight("H2O")
  println("H2O: \{mw1}")
  let mw2 = @bio.compound_molecular_weight("C6H12O6")
  println("C6H12O6: \{mw2}")
  let mw3 = @bio.compound_molecular_weight("")
  println("Empty: \{mw3}")
  println("")

  // Chemical reactions
  println("=== Chemical Reactions ===")
  let r1 = @bio.ChemicalReaction::new("R0001", "Glucose uptake")
  println("Reaction: \{r1.to_string()}")
  println("ID: \{r1.get_id()}")
  println("Name: \{r1.get_name()}")
  println("Direction: \{r1.get_direction()}")
  println("")

  let r2 = @bio.ChemicalReaction::new("R0001", "Test")
  let r2b = r2.add_substrate("A").add_substrate("B").add_product("C")
  println("Reaction equation: \{r2b.get_equation()}")
  println("Substrates: \{r2b.get_substrates()}")
  println("Products: \{r2b.get_products()}")
  println("")

  let r3 = @bio.ChemicalReaction::new("R0001", "Test")
  let r3b = r3.set_direction("reversible")
  println("Reversible direction: \{r3b.get_direction()}")
  println("Equation: \{r3b.get_equation()}")
  println("")

  let r4 = @bio.ChemicalReaction::new("R0001", "Test")
  let r4b = r4.set_direction("reverse")
  println("Reverse direction: \{r4b.get_direction()}")
  println("Equation: \{r4b.get_equation()}")
  println("")

  // Reaction enzymes
  println("=== Reaction Enzymes ===")
  let r5 = @bio.ChemicalReaction::new("R0001", "Test")
  let r5b = r5.add_enzyme("HK1")
  println("Enzymes: \{r5b.get_enzymes()}")
  println("")

  // Pathway map
  println("=== Pathway Map ===")
  let p1 = @bio.CompoundPathwayMap::new("Glycolysis")
  println("Pathway: \{p1.get_name()}")
  println("Compounds: \{p1.get_compounds().length()}")
  println("Reactions: \{p1.get_reactions().length()}")
  println("")

  let glucose = @bio.Compound::with_chemical("C0003", "Glucose", "C6H12O6", 0, "CCO")
  let r6 = @bio.ChemicalReaction::new("R0001", "HK1 reaction")
  let p2 = p1.add_compound(glucose).add_reaction(r6)
  println("After additions:")
  println("  Compounds: \{p2.get_compounds().length()}")
  println("  Reactions: \{p2.get_reactions().length()}")
  println("")

  // Find compound/reaction
  println("=== Find Compound/Reaction ===")
  let found_c = p2.find_compound("C0003")
  match found_c {
    Some(comp) => println("Found compound: \{comp.get_name()}")
    None => println("Compound not found")
  }

  let not_found_c = p2.find_compound("C9999")
  println("Non-existent compound found: \{not_found_c is Option::Some}")

  let found_r = p2.find_reaction("R0001")
  println("Reaction found: \{found_r is Option::Some}")
  let not_found_r = p2.find_reaction("R9999")
  println("Non-existent reaction found: \{not_found_r is Option::Some}")
  println("")

  // Pathway summary
  println("=== Pathway Summary ===")
  let summary = p2.summary()
  println(summary)

  // Compound templates
  println("=== Compound Templates ===")
  let t_glucose = @bio.compound_glucose()
  println("Glucose template: id=\{t_glucose.get_id()}, formula=\{t_glucose.get_formula()}")

  let t_fructose = @bio.compound_fructose()
  println("Fructose template: id=\{t_fructose.get_id()}")

  let t_pyruvate = @bio.compound_pyruvate()
  println("Pyruvate template: id=\{t_pyruvate.get_id()}, charge=\{t_pyruvate.get_charge()}")

  let t_citrate = @bio.compound_citrate()
  println("Citrate template: id=\{t_citrate.get_id()}, charge=\{t_citrate.get_charge()}")
  println("")

  println("=== Done ===")
}