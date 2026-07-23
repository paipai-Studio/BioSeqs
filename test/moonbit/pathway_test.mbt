///|
test "Species_new" {
  let species = @src.Species::new("glucose", "Glucose")
  assert_eq(species.id, "glucose")
  assert_eq(species.name, "Glucose")
}

///|
test "Reaction_new" {
  let reactants : Array[String] = Array::new()
  reactants.push("glucose")
  let products : Array[String] = Array::new()
  products.push("g6p")
  let reaction = @src.Reaction::new("r1", "Hexokinase", reactants, products, false)
  assert_eq(reaction.id, "r1")
  assert_eq(reaction.name, "Hexokinase")
  assert_eq(reaction.reactants.length(), 1)
  assert_eq(reaction.products.length(), 1)
  assert_eq(reaction.reversible, false)
}

///|
test "Pathway_new" {
  let species : Array[@src.Species] = Array::new()
  species.push(@src.Species::new("glucose", "Glucose"))
  
  let reactants : Array[String] = Array::new()
  reactants.push("glucose")
  let products : Array[String] = Array::new()
  products.push("g6p")
  let reactions : Array[@src.Reaction] = Array::new()
  reactions.push(@src.Reaction::new("r1", "Hexokinase", reactants, products, false))
  
  let pathway = @src.Pathway::new("test", "Test Pathway", species, reactions)
  assert_eq(pathway.id, "test")
  assert_eq(pathway.name, "Test Pathway")
  assert_eq(pathway.species.length(), 1)
  assert_eq(pathway.reactions.length(), 1)
}

///|
test "pathway_find_species" {
  let pathway = @src.create_example_pathway()
  let species = @src.pathway_find_species(pathway, "glucose")
  assert_true(species.is_some())
  assert_eq(species.unwrap().name, "Glucose")
}

///|
test "pathway_find_species_not_found" {
  let pathway = @src.create_example_pathway()
  let species = @src.pathway_find_species(pathway, "nonexistent")
  assert_true(species.is_none())
}

///|
test "pathway_find_reactions_with_species" {
  let pathway = @src.create_example_pathway()
  let reactions = @src.pathway_find_reactions_with_species(pathway, "g6p")
  assert_true(reactions.length() >= 1)
}

///|
test "pathway_count_species" {
  let pathway = @src.create_example_pathway()
  let count = @src.pathway_count_species(pathway)
  assert_eq(count, 6)
}

///|
test "pathway_count_reactions" {
  let pathway = @src.create_example_pathway()
  let count = @src.pathway_count_reactions(pathway)
  assert_eq(count, 5)
}

///|
test "create_example_pathway" {
  let pathway = @src.create_example_pathway()
  assert_eq(pathway.id, "glycolysis")
  assert_eq(pathway.name, "Glycolysis")
}

///|
test "Reaction_to_string" {
  let reactants : Array[String] = Array::new()
  reactants.push("A")
  let products : Array[String] = Array::new()
  products.push("B")
  let reaction = @src.Reaction::new("r1", "Test", reactants, products, true)
  let str = @src.Reaction::to_string(reaction)
  assert_true(str.length() > 0)
}