///|
/// Pathway demo.
///
/// Run: moon run examples/pathway_demo/main.mbt
fn main {
  println("=== Pathway Demo ===\n")

  println("1. Create example pathway (Glycolysis):")
  let pathway = @bio.create_example_pathway()
  println("   Pathway ID: \{pathway.id}")
  println("   Pathway name: \{pathway.name}")

  println("\n2. Count species and reactions:")
  let species_count = @bio.pathway_count_species(pathway)
  let reaction_count = @bio.pathway_count_reactions(pathway)
  println("   Species count: \{species_count}")
  println("   Reaction count: \{reaction_count}")

  println("\n3. Find species:")
  let glucose = @bio.pathway_find_species(pathway, "glucose")
  match glucose {
    Some(s) => println("   Found: \{s.id} - \{s.name}")
    None => println("   Not found")
  }

  let g6p = @bio.pathway_find_species(pathway, "g6p")
  match g6p {
    Some(s) => println("   Found: \{s.id} - \{s.name}")
    None => println("   Not found")
  }

  println("\n4. Find reactions with species 'g6p':")
  let reactions_with_g6p = @bio.pathway_find_reactions_with_species(pathway, "g6p")
  println("   Reactions involving g6p: \{reactions_with_g6p.length()}")
  let mut i = 0
  while i < reactions_with_g6p.length() {
    let r = reactions_with_g6p[i]
    println("   - \{r.id}: \{r.name}")
    i = i + 1
  }

  println("\n5. Species degree (number of reactions):")
  let glucose_degree = @bio.species_degree(pathway, "glucose")
  let g6p_degree = @bio.species_degree(pathway, "g6p")
  let g3p_degree = @bio.species_degree(pathway, "g3p")
  println("   glucose degree: \{glucose_degree}")
  println("   g6p degree: \{g6p_degree}")
  println("   g3p degree: \{g3p_degree}")

  println("\n6. Reaction string representation:")
  let reactants : Array[String] = Array::new()
  reactants.push("A")
  let products : Array[String] = Array::new()
  products.push("B")
  let r1 = @bio.Reaction::new("r1", "Test", reactants, products, true)
  let r2 = @bio.Reaction::new("r2", "Test2", reactants, products, false)
  println("   Reversible: \{@bio.Reaction::to_string(r1)}")
  println("   Irreversible: \{@bio.Reaction::to_string(r2)}")

  println("\n=== Demo Complete ===")
}