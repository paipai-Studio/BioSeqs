// Demo for Bio.Alphabet module

fn main {
  // Create IUPAC unambiguous DNA alphabet
  let dna = @bio.iupac_unambiguous_dna()
  println("=== IUPAC Unambiguous DNA Alphabet ===")
  println("Name: \{dna.name}")
  println("Is gapped: \{dna.is_gapped}")
  println("Letters count: \{dna.letters.length()}")
  
  // Validate characters
  println("\nValidating characters:")
  println("A valid: \{dna.is_valid("A")}")
  println("C valid: \{dna.is_valid("C")}")
  println("G valid: \{dna.is_valid("G")}")
  println("T valid: \{dna.is_valid("T")}")
  println("U valid: \{dna.is_valid("U")}")
  
  // Create IUPAC ambiguous DNA alphabet
  let ambiguous_dna = @bio.iupac_ambiguous_dna()
  println("\n=== IUPAC Ambiguous DNA Alphabet ===")
  println("Name: \{ambiguous_dna.name}")
  println("Letters count: \{ambiguous_dna.letters.length()}")
  println("R valid: \{ambiguous_dna.is_valid("R")}")
  println("Y valid: \{ambiguous_dna.is_valid("Y")}")
  println("N valid: \{ambiguous_dna.is_valid("N")}")
  
  // Create protein alphabet
  let protein = @bio.iupac_protein()
  println("\n=== IUPAC Protein Alphabet ===")
  println("Name: \{protein.name}")
  println("Letters count: \{protein.letters.length()}")
  println("A valid: \{protein.is_valid("A")}")
  println("W valid: \{protein.is_valid("W")}")
  println("X valid: \{protein.is_valid("X")}")
  
  // Create gapped DNA alphabet
  let gapped_dna = @bio.gapped_dna()
  println("\n=== Gapped DNA Alphabet ===")
  println("Name: \{gapped_dna.name}")
  println("Is gapped: \{gapped_dna.is_gapped}")
  println("- valid: \{gapped_dna.is_valid("-")}")
  
  // Get codon mapping
  let mapping = @bio.three_to_one()
  println("\n=== Three to One Letter Code ===")
  let ala = mapping.get("Ala")
  if ala is Some(v) {
    println("Ala -> \{v}")
  }
  let gly = mapping.get("Gly")
  if gly is Some(v) {
    println("Gly -> \{v}")
  }
  let ter = mapping.get("Ter")
  if ter is Some(v) {
    println("Ter -> \{v}")
  }
  
  println("\nDemo completed successfully!")
}
