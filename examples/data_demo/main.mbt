///|
/// Bio.Data module example - Biological data constants and utilities.
///
/// Run: moon run examples/data_demo/main.mbt
fn main {
  println("=== Bio.Data Module Demo ===")
  println("")

  // IUPAC ambiguous bases
  println("=== IUPAC Ambiguous DNA Bases ===")
  let dna_amb = @bio.iupac_ambiguous_dna_map()
  println("Total ambiguous bases: \{dna_amb.size()}")
  println("R (G/A): \{dna_amb.get("R").unwrap()}")
  println("Y (C/T): \{dna_amb.get("Y").unwrap()}")
  println("N (A/C/G/T): \{dna_amb.get("N").unwrap()}")
  println("")

  println("=== IUPAC Ambiguous RNA Bases ===")
  let rna_amb = @bio.iupac_ambiguous_rna_map()
  println("Total ambiguous bases: \{rna_amb.size()}")
  println("Y (C/U): \{rna_amb.get("Y").unwrap()}")
  println("")

  // Check if bases are ambiguous
  println("=== Ambiguity Checks ===")
  println("Is 'R' ambiguous DNA? \{@bio.is_ambiguous_dna("R")}")
  println("Is 'A' ambiguous DNA? \{@bio.is_ambiguous_dna("A")}")
  println("Is 'N' ambiguous DNA? \{@bio.is_ambiguous_dna("N")}")
  println("")

  // Expand ambiguous bases
  println("=== Expand Ambiguous Bases ===")
  println("Expand 'R': \{@bio.expand_ambiguous_dna("R")}")
  println("Expand 'N': \{@bio.expand_ambiguous_dna("N")}")
  println("Expand 'A': \{@bio.expand_ambiguous_dna("A")}")
  println("")

  // Amino acid mappings
  println("=== Amino Acid 3-to-1 Letter Code ===")
  let aa_map = @bio.amino_3_to_1()
  println("Total amino acids: \{aa_map.size()}")
  println("Ala -> \{aa_map.get("Ala").unwrap()}")
  println("Ter -> \{aa_map.get("Ter").unwrap()}")
  println("")

  println("=== Amino Acid 1-to-3 Letter Code ===")
  let aa1to3 = @bio.amino_1_to_3()
  println("Total: \{aa1to3.size()}")
  println("A -> \{aa1to3.get("A").unwrap()}")
  println("* -> \{aa1to3.get("*").unwrap()}")
  println("")

  // Amino acid lookup
  println("=== Direct Amino Acid Lookup ===")
  println("amino3_to_1('Ala'): \{@bio.amino3_to_1("Ala")}")
  println("amino3_to_1('Ter'): \{@bio.amino3_to_1("Ter")}")
  println("amino3_to_1('XXX'): \{@bio.amino3_to_1("XXX")}")
  println("amino1_to_3('A'): \{@bio.amino1_to_3("A")}")
  println("amino1_to_3('*'): \{@bio.amino1_to_3("*")}")
  println("")

  // Amino acid properties
  println("=== Amino Acid Properties ===")
  let ala = @bio.amino_acid_properties("A")
  println("Ala: property=\{ala.property}, charge=\{ala.charge}, hydropathy=\{ala.hydropathy}")
  let arg = @bio.amino_acid_properties("R")
  println("Arg: property=\{arg.property}, charge=\{arg.charge}, hydropathy=\{arg.hydropathy}")
  let unknown = @bio.amino_acid_properties("X")
  println("Unknown: property=\{unknown.property}")
  println("")

  // Codon table
  println("=== Standard Codon Table ===")
  let codons = @bio.standard_codon_table()
  println("Total codons: \{codons.size()}")
  println("ATG (start) -> \{codons.get("ATG").unwrap()}")
  println("TAA (stop) -> \{codons.get("TAA").unwrap()}")
  println("TAG (stop) -> \{codons.get("TAG").unwrap()}")
  println("TGA (stop) -> \{codons.get("TGA").unwrap()}")
  println("")

  // Codon lookup
  println("=== Codon Lookup ===")
  println("codon_lookup('ATG'): \{@bio.codon_lookup("ATG")}")
  println("codon_lookup('TAA'): \{@bio.codon_lookup("TAA")}")
  println("codon_lookup('atg'): \{@bio.codon_lookup("atg")}")
  println("codon_lookup('XXX'): \{@bio.codon_lookup("XXX")}")
  println("")

  // Reverse codon table
  println("=== Reverse Codon Table ===")
  let rev = @bio.reverse_codon_table()
  println("Total amino acids: \{rev.size()}")
  println("M (Met) codons: \{rev.get("M").unwrap()}")
  println("L (Leu) codons: \{rev.get("L").unwrap()} (count: \{rev.get("L").unwrap().length()})")
  println("* (Stop) codons: \{rev.get("*").unwrap()}")
  println("")

  // Synonymous codons
  println("=== Synonymous Codon Counts ===")
  println("M: \{@bio.num_synonymous_codons("M")}")
  println("L: \{@bio.num_synonymous_codons("L")}")
  println("*: \{@bio.num_synonymous_codons("*")}")
  println("X: \{@bio.num_synonymous_codons("X")}")
  println("")

  // Reverse complement with IUPAC
  println("=== IUPAC Reverse Complement ===")
  println("ATCG -> \{@bio.reverse_complement_iupac("ATCG")}")
  println("R -> \{@bio.reverse_complement_iupac("R")}")
  println("Y -> \{@bio.reverse_complement_iupac("Y")}")
  println("N -> \{@bio.reverse_complement_iupac("N")}")
  println("K -> \{@bio.reverse_complement_iupac("K")}")
  println("")

  println("=== Done ===")
}