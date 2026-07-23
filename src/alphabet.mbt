///|
/// Bio.Alphabet module - Sequence alphabet definitions.
/// Inspired by Biopython's Bio.Alphabet package.
/// Provides IUPAC alphabets, reduced alphabets, and gapped alphabets.

///|
/// Represents a sequence alphabet type.
pub struct Alphabet {
  name : String
  letters : Array[String]
  is_gapped : Bool
}

///|
/// Create a new alphabet.
pub fn Alphabet::new(name : String, letters : Array[String], is_gapped : Bool) -> Alphabet {
  Alphabet::{ name, letters, is_gapped }
}

///|
/// Get alphabet name.
pub fn get_name(self : Alphabet) -> String {
  self.name
}

///|
/// Get alphabet letters.
pub fn get_letters(self : Alphabet) -> Array[String] {
  self.letters
}

///|
/// Check if alphabet is gapped.
pub fn is_gapped(self : Alphabet) -> Bool {
  self.is_gapped
}

///|
/// Check if a character is valid in this alphabet.
pub fn is_valid(self : Alphabet, c : String) -> Bool {
  let mut i = 0
  while i < self.letters.length() {
    if self.letters[i] == c {
      return true
    }
    i = i + 1
  }
  false
}

///|
/// IUPAC unambiguous DNA alphabet (A, C, G, T).
pub fn iupac_unambiguous_dna() -> Alphabet {
  let letters = [
    "A", "C", "G", "T"
  ]
  Alphabet::new("IUPACUnambiguousDNA", letters, false)
}

///|
/// IUPAC unambiguous RNA alphabet (A, C, G, U).
pub fn iupac_unambiguous_rna() -> Alphabet {
  let letters = [
    "A", "C", "G", "U"
  ]
  Alphabet::new("IUPACUnambiguousRNA", letters, false)
}

///|
/// IUPAC ambiguous DNA alphabet (A, C, G, T, R, Y, S, W, K, M, B, D, H, V, N).
pub fn iupac_ambiguous_dna() -> Alphabet {
  let letters = [
    "A", "C", "G", "T",
    "R", "Y", "S", "W", "K", "M",
    "B", "D", "H", "V", "N"
  ]
  Alphabet::new("IUPACAmbiguousDNA", letters, false)
}

///|
/// IUPAC ambiguous RNA alphabet (A, C, G, U, R, Y, S, W, K, M, B, D, H, V, N).
pub fn iupac_ambiguous_rna() -> Alphabet {
  let letters = [
    "A", "C", "G", "U",
    "R", "Y", "S", "W", "K", "M",
    "B", "D", "H", "V", "N"
  ]
  Alphabet::new("IUPACAmbiguousRNA", letters, false)
}

///|
/// IUPAC protein alphabet (20 standard amino acids + B, Z, X).
pub fn iupac_protein() -> Alphabet {
  let letters = [
    "A", "C", "D", "E", "F", "G", "H", "I", "K", "L",
    "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y",
    "B", "Z", "X"
  ]
  Alphabet::new("IUPACProtein", letters, false)
}

///|
/// Extended IUPAC protein alphabet (20 standard + B, Z, X, U, O).
pub fn iupac_extended_protein() -> Alphabet {
  let letters = [
    "A", "C", "D", "E", "F", "G", "H", "I", "K", "L",
    "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y",
    "B", "Z", "X", "U", "O"
  ]
  Alphabet::new("IUPACExtendedProtein", letters, false)
}

///|
/// Gapped DNA alphabet (includes gap character '-').
pub fn gapped_dna() -> Alphabet {
  let letters = [
    "A", "C", "G", "T", "-"
  ]
  Alphabet::new("GappedDNA", letters, true)
}

///|
/// Gapped RNA alphabet (includes gap character '-').
pub fn gapped_rna() -> Alphabet {
  let letters = [
    "A", "C", "G", "U", "-"
  ]
  Alphabet::new("GappedRNA", letters, true)
}

///|
/// Gapped protein alphabet (includes gap character '-').
pub fn gapped_protein() -> Alphabet {
  let letters = [
    "A", "C", "D", "E", "F", "G", "H", "I", "K", "L",
    "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y",
    "-"
  ]
  Alphabet::new("GappedProtein", letters, true)
}

///|
/// Reduced alphabet classes for amino acids.
/// Grouped by chemical properties.
pub fn reduced_protein() -> Alphabet {
  let letters = [
    "A", "C", "D", "E", "F", "G", "H", "I", "K", "L",
    "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y"
  ]
  Alphabet::new("ReducedProtein", letters, false)
}

///|
/// Get all ambiguous codes for IUPAC DNA.
pub fn iupac_dna_ambiguities() -> Map[String, Array[String]] {
  let ambiguities : Map[String, Array[String]] = Map([], capacity=16)
  ambiguities.set("R", ["A", "G"])
  ambiguities.set("Y", ["C", "T"])
  ambiguities.set("S", ["G", "C"])
  ambiguities.set("W", ["A", "T"])
  ambiguities.set("K", ["G", "T"])
  ambiguities.set("M", ["A", "C"])
  ambiguities.set("B", ["C", "G", "T"])
  ambiguities.set("D", ["A", "G", "T"])
  ambiguities.set("H", ["A", "C", "T"])
  ambiguities.set("V", ["A", "C", "G"])
  ambiguities.set("N", ["A", "C", "G", "T"])
  ambiguities
}

///|
/// Get all ambiguous codes for IUPAC RNA.
pub fn iupac_rna_ambiguities() -> Map[String, Array[String]] {
  let ambiguities : Map[String, Array[String]] = Map([], capacity=16)
  ambiguities.set("R", ["A", "G"])
  ambiguities.set("Y", ["C", "U"])
  ambiguities.set("S", ["G", "C"])
  ambiguities.set("W", ["A", "U"])
  ambiguities.set("K", ["G", "U"])
  ambiguities.set("M", ["A", "C"])
  ambiguities.set("B", ["C", "G", "U"])
  ambiguities.set("D", ["A", "G", "U"])
  ambiguities.set("H", ["A", "C", "U"])
  ambiguities.set("V", ["A", "C", "G"])
  ambiguities.set("N", ["A", "C", "G", "U"])
  ambiguities
}

///|
/// Map from three-letter amino acid codes to one-letter codes.
pub fn three_to_one() -> Map[String, String] {
  let mapping : Map[String, String] = Map([], capacity=25)
  mapping.set("Ala", "A")
  mapping.set("Arg", "R")
  mapping.set("Asn", "N")
  mapping.set("Asp", "D")
  mapping.set("Cys", "C")
  mapping.set("Gln", "Q")
  mapping.set("Glu", "E")
  mapping.set("Gly", "G")
  mapping.set("His", "H")
  mapping.set("Ile", "I")
  mapping.set("Leu", "L")
  mapping.set("Lys", "K")
  mapping.set("Met", "M")
  mapping.set("Phe", "F")
  mapping.set("Pro", "P")
  mapping.set("Ser", "S")
  mapping.set("Thr", "T")
  mapping.set("Trp", "W")
  mapping.set("Tyr", "Y")
  mapping.set("Val", "V")
  mapping.set("Ter", "*")
  mapping
}

///|
/// Map from one-letter amino acid codes to three-letter codes.
pub fn one_to_three() -> Map[String, String] {
  let mapping : Map[String, String] = Map([], capacity=25)
  mapping.set("A", "Ala")
  mapping.set("R", "Arg")
  mapping.set("N", "Asn")
  mapping.set("D", "Asp")
  mapping.set("C", "Cys")
  mapping.set("Q", "Gln")
  mapping.set("E", "Glu")
  mapping.set("G", "Gly")
  mapping.set("H", "His")
  mapping.set("I", "Ile")
  mapping.set("L", "Leu")
  mapping.set("K", "Lys")
  mapping.set("M", "Met")
  mapping.set("F", "Phe")
  mapping.set("P", "Pro")
  mapping.set("S", "Ser")
  mapping.set("T", "Thr")
  mapping.set("W", "Trp")
  mapping.set("Y", "Tyr")
  mapping.set("V", "Val")
  mapping.set("*", "Ter")
  mapping
}
