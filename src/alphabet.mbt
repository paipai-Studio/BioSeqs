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
pub fn Alphabet::new(
  name : String,
  letters : Array[String],
  is_gapped : Bool,
) -> Alphabet {
  Alphabet::{ name, letters, is_gapped, }
}

///|
/// Get alphabet name.
pub fn Alphabet::get_name(self : Alphabet) -> String {
  self.name
}

///|
/// Get alphabet letters.
pub fn Alphabet::get_letters(self : Alphabet) -> Array[String] {
  self.letters
}

///|
/// Check if alphabet is gapped.
pub fn Alphabet::is_gapped(self : Alphabet) -> Bool {
  self.is_gapped
}

///|
/// Check if a character is valid in this alphabet.
pub fn Alphabet::is_valid(self : Alphabet, c : String) -> Bool {
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
  let letters = ["A", "C", "G", "T"]
  Alphabet::new("IUPACUnambiguousDNA", letters, false)
}

///|
/// IUPAC unambiguous RNA alphabet (A, C, G, U).
pub fn iupac_unambiguous_rna() -> Alphabet {
  let letters = ["A", "C", "G", "U"]
  Alphabet::new("IUPACUnambiguousRNA", letters, false)
}

///|
/// IUPAC ambiguous DNA alphabet (A, C, G, T, R, Y, S, W, K, M, B, D, H, V, N).
pub fn iupac_ambiguous_dna() -> Alphabet {
  let letters = [
    "A", "C", "G", "T", "R", "Y", "S", "W", "K", "M", "B", "D", "H", "V", "N",
  ]
  Alphabet::new("IUPACAmbiguousDNA", letters, false)
}

///|
/// IUPAC protein alphabet (20 standard amino acids + B, Z, X).
pub fn iupac_protein() -> Alphabet {
  let letters = [
    "A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S",
    "T", "V", "W", "Y", "B", "Z", "X",
  ]
  Alphabet::new("IUPACProtein", letters, false)
}

///|
/// Gapped DNA alphabet (includes gap character '-').
pub fn gapped_dna() -> Alphabet {
  let letters = ["A", "C", "G", "T", "-"]
  Alphabet::new("GappedDNA", letters, true)
}

///|
/// Reduced alphabet classes for amino acids.
/// Grouped by chemical properties.
pub fn reduced_protein() -> Alphabet {
  let letters = [
    "A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S",
    "T", "V", "W", "Y",
  ]
  Alphabet::new("ReducedProtein", letters, false)
}

///|
/// Map from three-letter amino acid codes to one-letter codes.
pub fn three_to_one() -> Map[String, String] {
  amino_3_to_1()
}

///|
/// Map from one-letter amino acid codes to three-letter codes.
pub fn one_to_three() -> Map[String, String] {
  amino_1_to_3()
}
