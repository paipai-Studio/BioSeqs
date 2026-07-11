// Standard genetic code (NCBI translation table 1, "Standard"), mirroring
// BioPython's Bio.Data.CodonTable.standard_dna_table.
//
// forward_table: codon (uppercase DNA, e.g. "AAA") -> amino acid letter.
// stop_codons:   the three stop codons TAA, TAG, TGA.
// start_codons:  alternative start codons ATG, CTG, TTG (translated as M when
//                cds=True).

///|
/// Nucleotide to 2-bit encoding: A=0, C=1, G=2, T/U=3.
fn nuc_to_2bit(c : Int) -> Int? {
  if c == 'A'.to_int() || c == 'a'.to_int() {
    Some(0)
  } else if c == 'C'.to_int() || c == 'c'.to_int() {
    Some(1)
  } else if c == 'G'.to_int() || c == 'g'.to_int() {
    Some(2)
  } else if c == 'T'.to_int() ||
    c == 't'.to_int() ||
    c == 'U'.to_int() ||
    c == 'u'.to_int() {
    Some(3)
  } else {
    None
  }
}

///|
/// Forward translation table keyed by Int encoding of the codon.
/// The key is computed as: c0 * 65536 + c1 * 256 + c2, where c0/c1/c2 are
/// the UInt16 code units of the three characters. This avoids allocating a
/// String for every codon lookup.
let codon_forward_table : Map[Int, Char] = {
  let m : Map[Int, Char] = Map([], capacity=64)
  fn reg(c0 : Char, c1 : Char, c2 : Char, aa : Char) -> Unit {
    m[c0.to_int() * 65536 + c1.to_int() * 256 + c2.to_int()] = aa
  }
  // Phenylalanine (F)
  reg('T', 'T', 'T', 'F')
  reg('T', 'T', 'C', 'F')
  // Leucine (L)
  reg('T', 'T', 'A', 'L')
  reg('T', 'T', 'G', 'L')
  reg('C', 'T', 'T', 'L')
  reg('C', 'T', 'C', 'L')
  reg('C', 'T', 'A', 'L')
  reg('C', 'T', 'G', 'L')
  // Isoleucine (I)
  reg('A', 'T', 'T', 'I')
  reg('A', 'T', 'C', 'I')
  reg('A', 'T', 'A', 'I')
  // Methionine / Start
  reg('A', 'T', 'G', 'M')
  // Valine (V)
  reg('G', 'T', 'T', 'V')
  reg('G', 'T', 'C', 'V')
  reg('G', 'T', 'A', 'V')
  reg('G', 'T', 'G', 'V')
  // Serine (S)
  reg('T', 'C', 'T', 'S')
  reg('T', 'C', 'C', 'S')
  reg('T', 'C', 'A', 'S')
  reg('T', 'C', 'G', 'S')
  // Proline (P)
  reg('C', 'C', 'T', 'P')
  reg('C', 'C', 'C', 'P')
  reg('C', 'C', 'A', 'P')
  reg('C', 'C', 'G', 'P')
  // Threonine (T)
  reg('A', 'C', 'T', 'T')
  reg('A', 'C', 'C', 'T')
  reg('A', 'C', 'A', 'T')
  reg('A', 'C', 'G', 'T')
  // Alanine (A)
  reg('G', 'C', 'T', 'A')
  reg('G', 'C', 'C', 'A')
  reg('G', 'C', 'A', 'A')
  reg('G', 'C', 'G', 'A')
  // Tyrosine (Y)
  reg('T', 'A', 'T', 'Y')
  reg('T', 'A', 'C', 'Y')
  // Histidine (H)
  reg('C', 'A', 'T', 'H')
  reg('C', 'A', 'C', 'H')
  // Glutamine (Q)
  reg('C', 'A', 'A', 'Q')
  reg('C', 'A', 'G', 'Q')
  // Asparagine (N)
  reg('A', 'A', 'T', 'N')
  reg('A', 'A', 'C', 'N')
  // Lysine (K)
  reg('A', 'A', 'A', 'K')
  reg('A', 'A', 'G', 'K')
  // Aspartate (D)
  reg('G', 'A', 'T', 'D')
  reg('G', 'A', 'C', 'D')
  // Glutamate (E)
  reg('G', 'A', 'A', 'E')
  reg('G', 'A', 'G', 'E')
  // Cysteine (C)
  reg('T', 'G', 'T', 'C')
  reg('T', 'G', 'C', 'C')
  // Tryptophan (W)
  reg('T', 'G', 'G', 'W')
  // Arginine (R)
  reg('C', 'G', 'T', 'R')
  reg('C', 'G', 'C', 'R')
  reg('C', 'G', 'A', 'R')
  reg('C', 'G', 'G', 'R')
  reg('A', 'G', 'A', 'R')
  reg('A', 'G', 'G', 'R')
  // Glycine (G)
  reg('G', 'G', 'T', 'G')
  reg('G', 'G', 'C', 'G')
  reg('G', 'G', 'A', 'G')
  reg('G', 'G', 'G', 'G')
  // Serine (AG-)
  reg('A', 'G', 'T', 'S')
  reg('A', 'G', 'C', 'S')
  m
}

///|
/// Lookup table for valid nucleotide letters (ambiguous DNA + RNA), indexed by
/// UInt16 code unit. Mirrors BioPython's IUPACData.ambiguous_dna_letters +
/// ambiguous_rna_letters. Using FixedArray for O(1) lookup.
let valid_letter_table : FixedArray[Bool] = {
  let t = FixedArray::make(128, false)
  for c in "ACGTUMRWSYKVHDBXN".iter() {
    t[c.to_int()] = true
  }
  t
}

///|
/// Check if a code unit (as Int) is a valid (ambiguous) nucleotide letter.
fn valid_letter(c : Int) -> Bool {
  if c < 128 {
    valid_letter_table[c]
  } else {
    false
  }
}
