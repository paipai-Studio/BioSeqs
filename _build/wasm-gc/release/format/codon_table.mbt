// Standard genetic code (NCBI translation table 1, "Standard"), mirroring
// BioPython's Bio.Data.CodonTable.standard_dna_table.
//
// forward_table: codon (uppercase DNA, e.g. "AAA") -> amino acid letter.
// stop_codons:   the three stop codons TAA, TAG, TGA.
// start_codons:  alternative start codons ATG, CTG, TTG (translated as M when
//                cds=True).

///|
/// Forward translation table keyed by Int encoding of the codon.
/// The key is computed as: c0 * 65536 + c1 * 256 + c2, where c0/c1/c2 are
/// the UInt16 code units of the three characters. This avoids allocating a
/// String for every codon lookup.
let codon_forward_table : Map[Int, Char] = {
  let m : Map[Int, Char] = Map([], capacity=64)
  // Helper to register a codon
  fn reg(c0 : Int, c1 : Int, c2 : Int, aa : Char) -> Unit {
    m[c0 * 65536 + c1 * 256 + c2] = aa
  }
  // Phenylalanine (F)
  reg('T'.to_int(), 'T'.to_int(), 'T'.to_int(), 'F')
  reg('T'.to_int(), 'T'.to_int(), 'C'.to_int(), 'F')
  // Leucine (L)
  reg('T'.to_int(), 'T'.to_int(), 'A'.to_int(), 'L')
  reg('T'.to_int(), 'T'.to_int(), 'G'.to_int(), 'L')
  reg('C'.to_int(), 'T'.to_int(), 'T'.to_int(), 'L')
  reg('C'.to_int(), 'T'.to_int(), 'C'.to_int(), 'L')
  reg('C'.to_int(), 'T'.to_int(), 'A'.to_int(), 'L')
  reg('C'.to_int(), 'T'.to_int(), 'G'.to_int(), 'L')
  // Isoleucine (I)
  reg('A'.to_int(), 'T'.to_int(), 'T'.to_int(), 'I')
  reg('A'.to_int(), 'T'.to_int(), 'C'.to_int(), 'I')
  reg('A'.to_int(), 'T'.to_int(), 'A'.to_int(), 'I')
  // Methionine / Start
  reg('A'.to_int(), 'T'.to_int(), 'G'.to_int(), 'M')
  // Valine (V)
  reg('G'.to_int(), 'T'.to_int(), 'T'.to_int(), 'V')
  reg('G'.to_int(), 'T'.to_int(), 'C'.to_int(), 'V')
  reg('G'.to_int(), 'T'.to_int(), 'A'.to_int(), 'V')
  reg('G'.to_int(), 'T'.to_int(), 'G'.to_int(), 'V')
  // Serine (S)
  reg('T'.to_int(), 'C'.to_int(), 'T'.to_int(), 'S')
  reg('T'.to_int(), 'C'.to_int(), 'C'.to_int(), 'S')
  reg('T'.to_int(), 'C'.to_int(), 'A'.to_int(), 'S')
  reg('T'.to_int(), 'C'.to_int(), 'G'.to_int(), 'S')
  // Proline (P)
  reg('C'.to_int(), 'C'.to_int(), 'T'.to_int(), 'P')
  reg('C'.to_int(), 'C'.to_int(), 'C'.to_int(), 'P')
  reg('C'.to_int(), 'C'.to_int(), 'A'.to_int(), 'P')
  reg('C'.to_int(), 'C'.to_int(), 'G'.to_int(), 'P')
  // Threonine (T)
  reg('A'.to_int(), 'C'.to_int(), 'T'.to_int(), 'T')
  reg('A'.to_int(), 'C'.to_int(), 'C'.to_int(), 'T')
  reg('A'.to_int(), 'C'.to_int(), 'A'.to_int(), 'T')
  reg('A'.to_int(), 'C'.to_int(), 'G'.to_int(), 'T')
  // Alanine (A)
  reg('G'.to_int(), 'C'.to_int(), 'T'.to_int(), 'A')
  reg('G'.to_int(), 'C'.to_int(), 'C'.to_int(), 'A')
  reg('G'.to_int(), 'C'.to_int(), 'A'.to_int(), 'A')
  reg('G'.to_int(), 'C'.to_int(), 'G'.to_int(), 'A')
  // Tyrosine (Y)
  reg('T'.to_int(), 'A'.to_int(), 'T'.to_int(), 'Y')
  reg('T'.to_int(), 'A'.to_int(), 'C'.to_int(), 'Y')
  // Histidine (H)
  reg('C'.to_int(), 'A'.to_int(), 'T'.to_int(), 'H')
  reg('C'.to_int(), 'A'.to_int(), 'C'.to_int(), 'H')
  // Glutamine (Q)
  reg('C'.to_int(), 'A'.to_int(), 'A'.to_int(), 'Q')
  reg('C'.to_int(), 'A'.to_int(), 'G'.to_int(), 'Q')
  // Asparagine (N)
  reg('A'.to_int(), 'A'.to_int(), 'T'.to_int(), 'N')
  reg('A'.to_int(), 'A'.to_int(), 'C'.to_int(), 'N')
  // Lysine (K)
  reg('A'.to_int(), 'A'.to_int(), 'A'.to_int(), 'K')
  reg('A'.to_int(), 'A'.to_int(), 'G'.to_int(), 'K')
  // Aspartate (D)
  reg('G'.to_int(), 'A'.to_int(), 'T'.to_int(), 'D')
  reg('G'.to_int(), 'A'.to_int(), 'C'.to_int(), 'D')
  // Glutamate (E)
  reg('G'.to_int(), 'A'.to_int(), 'A'.to_int(), 'E')
  reg('G'.to_int(), 'A'.to_int(), 'G'.to_int(), 'E')
  // Cysteine (C)
  reg('T'.to_int(), 'G'.to_int(), 'T'.to_int(), 'C')
  reg('T'.to_int(), 'G'.to_int(), 'C'.to_int(), 'C')
  // Tryptophan (W)
  reg('T'.to_int(), 'G'.to_int(), 'G'.to_int(), 'W')
  // Arginine (R)
  reg('C'.to_int(), 'G'.to_int(), 'T'.to_int(), 'R')
  reg('C'.to_int(), 'G'.to_int(), 'C'.to_int(), 'R')
  reg('C'.to_int(), 'G'.to_int(), 'A'.to_int(), 'R')
  reg('C'.to_int(), 'G'.to_int(), 'G'.to_int(), 'R')
  reg('A'.to_int(), 'G'.to_int(), 'A'.to_int(), 'R')
  reg('A'.to_int(), 'G'.to_int(), 'G'.to_int(), 'R')
  // Glycine (G)
  reg('G'.to_int(), 'G'.to_int(), 'T'.to_int(), 'G')
  reg('G'.to_int(), 'G'.to_int(), 'C'.to_int(), 'G')
  reg('G'.to_int(), 'G'.to_int(), 'A'.to_int(), 'G')
  reg('G'.to_int(), 'G'.to_int(), 'G'.to_int(), 'G')
  // Serine (AG-)
  reg('A'.to_int(), 'G'.to_int(), 'T'.to_int(), 'S')
  reg('A'.to_int(), 'G'.to_int(), 'C'.to_int(), 'S')
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
