// Standard genetic code (NCBI translation table 1, "Standard"), mirroring
// BioPython's Bio.Data.CodonTable.standard_dna_table.
//
// forward_table: codon (uppercase DNA, e.g. "AAA") -> amino acid letter.
// stop_codons:   the three stop codons TAA, TAG, TGA.
// start_codons:  alternative start codons ATG, CTG, TTG (translated as M when
//                cds=True).

///|
/// Nucleotide to 2-bit encoding lookup table.
/// A=0, C=1, G=2, T/U=3, 4 for ambiguous, -1 for invalid.
pub let nuc_to_2bit_table : FixedArray[Int] = {
  let t = FixedArray::make(256, -1)
  t[65] = 0
  t[97] = 0 // A, a
  t[67] = 1
  t[99] = 1 // C, c
  t[71] = 2
  t[103] = 2 // G, g
  t[84] = 3
  t[116] = 3 // T, t
  t[85] = 3
  t[117] = 3 // U, u
  t[78] = 4
  t[110] = 4 // N, n - ambiguous
  t[77] = 4
  t[109] = 4 // M, m
  t[82] = 4
  t[114] = 4 // R, r
  t[87] = 4
  t[119] = 4 // W, w
  t[83] = 4
  t[115] = 4 // S, s
  t[89] = 4
  t[121] = 4 // Y, y
  t[75] = 4
  t[107] = 4 // K, k
  t[86] = 4
  t[118] = 4 // V, v
  t[72] = 4
  t[104] = 4 // H, h
  t[68] = 4
  t[100] = 4 // D, d
  t[66] = 4
  t[98] = 4 // B, b
  t[88] = 4
  t[120] = 4 // X, x
  t
}

///|
/// Nucleotide to 2-bit encoding: A=0, C=1, G=2, T/U=3.
/// Handles both uppercase and lowercase. Returns -1 for invalid nucleotides.
pub fn nuc_to_2bit(c : UInt16) -> Int {
  nuc_to_2bit_table[c.to_int()]
}

///|
///|
/// FixedArray-based codon lookup table for fast translation.
/// Index: c0*16 + c1*4 + c2 (2-bit encoding: A=0, C=1, G=2, T/U=3)
/// Values: UInt16 code of amino acid, or special sentinel values:
///   - 0xFFFF: stop codon
///   - 0xFFFE: ambiguous/unmapped codon (needs fallback)
pub let codon_fixed_table : FixedArray[UInt16] = {
  let t : FixedArray[UInt16] = FixedArray::make(64, 0xFFFE)

  // Phenylalanine (F)
  t[3 * 16 + 3 * 4 + 3] = 'F'.to_int().to_uint16() // TTT
  t[3 * 16 + 3 * 4 + 1] = 'F'.to_int().to_uint16() // TTC
  // Leucine (L)
  t[3 * 16 + 3 * 4 + 0] = 'L'.to_int().to_uint16() // TTA
  t[3 * 16 + 3 * 4 + 2] = 'L'.to_int().to_uint16() // TTG
  t[1 * 16 + 3 * 4 + 3] = 'L'.to_int().to_uint16() // CTT
  t[1 * 16 + 3 * 4 + 1] = 'L'.to_int().to_uint16() // CTC
  t[1 * 16 + 3 * 4 + 0] = 'L'.to_int().to_uint16() // CTA
  t[1 * 16 + 3 * 4 + 2] = 'L'.to_int().to_uint16() // CTG
  // Isoleucine (I)
  t[0 * 16 + 3 * 4 + 3] = 'I'.to_int().to_uint16() // ATT
  t[0 * 16 + 3 * 4 + 1] = 'I'.to_int().to_uint16() // ATC
  t[0 * 16 + 3 * 4 + 0] = 'I'.to_int().to_uint16() // ATA
  // Methionine / Start (M)
  t[0 * 16 + 3 * 4 + 2] = 'M'.to_int().to_uint16() // ATG
  // Valine (V)
  t[2 * 16 + 3 * 4 + 3] = 'V'.to_int().to_uint16() // GTT
  t[2 * 16 + 3 * 4 + 1] = 'V'.to_int().to_uint16() // GTC
  t[2 * 16 + 3 * 4 + 0] = 'V'.to_int().to_uint16() // GTA
  t[2 * 16 + 3 * 4 + 2] = 'V'.to_int().to_uint16() // GTG
  // Serine (S)
  t[3 * 16 + 1 * 4 + 3] = 'S'.to_int().to_uint16() // TCT
  t[3 * 16 + 1 * 4 + 1] = 'S'.to_int().to_uint16() // TCC
  t[3 * 16 + 1 * 4 + 0] = 'S'.to_int().to_uint16() // TCA
  t[3 * 16 + 1 * 4 + 2] = 'S'.to_int().to_uint16() // TCG
  // Proline (P)
  t[1 * 16 + 1 * 4 + 3] = 'P'.to_int().to_uint16() // CCT
  t[1 * 16 + 1 * 4 + 1] = 'P'.to_int().to_uint16() // CCC
  t[1 * 16 + 1 * 4 + 0] = 'P'.to_int().to_uint16() // CCA
  t[1 * 16 + 1 * 4 + 2] = 'P'.to_int().to_uint16() // CCG
  // Threonine (T)
  t[0 * 16 + 1 * 4 + 3] = 'T'.to_int().to_uint16() // ACT
  t[0 * 16 + 1 * 4 + 1] = 'T'.to_int().to_uint16() // ACC
  t[0 * 16 + 1 * 4 + 0] = 'T'.to_int().to_uint16() // ACA
  t[0 * 16 + 1 * 4 + 2] = 'T'.to_int().to_uint16() // ACG
  // Alanine (A)
  t[2 * 16 + 1 * 4 + 3] = 'A'.to_int().to_uint16() // GCT
  t[2 * 16 + 1 * 4 + 1] = 'A'.to_int().to_uint16() // GCC
  t[2 * 16 + 1 * 4 + 0] = 'A'.to_int().to_uint16() // GCA
  t[2 * 16 + 1 * 4 + 2] = 'A'.to_int().to_uint16() // GCG
  // Tyrosine (Y)
  t[3 * 16 + 0 * 4 + 3] = 'Y'.to_int().to_uint16() // TAT
  t[3 * 16 + 0 * 4 + 1] = 'Y'.to_int().to_uint16() // TAC
  // Stop codons
  t[3 * 16 + 0 * 4 + 0] = 0xFFFF // TAA
  t[3 * 16 + 0 * 4 + 2] = 0xFFFF // TAG
  // Histidine (H)
  t[1 * 16 + 0 * 4 + 3] = 'H'.to_int().to_uint16() // CAT
  t[1 * 16 + 0 * 4 + 1] = 'H'.to_int().to_uint16() // CAC
  // Glutamine (Q)
  t[1 * 16 + 0 * 4 + 0] = 'Q'.to_int().to_uint16() // CAA
  t[1 * 16 + 0 * 4 + 2] = 'Q'.to_int().to_uint16() // CAG
  // Asparagine (N)
  t[0 * 16 + 0 * 4 + 3] = 'N'.to_int().to_uint16() // AAT
  t[0 * 16 + 0 * 4 + 1] = 'N'.to_int().to_uint16() // AAC
  // Lysine (K)
  t[0 * 16 + 0 * 4 + 0] = 'K'.to_int().to_uint16() // AAA
  t[0 * 16 + 0 * 4 + 2] = 'K'.to_int().to_uint16() // AAG
  // Aspartate (D)
  t[2 * 16 + 0 * 4 + 3] = 'D'.to_int().to_uint16() // GAT
  t[2 * 16 + 0 * 4 + 1] = 'D'.to_int().to_uint16() // GAC
  // Glutamate (E)
  t[2 * 16 + 0 * 4 + 0] = 'E'.to_int().to_uint16() // GAA
  t[2 * 16 + 0 * 4 + 2] = 'E'.to_int().to_uint16() // GAG
  // Cysteine (C)
  t[3 * 16 + 2 * 4 + 3] = 'C'.to_int().to_uint16() // TGT
  t[3 * 16 + 2 * 4 + 1] = 'C'.to_int().to_uint16() // TGC
  // Stop codon
  t[3 * 16 + 2 * 4 + 0] = 0xFFFF // TGA
  // Tryptophan (W)
  t[3 * 16 + 2 * 4 + 2] = 'W'.to_int().to_uint16() // TGG
  // Arginine (R)
  t[1 * 16 + 2 * 4 + 3] = 'R'.to_int().to_uint16() // CGT
  t[1 * 16 + 2 * 4 + 1] = 'R'.to_int().to_uint16() // CGC
  t[1 * 16 + 2 * 4 + 0] = 'R'.to_int().to_uint16() // CGA
  t[1 * 16 + 2 * 4 + 2] = 'R'.to_int().to_uint16() // CGG
  t[0 * 16 + 2 * 4 + 0] = 'R'.to_int().to_uint16() // AGA
  t[0 * 16 + 2 * 4 + 2] = 'R'.to_int().to_uint16() // AGG
  // Glycine (G)
  t[2 * 16 + 2 * 4 + 3] = 'G'.to_int().to_uint16() // GGT
  t[2 * 16 + 2 * 4 + 1] = 'G'.to_int().to_uint16() // GGC
  t[2 * 16 + 2 * 4 + 0] = 'G'.to_int().to_uint16() // GGA
  t[2 * 16 + 2 * 4 + 2] = 'G'.to_int().to_uint16() // GGG
  // Serine (AG-)
  t[0 * 16 + 2 * 4 + 3] = 'S'.to_int().to_uint16() // AGT
  t[0 * 16 + 2 * 4 + 1] = 'S'.to_int().to_uint16() // AGC

  t
}

///|
