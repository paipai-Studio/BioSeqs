// Standard genetic code (NCBI translation table 1, "Standard"), mirroring
// BioPython's Bio.Data.CodonTable.standard_dna_table.
//
// forward_table: codon (uppercase DNA, e.g. "AAA") -> amino acid letter.
// stop_codons:   the three stop codons TAA, TAG, TGA.
// start_codons:  alternative start codons ATG, CTG, TTG (translated as M when
//                cds=True).

///|
/// Forward translation table: codon -> amino acid (single-letter code).
/// Derived from CodonTable.standard_dna_table.forward_table.
let standard_forward_table : Map[String, Char] = {
  let m = Map::new()
  // Phenylalanine (F)
  m["TTT"] = 'F'
  m["TTC"] = 'F'
  // Leucine (L)
  m["TTA"] = 'L'
  m["TTG"] = 'L'
  m["CTT"] = 'L'
  m["CTC"] = 'L'
  m["CTA"] = 'L'
  m["CTG"] = 'L'
  // Isoleucine (I)
  m["ATT"] = 'I'
  m["ATC"] = 'I'
  m["ATA"] = 'I'
  // Methionine / Start
  m["ATG"] = 'M'
  // Valine (V)
  m["GTT"] = 'V'
  m["GTC"] = 'V'
  m["GTA"] = 'V'
  m["GTG"] = 'V'
  // Serine (S)
  m["TCT"] = 'S'
  m["TCC"] = 'S'
  m["TCA"] = 'S'
  m["TCG"] = 'S'
  // Proline (P)
  m["CCT"] = 'P'
  m["CCC"] = 'P'
  m["CCA"] = 'P'
  m["CCG"] = 'P'
  // Threonine (T)
  m["ACT"] = 'T'
  m["ACC"] = 'T'
  m["ACA"] = 'T'
  m["ACG"] = 'T'
  // Alanine (A)
  m["GCT"] = 'A'
  m["GCC"] = 'A'
  m["GCA"] = 'A'
  m["GCG"] = 'A'
  // Tyrosine (Y)
  m["TAT"] = 'Y'
  m["TAC"] = 'Y'
  // Histidine (H)
  m["CAT"] = 'H'
  m["CAC"] = 'H'
  // Glutamine (Q)
  m["CAA"] = 'Q'
  m["CAG"] = 'Q'
  // Asparagine (N)
  m["AAT"] = 'N'
  m["AAC"] = 'N'
  // Lysine (K)
  m["AAA"] = 'K'
  m["AAG"] = 'K'
  // Aspartate (D)
  m["GAT"] = 'D'
  m["GAC"] = 'D'
  // Glutamate (E)
  m["GAA"] = 'E'
  m["GAG"] = 'E'
  // Cysteine (C)
  m["TGT"] = 'C'
  m["TGC"] = 'C'
  // Tryptophan (W)
  m["TGG"] = 'W'
  // Arginine (R)
  m["CGT"] = 'R'
  m["CGC"] = 'R'
  m["CGA"] = 'R'
  m["CGG"] = 'R'
  m["AGA"] = 'R'
  m["AGG"] = 'R'
  // Glycine (G)
  m["GGT"] = 'G'
  m["GGC"] = 'G'
  m["GGA"] = 'G'
  m["GGG"] = 'G'
  // Serine (AG-)
  m["AGT"] = 'S'
  m["AGC"] = 'S'
  m
}

///|
/// Stop codons for the standard table.
let standard_stop_codons : Array[String] = ["TAA", "TAG", "TGA"]

///|
/// Alternative start codons for the standard table (translated as M when
/// cds=True).
let standard_start_codons : Array[String] = ["ATG", "CTG", "TTG"]

///|
/// Set of valid nucleotide letters (ambiguous DNA + RNA) used to detect
/// "possible stop" codons such as TAN or NNN. Mirrors BioPython's
/// IUPACData.ambiguous_dna_letters + ambiguous_rna_letters.
let valid_nucleotide_letters : Map[Char, Bool] = {
  let m = Map::new()
  for c in "ACGTUMRWSYKVHDBXN".iter() {
    m[c] = true
  }
  m
}
