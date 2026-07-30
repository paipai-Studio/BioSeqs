///|
/// Bio.Data module - Biological data constants and utilities.
/// Inspired by Biopython's Bio.Data package.
/// Provides IUPAC ambiguous bases, amino acid property tables, and codon mappings.

///|
/// IUPAC DNA ambiguous bases mapping.
/// R = G or A, Y = C or T, S = G or C, W = A or T, K = G or T,
/// M = A or C, B = C or G or T, D = A or G or T, H = A or C or T,
/// V = A or C or G, N = A or C or G or T
pub fn iupac_ambiguous_dna_map() -> Map[String, Array[String]] {
  let m : Map[String, Array[String]] = Map([], capacity=11)
  m["R"] = ["G", "A"]
  m["Y"] = ["C", "T"]
  m["S"] = ["G", "C"]
  m["W"] = ["A", "T"]
  m["K"] = ["G", "T"]
  m["M"] = ["A", "C"]
  m["B"] = ["C", "G", "T"]
  m["D"] = ["A", "G", "T"]
  m["H"] = ["A", "C", "T"]
  m["V"] = ["A", "C", "G"]
  m["N"] = ["A", "C", "G", "T"]
  m
}

///|
/// IUPAC RNA ambiguous bases mapping.
pub fn iupac_ambiguous_rna_map() -> Map[String, Array[String]] {
  let m : Map[String, Array[String]] = Map([], capacity=11)
  m["R"] = ["G", "A"]
  m["Y"] = ["C", "U"]
  m["S"] = ["G", "C"]
  m["W"] = ["A", "U"]
  m["K"] = ["G", "U"]
  m["M"] = ["A", "C"]
  m["B"] = ["C", "G", "U"]
  m["D"] = ["A", "G", "U"]
  m["H"] = ["A", "C", "U"]
  m["V"] = ["A", "C", "G"]
  m["N"] = ["A", "C", "G", "U"]
  m
}

///|
/// Check if a DNA base is ambiguous (IUPAC ambiguous).
pub fn is_ambiguous_dna(base : String) -> Bool {
  let amb = iupac_ambiguous_dna_map()
  amb.contains(base.to_upper())
}

///|
/// Check if an RNA base is ambiguous (IUPAC ambiguous).
pub fn is_ambiguous_rna(base : String) -> Bool {
  let amb = iupac_ambiguous_rna_map()
  amb.contains(base.to_upper())
}

///|
/// Expand an IUPAC ambiguous DNA base into its constituent bases.
/// Returns a single-element array for standard bases, or the expansion for ambiguous ones.
pub fn expand_ambiguous_dna(base : String) -> Array[String] {
  let amb = iupac_ambiguous_dna_map()
  let upper = base.to_upper()
  match amb.get(upper) {
    Some(bases) => bases
    None => [upper]
  }
}

///|
/// Expand an IUPAC ambiguous RNA base into its constituent bases.
pub fn expand_ambiguous_rna(base : String) -> Array[String] {
  let amb = iupac_ambiguous_rna_map()
  let upper = base.to_upper()
  match amb.get(upper) {
    Some(bases) => bases
    None => [upper]
  }
}

///|
/// Amino acid 3-letter to 1-letter code mapping (standard 20 + stop).
pub fn amino_3_to_1() -> Map[String, String] {
  let m : Map[String, String] = Map([], capacity=21)
  m["Ala"] = "A"
  m["Arg"] = "R"
  m["Asn"] = "N"
  m["Asp"] = "D"
  m["Cys"] = "C"
  m["Gln"] = "Q"
  m["Glu"] = "E"
  m["Gly"] = "G"
  m["His"] = "H"
  m["Ile"] = "I"
  m["Leu"] = "L"
  m["Lys"] = "K"
  m["Met"] = "M"
  m["Phe"] = "F"
  m["Pro"] = "P"
  m["Ser"] = "S"
  m["Thr"] = "T"
  m["Trp"] = "W"
  m["Tyr"] = "Y"
  m["Val"] = "V"
  m["Ter"] = "*"
  m
}

///|
/// Amino acid 1-letter to 3-letter code mapping.
pub fn amino_1_to_3() -> Map[String, String] {
  let m : Map[String, String] = Map([], capacity=21)
  m["A"] = "Ala"
  m["R"] = "Arg"
  m["N"] = "Asn"
  m["D"] = "Asp"
  m["C"] = "Cys"
  m["Q"] = "Gln"
  m["E"] = "Glu"
  m["G"] = "Gly"
  m["H"] = "His"
  m["I"] = "Ile"
  m["L"] = "Leu"
  m["K"] = "Lys"
  m["M"] = "Met"
  m["F"] = "Phe"
  m["P"] = "Pro"
  m["S"] = "Ser"
  m["T"] = "Thr"
  m["W"] = "Trp"
  m["Y"] = "Tyr"
  m["V"] = "Val"
  m["*"] = "Ter"
  m
}

///|
/// Convert 3-letter amino acid code to 1-letter code.
/// Returns "?" for unknown codes.
pub fn amino3_to_1(code : String) -> String {
  let m = amino_3_to_1()
  match m.get(code) {
    Some(one) => one
    None => "?"
  }
}

///|
/// Convert 1-letter amino acid code to 3-letter code.
/// Returns "???" for unknown codes.
pub fn amino1_to_3(code : String) -> String {
  let m = amino_1_to_3()
  match m.get(code) {
    Some(three) => three
    None => "???"
  }
}

///|
/// Standard amino acid chemical property groups.
pub struct AminoAcidProperty {
  amino_acid : String
  property : String
  polarity : String
  charge : Double
  hydropathy : Double
}

///|
/// Get amino acid property information.
pub fn amino_acid_properties(aa : String) -> AminoAcidProperty {
  let upper = aa.to_upper()
  match upper {
    "A" => AminoAcidProperty::{ amino_acid: "A", property: "hydrophobic", polarity: "nonpolar", charge: 0.0, hydropathy: 1.8 }
    "R" => AminoAcidProperty::{ amino_acid: "R", property: "positive", polarity: "polar", charge: 1.0, hydropathy: -4.5 }
    "N" => AminoAcidProperty::{ amino_acid: "N", property: "polar", polarity: "polar", charge: 0.0, hydropathy: -3.5 }
    "D" => AminoAcidProperty::{ amino_acid: "D", property: "negative", polarity: "polar", charge: -1.0, hydropathy: -3.5 }
    "C" => AminoAcidProperty::{ amino_acid: "C", property: "polar", polarity: "polar", charge: 0.0, hydropathy: 2.5 }
    "Q" => AminoAcidProperty::{ amino_acid: "Q", property: "polar", polarity: "polar", charge: 0.0, hydropathy: -3.5 }
    "E" => AminoAcidProperty::{ amino_acid: "E", property: "negative", polarity: "polar", charge: -1.0, hydropathy: -3.5 }
    "G" => AminoAcidProperty::{ amino_acid: "G", property: "hydrophobic", polarity: "nonpolar", charge: 0.0, hydropathy: -0.4 }
    "H" => AminoAcidProperty::{ amino_acid: "H", property: "positive", polarity: "polar", charge: 1.0, hydropathy: -3.2 }
    "I" => AminoAcidProperty::{ amino_acid: "I", property: "hydrophobic", polarity: "nonpolar", charge: 0.0, hydropathy: 4.5 }
    "L" => AminoAcidProperty::{ amino_acid: "L", property: "hydrophobic", polarity: "nonpolar", charge: 0.0, hydropathy: 3.8 }
    "K" => AminoAcidProperty::{ amino_acid: "K", property: "positive", polarity: "polar", charge: 1.0, hydropathy: -3.9 }
    "M" => AminoAcidProperty::{ amino_acid: "M", property: "hydrophobic", polarity: "nonpolar", charge: 0.0, hydropathy: 1.9 }
    "F" => AminoAcidProperty::{ amino_acid: "F", property: "hydrophobic", polarity: "nonpolar", charge: 0.0, hydropathy: 2.8 }
    "P" => AminoAcidProperty::{ amino_acid: "P", property: "hydrophobic", polarity: "nonpolar", charge: 0.0, hydropathy: -1.6 }
    "S" => AminoAcidProperty::{ amino_acid: "S", property: "polar", polarity: "polar", charge: 0.0, hydropathy: -0.8 }
    "T" => AminoAcidProperty::{ amino_acid: "T", property: "polar", polarity: "polar", charge: 0.0, hydropathy: -0.7 }
    "W" => AminoAcidProperty::{ amino_acid: "W", property: "hydrophobic", polarity: "nonpolar", charge: 0.0, hydropathy: -0.9 }
    "Y" => AminoAcidProperty::{ amino_acid: "Y", property: "polar", polarity: "polar", charge: 0.0, hydropathy: -1.3 }
    "V" => AminoAcidProperty::{ amino_acid: "V", property: "hydrophobic", polarity: "nonpolar", charge: 0.0, hydropathy: 4.2 }
    _ => AminoAcidProperty::{ amino_acid: "?", property: "unknown", polarity: "unknown", charge: 0.0, hydropathy: 0.0 }
  }
}

///|
/// Get codon table for standard genetic code (Table 1).
pub fn standard_codon_table() -> Map[String, String] {
  let m : Map[String, String] = Map([], capacity=64)
  m["TTT"] = "F"
  m["TTC"] = "F"
  m["TTA"] = "L"
  m["TTG"] = "L"
  m["CTT"] = "L"
  m["CTC"] = "L"
  m["CTA"] = "L"
  m["CTG"] = "L"
  m["ATT"] = "I"
  m["ATC"] = "I"
  m["ATA"] = "I"
  m["ATG"] = "M"
  m["GTT"] = "V"
  m["GTC"] = "V"
  m["GTA"] = "V"
  m["GTG"] = "V"
  m["TCT"] = "S"
  m["TCC"] = "S"
  m["TCA"] = "S"
  m["TCG"] = "S"
  m["CCT"] = "P"
  m["CCC"] = "P"
  m["CCA"] = "P"
  m["CCG"] = "P"
  m["ACT"] = "T"
  m["ACC"] = "T"
  m["ACA"] = "T"
  m["ACG"] = "T"
  m["GCT"] = "A"
  m["GCC"] = "A"
  m["GCA"] = "A"
  m["GCG"] = "A"
  m["TAT"] = "Y"
  m["TAC"] = "Y"
  m["TAA"] = "*"
  m["TAG"] = "*"
  m["CAT"] = "H"
  m["CAC"] = "H"
  m["CAA"] = "Q"
  m["CAG"] = "Q"
  m["AAT"] = "N"
  m["AAC"] = "N"
  m["AAA"] = "K"
  m["AAG"] = "K"
  m["GAT"] = "D"
  m["GAC"] = "D"
  m["GAA"] = "E"
  m["GAG"] = "E"
  m["TGT"] = "C"
  m["TGC"] = "C"
  m["TGA"] = "*"
  m["TGG"] = "W"
  m["CGT"] = "R"
  m["CGC"] = "R"
  m["CGA"] = "R"
  m["CGG"] = "R"
  m["AGT"] = "S"
  m["AGC"] = "S"
  m["AGA"] = "R"
  m["AGG"] = "R"
  m["GGT"] = "G"
  m["GGC"] = "G"
  m["GGA"] = "G"
  m["GGG"] = "G"
  m
}

///|
/// Look up amino acid for a codon using the standard genetic code.
pub fn codon_lookup(codon : String) -> String {
  let table = standard_codon_table()
  match table.get(codon.to_upper()) {
    Some(aa) => aa
    None => "?"
  }
}

///|
/// Reverse map: amino acid -> list of synonymous codons.
pub fn reverse_codon_table() -> Map[String, Array[String]] {
  let forward = standard_codon_table()
  let result : Map[String, Array[String]] = Map([], capacity=21)

  for codon in forward.keys() {
    match forward.get(codon) {
      Some(aa) => {
        match result.get(aa) {
          Some(arr) => {
            let new_arr = arr.copy()
            new_arr.push(codon)
            result[aa] = new_arr
          }
          None => {
            result[aa] = [codon]
          }
        }
      }
      None => ()
    }
  }

  result
}

///|
/// Calculate the number of synonymous codons for an amino acid.
pub fn num_synonymous_codons(aa : String) -> Int {
  let rev = reverse_codon_table()
  match rev.get(aa) {
    Some(arr) => arr.length()
    None => 0
  }
}

///|
/// Generate the reverse complement of a DNA sequence.
/// Handles IUPAC ambiguous bases correctly.
pub fn reverse_complement_iupac(seq : String) -> String {
  let comp : Map[String, String] = Map([], capacity=16)
  comp["A"] = "T"
  comp["T"] = "A"
  comp["C"] = "G"
  comp["G"] = "C"
  comp["U"] = "A"
  comp["R"] = "Y"
  comp["Y"] = "R"
  comp["S"] = "S"
  comp["W"] = "W"
  comp["K"] = "M"
  comp["M"] = "K"
  comp["B"] = "V"
  comp["V"] = "B"
  comp["D"] = "H"
  comp["H"] = "D"
  comp["N"] = "N"

  let mut result = ""
  let mut i = seq.length() - 1
  while i >= 0 {
    let base = seq.substring(start=i, end=i + 1)
    let upper = base.to_upper()
    let comp_base = match comp.get(upper) {
      Some(c) => c
      None => base
    }
    result = result + comp_base
    i = i - 1
  }
  result
}
