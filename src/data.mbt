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
  let m : Map[String, String] = Map([], capacity=23)
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
  m["Sec"] = "U"
  m["Pyl"] = "O"
  m["Ter"] = "*"
  m
}

///|
/// Amino acid 1-letter to 3-letter code mapping.
pub fn amino_1_to_3() -> Map[String, String] {
  let m : Map[String, String] = Map([], capacity=23)
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
  m["U"] = "Sec"
  m["O"] = "Pyl"
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
    "A" =>
      AminoAcidProperty::{
        amino_acid: "A",
        property: "hydrophobic",
        polarity: "nonpolar",
        charge: 0.0,
        hydropathy: 1.8,
      }
    "R" =>
      AminoAcidProperty::{
        amino_acid: "R",
        property: "positive",
        polarity: "polar",
        charge: 1.0,
        hydropathy: -4.5,
      }
    "N" =>
      AminoAcidProperty::{
        amino_acid: "N",
        property: "polar",
        polarity: "polar",
        charge: 0.0,
        hydropathy: -3.5,
      }
    "D" =>
      AminoAcidProperty::{
        amino_acid: "D",
        property: "negative",
        polarity: "polar",
        charge: -1.0,
        hydropathy: -3.5,
      }
    "C" =>
      AminoAcidProperty::{
        amino_acid: "C",
        property: "polar",
        polarity: "polar",
        charge: 0.0,
        hydropathy: 2.5,
      }
    "Q" =>
      AminoAcidProperty::{
        amino_acid: "Q",
        property: "polar",
        polarity: "polar",
        charge: 0.0,
        hydropathy: -3.5,
      }
    "E" =>
      AminoAcidProperty::{
        amino_acid: "E",
        property: "negative",
        polarity: "polar",
        charge: -1.0,
        hydropathy: -3.5,
      }
    "G" =>
      AminoAcidProperty::{
        amino_acid: "G",
        property: "hydrophobic",
        polarity: "nonpolar",
        charge: 0.0,
        hydropathy: -0.4,
      }
    "H" =>
      AminoAcidProperty::{
        amino_acid: "H",
        property: "positive",
        polarity: "polar",
        charge: 1.0,
        hydropathy: -3.2,
      }
    "I" =>
      AminoAcidProperty::{
        amino_acid: "I",
        property: "hydrophobic",
        polarity: "nonpolar",
        charge: 0.0,
        hydropathy: 4.5,
      }
    "L" =>
      AminoAcidProperty::{
        amino_acid: "L",
        property: "hydrophobic",
        polarity: "nonpolar",
        charge: 0.0,
        hydropathy: 3.8,
      }
    "K" =>
      AminoAcidProperty::{
        amino_acid: "K",
        property: "positive",
        polarity: "polar",
        charge: 1.0,
        hydropathy: -3.9,
      }
    "M" =>
      AminoAcidProperty::{
        amino_acid: "M",
        property: "hydrophobic",
        polarity: "nonpolar",
        charge: 0.0,
        hydropathy: 1.9,
      }
    "F" =>
      AminoAcidProperty::{
        amino_acid: "F",
        property: "hydrophobic",
        polarity: "nonpolar",
        charge: 0.0,
        hydropathy: 2.8,
      }
    "P" =>
      AminoAcidProperty::{
        amino_acid: "P",
        property: "hydrophobic",
        polarity: "nonpolar",
        charge: 0.0,
        hydropathy: -1.6,
      }
    "S" =>
      AminoAcidProperty::{
        amino_acid: "S",
        property: "polar",
        polarity: "polar",
        charge: 0.0,
        hydropathy: -0.8,
      }
    "T" =>
      AminoAcidProperty::{
        amino_acid: "T",
        property: "polar",
        polarity: "polar",
        charge: 0.0,
        hydropathy: -0.7,
      }
    "W" =>
      AminoAcidProperty::{
        amino_acid: "W",
        property: "hydrophobic",
        polarity: "nonpolar",
        charge: 0.0,
        hydropathy: -0.9,
      }
    "Y" =>
      AminoAcidProperty::{
        amino_acid: "Y",
        property: "polar",
        polarity: "polar",
        charge: 0.0,
        hydropathy: -1.3,
      }
    "V" =>
      AminoAcidProperty::{
        amino_acid: "V",
        property: "hydrophobic",
        polarity: "nonpolar",
        charge: 0.0,
        hydropathy: 4.2,
      }
    _ =>
      AminoAcidProperty::{
        amino_acid: "?",
        property: "unknown",
        polarity: "unknown",
        charge: 0.0,
        hydropathy: 0.0,
      }
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
      Some(aa) =>
        match result.get(aa) {
          Some(arr) => {
            let new_arr = arr.copy()
            new_arr.push(codon)
            result[aa] = new_arr
          }
          None => result[aa] = [codon]
        }
      None => ()
    }
  }

  result
}

///|
/// Calculate the number of synonymous codons for an amino acid.
pub fn num_synonymous_codons(aa : String) -> Int32 {
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
  Seq::new(seq.to_upper()).reverse_complement().to_string()
}

///|
/// Reverse a string character-wise (O(n)).
/// A single-pass implementation used by alignment/assembly/BWT modules.
pub fn reverse_string(s : String) -> String {
  let n = s.length()
  let buf : FixedArray[UInt16] = FixedArray::make(n, 0)
  for i = 0; i < n; i = i + 1 {
    buf[i] = s.unsafe_get(n - 1 - i)
  }
  unsafe_from_uint16_fixedarray(buf)
}

///|
/// Standard stop codons (NCBI table 1), DNA and RNA forms.
pub fn stop_codon_list() -> Array[String] {
  ["TAA", "TAG", "TGA", "UAA", "UAG", "UGA"]
}

///|
/// Check whether a 3-letter codon is a standard stop codon.
/// Case-insensitive, works with both DNA (T) and RNA (U) input.
pub fn is_stop_codon(codon : String) -> Bool {
  let upper = codon.to_upper()
  upper == "TAA" ||
  upper == "TAG" ||
  upper == "TGA" ||
  upper == "UAA" ||
  upper == "UAG" ||
  upper == "UGA"
}

///|
/// Extract a 3-letter codon at position `i` and uppercase it.
/// Returns an empty string if `i` is out of bounds.
pub fn codon_at(seq : String, i : Int32) -> String {
  let n = seq.length()
  if i < 0 || i + 3 > n {
    ""
  } else {
    seq[i:i + 3].to_owned().to_upper()
  }
}

// ============================================================
// Protein molecular weight infrastructure
// Unified source for 6 modules that previously had 6 separate
// hard-coded amino acid weight tables.
// ============================================================

///|
/// Water molecule average isotopic mass (Da).
/// Used by most molecular weight calculations.
pub const WATER_MASS_AVERAGE : Double = 18.01524

///|
/// Water molecule monoisotopic mass (Da).
/// Used for monoisotopic mass calculations.
pub const WATER_MASS_MONOISOTOPIC : Double = 18.01056

///|
/// Proton mass (Da), used in charge/mass calculations.
pub const PROTON_MASS : Double = 1.007276

///|
/// Amino acid **complete molecular weights** (average isotopic).
/// These are full amino acid molecule masses (residue + water),
/// e.g. Alanine = C3H7NO2 = 89.0935 Da.
/// Source: Bio.SeqUtils.MolWt (Biopython) / mol_wt.mbt.
/// 22 entries: standard 20 + Sec (U) + Pyl (O).
pub fn protein_complete_mass_average() -> Map[String, Double] {
  let m : Map[String, Double] = Map([], capacity=22)
  m["A"] = 89.0935
  m["R"] = 174.2027
  m["N"] = 132.1184
  m["D"] = 133.1032
  m["C"] = 121.1590
  m["E"] = 147.1299
  m["Q"] = 146.1451
  m["G"] = 75.0669
  m["H"] = 155.1552
  m["I"] = 131.1736
  m["L"] = 131.1736
  m["K"] = 146.1882
  m["M"] = 149.2124
  m["F"] = 165.1900
  m["P"] = 115.1310
  m["S"] = 105.0930
  m["T"] = 119.1197
  m["W"] = 204.2262
  m["Y"] = 181.1894
  m["V"] = 117.1469
  m["U"] = 150.0388 /// Selenocysteine
  m["O"] = 237.2984 /// Pyrrolysine
  m
}

///|
/// Amino acid **residue masses** (average isotopic).
/// Complete mass minus one water molecule (18.01524).
/// This is the mass contribution of an amino acid in a polypeptide chain.
/// 20 standard amino acids only (U/O rarely used in residue mass context).
pub fn protein_residue_mass_average() -> Map[String, Double] {
  let complete = protein_complete_mass_average()
  let m : Map[String, Double] = Map([], capacity=22)
  for (aa, mass) in complete {
    m[aa] = mass - WATER_MASS_AVERAGE
  }
  m
}

///|
/// Amino acid **residue masses** (monoisotopic).
/// Source: IUPAC / Bio.Data.IUPACData (Biopython),
/// same values used by proteomics.mbt and bp3_iupac_pdbdata.mbt.
/// 20 standard amino acids.
pub fn protein_residue_mass_monoisotopic() -> Map[String, Double] {
  let m : Map[String, Double] = Map([], capacity=20)
  m["A"] = 71.03711
  m["C"] = 103.00919
  m["D"] = 115.02694
  m["E"] = 129.04259
  m["F"] = 147.06841
  m["G"] = 57.02146
  m["H"] = 137.05891
  m["I"] = 113.08406
  m["K"] = 128.09496
  m["L"] = 113.08406
  m["M"] = 131.04049
  m["N"] = 114.04293
  m["P"] = 97.05276
  m["Q"] = 128.05858
  m["R"] = 156.10111
  m["S"] = 87.03203
  m["T"] = 101.04768
  m["V"] = 99.06841
  m["W"] = 186.07931
  m["Y"] = 163.06333
  m
}

///|
/// Look up average isotopic complete mass for one amino acid (1-letter code).
/// Case-insensitive. Returns 0.0 for unknown codes.
pub fn amino_complete_mass(aa : String) -> Double {
  let table = protein_complete_mass_average()
  table.get(aa.to_upper()).unwrap_or(0.0)
}

///|
/// Look up average isotopic residue mass for one amino acid.
/// Case-insensitive. Returns 0.0 for unknown codes.
pub fn amino_residue_mass(aa : String) -> Double {
  let table = protein_residue_mass_average()
  table.get(aa.to_upper()).unwrap_or(0.0)
}

///|
/// Look up monoisotopic residue mass for one amino acid.
/// Case-insensitive. Returns 0.0 for unknown codes.
pub fn amino_residue_mass_mono(aa : String) -> Double {
  let table = protein_residue_mass_monoisotopic()
  table.get(aa.to_upper()).unwrap_or(0.0)
}

///|
/// Unified protein molecular weight calculator.
/// Uses the standard formula: sum(complete masses) - (n-1) * water.
/// This is mathematically equivalent to: water + sum(residue masses).
///
/// By default uses average isotopic masses and WATER_MASS_AVERAGE.
/// When is_monoisotopic=true, uses monoisotopic residue masses + WATER_MASS_MONOISOTOPIC.
pub fn protein_molecular_weight(
  seq : String,
  is_monoisotopic~ : Bool = false,
) -> Double {
  let n = seq.length()
  if n == 0 {
    return 0.0
  }
  if is_monoisotopic {
    let residues = protein_residue_mass_monoisotopic()
    let mut mass = WATER_MASS_MONOISOTOPIC
    let mut i = 0
    while i < n {
      let aa = seq[i:i + 1].to_owned().to_upper()
      mass = mass + residues.get(aa).unwrap_or(0.0)
      i = i + 1
    }
    mass
  } else {
    let complete = protein_complete_mass_average()
    let mut total = 0.0
    let mut i = 0
    while i < n {
      let aa = seq[i:i + 1].to_owned().to_upper()
      total = total + complete.get(aa).unwrap_or(0.0)
      i = i + 1
    }
    total - (n - 1).to_double() * WATER_MASS_AVERAGE
  }
}

///|
/// DNA nucleotide molecular weights (average isotopic, per-nucleotide contribution).
/// Source: Bio.SeqUtils.MolWt (Biopython). Values are dNTP masses
/// (base + deoxyribose + triphosphate), used as per-residue weights.
/// NOTE: DNA weight tables vary by source (emboss uses 313.21 for A,
/// while mol_wt uses 335.2444). This table follows the mol_wt convention.
pub fn dna_nucleotide_mass() -> Map[String, Double] {
  let m : Map[String, Double] = Map([], capacity=16)
  m["A"] = 335.2444
  m["C"] = 311.2133
  m["G"] = 351.2407
  m["T"] = 326.2079
  m["U"] = 324.1968
  m["I"] = 350.2559 /// Inosine
  m["X"] = 335.2444 /// Unknown purine → A
  m["Y"] = 326.2079 /// Unknown pyrimidine → T
  m["N"] = 335.2444 /// Any → average
  m["R"] = 343.2426 /// A|G
  m["S"] = 331.2270 /// G|C
  m["W"] = 330.7262 /// A|T
  m["K"] = 338.7243 /// G|T
  m["M"] = 323.2289 /// A|C
  m["B"] = 318.7076 /// C|G|T
  m["D"] = 337.5837 /// A|G|T
  m
}
