///|
/// Bio.Variation module, inspired by Biopython Bio.Variation.
/// Provides functions for variation analysis, SNP detection, and mutation annotation.

///|
pub struct SNP {
  position : Int
  ref_allele : String
  alt_allele : String
  quality : Double
  filter : String
}

///|
pub fn SNP::new(
  position : Int,
  ref_allele : String,
  alt_allele : String,
  quality : Double,
  filter : String,
) -> SNP {
  SNP::{ position, ref_allele, alt_allele, quality, filter }
}

///|
pub struct Mutation {
  gene : String
  protein : String
  position : Int
  original_aa : String
  mutated_aa : String
  mutation_type : String
  impact : String
}

///|
pub fn Mutation::new(
  gene : String,
  protein : String,
  position : Int,
  original_aa : String,
  mutated_aa : String,
  mutation_type : String,
  impact : String,
) -> Mutation {
  Mutation::{ gene, protein, position, original_aa, mutated_aa, mutation_type, impact }
}

///|
pub struct VariationRecord {
  chromosome : String
  position : Int
  ref_allele : String
  alt_allele : String
  mutations : Array[Mutation]
}

///|
pub fn VariationRecord::new(
  chromosome : String,
  position : Int,
  ref_allele : String,
  alt_allele : String,
  mutations : Array[Mutation],
) -> VariationRecord {
  VariationRecord::{ chromosome, position, ref_allele, alt_allele, mutations }
}

///|
pub struct AminoAcidSubstitution {
  original : String
  mutated : String
  blosum62_score : Int
  grantham_distance : Double
}

///|
pub fn AminoAcidSubstitution::new(
  original : String,
  mutated : String,
  blosum62_score : Int,
  grantham_distance : Double,
) -> AminoAcidSubstitution {
  AminoAcidSubstitution::{ original, mutated, blosum62_score, grantham_distance }
}

///|
fn blosum62_matrix() -> Map[String, Map[String, Int]] {
  let matrix : Map[String, Map[String, Int]] = Map([], capacity=20)
  
  let amino_acids = ["A", "R", "N", "D", "C", "Q", "E", "G", "H", "I", "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V"]
  
  let scores = [
    [4, -1, -2, -2, 0, -1, -1, 0, -2, -1, -1, -1, -1, -2, -1, 1, 0, -3, -2, 0],
    [-1, 5, 0, -2, -3, 1, 0, -2, 0, -3, -2, 2, -1, -3, -2, -1, -1, -3, -2, -3],
    [-2, 0, 6, 1, -3, 0, 0, 0, 1, -3, -3, 0, -2, -3, -2, 1, 0, -4, -2, -3],
    [-2, -2, 1, 6, -3, 0, 2, -1, -1, -3, -4, -1, -3, -3, -1, 0, -1, -4, -3, -3],
    [0, -3, -3, -3, 9, -3, -4, -3, -3, -1, -1, -3, -1, -2, -3, -1, -1, -2, -2, -1],
    [-1, 1, 0, 0, -3, 5, 2, -2, 0, -3, -2, 1, 0, -3, -1, 0, -1, -2, -1, -2],
    [-1, 0, 0, 2, -4, 2, 5, -2, 0, -3, -3, 1, -2, -3, -1, 0, -1, -3, -2, -2],
    [0, -2, 0, -1, -3, -2, -2, 6, -2, -4, -4, -2, -3, -3, -2, 0, -2, -2, -3, -3],
    [-2, 0, 1, -1, -3, 0, 0, -2, 8, -3, -3, -1, -2, -1, -2, -1, -2, -2, 2, -3],
    [-1, -3, -3, -3, -1, -3, -3, -4, -3, 4, 2, -3, 1, 0, -3, -2, -1, -3, -1, 3],
    [-1, -2, -3, -4, -1, -2, -3, -4, -3, 2, 4, -2, 2, 0, -3, -2, -1, -2, -1, 1],
    [-1, 2, 0, -1, -3, 1, 1, -2, -1, -3, -2, 5, -1, -3, -1, 0, -1, -3, -2, -2],
    [-1, -1, -2, -3, -1, 0, -2, -3, -2, 1, 2, -1, 5, 0, -2, -1, -1, -1, -1, 1],
    [-2, -3, -3, -3, -2, -3, -3, -3, -1, 0, 0, -3, 0, 6, -4, -2, -2, 1, 3, -1],
    [-1, -2, -2, -1, -3, -1, -1, -2, -2, -3, -3, -1, -2, -4, 7, -1, -1, -4, -3, -2],
    [1, -1, 1, 0, -1, 0, 0, 0, -1, -2, -2, 0, -1, -2, -1, 4, 1, -3, -2, -2],
    [0, -1, 0, -1, -1, -1, -1, -2, -2, -1, -1, -1, -1, -2, -1, 1, 5, -2, -2, 0],
    [-3, -3, -4, -4, -2, -2, -3, -2, -2, -3, -2, -3, -1, 1, -4, -3, -2, 11, 2, -3],
    [-2, -2, -2, -3, -2, -1, -2, -3, 2, -1, -1, -2, -1, 3, -3, -2, -2, 2, 7, -1],
    [0, -3, -3, -3, -1, -2, -2, -3, -3, 3, 1, -2, 1, -1, -2, -2, 0, -3, -1, 4],
  ]
  
  let mut i = 0
  while i < amino_acids.length() {
    let row_map : Map[String, Int] = Map([], capacity=20)
    let mut j = 0
    while j < amino_acids.length() {
      row_map.set(amino_acids[j], scores[i][j])
      j = j + 1
    }
    matrix.set(amino_acids[i], row_map)
    i = i + 1
  }
  
  matrix
}

///|
fn grantham_matrix() -> Map[String, Map[String, Double]] {
  let matrix : Map[String, Map[String, Double]] = Map([], capacity=20)
  
  let amino_acids = ["A", "R", "N", "D", "C", "Q", "E", "G", "H", "I", "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V"]
  
  let distances = [
    [0.0, 112.0, 145.0, 145.0, 130.0, 119.0, 118.0, 75.0, 127.0, 96.0, 96.0, 125.0, 107.0, 130.0, 97.0, 56.0, 60.0, 157.0, 147.0, 90.0],
    [112.0, 0.0, 85.0, 85.0, 145.0, 60.0, 64.0, 125.0, 48.0, 142.0, 142.0, 26.0, 130.0, 165.0, 143.0, 97.0, 101.0, 153.0, 152.0, 136.0],
    [145.0, 85.0, 0.0, 20.0, 157.0, 38.0, 55.0, 110.0, 66.0, 171.0, 171.0, 70.0, 159.0, 189.0, 157.0, 82.0, 86.0, 180.0, 168.0, 165.0],
    [145.0, 85.0, 20.0, 0.0, 157.0, 58.0, 35.0, 110.0, 86.0, 171.0, 171.0, 70.0, 159.0, 189.0, 157.0, 82.0, 86.0, 180.0, 168.0, 165.0],
    [130.0, 145.0, 157.0, 157.0, 0.0, 148.0, 147.0, 95.0, 145.0, 108.0, 108.0, 153.0, 119.0, 126.0, 127.0, 78.0, 82.0, 144.0, 134.0, 102.0],
    [119.0, 60.0, 38.0, 58.0, 148.0, 0.0, 25.0, 104.0, 46.0, 151.0, 151.0, 45.0, 139.0, 178.0, 146.0, 73.0, 77.0, 169.0, 157.0, 145.0],
    [118.0, 64.0, 55.0, 35.0, 147.0, 25.0, 0.0, 103.0, 66.0, 150.0, 150.0, 49.0, 138.0, 177.0, 145.0, 72.0, 76.0, 168.0, 156.0, 144.0],
    [75.0, 125.0, 110.0, 110.0, 95.0, 104.0, 103.0, 0.0, 112.0, 91.0, 91.0, 130.0, 102.0, 125.0, 87.0, 60.0, 64.0, 152.0, 142.0, 85.0],
    [127.0, 48.0, 66.0, 86.0, 145.0, 46.0, 66.0, 112.0, 0.0, 149.0, 149.0, 58.0, 137.0, 175.0, 154.0, 89.0, 93.0, 166.0, 148.0, 143.0],
    [96.0, 142.0, 171.0, 171.0, 108.0, 151.0, 150.0, 91.0, 149.0, 0.0, 9.0, 159.0, 21.0, 41.0, 112.0, 100.0, 66.0, 84.0, 66.0, 28.0],
    [96.0, 142.0, 171.0, 171.0, 108.0, 151.0, 150.0, 91.0, 149.0, 9.0, 0.0, 159.0, 12.0, 32.0, 112.0, 100.0, 66.0, 75.0, 57.0, 19.0],
    [125.0, 26.0, 70.0, 70.0, 153.0, 45.0, 49.0, 130.0, 58.0, 159.0, 159.0, 0.0, 147.0, 182.0, 160.0, 104.0, 108.0, 167.0, 166.0, 153.0],
    [107.0, 130.0, 159.0, 159.0, 119.0, 139.0, 138.0, 102.0, 137.0, 21.0, 12.0, 147.0, 0.0, 53.0, 123.0, 112.0, 78.0, 96.0, 78.0, 40.0],
    [130.0, 165.0, 189.0, 189.0, 126.0, 178.0, 177.0, 125.0, 175.0, 41.0, 32.0, 182.0, 53.0, 0.0, 143.0, 133.0, 99.0, 67.0, 17.0, 63.0],
    [97.0, 143.0, 157.0, 157.0, 127.0, 146.0, 145.0, 87.0, 154.0, 112.0, 112.0, 160.0, 123.0, 143.0, 0.0, 63.0, 67.0, 174.0, 164.0, 106.0],
    [56.0, 97.0, 82.0, 82.0, 78.0, 73.0, 72.0, 60.0, 89.0, 100.0, 100.0, 104.0, 112.0, 133.0, 63.0, 0.0, 19.0, 145.0, 135.0, 94.0],
    [60.0, 101.0, 86.0, 86.0, 82.0, 77.0, 76.0, 64.0, 93.0, 66.0, 66.0, 108.0, 78.0, 99.0, 67.0, 19.0, 0.0, 149.0, 139.0, 60.0],
    [157.0, 153.0, 180.0, 180.0, 144.0, 169.0, 168.0, 152.0, 166.0, 84.0, 75.0, 167.0, 96.0, 67.0, 174.0, 145.0, 149.0, 0.0, 80.0, 106.0],
    [147.0, 152.0, 168.0, 168.0, 134.0, 157.0, 156.0, 142.0, 148.0, 66.0, 57.0, 166.0, 78.0, 17.0, 164.0, 135.0, 139.0, 80.0, 0.0, 88.0],
    [90.0, 136.0, 165.0, 165.0, 102.0, 145.0, 144.0, 85.0, 143.0, 28.0, 19.0, 153.0, 40.0, 63.0, 106.0, 94.0, 60.0, 106.0, 88.0, 0.0],
  ]
  
  let mut i = 0
  while i < amino_acids.length() {
    let row_map : Map[String, Double] = Map([], capacity=20)
    let mut j = 0
    while j < amino_acids.length() {
      row_map.set(amino_acids[j], distances[i][j])
      j = j + 1
    }
    matrix.set(amino_acids[i], row_map)
    i = i + 1
  }
  
  matrix
}

///|
pub fn analyze_snp(
  position : Int,
  ref_allele : String,
  alt_allele : String,
  quality : Double,
  filter : String,
) -> SNP {
  SNP::new(position, ref_allele, alt_allele, quality, filter)
}

///|
pub fn detect_mutations(
  gene : String,
  protein : String,
  dna_sequence : String,
  variant_position : Int,
  alt_allele : String,
) -> Array[Mutation] {
  let mutations : Array[Mutation] = Array::new()
  
  let codon_start = (variant_position - 1) / 3 * 3
  let codon_end = codon_start + 3
  
  if codon_end <= dna_sequence.length() {
    let original_codon = substring(dna_sequence, codon_start, 3)
    let new_sequence = replace_char_at(dna_sequence, variant_position - 1, alt_allele)
    let mutated_codon = substring(new_sequence, codon_start, 3)
    
    let original_aa = variation_translate_codon(original_codon)
    let mutated_aa = variation_translate_codon(mutated_codon)
    
    let mutation_type = classify_mutation(original_aa, mutated_aa)
    let impact = predict_impact(original_aa, mutated_aa, mutation_type)
    
    mutations.push(Mutation::new(
      gene, protein,
      codon_start / 3 + 1,
      original_aa, mutated_aa,
      mutation_type, impact,
    ))
  }
  
  mutations
}

///|
fn replace_char_at(s : String, idx : Int, replacement : String) -> String {
  if idx < 0 || idx >= s.length() || replacement.length() == 0 {
    return s
  }
  
  let prefix = substring(s, 0, idx)
  let suffix = substring(s, idx + 1, s.length() - idx - 1)
  
  prefix + substring(replacement, 0, 1) + suffix
}

///|
fn variation_translate_codon(codon : String) -> String {
  let codon_upper = variation_to_uppercase(codon)
  
  let codon_table = Map([
    ("TTT", "F"), ("TTC", "F"), ("TTA", "L"), ("TTG", "L"),
    ("CTT", "L"), ("CTC", "L"), ("CTA", "L"), ("CTG", "L"),
    ("ATT", "I"), ("ATC", "I"), ("ATA", "I"), ("ATG", "M"),
    ("GTT", "V"), ("GTC", "V"), ("GTA", "V"), ("GTG", "V"),
    ("TCT", "S"), ("TCC", "S"), ("TCA", "S"), ("TCG", "S"),
    ("CCT", "P"), ("CCC", "P"), ("CCA", "P"), ("CCG", "P"),
    ("ACT", "T"), ("ACC", "T"), ("ACA", "T"), ("ACG", "T"),
    ("GCT", "A"), ("GCC", "A"), ("GCA", "A"), ("GCG", "A"),
    ("TAT", "Y"), ("TAC", "Y"), ("TAA", "*"), ("TAG", "*"),
    ("CAT", "H"), ("CAC", "H"), ("CAA", "Q"), ("CAG", "Q"),
    ("AAT", "N"), ("AAC", "N"), ("AAA", "K"), ("AAG", "K"),
    ("GAT", "D"), ("GAC", "D"), ("GAA", "E"), ("GAG", "E"),
    ("TGT", "C"), ("TGC", "C"), ("TGA", "*"), ("TGG", "W"),
    ("CGT", "R"), ("CGC", "R"), ("CGA", "R"), ("CGG", "R"),
    ("AGT", "S"), ("AGC", "S"), ("AGA", "R"), ("AGG", "R"),
    ("GGT", "G"), ("GGC", "G"), ("GGA", "G"), ("GGG", "G"),
  ], capacity=64)
  
  codon_table.get(codon_upper).unwrap_or("?")
}

///|
fn variation_to_uppercase(s : String) -> String {
  let mut result = ""
  let mut i = 0
  while i < s.length() {
    let c = s.unsafe_get(i)
    let c_int = c.to_int()
    if c_int >= 'a'.to_int() && c_int <= 'z'.to_int() {
      result = result + (c_int - 'a'.to_int() + 'A'.to_int()).to_char().to_string()
    } else {
      result = result + c.to_string()
    }
    i = i + 1
  }
  result
}

///|
fn classify_mutation(original : String, mutated : String) -> String {
  if original == mutated {
    return "Silent"
  }
  if mutated == "*" {
    return "Nonsense"
  }
  if original == "*" {
    return "Frameshift"
  }
  
  let aa_groups = [
    ["G", "A", "V", "L", "I", "M"],
    ["F", "Y", "W"],
    ["S", "T", "N", "Q"],
    ["D", "E"],
    ["R", "K", "H"],
    ["C", "P"],
  ]
  
  let mut original_group = -1
  let mut mutated_group = -1
  
  let mut i = 0
  while i < aa_groups.length() {
    let mut j = 0
    while j < aa_groups[i].length() {
      if aa_groups[i][j] == original {
        original_group = i
      }
      if aa_groups[i][j] == mutated {
        mutated_group = i
      }
      j = j + 1
    }
    i = i + 1
  }
  
  if original_group == mutated_group {
    "Conservative"
  } else {
    "Non-conservative"
  }
}

///|
fn predict_impact(original : String, mutated : String, mutation_type : String) -> String {
  if mutation_type == "Nonsense" {
    return "High"
  }
  if mutation_type == "Frameshift" {
    return "High"
  }
  if mutation_type == "Non-conservative" {
    return "Medium"
  }
  if mutation_type == "Conservative" {
    return "Low"
  }
  "Low"
}

///|
pub fn analyze_amino_acid_substitution(original : String, mutated : String) -> AminoAcidSubstitution {
  let blosum62 = blosum62_matrix()
  let grantham = grantham_matrix()
  
  let blosum_val = blosum62.get(original)
  let blosum_score = if blosum_val.is_some() {
    blosum_val.unwrap().get(mutated).unwrap_or(0)
  } else {
    0
  }
  
  let grantham_val = grantham.get(original)
  let grantham_dist = if grantham_val.is_some() {
    grantham_val.unwrap().get(mutated).unwrap_or(0.0)
  } else {
    0.0
  }
  
  AminoAcidSubstitution::new(original, mutated, blosum_score, grantham_dist)
}

///|
pub fn predict_damaging(substitution : AminoAcidSubstitution) -> Bool {
  substitution.blosum62_score < 0 || substitution.grantham_distance > 100.0
}

///|
pub fn parse_vcf_line(line : String) -> VariationRecord? {
  let parts = split_by_tab(line)
  if parts.length() < 5 {
    return None
  }
  
  let chromosome = parts[0]
  let position_str = parts[1]
  let ref_allele = parts[3]
  let alt_allele = parts[4]
  
  let position = parse_int(position_str)
  
  Some(VariationRecord::new(chromosome, position, ref_allele, alt_allele, Array::new()))
}

///|
pub fn create_example_variation_data() -> Array[VariationRecord] {
  let records : Array[VariationRecord] = Array::new()
  
  let mutations1 = [
    Mutation::new("BRCA1", "BRCA1_HUMAN", 178, "G", "A", "Non-conservative", "Medium"),
    Mutation::new("BRCA1", "BRCA1_HUMAN", 356, "R", "*", "Nonsense", "High"),
  ]
  records.push(VariationRecord::new("17", 43044295, "C", "T", mutations1))
  
  let mutations2 = [
    Mutation::new("TP53", "P53_HUMAN", 273, "R", "H", "Conservative", "Low"),
  ]
  records.push(VariationRecord::new("17", 7577538, "G", "A", mutations2))
  
  records
}