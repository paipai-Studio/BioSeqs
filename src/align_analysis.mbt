///|
/// Bio.Align.analysis module - Evolutionary analysis from alignments.
/// Inspired by Biopython's Bio.Align.analysis module.
/// Provides dn/ds calculations, evolutionary distances, and selection pressure analysis.

///|
/// Represents a dn/ds analysis result.
pub struct AlnAnalysisResult {
  dn : Double
  ds : Double
  dn_ds_ratio : Double
  dn_ds_ratio_sem : Double
  n_synonymous : Int
  n_nonsynonymous : Int
}

///|
/// Create a new dn/ds result.
pub fn AlnAnalysisResult::new(
  dn : Double,
  ds : Double,
  dn_ds_ratio : Double,
  dn_ds_ratio_sem : Double,
  n_synonymous : Int,
  n_nonsynonymous : Int
) -> AlnAnalysisResult {
  AlnAnalysisResult::{ dn, ds, dn_ds_ratio, dn_ds_ratio_sem, n_synonymous, n_nonsynonymous }
}

///|
/// Get dn (non-synonymous substitution rate).
pub fn aln_get_dn(self : AlnAnalysisResult) -> Double {
  self.dn
}

///|
/// Get ds (synonymous substitution rate).
pub fn aln_get_ds(self : AlnAnalysisResult) -> Double {
  self.ds
}

///|
/// Get dn/ds ratio (omega).
pub fn aln_get_ratio(self : AlnAnalysisResult) -> Double {
  self.dn_ds_ratio
}

///|
/// Get standard error of mean for dn/ds ratio.
pub fn aln_get_ratio_sem(self : AlnAnalysisResult) -> Double {
  self.dn_ds_ratio_sem
}

///|
/// Get number of synonymous substitutions.
pub fn aln_get_n_syn(self : AlnAnalysisResult) -> Int {
  self.n_synonymous
}

///|
/// Get number of non-synonymous substitutions.
pub fn aln_get_n_nonsyn(self : AlnAnalysisResult) -> Int {
  self.n_nonsynonymous
}

///|
/// Check if the result indicates positive selection (omega > 1).
pub fn aln_has_positive_selection(self : AlnAnalysisResult) -> Bool {
  self.dn_ds_ratio > 1.0
}

///|
/// Check if the result indicates purifying selection (omega < 1).
pub fn aln_has_purifying_selection(self : AlnAnalysisResult) -> Bool {
  self.dn_ds_ratio < 1.0
}

///|
/// Simple genetic code table (standard code 1).
/// Returns amino acid for a codon.
pub fn aln_aa_from_codon(codon : String) -> String {
  if codon == "TTT" || codon == "TTC" { "F" }
  else if codon == "TTA" || codon == "TTG" { "L" }
  else if codon == "CTT" || codon == "CTC" || codon == "CTA" || codon == "CTG" { "L" }
  else if codon == "ATT" || codon == "ATC" || codon == "ATA" { "I" }
  else if codon == "ATG" { "M" }
  else if codon == "GTT" || codon == "GTC" || codon == "GTA" || codon == "GTG" { "V" }
  else if codon == "TCT" || codon == "TCC" || codon == "TCA" || codon == "TCG" { "S" }
  else if codon == "CCT" || codon == "CCC" || codon == "CCA" || codon == "CCG" { "P" }
  else if codon == "ACT" || codon == "ACC" || codon == "ACA" || codon == "ACG" { "T" }
  else if codon == "GCT" || codon == "GCC" || codon == "GCA" || codon == "GCG" { "A" }
  else if codon == "TAT" || codon == "TAC" { "Y" }
  else if codon == "TAA" || codon == "TAG" || codon == "TGA" { "*" }
  else if codon == "CAT" || codon == "CAC" { "H" }
  else if codon == "CAA" || codon == "CAG" { "Q" }
  else if codon == "AAT" || codon == "AAC" { "N" }
  else if codon == "AAA" || codon == "AAG" { "K" }
  else if codon == "GAT" || codon == "GAC" { "D" }
  else if codon == "GAA" || codon == "GAG" { "E" }
  else if codon == "TGT" || codon == "TGC" { "C" }
  else if codon == "TGG" { "W" }
  else if codon == "CGT" || codon == "CGC" || codon == "CGA" || codon == "CGG" { "R" }
  else if codon == "AGT" || codon == "AGC" { "S" }
  else if codon == "AGA" || codon == "AGG" { "R" }
  else if codon == "GGT" || codon == "GGC" || codon == "GGA" || codon == "GGG" { "G" }
  else { "?" }
}

///|
/// Classify a nucleotide substitution as synonymous or non-synonymous.
/// Returns true if the substitution is synonymous (doesn't change amino acid).
pub fn aln_is_synonymous(codon : String, pos : Int, new_base : String) -> Bool {
  let first_part = if pos > 0 { substring(codon, 0, pos) } else { "" }
  let second_part = if pos < 2 { substring(codon, pos + 1, 3 - pos - 1) } else { "" }
  let new_codon = first_part + new_base + second_part
  let old_aa = aln_aa_from_codon(codon)
  let new_aa = aln_aa_from_codon(new_codon)
  old_aa == new_aa
}

///|
/// Calculate the number of synonymous and non-synonymous sites.
pub fn aln_count_sites(sequence : String) -> Array[Int] {
  let seq_len = sequence.length()
  let n_codons = seq_len / 3
  
  let mut n_syn = 0
  let mut n_nonsyn = 0
  
  let bases = ["A", "C", "G", "T"]
  
  let mut codon_idx = 0
  while codon_idx < n_codons {
    let codon_start = codon_idx * 3
    let codon = substring(sequence, codon_start, 3)
    
    // Check each position in the codon
    let mut pos = 0
    while pos < 3 {
      let old_base = substring(codon, pos, 1)
      let mut b = 0
      while b < bases.length() {
        let new_base = bases[b]
        if new_base != old_base {
          if aln_is_synonymous(codon, pos, new_base) {
            n_syn = n_syn + 1
          } else {
            n_nonsyn = n_nonsyn + 1
          }
        }
        b = b + 1
      }
      pos = pos + 1
    }
    codon_idx = codon_idx + 1
  }
  
  // Each site has 3 possible changes (excluding the original), so divide by 3
  let sites = [n_syn / 3, n_nonsyn / 3]
  sites
}

///|
/// Calculate dn/ds ratio between two sequences using the NG86 method.
/// Nei-Gojobori method for estimating dn and ds.
pub fn aln_analyze_dn_ds(seq1 : String, seq2 : String) -> AlnAnalysisResult {
  let len = if seq1.length() < seq2.length() { seq1.length() } else { seq2.length() }
  let n_codons = len / 3
  
  let sites = aln_count_sites(seq1)
  let total_syn_sites = sites[0].to_double()
  let total_nonsyn_sites = sites[1].to_double()
  
  if total_syn_sites == 0.0 || total_nonsyn_sites == 0.0 {
    return AlnAnalysisResult::new(0.0, 0.0, 0.0, 0.0, 0, 0)
  }
  
  let mut n_syn_changes = 0.0
  let mut n_nonsyn_changes = 0.0
  
  let mut codon_idx = 0
  while codon_idx < n_codons {
    let codon_start = codon_idx * 3
    let codon1 = substring(seq1, codon_start, 3)
    let codon2 = substring(seq2, codon_start, 3)
    
    // Check each position in the codon
    let mut pos = 0
    while pos < 3 {
      let base1 = substring(codon1, pos, 1)
      let base2 = substring(codon2, pos, 1)
      
      if base1 != base2 {
        if aln_is_synonymous(codon1, pos, base2) {
          n_syn_changes = n_syn_changes + 1.0
        } else {
          n_nonsyn_changes = n_nonsyn_changes + 1.0
        }
      }
      pos = pos + 1
    }
    codon_idx = codon_idx + 1
  }
  
  // Calculate rates
  let ds = if total_syn_sites > 0.0 {
    n_syn_changes / total_syn_sites
  } else {
    0.0
  }
  
  let dn = if total_nonsyn_sites > 0.0 {
    n_nonsyn_changes / total_nonsyn_sites
  } else {
    0.0
  }
  
  // Calculate ratio
  let ratio = if ds > 0.0 { dn / ds } else { 0.0 }
  
  // Calculate SEM (simplified)
  let variance = if total_syn_sites > 0.0 && total_nonsyn_sites > 0.0 {
    (1.0 / total_syn_sites) + (1.0 / total_nonsyn_sites)
  } else {
    1.0
  }
  let sem = ratio * variance.sqrt()
  
  AlnAnalysisResult::new(
    dn,
    ds,
    ratio,
    sem,
    n_syn_changes.to_int(),
    n_nonsyn_changes.to_int()
  )
}

///|
/// Calculate dn/ds matrix for multiple sequence comparisons.
pub fn aln_calculate_dn_ds_matrix(sequences : Array[String]) -> Array[Array[AlnAnalysisResult]] {
  let n = sequences.length()
  let matrix : Array[Array[AlnAnalysisResult]] = Array::new()
  
  let mut i = 0
  while i < n {
    let row : Array[AlnAnalysisResult] = Array::new()
    let mut j = 0
    while j < n {
      if i == j {
        // Diagonal: self-comparison
        row.push(AlnAnalysisResult::new(0.0, 0.0, 0.0, 0.0, 0, 0))
      } else if j > i {
        let result = aln_analyze_dn_ds(sequences[i], sequences[j])
        row.push(result)
      } else {
        // Symmetric: reuse previous result
        row.push(AlnAnalysisResult::new(0.0, 0.0, 0.0, 0.0, 0, 0))
      }
      j = j + 1
    }
    matrix.push(row)
    i = i + 1
  }
  
  // Fill symmetric entries
  let mut x = 0
  while x < n {
    let mut y = x + 1
    while y < n {
      matrix[y].set(x, matrix[x][y])
      y = y + 1
    }
    x = x + 1
  }
  
  matrix
}

///|
/// Calculate Jukes-Cantor distance between two sequences.
pub fn aln_jukes_cantor_distance(seq1 : String, seq2 : String) -> Double {
  let len = if seq1.length() < seq2.length() { seq1.length() } else { seq2.length() }
  
  let mut diffs = 0
  let mut i = 0
  while i < len {
    if substring(seq1, i, 1) != substring(seq2, i, 1) {
      diffs = diffs + 1
    }
    i = i + 1
  }
  
  let p = diffs.to_double() / len.to_double()
  if p >= 0.75 {
    // Distance is too large for JC correction
    10.0
  } else {
    let arg = 1.0 - (4.0 / 3.0) * p
    if arg > 0.0 {
      -0.75 * @math.ln(arg)
    } else {
      10.0
    }
  }
}

///|
/// Calculate Kimura 2-parameter distance.
pub fn aln_kimura_2p_distance(seq1 : String, seq2 : String) -> Double {
  let len = if seq1.length() < seq2.length() { seq1.length() } else { seq2.length() }
  
  let mut transitions = 0
  let mut transversions = 0
  
  let mut i = 0
  while i < len {
    let base1 = substring(seq1, i, 1)
    let base2 = substring(seq2, i, 1)
    
    if base1 != base2 {
      let is_base1_purine = if base1 == "A" || base1 == "G" { true } else { false }
      let is_base2_purine = if base2 == "A" || base2 == "G" { true } else { false }
      
      if is_base1_purine == is_base2_purine {
        // Same type (both purine or both pyrimidine) = transition
        transitions = transitions + 1
      } else {
        // Different type = transversion
        transversions = transversions + 1
      }
    }
    i = i + 1
  }
  
  let p = transitions.to_double() / len.to_double()
  let q = transversions.to_double() / len.to_double()
  
  let arg1 = 1.0 - 2.0 * p - q
  let arg2 = 1.0 - 2.0 * q
  
  if arg1 <= 0.0 || arg2 <= 0.0 {
    10.0
  } else {
    -0.5 * @math.ln(arg1) - 0.25 * @math.ln(arg2)
  }
}

///|
/// Calculate distance matrix using Jukes-Cantor method.
pub fn aln_jukes_cantor_matrix(sequences : Array[String]) -> Array[Array[Double]] {
  let n = sequences.length()
  let matrix : Array[Array[Double]] = Array::new()
  
  let mut i = 0
  while i < n {
    let row : Array[Double] = Array::make(n, 0.0)
    let mut j = i + 1
    while j < n {
      let dist = aln_jukes_cantor_distance(sequences[i], sequences[j])
      row[j] = dist
      j = j + 1
    }
    matrix.push(row)
    i = i + 1
  }
  
  // Fill symmetric part
  let mut x = 0
  while x < n {
    let mut y = x + 1
    while y < n {
      matrix[y].set(x, matrix[x][y])
      y = y + 1
    }
    x = x + 1
  }
  
  matrix
}

///|
/// Calculate distance matrix using Kimura 2-parameter method.
pub fn aln_kimura_2p_matrix(sequences : Array[String]) -> Array[Array[Double]] {
  let n = sequences.length()
  let matrix : Array[Array[Double]] = Array::new()
  
  let mut i = 0
  while i < n {
    let row : Array[Double] = Array::make(n, 0.0)
    let mut j = i + 1
    while j < n {
      let dist = aln_kimura_2p_distance(sequences[i], sequences[j])
      row[j] = dist
      j = j + 1
    }
    matrix.push(row)
    i = i + 1
  }
  
  // Fill symmetric part
  let mut x = 0
  while x < n {
    let mut y = x + 1
    while y < n {
      matrix[y].set(x, matrix[x][y])
      y = y + 1
    }
    x = x + 1
  }
  
  matrix
}

///|
/// Calculate number of substitutions per site (simple p-distance).
pub fn aln_p_distance(seq1 : String, seq2 : String) -> Double {
  let len = if seq1.length() < seq2.length() { seq1.length() } else { seq2.length() }
  
  let mut diffs = 0
  let mut i = 0
  while i < len {
    if substring(seq1, i, 1) != substring(seq2, i, 1) {
      diffs = diffs + 1
    }
    i = i + 1
  }
  
  diffs.to_double() / len.to_double()
}
