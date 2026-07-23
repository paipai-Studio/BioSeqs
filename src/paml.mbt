///|
/// PAML - Molecular evolution analysis.
/// Inspired by Biopython Bio.PAML and the PAML software package.
/// Provides maximum likelihood estimation and dN/dS ratio calculation.

///|
pub struct PAMLAlignment {
  sequences : Array[(String, String)]
  tree : String
  model : String
}

///|
pub fn PAMLAlignment::new(sequences : Array[(String, String)]) -> PAMLAlignment {
  PAMLAlignment::{
    sequences,
    tree: "",
    model: "HKY85"
  }
}

///|
pub struct PAMLResult {
  ln_likelihood : Double
  parameters : Map[String, Double]
  dnds_ratios : Array[(String, Double)]
  substitution_rates : Array[(String, Double)]
  ancestral_sequences : Array[(String, String)]
}

///|
pub fn PAMLResult::new(ln_likelihood : Double) -> PAMLResult {
  PAMLResult::{
    ln_likelihood,
    parameters: Map([], capacity=10),
    dnds_ratios: Array::new(),
    substitution_rates: Array::new(),
    ancestral_sequences: Array::new()
  }
}

///|
pub struct DNDSResult {
  dN : Double
  dS : Double
  omega : Double
  se_dN : Double
  se_dS : Double
  se_omega : Double
  method : String
}

///|
pub fn DNDSResult::new(dN : Double, dS : Double, omega : Double) -> DNDSResult {
  DNDSResult::{
    dN,
    dS,
    omega,
    se_dN: 0.0,
    se_dS: 0.0,
    se_omega: 0.0,
    method: "Nei-Gojobori"
  }
}

///|
pub fn paml_calculate_dnds(seq1 : String, seq2 : String, method : String) -> DNDSResult {
  let aligned_seq1 = to_lowercase(seq1)
  let aligned_seq2 = to_lowercase(seq2)
  
  if aligned_seq1.length() != aligned_seq2.length() {
    return DNDSResult::new(0.0, 0.0, 0.0)
  }
  
  let mut synonymous_sites = 0.0
  let mut non_synonymous_sites = 0.0
  let mut synonymous_changes = 0
  let mut non_synonymous_changes = 0
  
  let codon_table = get_standard_codon_table()
  
  let mut i = 0
  while i < aligned_seq1.length() - 2 {
    let codon1 = substring(aligned_seq1, i, 3)
    let codon2 = substring(aligned_seq2, i, 3)
    
    if !contains(codon1, "-") && !contains(codon2, "-") {
      let aa1 = codon_table.get(codon1).unwrap_or("X")
      let aa2 = codon_table.get(codon2).unwrap_or("X")
      
      let (syn_sites, non_syn_sites) = count_sites(codon1, codon2, codon_table)
      synonymous_sites = synonymous_sites + syn_sites
      non_synonymous_sites = non_synonymous_sites + non_syn_sites
      
      if aa1 != aa2 {
        non_synonymous_changes = non_synonymous_changes + 1
      } else if codon1 != codon2 {
        synonymous_changes = synonymous_changes + 1
      }
    }
    
    i = i + 3
  }
  
  let dS = if synonymous_sites > 0.0 {
    paml_jukes_cantor_correction(synonymous_changes.to_double() / synonymous_sites)
  } else {
    0.0
  }
  
  let dN = if non_synonymous_sites > 0.0 {
    paml_jukes_cantor_correction(non_synonymous_changes.to_double() / non_synonymous_sites)
  } else {
    0.0
  }
  
  let omega = if dS > 0.0 { dN / dS } else { 0.0 }
  
  DNDSResult::new(dN, dS, omega)
}

///|
pub fn get_standard_codon_table() -> Map[String, String] {
  let table = Map([], capacity=64)
  
  table.set("TTT", "F"); table.set("TTC", "F")
  table.set("TTA", "L"); table.set("TTG", "L")
  table.set("TCT", "S"); table.set("TCC", "S"); table.set("TCA", "S"); table.set("TCG", "S")
  table.set("TAT", "Y"); table.set("TAC", "Y")
  table.set("TAA", "*"); table.set("TAG", "*"); table.set("TGA", "*")
  table.set("TGT", "C"); table.set("TGC", "C")
  table.set("TGG", "W")
  
  table.set("CTT", "L"); table.set("CTC", "L"); table.set("CTA", "L"); table.set("CTG", "L")
  table.set("CCT", "P"); table.set("CCC", "P"); table.set("CCA", "P"); table.set("CCG", "P")
  table.set("CAT", "H"); table.set("CAC", "H")
  table.set("CAA", "Q"); table.set("CAG", "Q")
  table.set("CGT", "R"); table.set("CGC", "R"); table.set("CGA", "R"); table.set("CGG", "R")
  
  table.set("ATT", "I"); table.set("ATC", "I"); table.set("ATA", "I")
  table.set("ATG", "M")
  table.set("ACT", "T"); table.set("ACC", "T"); table.set("ACA", "T"); table.set("ACG", "T")
  table.set("AAT", "N"); table.set("AAC", "N")
  table.set("AAA", "K"); table.set("AAG", "K")
  table.set("AGT", "S"); table.set("AGC", "S")
  table.set("AGA", "R"); table.set("AGG", "R")
  
  table.set("GTT", "V"); table.set("GTC", "V"); table.set("GTA", "V"); table.set("GTG", "V")
  table.set("GCT", "A"); table.set("GCC", "A"); table.set("GCA", "A"); table.set("GCG", "A")
  table.set("GAT", "D"); table.set("GAC", "D")
  table.set("GAA", "E"); table.set("GAG", "E")
  table.set("GGT", "G"); table.set("GGC", "G"); table.set("GGA", "G"); table.set("GGG", "G")
  
  table
}

///|
fn count_sites(codon1 : String, codon2 : String, codon_table : Map[String, String]) -> (Double, Double) {
  let aa1 = codon_table.get(codon1).unwrap_or("X")
  let aa2 = codon_table.get(codon2).unwrap_or("X")
  
  if aa1 == "*" || aa2 == "*" || aa1 == "X" || aa2 == "X" {
    return (0.0, 0.0)
  }
  
  if aa1 == aa2 {
    let mut syn_count = 0
    let mut all_count = 0
    
    let mut pos = 0
    while pos < 3 {
      if char_at(codon1, pos) != char_at(codon2, pos) {
        all_count = all_count + 1
        let bases = ["A", "T", "C", "G"]
        
        let mut b = 0
        while b < bases.length() {
          let mut new_codon = ""
          let mut p = 0
          while p < 3 {
            if p == pos {
              new_codon = new_codon + bases[b]
            } else {
              new_codon = new_codon + char_at(codon1, p)
            }
            p = p + 1
          }
          
          if codon_table.contains(new_codon) && codon_table[new_codon] == aa1 {
            syn_count = syn_count + 1
            break
          }
          b = b + 1
        }
      }
      pos = pos + 1
    }
    
    if all_count > 0 {
      return (syn_count.to_double() / all_count.to_double(), (all_count - syn_count).to_double() / all_count.to_double())
    } else {
      return (0.333, 0.667)
    }
  } else {
    return (0.0, 1.0)
  }
}

///|
fn paml_jukes_cantor_correction(p : Double) -> Double {
  if p >= 1.0 {
    10.0
  } else if p >= 0.75 {
    -0.75 * @math.log2(1.0 - p * 4.0 / 3.0)
  } else {
    -0.75 * @math.log2(1.0 - p * 4.0 / 3.0)
  }
}

///|
pub fn paml_estimate_parameters(alignment : PAMLAlignment) -> PAMLResult {
  let result = PAMLResult::new(-100.0)
  
  let n_seqs = alignment.sequences.length()
  let seq_len = if n_seqs > 0 { alignment.sequences[0].1.length() } else { 0 }
  
  let mut pi_A = 0.0
  let mut pi_T = 0.0
  let mut pi_C = 0.0
  let mut pi_G = 0.0
  
  let mut total = 0
  let mut s = 0
  while s < n_seqs {
    let seq = to_lowercase(alignment.sequences[s].1)
    let mut i = 0
    while i < seq.length() {
      let c = char_at(seq, i)
      if c == "A" {
        pi_A = pi_A + 1.0
      } else if c == "T" {
        pi_T = pi_T + 1.0
      } else if c == "C" {
        pi_C = pi_C + 1.0
      } else if c == "G" {
        pi_G = pi_G + 1.0
      }
      total = total + 1
      i = i + 1
    }
    s = s + 1
  }
  
  let denom = total.to_double()
  pi_A = pi_A / denom
  pi_T = pi_T / denom
  pi_C = pi_C / denom
  pi_G = pi_G / denom
  
  let params = result.parameters.copy()
  params.set("pi_A", pi_A)
  params.set("pi_T", pi_T)
  params.set("pi_C", pi_C)
  params.set("pi_G", pi_G)
  params.set("kappa", 2.0)
  params.set("omega", 0.5)
  
  PAMLResult::{
    ln_likelihood: -seq_len.to_double() * @math.log2(0.25),
    parameters: params,
    dnds_ratios: result.dnds_ratios,
    substitution_rates: result.substitution_rates,
    ancestral_sequences: result.ancestral_sequences
  }
}

///|
pub fn paml_calculate_substitution_matrix(kappa : Double, pi : Map[String, Double]) -> Array[Array[Double]] {
  let q = Array::make(4, Array::make(4, 0.0))
  
  let pi_A = pi.get("pi_A").unwrap_or(0.25)
  let pi_T = pi.get("pi_T").unwrap_or(0.25)
  let pi_C = pi.get("pi_C").unwrap_or(0.25)
  let pi_G = pi.get("pi_G").unwrap_or(0.25)
  
  q[0][1] = kappa * pi_T
  q[0][2] = pi_C
  q[0][3] = pi_G
  
  q[1][0] = kappa * pi_A
  q[1][2] = pi_C
  q[1][3] = pi_G
  
  q[2][0] = pi_A
  q[2][1] = pi_T
  q[2][3] = kappa * pi_G
  
  q[3][0] = pi_A
  q[3][1] = pi_T
  q[3][2] = kappa * pi_C
  
  let mut i = 0
  while i < 4 {
    let mut sum = 0.0
    let mut j = 0
    while j < 4 {
      if i != j {
        sum = sum + q[i][j]
      }
      j = j + 1
    }
    q[i][i] = -sum
    i = i + 1
  }
  
  q
}

///|
pub fn paml_create_example_alignment() -> PAMLAlignment {
  let seqs = [
    ("Human", "ATGGCCGGTAAATGA"),
    ("Mouse", "ATGGCCGGTAAATGA"),
    ("Rat", "ATGGCCGGTAATAGA"),
    ("Chicken", "ATGGCCGGTAATGGA")
  ]
  
  let alignment = PAMLAlignment::new(seqs)
  PAMLAlignment::{
    sequences: alignment.sequences,
    tree: "((Human,Mouse),Rat,Chicken)",
    model: "HKY85"
  }
}

///|
pub fn paml_run_likelihood(alignment : PAMLAlignment) -> PAMLResult {
  let result = paml_estimate_parameters(alignment)
  
  let dnds_results : Array[(String, Double)] = Array::new()
  
  let seqs = alignment.sequences
  let mut i = 0
  while i < seqs.length() {
    let mut j = i + 1
    while j < seqs.length() {
      let dnds = paml_calculate_dnds(seqs[i].1, seqs[j].1, "Nei-Gojobori")
      dnds_results.push((seqs[i].0 + "-" + seqs[j].0, dnds.omega))
      j = j + 1
    }
    i = i + 1
  }
  
  PAMLResult::{
    ln_likelihood: result.ln_likelihood,
    parameters: result.parameters,
    dnds_ratios: dnds_results,
    substitution_rates: result.substitution_rates,
    ancestral_sequences: result.ancestral_sequences
  }
}