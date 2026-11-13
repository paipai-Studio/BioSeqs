/// biopython_paml_codon.mbt — PAML codon model utilities
/// Inspired by Biopython Bio.PAML: codon substitution models, dN/dS

///|
/// Codon frequency model types.
pub enum BiopythonCodonFreqModel {
  FEqual // equal frequencies (1/61)
  F1x4 //  nucleotide frequencies product
  F3x4 //  nucleotide frequencies by position
  F61 //  codon frequencies
}

///|
/// Count synonymous and nonsynonymous sites in a codon.
/// Returns (S_sites, N_sites) where S + N = 3.
fn biopython_paml_count_sites(codon : String) -> (Double, Double) {
  let aa = biopython_paml_codon_to_aa(codon)
  let mut s_sites = 0.0
  let mut n_sites = 0.0
  let bases = "ACGT"
  for pos = 0; pos < 3; pos = pos + 1 {
    for bi = 0; bi < 4; bi = bi + 1 {
      let new_base = bases[bi:bi + 1].to_owned()
      let orig_base = codon[pos:pos + 1].to_owned()
      if new_base == orig_base {
        continue
      }
      // Create mutant codon
      let mutant = if pos == 0 {
        new_base + codon[1:3].to_owned()
      } else if pos == 1 {
        codon[0:1].to_owned() + new_base + codon[2:3].to_owned()
      } else {
        codon[0:2].to_owned() + new_base
      }
      let mut_aa = biopython_paml_codon_to_aa(mutant)
      if mut_aa == aa {
        s_sites = s_sites + 1.0 / 3.0
      } else {
        n_sites = n_sites + 1.0 / 3.0
      }
    }
  }
  (s_sites, n_sites)
}

///|
/// Standard genetic code.
fn biopython_paml_codon_to_aa(codon : String) -> String {
  let aa = codon_lookup(codon)
  if aa == "?" {
    "X"
  } else {
    aa
  }
}

///|
/// Count synonymous (Sd) and nonsynonymous (Nd) differences between two codons.
fn biopython_paml_count_diffs(
  codon1 : String,
  codon2 : String,
) -> (Int32, Int32) {
  let mut n_diff = 0
  let positions : Array[Int32] = Array(capacity=0)
  for i = 0; i < 3; i = i + 1 {
    if codon1[i] != codon2[i] {
      n_diff = n_diff + 1
      positions.push(i)
    }
  }
  if n_diff == 0 {
    return (0, 0)
  }
  let aa1 = biopython_paml_codon_to_aa(codon1)
  let aa2 = biopython_paml_codon_to_aa(codon2)
  if n_diff == 1 {
    if aa1 == aa2 {
      (1, 0) // synonymous
    } else {
      (0, 1) // nonsynonymous
    }
  } else if n_diff == 2 {
    // Need to check both pathways
    let p1 = biopython_paml_mutate_codon(codon1, positions[0], codon2)
    let p1_aa = biopython_paml_codon_to_aa(p1)
    let p2 = biopython_paml_mutate_codon(codon1, positions[1], codon2)
    let p2_aa = biopython_paml_codon_to_aa(p2)
    // Pathway 1: codon1 -> p1 -> codon2
    let (s1a, n1a) = if p1_aa == aa1 { (1, 0) } else { (0, 1) }
    let (s1b, n1b) = if aa2 == p1_aa { (1, 0) } else { (0, 1) }
    // Pathway 2: codon1 -> p2 -> codon2
    let (s2a, n2a) = if p2_aa == aa1 { (1, 0) } else { (0, 1) }
    let (s2b, n2b) = if aa2 == p2_aa { (1, 0) } else { (0, 1) }
    ((s1a + s1b + s2a + s2b) / 2, (n1a + n1b + n2a + n2b) / 2)
    // 3 differences: approximate
  } else if aa1 == aa2 {
    (3, 0)
  } else {
    (0, 3)
  }
}

///|
/// Mutate a codon at a specific position.
fn biopython_paml_mutate_codon(
  codon : String,
  pos : Int32,
  target : String,
) -> String {
  if pos == 0 {
    target[0:1].to_owned() + codon[1:3].to_owned()
  } else if pos == 1 {
    codon[0:1].to_owned() + target[1:2].to_owned() + codon[2:3].to_owned()
  } else {
    codon[0:2].to_owned() + target[2:3].to_owned()
  }
}

///|
/// Nei-Gojobori method for dN/dS estimation.
/// Returns (dN, dS, dN_dS_ratio).
pub fn biopython_nei_gojobori(
  seq1 : String,
  seq2 : String,
) -> (Double, Double, Double) {
  let n = seq1.length()
  if n < 3 || n != seq2.length() {
    return (0.0, 0.0, 0.0)
  }
  let mut total_sd = 0.0 // synonymous differences
  let mut total_nd = 0.0 // nonsynonymous differences
  let mut total_s = 0.0 // synonymous sites
  let mut total_n = 0.0 // nonsynonymous sites
  let mut n_codons = 0
  let mut i = 0
  while i + 3 <= n {
    let codon1 = seq1[i:i + 3].to_owned()
    let codon2 = seq2[i:i + 3].to_owned()
    let (s, ns) = biopython_paml_count_sites(codon1)
    total_s = total_s + s
    total_n = total_n + ns
    let (sd, nd) = biopython_paml_count_diffs(codon1, codon2)
    total_sd = total_sd + sd.to_double()
    total_nd = total_nd + nd.to_double()
    n_codons = n_codons + 1
    i = i + 3
  }
  if n_codons == 0 || total_s == 0.0 || total_n == 0.0 {
    return (0.0, 0.0, 0.0)
  }
  let p_s = total_sd / total_s
  let p_n = total_nd / total_n
  // Jukes-Cantor correction
  let d_s = if p_s < 0.75 {
    -0.75 * @math.ln(1.0 - 4.0 * p_s / 3.0)
  } else {
    0.0
  }
  let d_n = if p_n < 0.75 {
    -0.75 * @math.ln(1.0 - 4.0 * p_n / 3.0)
  } else {
    0.0
  }
  let ratio = if d_s > 0.0 { d_n / d_s } else { 0.0 }
  (d_n, d_s, ratio)
}

///|
/// Dummy function to suppress unused constructor warnings.
pub fn biopython_paml_dummy() -> Unit {
  ignore(BiopythonCodonFreqModel::FEqual)
  ignore(BiopythonCodonFreqModel::F1x4)
  ignore(BiopythonCodonFreqModel::F3x4)
  ignore(BiopythonCodonFreqModel::F61)
}
