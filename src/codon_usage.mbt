///|
/// Codon Adaptation Index (CAI) and codon usage analysis.
///
/// Biotite provides codon usage analysis. This module implements CAI
/// calculation and relative synonymous codon usage (RSCU).

///|
/// Standard genetic code: codon -> amino acid.
fn cu_codon_to_aa(codon : String) -> Char {
  let aa = codon_lookup(codon)
  let s = if aa == "?" { "X" } else { aa }
  s.unsafe_get(0).unsafe_to_char()
}

///|
/// Calculate codon frequencies from a coding sequence.
///
/// Returns a map from codon to count.
pub fn codon_counts(coding_sequence : String) -> Map[String, Int32] {
  let counts : Map[String, Int32] = Map([], capacity=64)
  let n = coding_sequence.length()
  let mut i = 0
  while i + 3 <= n {
    let codon = sc_extract(coding_sequence, i, 3)
    counts[codon] = counts.get(codon).unwrap_or(0) + 1
    i = i + 3
  }
  counts
}

///|
/// Extract substring helper.
fn sc_extract(s : String, start : Int32, length : Int32) -> String {
  let mut result = ""
  for i in 0..<length {
    if start + i < s.length() {
      result = result + s.unsafe_get(start + i).unsafe_to_char().to_string()
    }
  }
  result
}

///|
/// Calculate Relative Synonymous Codon Usage (RSCU).
///
/// RSCU = (count of codon) / (expected count if all synonymous codons
/// were used equally). RSCU > 1 means the codon is used more than expected.
pub fn rscu_table(coding_sequence : String) -> Map[String, Double] {
  let counts = codon_counts(coding_sequence)
  let rscu : Map[String, Double] = Map([], capacity=64)
  // Group codons by amino acid
  let aa_codons : Map[Int32, Array[String]] = Map([], capacity=21)
  let all_codons = [
    "TTT", "TTC", "TTA", "TTG", "CTT", "CTC", "CTA", "CTG", "ATT", "ATC", "ATA",
    "ATG", "GTT", "GTC", "GTA", "GTG", "TCT", "TCC", "TCA", "TCG", "AGT", "AGC",
    "CCT", "CCC", "CCA", "CCG", "ACT", "ACC", "ACA", "ACG", "GCT", "GCC", "GCA",
    "GCG", "TAT", "TAC", "TAA", "TAG", "CAT", "CAC", "CAA", "CAG", "AAT", "AAC",
    "AAA", "AAG", "GAT", "GAC", "GAA", "GAG", "TGT", "TGC", "TGA", "TGG", "CGT",
    "CGC", "CGA", "CGG", "AGA", "AGG", "GGT", "GGC", "GGA", "GGG",
  ]
  for codon in all_codons {
    let aa = cu_codon_to_aa(codon).to_int()
    let existing = aa_codons.get(aa).unwrap_or([])
    let new_arr : Array[String] = []
    for c in existing {
      new_arr.push(c)
    }
    new_arr.push(codon)
    aa_codons[aa] = new_arr
  }
  // Calculate RSCU for each codon
  let aa_keys = aa_codons.keys().collect()
  for aa in aa_keys {
    let codons = aa_codons[aa]
    let mut total = 0
    for c in codons {
      total = total + counts.get(c).unwrap_or(0)
    }
    for c in codons {
      let count = counts.get(c).unwrap_or(0)
      if total > 0 {
        rscu[c] = count.to_double() *
          codons.length().to_double() /
          total.to_double()
      } else {
        rscu[c] = 0.0
      }
    }
  }
  rscu
}

///|
/// Calculate Codon Adaptation Index (CAI).
///
/// CAI = geometric mean of relative adaptiveness (w) values.
/// w_i = RSCU_i / RSCU_max for the amino acid encoded by codon i.
///
/// Parameters:
/// - coding_sequence: the coding DNA sequence
/// - rscu_table: pre-computed RSCU table (from reference set)
pub fn cai_score(
  coding_sequence : String,
  rscu_table : Map[String, Double],
) -> Double {
  let n = coding_sequence.length()
  let mut log_sum = 0.0
  let mut codon_count = 0
  let mut i = 0
  while i + 3 <= n {
    let codon = sc_extract(coding_sequence, i, 3)
    let aa = cu_codon_to_aa(codon).to_int()
    if aa != '*' && aa != 'X' {
      let rscu = rscu_table.get(codon).unwrap_or(0.0)
      if rscu > 0.0 {
        log_sum = log_sum + @math.ln(rscu)
        codon_count = codon_count + 1
      }
    }
    i = i + 3
  }
  if codon_count == 0 {
    return 0.0
  }
  // CAI = exp(mean(ln(w)))
  @math.exp(log_sum / codon_count.to_double())
}

///|
/// Calculate GC content of a coding sequence at codon positions.
///
/// Returns (GC1, GC2, GC3) — GC content at each codon position.
pub fn gc_by_position(coding_sequence : String) -> (Double, Double, Double) {
  let n = coding_sequence.length()
  let mut gc1 = 0
  let mut gc2 = 0
  let mut gc3 = 0
  let mut total_codons = 0
  let mut i = 0
  while i + 3 <= n {
    for pos in 0..<3 {
      let c = coding_sequence.unsafe_get(i + pos)
      // 'G' = 71, 'C' = 67
      if c == 71 || c == 67 {
        match pos {
          0 => gc1 = gc1 + 1
          1 => gc2 = gc2 + 1
          _ => gc3 = gc3 + 1
        }
      }
    }
    total_codons = total_codons + 1
    i = i + 3
  }
  if total_codons == 0 {
    return (0.0, 0.0, 0.0)
  }
  let tot = total_codons.to_double()
  (gc1.to_double() / tot, gc2.to_double() / tot, gc3.to_double() / tot)
}

///|
/// Result struct for complete codon usage analysis.
pub struct CodonUsageResult {
  codon_counts : Map[String, Int32]
  rscu : Map[String, Double]
  cai : Double
  enc : Double
  gc3 : Double
  optimal_codons : Array[String]
}

///|
/// Count codons in a coding sequence. Returns Map[String, Int].
pub fn count_codons(seq : String) -> Map[String, Int32] {
  let counts : Map[String, Int32] = Map([], capacity=64)
  let n = seq.length()
  let mut i = 0
  while i + 3 <= n {
    let codon = cu_extract3(seq, i)
    counts[codon] = counts.get(codon).unwrap_or(0) + 1
    i = i + 3
  }
  counts
}

///|
/// Extract 3 characters starting at position i.
fn cu_extract3(s : String, start : Int32) -> String {
  let mut result = ""
  for j in 0..<3 {
    if start + j < s.length() {
      result = result + s.unsafe_get(start + j).unsafe_to_char().to_string()
    }
  }
  result
}

///|
/// Calculate RSCU from codon counts.
pub fn codon_rscu(counts : Map[String, Int32]) -> Map[String, Double] {
  let rscu : Map[String, Double] = Map([], capacity=64)
  let aa_codons : Map[Int32, Array[String]] = Map([], capacity=21)
  let all_codons = [
    "TTT", "TTC", "TTA", "TTG", "CTT", "CTC", "CTA", "CTG", "ATT", "ATC", "ATA",
    "ATG", "GTT", "GTC", "GTA", "GTG", "TCT", "TCC", "TCA", "TCG", "AGT", "AGC",
    "CCT", "CCC", "CCA", "CCG", "ACT", "ACC", "ACA", "ACG", "GCT", "GCC", "GCA",
    "GCG", "TAT", "TAC", "TAA", "TAG", "CAT", "CAC", "CAA", "CAG", "AAT", "AAC",
    "AAA", "AAG", "GAT", "GAC", "GAA", "GAG", "TGT", "TGC", "TGA", "TGG", "CGT",
    "CGC", "CGA", "CGG", "AGA", "AGG", "GGT", "GGC", "GGA", "GGG",
  ]
  for codon in all_codons {
    let aa = cu_codon_to_aa(codon).to_int()
    let existing = aa_codons.get(aa).unwrap_or([])
    let new_arr : Array[String] = []
    for c in existing {
      new_arr.push(c)
    }
    new_arr.push(codon)
    aa_codons[aa] = new_arr
  }
  let aa_keys = aa_codons.keys().collect()
  for aa in aa_keys {
    let codons = aa_codons[aa]
    let mut total = 0
    for c in codons {
      total = total + counts.get(c).unwrap_or(0)
    }
    for c in codons {
      let count = counts.get(c).unwrap_or(0)
      if total > 0 {
        rscu[c] = count.to_double() *
          codons.length().to_double() /
          total.to_double()
      } else {
        rscu[c] = 0.0
      }
    }
  }
  rscu
}

///|
/// Calculate CAI from a sequence and optimal codon weights.
pub fn codon_cai(seq : String, optimal : Map[String, Double]) -> Double {
  let n = seq.length()
  let mut log_sum = 0.0
  let mut codon_count = 0
  let mut i = 0
  while i + 3 <= n {
    let codon = cu_extract3(seq, i)
    let aa = cu_codon_to_aa(codon)
    if aa != '*' && aa != 'X' {
      let w = optimal.get(codon).unwrap_or(0.0)
      if w > 0.0 {
        log_sum = log_sum + @math.ln(w)
        codon_count = codon_count + 1
      }
    }
    i = i + 3
  }
  if codon_count == 0 {
    return 0.0
  }
  @math.exp(log_sum / codon_count.to_double())
}

///|
/// Effective Number of Codons (ENC).
pub fn codon_enc(counts : Map[String, Int32]) -> Double {
  // ENC = sum over amino acids of F_aa
  // where F_aa = (n * sum(p_i^2))^(-1) for each synonymous family
  // Simplified: use Wright's method
  let aa_codons : Map[Int32, Array[String]] = Map([], capacity=21)
  let all_codons = [
    "TTT", "TTC", "TTA", "TTG", "CTT", "CTC", "CTA", "CTG", "ATT", "ATC", "ATA",
    "ATG", "GTT", "GTC", "GTA", "GTG", "TCT", "TCC", "TCA", "TCG", "AGT", "AGC",
    "CCT", "CCC", "CCA", "CCG", "ACT", "ACC", "ACA", "ACG", "GCT", "GCC", "GCA",
    "GCG", "TAT", "TAC", "TAA", "TAG", "CAT", "CAC", "CAA", "CAG", "AAT", "AAC",
    "AAA", "AAG", "GAT", "GAC", "GAA", "GAG", "TGT", "TGC", "TGA", "TGG", "CGT",
    "CGC", "CGA", "CGG", "AGA", "AGG", "GGT", "GGC", "GGA", "GGG",
  ]
  for codon in all_codons {
    let aa = cu_codon_to_aa(codon).to_int()
    let existing = aa_codons.get(aa).unwrap_or([])
    let new_arr : Array[String] = []
    for c in existing {
      new_arr.push(c)
    }
    new_arr.push(codon)
    aa_codons[aa] = new_arr
  }
  let mut enc = 0.0
  let aa_keys = aa_codons.keys().collect()
  for aa in aa_keys {
    // Skip stop codons
    if aa == '*' {
      continue
    }
    let codons = aa_codons[aa]
    let mut total = 0
    for c in codons {
      total = total + counts.get(c).unwrap_or(0)
    }
    if total > 0 {
      // Amino acid present: F = 1 / sum(p_i^2)
      let mut sum_sq = 0.0
      for c in codons {
        let count = counts.get(c).unwrap_or(0)
        let p = count.to_double() / total.to_double()
        sum_sq = sum_sq + p * p
      }
      if sum_sq > 0.0 {
        enc = enc + 1.0 / sum_sq
      } else {
        enc = enc + codons.length().to_double()
      }
    } else {
      // Amino acid absent: contribute degeneracy (all equally used)
      enc = enc + codons.length().to_double()
    }
  }
  enc
}

///|
/// GC content at third codon position.
pub fn codon_gc3(seq : String) -> Double {
  let n = seq.length()
  let mut gc3 = 0
  let mut total = 0
  let mut i = 0
  while i + 3 <= n {
    let c = seq.unsafe_get(i + 2)
    if c == 71 || c == 67 {
      gc3 = gc3 + 1
    }
    total = total + 1
    i = i + 3
  }
  if total == 0 {
    return 0.0
  }
  gc3.to_double() / total.to_double()
}

///|
/// Translate a single codon to amino acid letter.
pub fn codon_translate(codon : String) -> String {
  cu_codon_to_aa(codon).to_string()
}

///|
/// Reverse complement of a codon.
pub fn codon_reverse_complement(codon : String) -> String {
  let mut result = ""
  for i in 0..<codon.length() {
    let ch = codon.unsafe_get(codon.length() - 1 - i).unsafe_to_char()
    let comp = match ch {
      'A' => 'T'
      'T' => 'A'
      'G' => 'C'
      'C' => 'G'
      _ => ch
    }
    result = result + comp.to_string()
  }
  result
}

///|
/// Find optimal codons (highest RSCU for each amino acid).
pub fn codon_optimal(rscu : Map[String, Double]) -> Array[String] {
  let aa_codons : Map[Int32, Array[String]] = Map([], capacity=21)
  let all_codons = [
    "TTT", "TTC", "TTA", "TTG", "CTT", "CTC", "CTA", "CTG", "ATT", "ATC", "ATA",
    "ATG", "GTT", "GTC", "GTA", "GTG", "TCT", "TCC", "TCA", "TCG", "AGT", "AGC",
    "CCT", "CCC", "CCA", "CCG", "ACT", "ACC", "ACA", "ACG", "GCT", "GCC", "GCA",
    "GCG", "TAT", "TAC", "TAA", "TAG", "CAT", "CAC", "CAA", "CAG", "AAT", "AAC",
    "AAA", "AAG", "GAT", "GAC", "GAA", "GAG", "TGT", "TGC", "TGA", "TGG", "CGT",
    "CGC", "CGA", "CGG", "AGA", "AGG", "GGT", "GGC", "GGA", "GGG",
  ]
  for codon in all_codons {
    let aa = cu_codon_to_aa(codon).to_int()
    let existing = aa_codons.get(aa).unwrap_or([])
    let new_arr : Array[String] = []
    for c in existing {
      new_arr.push(c)
    }
    new_arr.push(codon)
    aa_codons[aa] = new_arr
  }
  let optimal : Array[String] = []
  let aa_keys = aa_codons.keys().collect()
  for aa in aa_keys {
    let codons = aa_codons[aa]
    let mut best_codon = ""
    let mut best_rscu = 0.0
    for c in codons {
      let r = rscu.get(c).unwrap_or(0.0)
      if r > best_rscu {
        best_rscu = r
        best_codon = c
      }
    }
    if best_codon.length() > 0 {
      optimal.push(best_codon)
    }
  }
  optimal
}

///|
/// Codon Bias Index (CBI).
pub fn codon_cbi(seq : String, optimal : Array[String]) -> Double {
  let n = seq.length()
  let mut num = 0.0
  let mut den = 0.0
  let mut i = 0
  while i + 3 <= n {
    let codon = cu_extract3(seq, i)
    let aa = cu_codon_to_aa(codon)
    if aa != '*' && aa != 'X' {
      let is_optimal = optimal.contains(codon)
      // Count optimal and synonymous
      let mut opt_count = 0
      let mut total_count = 0
      let mut j = 0
      while j + 3 <= n {
        let c2 = cu_extract3(seq, j)
        if cu_codon_to_aa(c2) == aa {
          total_count = total_count + 1
          if optimal.contains(c2) {
            opt_count = opt_count + 1
          }
        }
        j = j + 3
      }
      if is_optimal {
        num = num + 1.0
      }
      den = den + 1.0
    }
    i = i + 3
  }
  if den == 0.0 {
    return 0.0
  }
  (num - 1.0) / (den - 1.0)
}

///|
/// Frequency of Optimal Codons (Fop).
pub fn codon_fop(seq : String, optimal : Array[String]) -> Double {
  let n = seq.length()
  let mut opt_count = 0
  let mut total = 0
  let mut i = 0
  while i + 3 <= n {
    let codon = cu_extract3(seq, i)
    let aa = cu_codon_to_aa(codon)
    if aa != '*' && aa != 'X' {
      total = total + 1
      if optimal.contains(codon) {
        opt_count = opt_count + 1
      }
    }
    i = i + 3
  }
  if total == 0 {
    return 0.0
  }
  opt_count.to_double() / total.to_double()
}

///|
/// Complete codon usage analysis.
pub fn codon_usage_analysis(seq : String) -> CodonUsageResult {
  let counts = count_codons(seq)
  let rscu = codon_rscu(counts)
  let optimal = codon_optimal(rscu)
  // Build optimal weights map for CAI
  let opt_weights : Map[String, Double] = Map([], capacity=64)
  for c in optimal {
    opt_weights[c] = rscu.get(c).unwrap_or(1.0)
  }
  let cai = codon_cai(seq, opt_weights)
  let enc = codon_enc(counts)
  let gc3 = codon_gc3(seq)
  { codon_counts: counts, rscu, cai, enc, gc3, optimal_codons: optimal, }
}

///|
/// Codon usage analysis with reference sequence.
pub fn codon_usage_with_ref(seq : String, ref_seq : String) -> CodonUsageResult {
  let ref_counts = count_codons(ref_seq)
  let ref_rscu = codon_rscu(ref_counts)
  let optimal = codon_optimal(ref_rscu)
  let opt_weights : Map[String, Double] = Map([], capacity=64)
  for c in optimal {
    opt_weights[c] = ref_rscu.get(c).unwrap_or(1.0)
  }
  let cai = codon_cai(seq, opt_weights)
  let counts = count_codons(seq)
  let enc = codon_enc(counts)
  let gc3 = codon_gc3(seq)
  {
    codon_counts: counts,
    rscu: ref_rscu,
    cai,
    enc,
    gc3,
    optimal_codons: optimal,
  }
}
