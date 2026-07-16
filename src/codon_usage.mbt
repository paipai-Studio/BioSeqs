///|
/// Codon usage analysis module, inspired by Biopython Bio.codonalign and Bio.CodonUsage.
///
/// Provides functions for analyzing codon usage patterns in DNA sequences,
/// including Codon Adaptation Index (CAI), Effective Number of Codons (ENC),
/// and optimal codon detection.

///|
fn codon_log_nat(x : Double) -> Double {
  if x <= 0.0 {
    0.0
  } else if x == 1.0 {
    0.0
  } else {
    let y = (x - 1.0) / (x + 1.0)
    let y2 = y * y
    let mut sum = y
    let mut term = y
    let mut i = 1
    while i <= 20 {
      term = term * y2
      sum = sum + term / (2 * i + 1).to_double()
      i = i + 1
    }
    2.0 * sum
  }
}

///|
fn exp_val(x : Double) -> Double {
  let mut sum = 1.0
  let mut term = 1.0
  let mut i = 1
  while i <= 20 {
    term = term * x / i.to_double()
    sum = sum + term
    i = i + 1
  }
  sum
}

///|
/// Get standard genetic code codon table.
/// Returns a Map mapping amino acid to list of codons.
pub fn get_standard_code() -> Map[String, Array[String]] {
  let code = Map([], capacity=21)
  code.set("A", ["GCA", "GCC", "GCG", "GCT"])
  code.set("R", ["AGA", "AGG", "CGA", "CGC", "CGG", "CGT"])
  code.set("N", ["AAC", "AAT"])
  code.set("D", ["GAC", "GAT"])
  code.set("C", ["TGC", "TGT"])
  code.set("E", ["GAA", "GAG"])
  code.set("Q", ["CAA", "CAG"])
  code.set("G", ["GGA", "GGC", "GGG", "GGT"])
  code.set("H", ["CAC", "CAT"])
  code.set("I", ["ATA", "ATC", "ATT"])
  code.set("L", ["CTA", "CTC", "CTG", "CTT", "TTA", "TTG"])
  code.set("K", ["AAA", "AAG"])
  code.set("M", ["ATG"])
  code.set("F", ["TTC", "TTT"])
  code.set("P", ["CCA", "CCC", "CCG", "CCT"])
  code.set("S", ["AGC", "AGT", "TCA", "TCC", "TCG", "TCT"])
  code.set("T", ["ACA", "ACC", "ACG", "ACT"])
  code.set("W", ["TGG"])
  code.set("Y", ["TAC", "TAT"])
  code.set("V", ["GTA", "GTC", "GTG", "GTT"])
  code.set("*", ["TAA", "TAG", "TGA"])
  code
}

///|
/// Count codon frequencies in a DNA sequence.
/// Returns a Map mapping codon to count.
pub fn count_codons(seq : String) -> Map[String, Int] {
  let counts = Map([], capacity=64)
  let seq_upper = seq.to_upper()
  let len = seq_upper.length()

  let mut i = 0
  while i + 3 <= len {
    let codon = seq_upper[i:i + 3].to_owned()
    counts.set(codon, counts.get(codon).unwrap_or(0) + 1)
    i = i + 3
  }

  counts
}

///|
/// Calculate Relative Synonymous Codon Usage (RSCU).
/// RSCU = observed codon frequency / expected frequency (1/number of synonyms)
pub fn codon_rscu(counts : Map[String, Int]) -> Map[String, Double] {
  let rscu_map = Map([], capacity=64)
  let code = get_standard_code()
  let aa_list = code.keys().collect()

  let mut i = 0
  while i < aa_list.length() {
    let aa = aa_list[i]
    let codons = code[aa]
    let n_syn = codons.length().to_double()

    let mut total = 0
    let mut j = 0
    while j < codons.length() {
      total = total + counts.get(codons[j]).unwrap_or(0)
      j = j + 1
    }

    if total > 0 {
      let mut k = 0
      while k < codons.length() {
        let codon = codons[k]
        let count = counts.get(codon).unwrap_or(0).to_double()
        let expected = total.to_double() / n_syn
        let rscu_val = if expected > 0.0 { count / expected } else { 0.0 }
        rscu_map.set(codon, rscu_val)
        k = k + 1
      }
    }
    i = i + 1
  }

  rscu_map
}

///|
/// Calculate Codon Adaptation Index (CAI).
/// CAI measures the adaptation of a gene's codon usage to the optimal codon usage
/// of a reference set. Range: 0-1, higher = better adaptation.
pub fn codon_cai(seq : String, optimal_codons : Map[String, Double]) -> Double {
  let seq_upper = seq.to_upper()
  let len = seq_upper.length()

  let mut total_log = 0.0
  let mut count = 0

  let mut i = 0
  while i + 3 <= len {
    let codon = seq_upper[i:i + 3].to_owned()
    if optimal_codons.contains(codon) {
      let val = optimal_codons.get(codon).unwrap_or(0.0)
      if val > 0.0 {
        total_log = total_log + codon_log_nat(val)
        count = count + 1
      }
    }
    i = i + 3
  }

  if count == 0 {
    return 0.0
  }

  exp_val(total_log / count.to_double())
}

///|
/// Calculate Effective Number of Codons (ENC).
/// ENC measures codon usage bias. Range: 20 (extreme bias) to 61 (no bias).
pub fn codon_enc(counts : Map[String, Int]) -> Double {
  let mut sum_f2 = 0.0
  let code = get_standard_code()
  let aa_list = code.keys().collect()

  let mut i = 0
  while i < aa_list.length() {
    let aa = aa_list[i]
    let codons = code[aa]

    let mut total = 0
    let mut j = 0
    while j < codons.length() {
      total = total + counts.get(codons[j]).unwrap_or(0)
      j = j + 1
    }

    if total > 0 {
      let mut f2_sum = 0.0
      let mut k = 0
      while k < codons.length() {
        let f = counts.get(codons[k]).unwrap_or(0).to_double() /
          total.to_double()
        f2_sum = f2_sum + f * f
        k = k + 1
      }

      let n_syn = codons.length().to_double()
      let term = if f2_sum > 0.0 { 1.0 / f2_sum } else { n_syn }
      sum_f2 = sum_f2 + term
    }
    i = i + 1
  }

  20.0 + 41.0 / sum_f2
}

///|
/// Detect optimal codons based on RSCU values.
/// Optimal codons are those with RSCU > 1.0 for each amino acid.
pub fn codon_optimal(rscu_map : Map[String, Double]) -> Map[String, Double] {
  let optimal = Map([], capacity=21)
  let code = get_standard_code()
  let aa_list = code.keys().collect()

  let mut i = 0
  while i < aa_list.length() {
    let aa = aa_list[i]
    let codons = code[aa]

    let mut max_rscu = 0.0
    let mut best_codon = ""

    let mut j = 0
    while j < codons.length() {
      let codon = codons[j]
      let r = rscu_map.get(codon).unwrap_or(0.0)
      if r > max_rscu {
        max_rscu = r
        best_codon = codon
      }
      j = j + 1
    }

    if max_rscu > 0.0 {
      optimal.set(best_codon, max_rscu)
    }
    i = i + 1
  }

  optimal
}

///|
/// Calculate GC content at third codon position (GC3).
pub fn codon_gc3(seq : String) -> Double {
  let seq_upper = seq.to_upper()
  let len = seq_upper.length()

  let mut gc_count = 0
  let mut total = 0

  let mut i = 0
  while i + 3 <= len {
    let third = seq_upper[i + 2]
    if third == 'G' || third == 'C' {
      gc_count = gc_count + 1
    }
    total = total + 1
    i = i + 3
  }

  if total == 0 {
    return 0.0
  }

  gc_count.to_double() / total.to_double()
}

///|
/// Calculate codon bias index (CBI).
/// CBI measures the extent to which the codon usage of a gene departs from
/// the average codon usage.
pub fn codon_cbi(seq : String, optimal_codons : Array[String]) -> Double {
  let seq_upper = seq.to_upper()
  let len = seq_upper.length()

  let mut f_opt = 0
  let mut f_nonopt = 0

  let mut i = 0
  while i + 3 <= len {
    let codon = seq_upper[i:i + 3].to_owned()
    let is_optimal = optimal_codons.contains(codon)

    let aa = codon_translate(codon)
    if aa != "*" {
      if is_optimal {
        f_opt = f_opt + 1
      } else {
        f_nonopt = f_nonopt + 1
      }
    }
    i = i + 3
  }

  if f_opt + f_nonopt == 0 {
    return 0.0
  }

  (f_opt.to_double() - f_nonopt.to_double()) /
  (f_opt.to_double() + f_nonopt.to_double())
}

///|
/// Translate a single codon to amino acid.
pub fn codon_translate(codon : String) -> String {
  let codon_upper = codon.to_upper()
  let code = get_standard_code()
  let aa_list = code.keys().collect()

  let mut i = 0
  while i < aa_list.length() {
    let aa = aa_list[i]
    let codons = code[aa]
    if codons.contains(codon_upper) {
      return aa
    }
    i = i + 1
  }

  "X"
}

///|
/// Calculate frequency of optimal codons (Fop).
pub fn codon_fop(seq : String, optimal_codons : Array[String]) -> Double {
  let seq_upper = seq.to_upper()
  let len = seq_upper.length()

  let mut opt_count = 0
  let mut total = 0

  let mut i = 0
  while i + 3 <= len {
    let codon = seq_upper[i:i + 3].to_owned()
    let aa = codon_translate(codon)

    if aa != "*" && aa != "X" {
      if optimal_codons.contains(codon) {
        opt_count = opt_count + 1
      }
      total = total + 1
    }
    i = i + 3
  }

  if total == 0 {
    return 0.0
  }

  opt_count.to_double() / total.to_double()
}

///|
/// Calculate reverse complement of a codon.
pub fn codon_reverse_complement(codon : String) -> String {
  if codon.length() != 3 {
    return codon
  }

  let mut result = ""
  let mut i = 0
  while i < 3 {
    let c = codon[2 - i:3 - i]
    if c == "A" {
      result = result + "T"
    } else if c == "T" {
      result = result + "A"
    } else if c == "C" {
      result = result + "G"
    } else if c == "G" {
      result = result + "C"
    } else {
      result = result + c.to_owned()
    }
    i = i + 1
  }
  result
}

///|
/// CodonUsageResult struct to hold all codon usage statistics.
pub struct CodonUsageResult {
  codon_counts : Map[String, Int]
  rscu : Map[String, Double]
  cai : Double
  enc : Double
  gc3 : Double
  optimal_codons : Map[String, Double]
}

///|
/// Create a CodonUsageResult.
pub fn CodonUsageResult::new(
  codon_counts : Map[String, Int],
  rscu : Map[String, Double],
  cai : Double,
  enc : Double,
  gc3 : Double,
  optimal_codons : Map[String, Double],
) -> CodonUsageResult {
  CodonUsageResult::{ codon_counts, rscu, cai, enc, gc3, optimal_codons }
}

///|
/// Compute all codon usage statistics for a sequence.
pub fn codon_usage_analysis(seq : String) -> CodonUsageResult {
  let counts = count_codons(seq)
  let rscu_map = codon_rscu(counts.copy())
  let optimal = codon_optimal(rscu_map.copy())
  let cai_val = codon_cai(seq, optimal.copy())
  let enc_val = codon_enc(counts.copy())
  let gc3_val = codon_gc3(seq)

  CodonUsageResult::new(counts, rscu_map, cai_val, enc_val, gc3_val, optimal)
}

///|
/// Compute codon usage statistics with a reference set for CAI.
pub fn codon_usage_with_ref(seq : String, ref_seq : String) -> CodonUsageResult {
  let counts = count_codons(seq)
  let rscu_map = codon_rscu(counts.copy())

  let ref_counts = count_codons(ref_seq)
  let ref_rscu = codon_rscu(ref_counts)
  let optimal = codon_optimal(ref_rscu)

  let cai_val = codon_cai(seq, optimal.copy())
  let enc_val = codon_enc(counts.copy())
  let gc3_val = codon_gc3(seq)

  CodonUsageResult::new(counts, rscu_map, cai_val, enc_val, gc3_val, optimal)
}
