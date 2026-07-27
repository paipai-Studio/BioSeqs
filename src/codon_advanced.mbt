///|
/// Advanced codon analysis module, inspired by Biopython's CodonUsage.
/// Provides comprehensive codon usage analysis with structured data types,
/// organism-specific reference tables, and support for multiple sequence
/// inputs.
///
/// Key features:
/// - CodonUsageTable: Maps codons to observed frequencies and counts
/// - CAIResult: Codon Adaptation Index result with per-codon scores
/// - RSCUResult: Relative Synonymous Codon Usage values
/// - CodonBiasResult: Comprehensive codon bias analysis
/// - Support for high GC, low GC, and general reference tables
/// - Optimal and rare codon identification

///|
/// Natural logarithm approximation using series expansion.
fn codon_ln(x : Double) -> Double {
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
/// Reverse codon lookup table: codon (e.g. "ATG") -> amino acid letter (e.g. "M").
/// Stop codons map to "*".
pub let codon_to_aa_map : Map[String, String] = {
  let m = Map([], capacity=64)
  m.set("TTT", "F")
  m.set("TTC", "F")
  m.set("TTA", "L")
  m.set("TTG", "L")
  m.set("CTT", "L")
  m.set("CTC", "L")
  m.set("CTA", "L")
  m.set("CTG", "L")
  m.set("ATT", "I")
  m.set("ATC", "I")
  m.set("ATA", "I")
  m.set("ATG", "M")
  m.set("GTT", "V")
  m.set("GTC", "V")
  m.set("GTA", "V")
  m.set("GTG", "V")
  m.set("TCT", "S")
  m.set("TCC", "S")
  m.set("TCA", "S")
  m.set("TCG", "S")
  m.set("CCT", "P")
  m.set("CCC", "P")
  m.set("CCA", "P")
  m.set("CCG", "P")
  m.set("ACT", "T")
  m.set("ACC", "T")
  m.set("ACA", "T")
  m.set("ACG", "T")
  m.set("GCT", "A")
  m.set("GCC", "A")
  m.set("GCA", "A")
  m.set("GCG", "A")
  m.set("TAT", "Y")
  m.set("TAC", "Y")
  m.set("CAT", "H")
  m.set("CAC", "H")
  m.set("CAA", "Q")
  m.set("CAG", "Q")
  m.set("AAT", "N")
  m.set("AAC", "N")
  m.set("AAA", "K")
  m.set("AAG", "K")
  m.set("GAT", "D")
  m.set("GAC", "D")
  m.set("GAA", "E")
  m.set("GAG", "E")
  m.set("TGT", "C")
  m.set("TGC", "C")
  m.set("TGG", "W")
  m.set("CGT", "R")
  m.set("CGC", "R")
  m.set("CGA", "R")
  m.set("CGG", "R")
  m.set("AGT", "S")
  m.set("AGC", "S")
  m.set("GGT", "G")
  m.set("GGC", "G")
  m.set("GGA", "G")
  m.set("GGG", "G")
  m.set("AGA", "R")
  m.set("AGG", "R")
  m.set("TAA", "*")
  m.set("TAG", "*")
  m.set("TGA", "*")
  m
}

///|
/// Amino acid to synonymous codons mapping (standard genetic code).
pub let aa_to_codons : Map[String, Array[String]] = {
  let m = Map([], capacity=21)
  m.set("A", ["GCA", "GCC", "GCG", "GCT"])
  m.set("R", ["AGA", "AGG", "CGA", "CGC", "CGG", "CGT"])
  m.set("N", ["AAC", "AAT"])
  m.set("D", ["GAC", "GAT"])
  m.set("C", ["TGC", "TGT"])
  m.set("E", ["GAA", "GAG"])
  m.set("Q", ["CAA", "CAG"])
  m.set("G", ["GGA", "GGC", "GGG", "GGT"])
  m.set("H", ["CAC", "CAT"])
  m.set("I", ["ATA", "ATC", "ATT"])
  m.set("L", ["CTA", "CTC", "CTG", "CTT", "TTA", "TTG"])
  m.set("K", ["AAA", "AAG"])
  m.set("M", ["ATG"])
  m.set("F", ["TTC", "TTT"])
  m.set("P", ["CCA", "CCC", "CCG", "CCT"])
  m.set("S", ["AGC", "AGT", "TCA", "TCC", "TCG", "TCT"])
  m.set("T", ["ACA", "ACC", "ACG", "ACT"])
  m.set("W", ["TGG"])
  m.set("Y", ["TAC", "TAT"])
  m.set("V", ["GTA", "GTC", "GTG", "GTT"])
  m.set("*", ["TAA", "TAG", "TGA"])
  m
}

///|
/// Compute codon frequencies from counts.
/// Returns a Map mapping each codon to its frequency (count / total).
fn compute_frequencies(counts : Map[String, Int]) -> Map[String, Double] {
  let freq = Map([], capacity=64)
  let keys = counts.keys().collect()
  let mut total = 0
  let mut i = 0
  while i < keys.length() {
    total = total + counts.get(keys[i]).unwrap_or(0)
    i = i + 1
  }
  if total == 0 {
    return freq
  }
  let mut k = 0
  while k < keys.length() {
    let codon = keys[k]
    let count = counts.get(codon).unwrap_or(0).to_double()
    freq.set(codon, count / total.to_double())
    k = k + 1
  }
  freq
}

///|
/// CodonUsageTable: Maps codons to observed frequencies and counts.
pub struct CodonUsageTable {
  codon_counts : Map[String, Int]
  total_codons : Int
  frequencies : Map[String, Double]
}

///|
/// Create a new CodonUsageTable from raw counts.
pub fn CodonUsageTable::new(
  codon_counts : Map[String, Int],
) -> CodonUsageTable {
  let frequencies = compute_frequencies(codon_counts.copy())
  let mut total = 0
  let keys = codon_counts.keys().collect()
  let mut i = 0
  while i < keys.length() {
    total = total + codon_counts.get(keys[i]).unwrap_or(0)
    i = i + 1
  }
  CodonUsageTable::{ codon_counts, total_codons: total, frequencies }
}

///|
/// Get the count for a specific codon.
pub fn CodonUsageTable::get_count(
  usage : CodonUsageTable,
  codon : String,
) -> Int {
  usage.codon_counts.get(codon).unwrap_or(0)
}

///|
/// Get the frequency for a specific codon.
pub fn CodonUsageTable::get_frequency(
  usage : CodonUsageTable,
  codon : String,
) -> Double {
  usage.frequencies.get(codon).unwrap_or(0.0)
}

///|
/// Get the total number of codons analyzed.
pub fn CodonUsageTable::total(usage : CodonUsageTable) -> Int {
  usage.total_codons
}

///|
/// Get all codon counts as a Map.
pub fn CodonUsageTable::counts(usage : CodonUsageTable) -> Map[String, Int] {
  usage.codon_counts.copy()
}

///|
/// Get all codon frequencies as a Map.
pub fn CodonUsageTable::frequencies(usage : CodonUsageTable) -> Map[String, Double] {
  usage.frequencies.copy()
}

///|
/// CAIResult: Codon Adaptation Index result with scores.
pub struct CAIResult {
  score : Double
  codon_scores : Map[String, Double]
  sequence_length : Int
}

///|
/// Create a new CAIResult.
pub fn CAIResult::new(
  score : Double,
  codon_scores : Map[String, Double],
  sequence_length : Int,
) -> CAIResult {
  CAIResult::{ score, codon_scores, sequence_length }
}

///|
/// Get the CAI score.
pub fn CAIResult::score(result : CAIResult) -> Double {
  result.score
}

///|
/// Get per-codon CAI scores.
pub fn CAIResult::codon_scores(result : CAIResult) -> Map[String, Double] {
  result.codon_scores.copy()
}

///|
/// RSCUResult: Relative Synonymous Codon Usage values.
pub struct RSCUResult {
  values : Map[String, Double]
  total_counts : Map[String, Int]
}

///|
/// Create a new RSCUResult.
pub fn RSCUResult::new(
  values : Map[String, Double],
  total_counts : Map[String, Int],
) -> RSCUResult {
  RSCUResult::{ values, total_counts }
}

///|
/// Get the RSCU value for a specific codon.
pub fn RSCUResult::get_value(result : RSCUResult, codon : String) -> Double {
  result.values.get(codon).unwrap_or(0.0)
}

///|
/// Get all RSCU values.
pub fn RSCUResult::values(result : RSCUResult) -> Map[String, Double] {
  result.values.copy()
}

///|
/// CodonBiasResult: Comprehensive codon bias analysis.
pub struct CodonBiasResult {
  enc : Double
  cai : Double
  gc3 : Double
  bias_score : Double
  optimal_fraction : Double
  rare_fraction : Double
  total_codons_analyzed : Int
}

///|
/// Create a new CodonBiasResult.
pub fn CodonBiasResult::new(
  enc : Double,
  cai : Double,
  gc3 : Double,
  bias_score : Double,
  optimal_fraction : Double,
  rare_fraction : Double,
  total_codons_analyzed : Int,
) -> CodonBiasResult {
  CodonBiasResult::{
    enc,
    cai,
    gc3,
    bias_score,
    optimal_fraction,
    rare_fraction,
    total_codons_analyzed,
  }
}

///|
/// Build a CodonUsageTable from a single DNA sequence.
/// Counts codons in-frame (starting from position 0, every 3 bases).
pub fn codon_usage_from_sequence(seq : String) -> CodonUsageTable {
  let seq_upper = seq.to_upper()
  let len = seq_upper.length()
  let counts = Map([], capacity=64)

  let mut i = 0
  while i + 3 <= len {
    let codon = seq_upper[i:i + 3].to_owned()
    if codon_to_aa_map.contains(codon) {
      counts.set(codon, counts.get(codon).unwrap_or(0) + 1)
    }
    i = i + 3
  }

  CodonUsageTable::new(counts)
}

///|
/// Build a CodonUsageTable from multiple DNA sequences.
/// Useful for building reference tables from a set of genes.
pub fn codon_usage_from_sequences(seqs : Array[String]) -> CodonUsageTable {
  let counts = Map([], capacity=64)

  let mut i = 0
  while i < seqs.length() {
    let seq_upper = seqs[i].to_upper()
    let len = seq_upper.length()
    let mut j = 0
    while j + 3 <= len {
      let codon = seq_upper[j:j + 3].to_owned()
      if codon_to_aa_map.contains(codon) {
        counts.set(codon, counts.get(codon).unwrap_or(0) + 1)
      }
      j = j + 3
    }
    i = i + 1
  }

  CodonUsageTable::new(counts)
}

///|
/// Get the amino acid for a given codon.
/// Returns "X" for unknown codons.
pub fn get_amino_acid(codon : String) -> String {
  let codon_upper = codon.to_upper()
  codon_to_aa_map.get(codon_upper).unwrap_or("X")
}

///|
/// Get all synonymous codons for an amino acid.
pub fn get_codons_for_amino_acid(aa : String) -> Array[String] {
  aa_to_codons.get(aa).unwrap_or([])
}

///|
/// Check if a codon is a start codon (ATG).
pub fn is_start_codon(codon : String) -> Bool {
  codon.to_upper() == "ATG"
}

///|
/// Calculate the relative adaptiveness (w) for each codon.
/// w_ij = x_ij / x_i_max, where x_ij is the frequency of codon j
/// for amino acid i, and x_i_max is the maximum frequency among
/// synonymous codons for that amino acid.
pub fn compute_codon_adaptiveness(
  usage : CodonUsageTable,
) -> Map[String, Double] {
  let w = Map([], capacity=64)
  let aa_list = aa_to_codons.keys().collect()

  let mut i = 0
  while i < aa_list.length() {
    let aa = aa_list[i]
    let codons = aa_to_codons[aa]

    let mut max_freq = 0.0
    let mut j = 0
    while j < codons.length() {
      let freq = CodonUsageTable::get_frequency(usage, codons[j])
      if freq > max_freq {
        max_freq = freq
      }
      j = j + 1
    }

    if max_freq > 0.0 {
      let mut k = 0
      while k < codons.length() {
        let codon = codons[k]
        let freq = CodonUsageTable::get_frequency(usage, codon)
        w.set(codon, freq / max_freq)
        k = k + 1
      }
    } else {
      let mut k = 0
      while k < codons.length() {
        w.set(codons[k], 0.0)
        k = k + 1
      }
    }
    i = i + 1
  }

  w
}

///|
/// Calculate Codon Adaptation Index (CAI) for a sequence
/// given a reference CodonUsageTable.
/// CAI = exp( (1/L) * sum(ln(w_i)) ) for each codon,
/// where w_i is the relative adaptiveness of codon i.
/// Returns a value between 0 and 1, where 1 indicates optimal codon usage.
pub fn calculate_cai(
  seq : String,
  reference_table : CodonUsageTable,
) -> Double {
  let seq_upper = seq.to_upper()
  let len = seq_upper.length()
  let adaptiveness = compute_codon_adaptiveness(reference_table)

  let mut total_log = 0.0
  let mut count = 0

  let mut i = 0
  while i + 3 <= len {
    let codon = seq_upper[i:i + 3].to_owned()
    let aa = get_amino_acid(codon)
    if aa != "*" && aa != "X" {
      let w = adaptiveness.get(codon).unwrap_or(0.0)
      if w > 0.0 {
        total_log = total_log + codon_ln(w)
        count = count + 1
      }
    }
    i = i + 3
  }

  if count == 0 {
    return 0.0
  }

  @math.exp(total_log / count.to_double())
}

///|
/// Calculate CAI with detailed per-codon scores.
pub fn calculate_cai_detailed(
  seq : String,
  reference_table : CodonUsageTable,
) -> CAIResult {
  let seq_upper = seq.to_upper()
  let len = seq_upper.length()
  let adaptiveness = compute_codon_adaptiveness(reference_table)
  let codon_scores = Map([], capacity=64)

  let mut total_log = 0.0
  let mut count = 0
  let mut analyzed = 0

  let mut i = 0
  while i + 3 <= len {
    let codon = seq_upper[i:i + 3].to_owned()
    let aa = get_amino_acid(codon)
    if aa != "*" && aa != "X" {
      let w = adaptiveness.get(codon).unwrap_or(0.0)
      codon_scores.set(codon, w)
      if w > 0.0 {
        total_log = total_log + codon_ln(w)
        count = count + 1
      }
      analyzed = analyzed + 1
    }
    i = i + 3
  }

  let score = if count == 0 {
    0.0
  } else {
    @math.exp(total_log / count.to_double())
  }

  CAIResult::new(score, codon_scores, analyzed)
}

///|
/// Calculate Relative Synonymous Codon Usage (RSCU) values.
/// RSCU = observed frequency / expected frequency under equal usage.
/// RSCU > 1 indicates preferred codon, RSCU < 1 indicates avoided codon.
pub fn calculate_rscu(usage : CodonUsageTable) -> RSCUResult {
  let values = Map([], capacity=64)
  let total_counts = Map([], capacity=64)
  let aa_list = aa_to_codons.keys().collect()

  let mut i = 0
  while i < aa_list.length() {
    let aa = aa_list[i]
    let codons = aa_to_codons[aa]
    let n_syn = codons.length().to_double()

    let mut total = 0
    let mut j = 0
    while j < codons.length() {
      total = total + CodonUsageTable::get_count(usage, codons[j])
      j = j + 1
    }

    if total > 0 {
      let mut k = 0
      while k < codons.length() {
        let codon = codons[k]
        let count = CodonUsageTable::get_count(usage, codon).to_double()
        let expected = total.to_double() / n_syn
        let rscu_val = if expected > 0.0 { count / expected } else { 0.0 }
        values.set(codon, rscu_val)
        total_counts.set(codon, count.to_int())
        k = k + 1
      }
    } else {
      let mut k = 0
      while k < codons.length() {
        values.set(codons[k], 0.0)
        total_counts.set(codons[k], 0)
        k = k + 1
      }
    }
    i = i + 1
  }

  RSCUResult::new(values, total_counts)
}

///|
/// Calculate Effective Number of Codons (ENC) from usage table.
/// ENC ranges from 20 (extreme bias) to 61 (no bias).
/// Uses Wright (1990) formula: ENC = 2 + 9*F2 + 1*F3 + 5*F4 + 3*F6
/// where Fk is the average effective codon count for k-fold families.
pub fn calculate_enc(usage : CodonUsageTable) -> Double {
  let mut f2_contrib = 0.0
  let mut f2_count = 0
  let mut f3_contrib = 0.0
  let mut f3_count = 0
  let mut f4_contrib = 0.0
  let mut f4_count = 0
  let mut f6_contrib = 0.0
  let mut f6_count = 0

  let aa_list = aa_to_codons.keys().collect()

  let mut i = 0
  while i < aa_list.length() {
    let aa = aa_list[i]
    if aa == "*" || aa == "M" || aa == "W" {
      i = i + 1
      continue
    }
    let codons = aa_to_codons[aa]
    let k = codons.length()

    let mut total = 0
    let mut j = 0
    while j < codons.length() {
      total = total + CodonUsageTable::get_count(usage, codons[j])
      j = j + 1
    }

    if total > 1 {
      let mut sum_pi2 = 0.0
      let mut kk = 0
      while kk < codons.length() {
        let pi = CodonUsageTable::get_count(usage, codons[kk]).to_double() /
          total.to_double()
        sum_pi2 = sum_pi2 + pi * pi
        kk = kk + 1
      }
      let f_hat = (total.to_double() * sum_pi2 - 1.0) /
        (total.to_double() - 1.0)
      let contrib = if f_hat > 0.0 { 1.0 / f_hat } else { k.to_double() }

      if k == 2 {
        f2_contrib = f2_contrib + contrib
        f2_count = f2_count + 1
      } else if k == 3 {
        f3_contrib = f3_contrib + contrib
        f3_count = f3_count + 1
      } else if k == 4 {
        f4_contrib = f4_contrib + contrib
        f4_count = f4_count + 1
      } else if k == 6 {
        f6_contrib = f6_contrib + contrib
        f6_count = f6_count + 1
      }
    } else if total == 1 {
      let contrib = k.to_double()
      if k == 2 {
        f2_contrib = f2_contrib + contrib
        f2_count = f2_count + 1
      } else if k == 3 {
        f3_contrib = f3_contrib + contrib
        f3_count = f3_count + 1
      } else if k == 4 {
        f4_contrib = f4_contrib + contrib
        f4_count = f4_count + 1
      } else if k == 6 {
        f6_contrib = f6_contrib + contrib
        f6_count = f6_count + 1
      }
    }
    i = i + 1
  }

  let mut enc = 2.0

  if f2_count > 0 {
    enc = enc + 9.0 * (f2_contrib / f2_count.to_double())
  }
  if f3_count > 0 {
    enc = enc + 1.0 * (f3_contrib / f3_count.to_double())
  }
  if f4_count > 0 {
    enc = enc + 5.0 * (f4_contrib / f4_count.to_double())
  }
  if f6_count > 0 {
    enc = enc + 3.0 * (f6_contrib / f6_count.to_double())
  }

  if enc > 61.0 { 61.0 } else if enc < 20.0 { 20.0 } else { enc }
}

///|
/// Calculate GC content at the third codon position (GC3).
pub fn calculate_gc3(seq : String) -> Double {
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
/// Calculate a composite bias score (0-1 scale).
/// Higher values indicate stronger codon usage bias.
fn calculate_bias_score(enc : Double, gc3 : Double, optimal_frac : Double) -> Double {
  let enc_score = (61.0 - enc) / 41.0
  let gc3_score = if gc3 > 0.5 { gc3 } else { 1.0 - gc3 }
  let opt_score = optimal_frac
  (enc_score + gc3_score + opt_score) / 3.0
}

///|
/// Count stop codons in a sequence.
fn count_stop_codons(seq : String) -> Int {
  let seq_upper = seq.to_upper()
  let len = seq_upper.length()
  let mut count = 0
  let mut i = 0
  while i + 3 <= len {
    let codon = seq_upper[i:i + 3].to_owned()
    if is_stop_codon(codon) {
      count = count + 1
    }
    i = i + 3
  }
  count
}

///|
/// Identify optimal codons for each amino acid.
/// Returns an array of (amino_acid, codon) tuples where the codon
/// has the highest frequency for that amino acid.
pub fn optimal_codons(
  usage : CodonUsageTable,
) -> Array[(String, String)] {
  let result = Array::make(0, ("", ""))
  let aa_list = aa_to_codons.keys().collect()

  let mut i = 0
  while i < aa_list.length() {
    let aa = aa_list[i]
    if aa == "*" {
      i = i + 1
      continue
    }
    let codons = aa_to_codons[aa]

    let mut max_freq = -1.0
    let mut best_codon = ""
    let mut j = 0
    while j < codons.length() {
      let freq = CodonUsageTable::get_frequency(usage, codons[j])
      if freq > max_freq {
        max_freq = freq
        best_codon = codons[j]
      }
      j = j + 1
    }

    if best_codon != "" && max_freq > 0.0 {
      result.push((aa, best_codon))
    }
    i = i + 1
  }

  result
}

///|
/// Identify rare codons in a usage table.
/// A codon is considered rare if its frequency is less than
/// 10% of the maximum frequency among synonymous codons
/// for the same amino acid.
pub fn rare_codons(usage : CodonUsageTable) -> Array[String] {
  let result = Array::make(0, "")
  let aa_list = aa_to_codons.keys().collect()

  let mut i = 0
  while i < aa_list.length() {
    let aa = aa_list[i]
    if aa == "*" {
      i = i + 1
      continue
    }
    let codons = aa_to_codons[aa]

    let mut max_freq = 0.0
    let mut j = 0
    while j < codons.length() {
      let freq = CodonUsageTable::get_frequency(usage, codons[j])
      if freq > max_freq {
        max_freq = freq
      }
      j = j + 1
    }

    let threshold = max_freq * 0.1
    let mut k = 0
    while k < codons.length() {
      let codon = codons[k]
      let freq = CodonUsageTable::get_frequency(usage, codon)
      if max_freq > 0.0 && freq < threshold && codons.length() > 1 {
        result.push(codon)
      }
      k = k + 1
    }
    i = i + 1
  }

  result
}

///|
/// Comprehensive codon bias analysis for a DNA sequence.
/// Returns ENC, CAI, GC3, bias score, optimal/rare codon fractions.
pub fn codon_bias_analysis(seq : String) -> CodonBiasResult {
  let usage = codon_usage_from_sequence(seq)
  let enc_val = calculate_enc(usage)
  let cai_val = calculate_cai(seq, usage)
  let gc3_val = calculate_gc3(seq)

  let total_codons = CodonUsageTable::total(usage)

  let stop_codon_count = count_stop_codons(seq)
  let non_stop_total = total_codons - stop_codon_count

  let optimal = optimal_codons(usage)
  let rare = rare_codons(usage)

  let optimal_set = Map([], capacity=optimal.length())
  let mut i = 0
  while i < optimal.length() {
    optimal_set.set(optimal[i].1, true)
    i = i + 1
  }

  let rare_set = Map([], capacity=rare.length())
  let mut r = 0
  while r < rare.length() {
    rare_set.set(rare[r], true)
    r = r + 1
  }

  let mut optimal_count = 0
  let mut rare_count = 0
  let seq_upper = seq.to_upper()
  let len = seq_upper.length()
  let mut j = 0
  while j + 3 <= len {
    let codon = seq_upper[j:j + 3].to_owned()
    let aa = get_amino_acid(codon)
    if aa != "*" && aa != "X" {
      if optimal_set.contains(codon) {
        optimal_count = optimal_count + 1
      }
      if rare_set.contains(codon) {
        rare_count = rare_count + 1
      }
    }
    j = j + 3
  }

  let optimal_frac = if non_stop_total > 0 {
    optimal_count.to_double() / non_stop_total.to_double()
  } else {
    0.0
  }

  let rare_frac = if non_stop_total > 0 {
    rare_count.to_double() / non_stop_total.to_double()
  } else {
    0.0
  }

  let bias_score = calculate_bias_score(enc_val, gc3_val, optimal_frac)

  CodonBiasResult::new(
    enc_val,
    cai_val,
    gc3_val,
    bias_score,
    optimal_frac,
    rare_frac,
    non_stop_total,
  )
}

///|
/// Get a high GC content reference table.
/// Returns adaptiveness values biased toward G/C ending codons.
pub fn get_high_gc_reference() -> Map[String, Double] {
  let ref_table = Map([], capacity=64)
  let gc_codons = [
    "GCA", "GCC", "GCG", "GCT",
    "CGA", "CGC", "CGG", "CGT",
    "GAC", "GAT", "GAA", "GAG",
    "GGA", "GGC", "GGG", "GGT",
    "CAC", "CAT", "CAA", "CAG",
    "ATA", "ATC", "ATT", "ATG",
    "CTA", "CTC", "CTG", "CTT",
    "AAA", "AAG", "TTA", "TTG",
    "TTC", "TTT", "CCA", "CCC",
    "CCG", "CCT", "AGC", "AGT",
    "TCA", "TCC", "TCG", "TCT",
    "ACA", "ACC", "ACG", "ACT",
    "GTA", "GTC", "GTG", "GTT",
    "TGC", "TGT", "TGG",
  ]
  let gc_bias = [
    0.7, 1.0, 0.9, 0.6,
    0.8, 1.0, 0.9, 0.6,
    0.8, 0.5, 0.9, 1.0,
    0.6, 1.0, 0.8, 0.5,
    0.8, 0.5, 0.7, 1.0,
    0.4, 0.8, 0.5, 1.0,
    0.6, 0.9, 1.0, 0.5,
    0.5, 1.0, 0.4, 0.3,
    0.5, 0.3, 0.7, 1.0,
    0.9, 0.6, 0.8, 0.5,
    0.4, 0.8, 1.0, 0.6,
    0.7, 1.0, 0.9, 0.5,
    0.5, 0.9, 1.0, 0.6,
    0.5, 0.6, 1.0,
  ]
  let mut i = 0
  while i < gc_codons.length() {
    ref_table.set(gc_codons[i], gc_bias[i])
    i = i + 1
  }
  ref_table
}

///|
/// Get a low GC content reference table.
/// Returns adaptiveness values biased toward A/T ending codons.
pub fn get_low_gc_reference() -> Map[String, Double] {
  let ref_table = Map([], capacity=64)
  let at_codons = [
    "GCT", "GCA", "GCC", "GCG",
    "CGT", "CGA", "CGC", "CGG",
    "GAT", "GAC", "GAA", "GAG",
    "GTT", "GTA", "GTC", "GTG",
    "CAT", "CAC", "CAA", "CAG",
    "ATT", "ATA", "ATC", "ATG",
    "CTT", "CTA", "CTC", "CTG",
    "AAA", "AAG", "TTA", "TTG",
    "TTT", "TTC", "CCT", "CCA",
    "CCC", "CCG", "AGT", "AGC",
    "TCT", "TCA", "TCC", "TCG",
    "ACT", "ACA", "ACC", "ACG",
    "GGT", "GGA", "GGC", "GGG",
    "TGT", "TGC", "TGG",
  ]
  let at_bias = [
    1.0, 0.5, 0.4, 0.3,
    1.0, 0.5, 0.4, 0.3,
    1.0, 0.5, 0.6, 0.5,
    1.0, 0.8, 0.5, 0.4,
    1.0, 0.4, 0.8, 0.6,
    1.0, 0.7, 0.4, 1.0,
    1.0, 0.8, 0.5, 0.4,
    1.0, 0.6, 0.9, 0.7,
    1.0, 0.6, 0.5, 0.3,
    1.0, 0.8, 0.7, 0.5,
    1.0, 0.8, 0.5, 0.4,
    1.0, 0.9, 0.6, 0.4,
    1.0, 0.9, 0.6, 0.4,
    1.0, 0.5, 1.0,
  ]
  let mut i = 0
  while i < at_codons.length() {
    ref_table.set(at_codons[i], at_bias[i])
    i = i + 1
  }
  ref_table
}

///|
/// Get a reference codon adaptiveness table for a given organism type.
/// Supported types: "high_gc", "low_gc", "general".
pub fn get_reference_table(org : String) -> Map[String, Double] {
  let org_lower = org.to_lower()
  if org_lower == "high_gc" {
    get_high_gc_reference()
  } else if org_lower == "low_gc" {
    get_low_gc_reference()
  } else {
    let general = Map([], capacity=64)
    let aa_list = aa_to_codons.keys().collect()
    let mut i = 0
    while i < aa_list.length() {
      let aa = aa_list[i]
      let codons = aa_to_codons[aa]
      if codons.length() > 0 {
        let mut j = 0
        while j < codons.length() {
          general.set(codons[j], 1.0)
          j = j + 1
        }
      }
      i = i + 1
    }
    general
  }
}

///|
/// Calculate CAI using a named organism reference table.
pub fn calculate_cai_with_org(
  seq : String,
  org : String,
) -> Double {
  let seq_upper = seq.to_upper()
  let len = seq_upper.length()
  let ref_table = get_reference_table(org)

  let mut total_log = 0.0
  let mut count = 0

  let mut i = 0
  while i + 3 <= len {
    let codon = seq_upper[i:i + 3].to_owned()
    let aa = get_amino_acid(codon)
    if aa != "*" && aa != "X" {
      let w = ref_table.get(codon).unwrap_or(0.0)
      if w > 0.0 {
        total_log = total_log + codon_ln(w)
        count = count + 1
      }
    }
    i = i + 3
  }

  if count == 0 {
    return 0.0
  }

  @math.exp(total_log / count.to_double())
}

///|
/// Get amino acid composition from a usage table.
/// Returns Map of amino acid -> fraction of total codons.
pub fn get_amino_acid_composition(
  usage : CodonUsageTable,
) -> Map[String, Double] {
  let comp = Map([], capacity=21)
  let aa_list = aa_to_codons.keys().collect()
  let total = CodonUsageTable::total(usage)

  if total == 0 {
    return comp
  }

  let mut i = 0
  while i < aa_list.length() {
    let aa = aa_list[i]
    if aa == "*" {
      i = i + 1
      continue
    }
    let codons = aa_to_codons[aa]

    let mut aa_count = 0
    let mut j = 0
    while j < codons.length() {
      aa_count = aa_count + CodonUsageTable::get_count(usage, codons[j])
      j = j + 1
    }

    comp.set(aa, aa_count.to_double() / total.to_double())
    i = i + 1
  }

  comp
}

///|
/// Format a CodonUsageTable as a human-readable string.
pub fn codon_usage_to_string(usage : CodonUsageTable) -> String {
  let mut result = "CodonUsageTable (total: " +
    CodonUsageTable::total(usage).to_string() +
    ")\n"
  let aa_list = aa_to_codons.keys().collect()

  let mut i = 0
  while i < aa_list.length() {
    let aa = aa_list[i]
    let codons = aa_to_codons[aa]
    result = result + aa + ": "

    let mut j = 0
    while j < codons.length() {
      let codon = codons[j]
      let count = CodonUsageTable::get_count(usage, codon)
      let freq = CodonUsageTable::get_frequency(usage, codon)
      result = result + codon + "(" + count.to_string() +
        "/" + freq.to_string() + ")"
      if j < codons.length() - 1 {
        result = result + ", "
      }
      j = j + 1
    }
    result = result + "\n"
    i = i + 1
  }

  result
}

///|
/// Count codons of a specific amino acid in a sequence.
pub fn count_amino_acid_codons(
  seq : String,
  aa : String,
) -> Int {
  let seq_upper = seq.to_upper()
  let len = seq_upper.length()
  let codons = aa_to_codons.get(aa).unwrap_or([])
  let mut count = 0

  let mut i = 0
  while i + 3 <= len {
    let codon = seq_upper[i:i + 3].to_owned()
    let mut j = 0
    while j < codons.length() {
      if codon == codons[j] {
        count = count + 1
        break
      }
      j = j + 1
    }
    i = i + 3
  }

  count
}

///|
/// Calculate the fraction of codons that are AT-rich (A or T at third position).
pub fn at_third_fraction(seq : String) -> Double {
  let seq_upper = seq.to_upper()
  let len = seq_upper.length()

  let mut at_count = 0
  let mut total = 0

  let mut i = 0
  while i + 3 <= len {
    let third = seq_upper[i + 2]
    if third == 'A' || third == 'T' {
      at_count = at_count + 1
    }
    total = total + 1
    i = i + 3
  }

  if total == 0 {
    return 0.0
  }

  at_count.to_double() / total.to_double()
}

///|
/// Calculate the fraction of codons with G or C at first position.
pub fn gc_first_fraction(seq : String) -> Double {
  let seq_upper = seq.to_upper()
  let len = seq_upper.length()

  let mut gc_count = 0
  let mut total = 0

  let mut i = 0
  while i + 3 <= len {
    let first = seq_upper[i]
    if first == 'G' || first == 'C' {
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
/// Calculate the fraction of codons with G or C at second position.
pub fn gc_second_fraction(seq : String) -> Double {
  let seq_upper = seq.to_upper()
  let len = seq_upper.length()

  let mut gc_count = 0
  let mut total = 0

  let mut i = 0
  while i + 3 <= len {
    let second = seq_upper[i + 1]
    if second == 'G' || second == 'C' {
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