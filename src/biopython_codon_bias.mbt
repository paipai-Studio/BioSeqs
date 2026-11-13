/// biopython_codon_bias.mbt — Codon usage bias metrics
/// Inspired by Biopython Bio.SeqUtils.CodonUsage: CAI, Nc, RSCU

///|
/// Standard genetic code: maps codon to amino acid.
fn biopython_cb_codon_to_aa(codon : String) -> String {
  let aa = codon_lookup(codon)
  if aa == "?" {
    "X"
  } else {
    aa
  }
}

///|
/// Count codon frequencies in a coding sequence.
pub fn biopython_codon_counts(sequence : String) -> Array[(String, Int32)] {
  let seq = sequence.to_upper()
  let n = seq.length()
  let counts : Array[(String, Int32)] = Array(capacity=0)
  let mut i = 0
  while i + 3 <= n {
    let codon = seq[i:i + 3].to_owned()
    let mut found = false
    for j = 0; j < counts.length(); j = j + 1 {
      if counts[j].0 == codon {
        let (c, cnt) = counts[j]
        counts[j] = (c, cnt + 1)
        found = true
        break
      }
    }
    if !found {
      counts.push((codon, 1))
    }
    i = i + 3
  }
  counts
}

///|
/// Relative Synonymous Codon Usage (RSCU).
/// RSCU = count(codon) / (sum of counts for synonymous codons / number of synonymous codons)
pub fn biopython_rscu(sequence : String) -> Array[(String, Double)] {
  let counts = biopython_codon_counts(sequence)
  // Group codons by amino acid
  let aas : Array[String] = Array(capacity=0)
  let aa_codons : Array[Array[String]] = Array(capacity=0)
  let aa_counts : Array[Array[Int32]] = Array(capacity=0)
  for i = 0; i < counts.length(); i = i + 1 {
    let (codon, count) = counts[i]
    let aa = biopython_cb_codon_to_aa(codon)
    let mut found = false
    for j = 0; j < aas.length(); j = j + 1 {
      if aas[j] == aa {
        aa_codons[j].push(codon)
        aa_counts[j].push(count)
        found = true
        break
      }
    }
    if !found {
      aas.push(aa)
      aa_codons.push([codon])
      aa_counts.push([count])
    }
  }
  // Calculate RSCU for each codon
  let result : Array[(String, Double)] = Array(capacity=0)
  for j = 0; j < aas.length(); j = j + 1 {
    let n_syn = aa_codons[j].length()
    let mut total = 0
    for k = 0; k < n_syn; k = k + 1 {
      total = total + aa_counts[j][k]
    }
    if total > 0 && n_syn > 0 {
      for k = 0; k < n_syn; k = k + 1 {
        let rscu_val = aa_counts[j][k].to_double() *
          n_syn.to_double() /
          total.to_double()
        result.push((aa_codons[j][k], rscu_val))
      }
    }
  }
  result
}

///|
/// Effective Number of Codons (Nc).
/// Measures codon usage bias: Nc ranges from 20 (one codon per aa) to 61 (uniform).
pub fn biopython_effective_number_of_codons(sequence : String) -> Double {
  let rscu = biopython_rscu(sequence)
  // Group RSCU by amino acid
  let aas : Array[String] = Array(capacity=0)
  let aa_rscus : Array[Array[Double]] = Array(capacity=0)
  for i = 0; i < rscu.length(); i = i + 1 {
    let (codon, rscu_val) = rscu[i]
    let aa = biopython_cb_codon_to_aa(codon)
    let mut found = false
    for j = 0; j < aas.length(); j = j + 1 {
      if aas[j] == aa {
        aa_rscus[j].push(rscu_val)
        found = true
        break
      }
    }
    if !found {
      aas.push(aa)
      aa_rscus.push([rscu_val])
    }
  }
  // Nc = sum over amino acids of F_aa
  // F_aa = n / sum(RSCU^2)
  let mut nc = 0.0
  for j = 0; j < aas.length(); j = j + 1 {
    let n = aa_rscus[j].length()
    let mut sum_sq = 0.0
    for k = 0; k < n; k = k + 1 {
      sum_sq = sum_sq + aa_rscus[j][k] * aa_rscus[j][k]
    }
    if sum_sq > 0.0 {
      nc = nc + n.to_double() * n.to_double() / sum_sq
    }
  }
  nc
}

///|
/// Codon Adaptation Index (CAI).
/// CAI = geometric mean of relative adaptiveness (w) values.
/// w_codon = RSCU_codon / RSCU_max(for that amino acid)
pub fn biopython_codon_adaptation_index(
  sequence : String,
  reference_rscu~ : Array[(String, Double)] = [],
) -> Double {
  let rscu = if reference_rscu.length() > 0 {
    reference_rscu
  } else {
    biopython_rscu(sequence)
  }
  // Build max RSCU per amino acid
  let aa_max : Array[(String, Double)] = Array(capacity=0)
  for i = 0; i < rscu.length(); i = i + 1 {
    let (codon, rscu_val) = rscu[i]
    let aa = biopython_cb_codon_to_aa(codon)
    let mut found = false
    for j = 0; j < aa_max.length(); j = j + 1 {
      if aa_max[j].0 == aa {
        if rscu_val > aa_max[j].1 {
          aa_max[j] = (aa, rscu_val)
        }
        found = true
        break
      }
    }
    if !found {
      aa_max.push((aa, rscu_val))
    }
  }
  // Calculate CAI
  let counts = biopython_codon_counts(sequence)
  let mut sum_ln_w = 0.0
  let mut n_codons = 0
  for i = 0; i < counts.length(); i = i + 1 {
    let (codon, count) = counts[i]
    let aa = biopython_cb_codon_to_aa(codon)
    // Find RSCU for this codon
    let mut codon_rscu = 0.0
    for j = 0; j < rscu.length(); j = j + 1 {
      if rscu[j].0 == codon {
        codon_rscu = rscu[j].1
        break
      }
    }
    // Find max RSCU for this aa
    let mut max_rscu = 0.0
    for j = 0; j < aa_max.length(); j = j + 1 {
      if aa_max[j].0 == aa {
        max_rscu = aa_max[j].1
        break
      }
    }
    if max_rscu > 0.0 && codon_rscu > 0.0 {
      let w = codon_rscu / max_rscu
      sum_ln_w = sum_ln_w + count.to_double() * @math.ln(w)
      n_codons = n_codons + count
    }
  }
  if n_codons == 0 {
    0.0
  } else {
    @math.exp(sum_ln_w / n_codons.to_double())
  }
}
