/// skbio_codon_genetic.mbt — Codon usage bias and genetic code tables
/// Items 8-9: skbio_codon_usage (CAI/RSCU/ENC), skbio_genetic_code

///|
/// Genetic code table.
pub struct SkbioGeneticCode {
  id : Int32
  name : String
  table : Array[(String, String)] // (codon, amino_acid) pairs
}

///|
/// Get standard genetic code (NCBI translation table 1).
pub fn skbio_genetic_code_standard() -> SkbioGeneticCode {
  let aa = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
  let bases = "TCAG"
  let table : Array[(String, String)] = Array(capacity=64)
  for i = 0; i < 4; i = i + 1 {
    for j = 0; j < 4; j = j + 1 {
      for k = 0; k < 4; k = k + 1 {
        let codon = StringBuilder()
        let _ = codon.write_char(bases[i].to_int().unsafe_to_char())
        let _ = codon.write_char(bases[j].to_int().unsafe_to_char())
        let _ = codon.write_char(bases[k].to_int().unsafe_to_char())
        let idx = i * 16 + j * 4 + k
        let amino = aa[idx].to_int().unsafe_to_char().to_string()
        table.push((codon.to_string(), amino))
      }
    }
  }
  SkbioGeneticCode::{ id: 1, name: "Standard", table, }
}

///|
/// Get vertebrate mitochondrial genetic code (NCBI table 2).
pub fn skbio_genetic_code_mito() -> SkbioGeneticCode {
  let code = skbio_genetic_code_standard()
  let overrides : Array[(String, String)] = [
    ("TGA", "W"), // TGA = W (not stop)
    ("AGA", "R"), // AGA = R
    ("AGG", "R"), // AGG = R
    ("ATA", "M"), // ATA = M
  ]
  let table : Array[(String, String)] = Array(capacity=64)
  for i = 0; i < code.table.length(); i = i + 1 {
    let (codon, aa) = code.table[i]
    let mut found = aa
    for j = 0; j < overrides.length(); j = j + 1 {
      if overrides[j].0 == codon {
        found = overrides[j].1
      }
    }
    table.push((codon, found))
  }
  SkbioGeneticCode::{ id: 2, name: "Vertebrate Mitochondrial", table, }
}

///|
/// Translate a DNA sequence using a genetic code table.
/// Returns amino acid sequence (stops as '*').
pub fn skbio_translate_seq(
  dna : String,
  code~ : SkbioGeneticCode = skbio_genetic_code_standard(),
) -> String {
  let sb = StringBuilder()
  let len = dna.length()
  let mut i = 0
  while i + 2 < len {
    let codon = dna[i:i + 3].to_owned().to_upper()
    let mut aa = "X"
    for j = 0; j < code.table.length(); j = j + 1 {
      if code.table[j].0 == codon {
        aa = code.table[j].1
      }
    }
    let _ = sb.write(aa)
    i = i + 3
  }
  sb.to_string()
}

///|
/// RSCU result: relative synonymous codon usage per codon.
pub struct SkbioCodonUsage {
  codon : String
  aa : String
  count : Int32
  rscu : Double // relative synonymous codon usage
}

///|
/// Compute codon usage statistics from a DNA coding sequence.
/// Returns per-codon: count, RSCU.
pub fn skbio_codon_usage(
  dna : String,
  code~ : SkbioGeneticCode = skbio_genetic_code_standard(),
) -> Array[SkbioCodonUsage] {
  // Count codons
  let counts : Array[Int32] = Array(capacity=64)
  for i = 0; i < 64; i = i + 1 {
    counts.push(0)
  }
  let len = dna.length()
  let mut pos = 0
  while pos + 2 < len {
    let codon = dna[pos:pos + 3].to_owned().to_upper()
    for j = 0; j < code.table.length(); j = j + 1 {
      if code.table[j].0 == codon {
        counts[j] = counts[j] + 1
      }
    }
    pos = pos + 3
  }
  // Compute total count per amino acid
  let aa_totals : Array[Int32] = Array(capacity=64)
  for i = 0; i < 64; i = i + 1 {
    aa_totals.push(0)
  }
  for i = 0; i < 64; i = i + 1 {
    let aa = code.table[i].1
    for j = 0; j < 64; j = j + 1 {
      if code.table[j].1 == aa {
        aa_totals[i] = aa_totals[i] + counts[j]
      }
    }
  }
  // Compute RSCU and build results
  let results : Array[SkbioCodonUsage] = Array(capacity=64)
  for i = 0; i < 64; i = i + 1 {
    let (codon, aa) = code.table[i]
    let rscu = if aa_totals[i] > 0 {
      counts[i].to_double() *
      skbio_count_synonymous(aa, code).to_double() /
      aa_totals[i].to_double()
    } else {
      0.0
    }
    results.push(SkbioCodonUsage::{ codon, aa, count: counts[i], rscu, })
  }
  results
}

///|
fn skbio_count_synonymous(aa : String, code : SkbioGeneticCode) -> Int32 {
  let mut count = 0
  for i = 0; i < code.table.length(); i = i + 1 {
    if code.table[i].1 == aa {
      count = count + 1
    }
  }
  count
}

///|
/// Codon Adaptation Index (CAI).
/// Given a coding sequence and a reference set of RSCU values,
/// CAI = geometric mean of (RSCU / RSCU_max) for each codon.
pub fn skbio_cai(
  dna : String,
  reference_rscu : Array[Double], // RSCU values for all 64 codons
  code~ : SkbioGeneticCode = skbio_genetic_code_standard(),
) -> Double {
  let len = dna.length()
  let mut sum_log = 0.0
  let mut n_codons = 0
  let mut pos = 0
  while pos + 2 < len {
    let codon = dna[pos:pos + 3].to_owned().to_upper()
    // Find codon index
    let mut idx = -1
    for j = 0; j < code.table.length(); j = j + 1 {
      if code.table[j].0 == codon {
        idx = j
      }
    }
    if idx >= 0 && idx < reference_rscu.length() {
      let aa = code.table[idx].1
      if aa != "*" && aa != "X" && reference_rscu[idx] > 0.0 {
        // Find max RSCU for this amino acid
        let mut max_rscu = 0.0
        for j = 0; j < code.table.length(); j = j + 1 {
          if code.table[j].1 == aa && j < reference_rscu.length() {
            if reference_rscu[j] > max_rscu {
              max_rscu = reference_rscu[j]
            }
          }
        }
        if max_rscu > 0.0 {
          sum_log = sum_log + @math.ln(reference_rscu[idx] / max_rscu)
          n_codons = n_codons + 1
        }
      }
    }
    pos = pos + 3
  }
  if n_codons == 0 {
    return 0.0
  }
  @math.pow(2.718281828459045, sum_log / n_codons.to_double())
}

///|
/// Effective Number of Codons (ENC).
/// Measures codon usage bias: 20 (max bias) to 61 (no bias).
pub fn skbio_enc(
  dna : String,
  code~ : SkbioGeneticCode = skbio_genetic_code_standard(),
) -> Double {
  let usage = skbio_codon_usage(dna, code~)
  // Group by amino acid and compute F values
  let aa_list : Array[String] = Array(capacity=0)
  for i = 0; i < usage.length(); i = i + 1 {
    let aa = usage[i].aa
    let mut found = false
    for j = 0; j < aa_list.length(); j = j + 1 {
      if aa_list[j] == aa {
        found = true
      }
    }
    if !found && aa != "*" {
      aa_list.push(aa)
    }
  }
  let mut enc = 0.0
  for i = 0; i < aa_list.length(); i = i + 1 {
    let aa = aa_list[i]
    // Sum of squared frequencies for this AA
    let mut total = 0
    let mut sum_sq = 0.0
    for j = 0; j < usage.length(); j = j + 1 {
      if usage[j].aa == aa {
        total = total + usage[j].count
      }
    }
    if total > 0 {
      for j = 0; j < usage.length(); j = j + 1 {
        if usage[j].aa == aa {
          let p = usage[j].count.to_double() / total.to_double()
          sum_sq = sum_sq + p * p
        }
      }
      enc = enc + 1.0 / sum_sq
    } else {
      enc = enc + 1.0
    }
  }
  enc
}
