///|
/// Bio.codonalign - Codon alignment and selection pressure analysis.
/// Inspired by Biopython's Bio.codonalign module and Ka/Ks calculators.

///|
/// Codon substitution type.
pub enum CodonSubstitution {
  Synonymous // dS - no amino acid change
  NonSynonymous // dN - amino acid change
  Unknown // ambiguous or stop codon
} derive(Eq, Show)

///|
/// Construct Synonymous CodonSubstitution.
pub fn codon_sub_synonymous() -> CodonSubstitution {
  Synonymous
}

///|
/// Construct NonSynonymous CodonSubstitution.
pub fn codon_sub_nonsynonymous() -> CodonSubstitution {
  NonSynonymous
}

///|
/// Construct Unknown CodonSubstitution.
pub fn codon_sub_unknown() -> CodonSubstitution {
  Unknown
}

///|
/// Codon substitution counts.
pub struct CodonSubCounts {
  synonymous : Int
  nonsynonymous : Int
  total : Int
}

///|
/// Construct CodonSubCounts.
pub fn CodonSubCounts::new(
  synonymous~ : Int,
  nonsynonymous~ : Int,
  total~ : Int,
) -> CodonSubCounts {
  CodonSubCounts::{ synonymous, nonsynonymous, total }
}

///|
/// dN/dS ratio result.
pub struct DnDsResult {
  dn : Double // non-synonymous substitution rate
  ds : Double // synonymous substitution rate
  dnds : Double // dn/ds ratio
  dn_std : Double? // standard deviation of dn
  ds_std : Double? // standard deviation of ds
}

///|
/// Construct DnDsResult.
pub fn DnDsResult::new(
  dn~ : Double,
  ds~ : Double,
  dnds~ : Double,
  dn_std? : Double? = None,
  ds_std? : Double? = None,
) -> DnDsResult {
  DnDsResult::{ dn, ds, dnds, dn_std, ds_std }
}

///|
/// Translate a codon to amino acid (NCBI table 1).
pub fn codon_to_aa(codon : String) -> String {
  let c = codon
  if c == "AAA" || c == "aaa" {
    "K"
  } else if c == "AAC" || c == "aac" {
    "N"
  } else if c == "AAG" || c == "aag" {
    "K"
  } else if c == "AAT" || c == "aat" {
    "N"
  } else if c == "ACA" || c == "aca" {
    "T"
  } else if c == "ACC" || c == "acc" {
    "T"
  } else if c == "ACG" || c == "acg" {
    "T"
  } else if c == "ACT" || c == "act" {
    "T"
  } else if c == "AGA" || c == "aga" {
    "R"
  } else if c == "AGC" || c == "agc" {
    "S"
  } else if c == "AGG" || c == "agg" {
    "R"
  } else if c == "AGT" || c == "agt" {
    "S"
  } else if c == "ATA" || c == "ata" {
    "I"
  } else if c == "ATC" || c == "atc" {
    "I"
  } else if c == "ATG" || c == "atg" {
    "M"
  } else if c == "ATT" || c == "att" {
    "I"
  } else if c == "CAA" || c == "caa" {
    "Q"
  } else if c == "CAC" || c == "cac" {
    "H"
  } else if c == "CAG" || c == "cag" {
    "Q"
  } else if c == "CAT" || c == "cat" {
    "H"
  } else if c == "CCA" || c == "cca" {
    "P"
  } else if c == "CCC" || c == "ccc" {
    "P"
  } else if c == "CCG" || c == "ccg" {
    "P"
  } else if c == "CCT" || c == "cct" {
    "P"
  } else if c == "CGA" || c == "cga" {
    "R"
  } else if c == "CGC" || c == "cgc" {
    "R"
  } else if c == "CGG" || c == "cgg" {
    "R"
  } else if c == "CGT" || c == "cgt" {
    "R"
  } else if c == "CTA" || c == "cta" {
    "L"
  } else if c == "CTC" || c == "ctc" {
    "L"
  } else if c == "CTG" || c == "ctg" {
    "L"
  } else if c == "CTT" || c == "ctt" {
    "L"
  } else if c == "GAA" || c == "gaa" {
    "E"
  } else if c == "GAC" || c == "gac" {
    "D"
  } else if c == "GAG" || c == "gag" {
    "E"
  } else if c == "GAT" || c == "gat" {
    "D"
  } else if c == "GCA" || c == "gca" {
    "A"
  } else if c == "GCC" || c == "gcc" {
    "A"
  } else if c == "GCG" || c == "gcg" {
    "A"
  } else if c == "GCT" || c == "gct" {
    "A"
  } else if c == "GGA" || c == "gga" {
    "G"
  } else if c == "GGC" || c == "ggc" {
    "G"
  } else if c == "GGG" || c == "ggg" {
    "G"
  } else if c == "GGT" || c == "ggt" {
    "G"
  } else if c == "GTA" || c == "gta" {
    "V"
  } else if c == "GTC" || c == "gtc" {
    "V"
  } else if c == "GTG" || c == "gtg" {
    "V"
  } else if c == "GTT" || c == "gtt" {
    "V"
  } else if c == "TAA" || c == "taa" {
    "*"
  } else if c == "TAC" || c == "tac" {
    "Y"
  } else if c == "TAG" || c == "tag" {
    "*"
  } else if c == "TAT" || c == "tat" {
    "Y"
  } else if c == "TCA" || c == "tca" {
    "S"
  } else if c == "TCC" || c == "tcc" {
    "S"
  } else if c == "TCG" || c == "tcg" {
    "S"
  } else if c == "TCT" || c == "tct" {
    "S"
  } else if c == "TGA" || c == "tga" {
    "*"
  } else if c == "TGC" || c == "tgc" {
    "C"
  } else if c == "TGG" || c == "tgg" {
    "W"
  } else if c == "TGT" || c == "tgt" {
    "C"
  } else if c == "TTA" || c == "tta" {
    "L"
  } else if c == "TTC" || c == "ttc" {
    "F"
  } else if c == "TTG" || c == "ttg" {
    "L"
  } else if c == "TTT" || c == "ttt" {
    "F"
  } else {
    "?"
  }
}

///|
/// Check if a codon is a stop codon.
pub fn is_stop_codon(codon : String) -> Bool {
  let aa = codon_to_aa(codon)
  aa == "*" || aa == "?"
}

///|
/// Classify a codon substitution.
pub fn classify_codon_substitution(
  codon1 : String,
  codon2 : String,
) -> CodonSubstitution {
  if codon1.length() != 3 || codon2.length() != 3 {
    return Unknown
  }
  let aa1 = codon_to_aa(codon1)
  let aa2 = codon_to_aa(codon2)
  if aa1 == "*" || aa2 == "*" || aa1 == "?" || aa2 == "?" {
    Unknown
  } else if aa1 == aa2 {
    Synonymous
  } else {
    NonSynonymous
  }
}

///|
/// Count synonymous and non-synonymous substitutions in aligned sequences.
pub fn count_codon_substitutions(
  seq1 : String,
  seq2 : String,
) -> CodonSubCounts {
  let n = if seq1.length() < seq2.length() {
    seq1.length()
  } else {
    seq2.length()
  }
  let codon_count = n / 3
  let mut synonymous = 0
  let mut nonsynonymous = 0
  let mut i = 0
  while i < codon_count {
    let start = i * 3
    let codon1 = seq1[start:start + 3].to_owned()
    let codon2 = seq2[start:start + 3].to_owned()
    match classify_codon_substitution(codon1, codon2) {
      Synonymous => synonymous = synonymous + 1
      NonSynonymous => nonsynonymous = nonsynonymous + 1
      Unknown => ()
    }
    i = i + 1
  }
  CodonSubCounts::new(
    synonymous~,
    nonsynonymous~,
    total=synonymous + nonsynonymous,
  )
}

///|
/// Calculate the number of synonymous and non-synonymous sites (Nei-Gojobori).
/// Returns (synonymous_sites, nonsynonymous_sites).
pub fn count_codon_sites(codon : String) -> (Double, Double) {
  if codon.length() != 3 || is_stop_codon(codon) {
    return (0.0, 0.0)
  }
  let mut s_sites = 0.0
  let mut ns_sites = 0.0
  let target_aa = codon_to_aa(codon)
  let possible_codons = get_possible_codons(codon)
  let n_possible = possible_codons.length().to_double()
  if n_possible == 0.0 {
    return (0.0, 0.0)
  }
  for c in possible_codons {
    let aa = codon_to_aa(c)
    if aa == target_aa {
      s_sites = s_sites + 1.0
    } else {
      ns_sites = ns_sites + 1.0
    }
  }
  s_sites = s_sites / n_possible * 3.0
  ns_sites = ns_sites / n_possible * 3.0
  (s_sites, ns_sites)
}

///|
/// Get all possible codons reachable via single nucleotide changes.
fn get_possible_codons(codon : String) -> Array[String] {
  let result : Array[String] = Array::new()
  let bases = ['A', 'T', 'C', 'G']
  let c = codon
  let mut pos = 0
  while pos < 3 {
    let mut base_idx = 0
    while base_idx < bases.length() {
      let base = bases[base_idx]
      let current_str = c[pos:pos + 1].to_owned()
      if base.to_string() != current_str {
        let mut new_codon = ""
        let mut i = 0
        while i < 3 {
          if i == pos {
            new_codon = new_codon + base.to_string()
          } else {
            new_codon = new_codon + c[i:i + 1].to_owned()
          }
          i = i + 1
        }
        result.push(new_codon)
      }
      base_idx = base_idx + 1
    }
    pos = pos + 1
  }
  result
}

///|
/// Calculate total synonymous and non-synonymous sites for a sequence.
pub fn total_codon_sites(seq : String) -> (Double, Double) {
  let n = seq.length() / 3
  let mut s_total = 0.0
  let mut ns_total = 0.0
  let mut i = 0
  while i < n {
    let start = i * 3
    let codon = seq[start:start + 3].to_owned()
    let (s, ns) = count_codon_sites(codon)
    s_total = s_total + s
    ns_total = ns_total + ns
    i = i + 1
  }
  (s_total, ns_total)
}

///|
/// Jukes-Cantor correction for multiple hits.
pub fn jukes_cantor_correction(p : Double) -> Double {
  if p >= 1.0 {
    return 0.0
  }
  let corrected = -3.0 / 4.0 * @math.ln(1.0 - 4.0 * p / 3.0)
  if corrected.is_nan() {
    0.0
  } else {
    corrected
  }
}

///|
/// Calculate dN/dS ratio using Nei-Gojobori method.
pub fn calc_dnds_nei_gojobori(seq1 : String, seq2 : String) -> DnDsResult {
  let n_codons = if seq1.length() < seq2.length() {
    seq1.length() / 3
  } else {
    seq2.length() / 3
  }
  let mut s_sites = 0.0
  let mut ns_sites = 0.0
  let mut s_subs = 0
  let mut ns_subs = 0
  let mut i = 0
  while i < n_codons {
    let start = i * 3
    let codon1 = seq1[start:start + 3].to_owned()
    let codon2 = seq2[start:start + 3].to_owned()
    if !is_stop_codon(codon1) && !is_stop_codon(codon2) {
      let (s1, ns1) = count_codon_sites(codon1)
      let (s2, ns2) = count_codon_sites(codon2)
      s_sites = s_sites + (s1 + s2) / 2.0
      ns_sites = ns_sites + (ns1 + ns2) / 2.0
      match classify_codon_substitution(codon1, codon2) {
        Synonymous => s_subs = s_subs + 1
        NonSynonymous => ns_subs = ns_subs + 1
        Unknown => ()
      }
    }
    i = i + 1
  }
  let p_s = if s_sites > 0.0 { s_subs.to_double() / s_sites } else { 0.0 }
  let p_ns = if ns_sites > 0.0 { ns_subs.to_double() / ns_sites } else { 0.0 }
  let d_s = jukes_cantor_correction(p_s)
  let d_n = jukes_cantor_correction(p_ns)
  let dnds = if d_s == 0.0 {
    if d_n > 0.0 {
      999.0
    } else {
      0.0
    }
  } else {
    d_n / d_s
  }
  DnDsResult::new(dn=d_n, ds=d_s, dnds~)
}

///|
/// Calculate dN/dS for all pairs of sequences.
pub fn pairwise_dnds(
  sequences : Array[String],
) -> Array[(String, String, DnDsResult)] {
  let result : Array[(String, String, DnDsResult)] = Array::new()
  let n = sequences.length()
  let mut i = 0
  while i < n {
    let mut j = i + 1
    while j < n {
      let dnds = calc_dnds_nei_gojobori(sequences[i], sequences[j])
      result.push(("seq" + i.to_string(), "seq" + j.to_string(), dnds))
      j = j + 1
    }
    i = i + 1
  }
  result
}

///|
/// Calculate codon usage bias for a sequence.
/// Returns relative synonymous codon usage (RSCU) values.
pub fn codon_usage_bias(seq : String) -> Map[String, Double] {
  let result : Map[String, Double] = Map([], capacity=64)
  let n = seq.length() / 3
  let aa_counts : Map[String, Int] = Map([], capacity=20)
  let codon_counts : Map[String, Int] = Map([], capacity=64)
  let mut i = 0
  while i < n {
    let start = i * 3
    let codon = seq[start:start + 3].to_owned()
    if !is_stop_codon(codon) {
      match codon_counts.get(codon) {
        Some(count) => codon_counts[codon] = count + 1
        None => codon_counts[codon] = 1
      }
      let aa = codon_to_aa(codon)
      match aa_counts.get(aa) {
        Some(count) => aa_counts[aa] = count + 1
        None => aa_counts[aa] = 1
      }
    }
    i = i + 1
  }
  for kv in codon_counts {
    let codon = kv.0
    let count = kv.1
    let aa = codon_to_aa(codon)
    let aa_total = aa_counts.get(aa).unwrap_or(0)
    let synonymous_codons = get_synonymous_codons(aa)
    let n_syn = synonymous_codons.length().to_double()
    if aa_total > 0 && n_syn > 0.0 {
      let expected = aa_total.to_double() / n_syn
      result[codon] = count.to_double() / expected
    }
  }
  result
}

///|
/// Get all synonymous codons for a given amino acid.
fn get_synonymous_codons(aa : String) -> Array[String] {
  let result : Array[String] = Array::new()
  let codons = [
    "AAA", "AAC", "AAG", "AAT", "ACA", "ACC", "ACG", "ACT", "AGA", "AGC", "AGG",
    "AGT", "ATA", "ATC", "ATG", "ATT", "CAA", "CAC", "CAG", "CAT", "CCA", "CCC",
    "CCG", "CCT", "CGA", "CGC", "CGG", "CGT", "CTA", "CTC", "CTG", "CTT", "GAA",
    "GAC", "GAG", "GAT", "GCA", "GCC", "GCG", "GCT", "GGA", "GGC", "GGG", "GGT",
    "GTA", "GTC", "GTG", "GTT", "TAA", "TAC", "TAG", "TAT", "TCA", "TCC", "TCG",
    "TCT", "TGA", "TGC", "TGG", "TGT", "TTA", "TTC", "TTG", "TTT",
  ]
  for codon in codons {
    if codon_to_aa(codon) == aa {
      result.push(codon)
    }
  }
  result
}

///|
/// Calculate the effective number of codons (ENC).
/// Lower ENC indicates stronger codon usage bias.
pub fn effective_number_of_codons(seq : String) -> Double {
  let n = seq.length() / 3
  if n == 0 {
    return 61.0
  }
  let aa_groups = [
    ["F", "L", "I", "M", "V"],
    ["S", "P", "T", "A", "Y"],
    ["H", "Q", "N", "K", "D", "E"],
    ["C", "W", "R", "G"],
  ]
  let mut enc = 0.0
  let mut n_groups = 0
  for group in aa_groups {
    let mut total_p2 = 0.0
    let mut total_n = 0
    for aa in group {
      let codons = get_synonymous_codons(aa)
      let n_syn = codons.length()
      if n_syn > 1 {
        let codon_freqs : Array[Double] = Array::new()
        let mut total = 0
        for codon in codons {
          let count = count_codon_in_sequence(seq, codon)
          codon_freqs.push(count.to_double())
          total = total + count
        }
        if total > 0 {
          let mut p2 = 0.0
          for freq in codon_freqs {
            let p = freq / total.to_double()
            p2 = p2 + p * p
          }
          total_p2 = total_p2 + 1.0 / (n_syn.to_double() * p2)
          total_n = total_n + 1
        }
      }
    }
    if total_n > 0 {
      enc = enc + total_p2 / total_n.to_double()
      n_groups = n_groups + 1
    }
  }
  if n_groups == 0 {
    61.0
  } else {
    enc * n_groups.to_double()
  }
}

///|
/// Count occurrences of a specific codon in a sequence.
fn count_codon_in_sequence(seq : String, codon : String) -> Int {
  let n = seq.length() / 3
  let mut count = 0
  let mut i = 0
  while i < n {
    let start = i * 3
    let c = seq[start:start + 3].to_owned()
    if c == codon {
      count = count + 1
    }
    i = i + 1
  }
  count
}

///|
/// Create a demo codon alignment.
pub fn create_demo_codon_alignment() -> Array[String] {
  ["ATGTCTGAAATG", "ATGTCTGAAATG", "ATGTCTGAGATG", "ATGTCCGAAATG"]
}

///|
/// Create a demo alignment with positive selection.
pub fn create_demo_positive_selection() -> Array[String] {
  [
    "ATGTCTGAATAA", "ATGTCTGAATAA", "ATGTCTGAGTAA", "ATGTCTGAGTAA", "ATGTCCGAGTAA",
    "ATGTCCGAGTAA",
  ]
}
