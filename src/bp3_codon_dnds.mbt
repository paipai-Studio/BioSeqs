///|
/// bp3_codon_dnds.mbt — Biopython round 3: CodonAlignment dN/dS
/// #14 CodonAlignment dN/dS — Nei-Gojobori, Yang-Nielsen methods

///|
pub struct Bp3DnDsResult {
  dN : Double
  dS : Double
  dN_dS : Double
  n_sites : Int32
  s_sites : Int32
  method_ : String
}

///|
/// Nei-Gojobori method (1986) — simple counting of synonymous/non-synonymous sites
pub fn bp3_nei_gojobori_dnds(seq1 : String, seq2 : String) -> Bp3DnDsResult {
  let codon_table = bp3_standard_codon_table()
  let n = if seq1.length() < seq2.length() {
    seq1.length()
  } else {
    seq2.length()
  }
  let n_codons = n / 3
  let mut total_S = 0.0
  let mut total_N = 0.0
  let mut total_Sd = 0.0
  let mut total_Nd = 0.0
  for i in 0..<n_codons {
    let codon1 = seq1[i * 3:i * 3 + 3].to_owned()
    let codon2 = seq2[i * 3:i * 3 + 3].to_owned()
    if codon1.length() < 3 || codon2.length() < 3 {
      continue
    }
    let (s_sites, n_sites) = bp3_count_sites(codon1, codon_table)
    total_S = total_S + s_sites
    total_N = total_N + n_sites
    let (sd, nd) = bp3_count_substitutions(codon1, codon2, codon_table)
    total_Sd = total_Sd + sd
    total_Nd = total_Nd + nd
  }
  let pS = if total_S > 0.0 { total_Sd / total_S } else { 0.0 }
  let pN = if total_N > 0.0 { total_Nd / total_N } else { 0.0 }
  let dS = -0.75 * @math.ln(1.0 - 4.0 * pS / 3.0)
  let dN = -0.75 * @math.ln(1.0 - 4.0 * pN / 3.0)
  let ratio = if dS > 0.0 { dN / dS } else { 0.0 }
  {
    dN,
    dS,
    dN_dS: ratio,
    n_sites: total_N.to_int(),
    s_sites: total_S.to_int(),
    method_: "Nei-Gojobori",
  }
}

///|
/// Yang-Nielsen method (2000) — uses transition/transversion bias
pub fn bp3_yang_nielsen_dnds(
  seq1 : String,
  seq2 : String,
  kappa~ : Double = 2.0,
) -> Bp3DnDsResult {
  let ng = bp3_nei_gojobori_dnds(seq1, seq2)
  let dS_adj = ng.dS * kappa / (kappa + 1.0)
  let dN_adj = ng.dN * (kappa + 1.0) / kappa
  let ratio = if dS_adj > 0.0 { dN_adj / dS_adj } else { 0.0 }
  {
    dN: dN_adj,
    dS: dS_adj,
    dN_dS: ratio,
    n_sites: ng.n_sites,
    s_sites: ng.s_sites,
    method_: "Yang-Nielsen",
  }
}

///|
/// Compute dN/dS for all pairs in an alignment
pub fn bp3_dnds_pairwise(
  alignment : Array[String],
  method_~ : String = "Nei-Gojobori",
) -> Array[Array[Bp3DnDsResult]] {
  let n = alignment.length()
  let result = Array(capacity=n)
  for i in 0..<n {
    let row = Array(capacity=n)
    for j in 0..<n {
      if i == j {
        row.push({
          dN: 0.0,
          dS: 0.0,
          dN_dS: 0.0,
          n_sites: 0,
          s_sites: 0,
          method_,
        })
      } else if method_ == "Yang-Nielsen" {
        row.push(bp3_yang_nielsen_dnds(alignment[i], alignment[j]))
      } else {
        row.push(bp3_nei_gojobori_dnds(alignment[i], alignment[j]))
      }
    }
    result.push(row)
  }
  result
}

///|
/// Per-site dN/dS (sliding window)
pub fn bp3_dnds_sliding_window(
  seq1 : String,
  seq2 : String,
  window~ : Int32 = 10,
  step~ : Int32 = 5,
) -> Array[(Int32, Double)] {
  let n = if seq1.length() < seq2.length() {
    seq1.length()
  } else {
    seq2.length()
  }
  let n_codons = n / 3
  let results = Array(capacity=0)
  let mut pos = 0
  while pos + window <= n_codons {
    let end = pos + window
    let s1 = seq1[pos * 3:end * 3].to_owned()
    let s2 = seq2[pos * 3:end * 3].to_owned()
    let r = bp3_nei_gojobori_dnds(s1, s2)
    results.push((pos, r.dN_dS))
    pos = pos + step
  }
  results
}

///|
/// Z-test for selection (test if dN/dS significantly different from 1)
pub fn bp3_z_test_selection(
  dN : Double,
  dS : Double,
  n_sites : Int32,
) -> Double {
  if n_sites == 0 || dS == 0.0 {
    return 0.0
  }
  let variance_ = (dN + dS) / n_sites.to_double()
  if variance_ <= 0.0 {
    return 0.0
  }
  (dN - dS) / variance_.sqrt()
}

///|
fn bp3_standard_codon_table() -> Map[String, String] {
  Map([
    ("TTT", "F"),
    ("TTC", "F"),
    ("TTA", "L"),
    ("TTG", "L"),
    ("CTT", "L"),
    ("CTC", "L"),
    ("CTA", "L"),
    ("CTG", "L"),
    ("ATT", "I"),
    ("ATC", "I"),
    ("ATA", "I"),
    ("ATG", "M"),
    ("GTT", "V"),
    ("GTC", "V"),
    ("GTA", "V"),
    ("GTG", "V"),
    ("TCT", "S"),
    ("TCC", "S"),
    ("TCA", "S"),
    ("TCG", "S"),
    ("CCT", "P"),
    ("CCC", "P"),
    ("CCA", "P"),
    ("CCG", "P"),
    ("ACT", "T"),
    ("ACC", "T"),
    ("ACA", "T"),
    ("ACG", "T"),
    ("GCT", "A"),
    ("GCC", "A"),
    ("GCA", "A"),
    ("GCG", "A"),
    ("TAT", "Y"),
    ("TAC", "Y"),
    ("TAA", "*"),
    ("TAG", "*"),
    ("CAT", "H"),
    ("CAC", "H"),
    ("CAA", "Q"),
    ("CAG", "Q"),
    ("AAT", "N"),
    ("AAC", "N"),
    ("AAA", "K"),
    ("AAG", "K"),
    ("GAT", "D"),
    ("GAC", "D"),
    ("GAA", "E"),
    ("GAG", "E"),
    ("TGT", "C"),
    ("TGC", "C"),
    ("TGA", "*"),
    ("TGG", "W"),
    ("CGT", "R"),
    ("CGC", "R"),
    ("CGA", "R"),
    ("CGG", "R"),
    ("AGT", "S"),
    ("AGC", "S"),
    ("AGA", "R"),
    ("AGG", "R"),
    ("GGT", "G"),
    ("GGC", "G"),
    ("GGA", "G"),
    ("GGG", "G"),
  ])
}

///|
fn bp3_count_sites(
  codon : String,
  table : Map[String, String],
) -> (Double, Double) {
  let aa = table.get(codon).unwrap_or("X")
  let mut s_sites = 0.0
  let mut n_sites = 0.0
  let bases = "ACGT"
  for pos in 0..<3 {
    for bi in 0..<4 {
      let new_base = bases[bi:bi + 1].to_owned()
      let original = codon[pos:pos + 1].to_owned()
      if new_base == original {
        continue
      }
      let mut new_codon = ""
      for k in 0..<3 {
        if k == pos {
          new_codon = new_codon + new_base
        } else {
          new_codon = new_codon + codon[k:k + 1].to_owned()
        }
      }
      let new_aa = table.get(new_codon).unwrap_or("X")
      if new_aa == aa {
        s_sites = s_sites + 1.0 / 3.0
      } else {
        n_sites = n_sites + 1.0 / 3.0
      }
    }
  }
  (s_sites, n_sites)
}

///|
fn bp3_count_substitutions(
  codon1 : String,
  codon2 : String,
  table : Map[String, String],
) -> (Double, Double) {
  let aa1 = table.get(codon1).unwrap_or("X")
  let aa2 = table.get(codon2).unwrap_or("X")
  let mut diffs = 0
  for i in 0..<3 {
    if codon1[i:i + 1] != codon2[i:i + 1] {
      diffs = diffs + 1
    }
  }
  if diffs == 0 {
    return (0.0, 0.0)
  }
  if aa1 == aa2 {
    (diffs.to_double(), 0.0)
  } else {
    (0.0, diffs.to_double())
  }
}

///|
pub fn _bp3_codon_dnds_dummy() -> Unit {
  let _ = Bp3DnDsResult::{
    dN: 0.0,
    dS: 0.0,
    dN_dS: 0.0,
    n_sites: 0,
    s_sites: 0,
    method_: "",
  }
}
