// Scikit-bio third-round: population genetics supplements (items 27-29)
// four_gamete_test, haplotype_diversity, linkage_disequilibrium

///|
pub struct Skbio3FourGameteResult {
  site1 : Int
  site2 : Int
  violates : Bool
  gametes : Array[String]  // observed gamete types
}

///|
pub fn skbio3_four_gamete_test(
  sequences : Array[String],
) -> Array[Skbio3FourGameteResult] {
  // Four-gamete test: for each pair of sites, check if all 4 gametes are present
  // If all 4 gametes present → recombination detected (violates infinite-sites no-recombination)
  let n_seqs = sequences.length()
  if n_seqs == 0 {
    return []
  }
  let n_sites = sequences[0].length()
  let results : Array[Skbio3FourGameteResult] = []
  for i in 0..<n_sites {
    let _next = i + 1
    for j in _next..<n_sites {
      // Collect gamete types at sites (i, j)
      let gametes : Array[String] = []
      for s in 0..<n_seqs {
        let c1 = sequences[s][i]
        let c2 = sequences[s][j]
        let gamete = {
          let sb = StringBuilder()
          let _ = sb.write_char(Int::unsafe_to_char(c1.to_int()))
          let _ = sb.write_char(Int::unsafe_to_char(c2.to_int()))
          sb.to_string()
        }
        let mut found = false
        for g in 0..<gametes.length() {
          if gametes[g] == gamete {
            found = true
            break
          }
        }
        if not(found) {
          gametes.push(gamete)
        }
      }
      let violates = gametes.length() >= 4
      results.push({ site1: i, site2: j, violates, gametes })
    }
  }
  results
}

///|
pub fn skbio3_haplotype_diversity(
  haplotype_counts : Array[Int],
) -> Double {
  // Haplotype diversity: h = n/(n-1) * (1 - sum(f_i^2))
  // where f_i = frequency of haplotype i, n = total samples
  let mut n = 0
  for i in 0..<haplotype_counts.length() {
    n = n + haplotype_counts[i]
  }
  if n <= 1 {
    return 0.0
  }
  let n_d = n.to_double()
  let mut sum_sq = 0.0
  for i in 0..<haplotype_counts.length() {
    let f = haplotype_counts[i].to_double() / n_d
    sum_sq = sum_sq + f * f
  }
  n_d / (n_d - 1.0) * (1.0 - sum_sq)
}

///|
pub struct Skbio3LDResult {
  locus_a : String
  locus_b : String
  d : Double          // raw D = f_AB - f_A * f_B
  d_prime : Double    // normalized D'
  r_squared : Double  // correlation coefficient squared
}

///|
pub fn skbio3_linkage_disequilibrium(
  locus_a_name : String,
  locus_b_name : String,
  count_AB : Int,
  count_Ab : Int,
  count_aB : Int,
  count_ab : Int,
) -> Skbio3LDResult {
  // Linkage disequilibrium between two loci
  // Gamete frequencies: AB, Ab, aB, ab
  let total = count_AB + count_Ab + count_aB + count_ab
  if total == 0 {
    return { locus_a: locus_a_name, locus_b: locus_b_name, d: 0.0, d_prime: 0.0, r_squared: 0.0 }
  }
  let t = total.to_double()
  let f_AB = count_AB.to_double() / t
  let f_Ab = count_Ab.to_double() / t
  let f_aB = count_aB.to_double() / t
  let f_ab = count_ab.to_double() / t
  // Allele frequencies
  let f_A = f_AB + f_Ab
  let f_B = f_AB + f_aB
  let f_a = f_aB + f_ab
  let f_b = f_Ab + f_ab
  // D = f_AB - f_A * f_B
  let d = f_AB - f_A * f_B
  // D' = D / D_max
  // D_max = min(f_A * f_b, f_a * f_B) if D > 0
  // D_max = max(-f_A * f_B, -f_a * f_b) if D < 0
  let d_max = if d > 0.0 {
    let v1 = f_A * f_b
    let v2 = f_a * f_B
    if v1 < v2 { v1 } else { v2 }
  } else if d < 0.0 {
    let v1 = 0.0 - f_A * f_B
    let v2 = 0.0 - f_a * f_b
    if v1 > v2 { v1 } else { v2 }
  } else {
    1.0
  }
  let d_prime = if d_max != 0.0 { d / d_max } else { 0.0 }
  // r^2 = D^2 / (f_A * f_a * f_B * f_b)
  let denom = f_A * f_a * f_B * f_b
  let r_squared = if denom > 0.0 { d * d / denom } else { 0.0 }
  { locus_a: locus_a_name, locus_b: locus_b_name, d, d_prime, r_squared }
}
