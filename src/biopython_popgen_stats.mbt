/// biopython_popgen_stats.mbt — Population genetics statistics
/// Inspired by Biopython Bio.PopGen: nucleotide diversity, Tajima's D, segregating sites

///|
/// Count segregating (polymorphic) sites in an alignment.
/// A site is segregating if at least two sequences differ at that position.
pub fn biopython_segregating_sites(sequences : Array[String]) -> Int32 {
  if sequences.length() < 2 {
    return 0
  }
  let aln_len = sequences[0].length()
  let mut count = 0
  for pos = 0; pos < aln_len; pos = pos + 1 {
    let first = sequences[0][pos:pos + 1].to_owned()
    let mut seg = false
    for i = 1; i < sequences.length(); i = i + 1 {
      if sequences[i][pos:pos + 1].to_owned() != first {
        seg = true
        break
      }
    }
    if seg {
      count = count + 1
    }
  }
  count
}

///|
/// Nucleotide diversity (pi): average pairwise differences.
/// pi = sum(d_ij) / (n choose 2)
pub fn biopython_nucleotide_diversity(sequences : Array[String]) -> Double {
  let n = sequences.length()
  if n < 2 {
    return 0.0
  }
  let mut total_diff = 0
  for i = 0; i < n; i = i + 1 {
    for j = i + 1; j < n; j = j + 1 {
      let aln_len = sequences[i].length()
      for pos = 0; pos < aln_len; pos = pos + 1 {
        if sequences[i][pos] != sequences[j][pos] {
          total_diff = total_diff + 1
        }
      }
    }
  }
  let n_pairs = n * (n - 1) / 2
  total_diff.to_double() / n_pairs.to_double()
}

///|
/// Tajima's D statistic.
/// D = (pi - S/a1) / sqrt(e1*S + e2*S*(S-1))
/// where S = segregating sites, a1 = sum(1/i), e1, e2 are corrections.
pub fn biopython_tajima_d(sequences : Array[String]) -> Double {
  let n = sequences.length()
  if n < 2 {
    return 0.0
  }
  let s = biopython_segregating_sites(sequences)
  if s == 0 {
    return 0.0
  }
  let pi = biopython_nucleotide_diversity(sequences)
  // a1 = sum(1/i for i=1..n-1)
  let mut a1 = 0.0
  let mut a2 = 0.0
  for i = 1; i < n; i = i + 1 {
    a1 = a1 + 1.0 / i.to_double()
    a2 = a2 + 1.0 / (i.to_double() * i.to_double())
  }
  let b1 = (n + 1).to_double() / (3.0 * (n - 1).to_double())
  let b2 = 2.0 *
    (n * n + n + 3).to_double() /
    (9.0 * (n - 1).to_double() * n.to_double())
  let c1 = b1 - 1.0 / a1
  let c2 = b2 - (n + 2).to_double() / (a1 * n.to_double()) + a2 / (a1 * a1)
  let e1 = c1 / a1
  let e2 = c2 / (a1 * a1 + a2)
  let numerator = pi - s.to_double() / a1
  let denominator = biopython_popgen_sqrt(
    e1 * s.to_double() + e2 * s.to_double() * (s - 1).to_double(),
  )
  if denominator == 0.0 {
    0.0
  } else {
    numerator / denominator
  }
}

///|
/// Watterson's theta (estimator of 4Nu).
/// theta = S / a1
pub fn biopython_watterson_theta(sequences : Array[String]) -> Double {
  let n = sequences.length()
  if n < 2 {
    return 0.0
  }
  let s = biopython_segregating_sites(sequences)
  let mut a1 = 0.0
  for i = 1; i < n; i = i + 1 {
    a1 = a1 + 1.0 / i.to_double()
  }
  if a1 == 0.0 {
    0.0
  } else {
    s.to_double() / a1
  }
}

///|
/// Fu and Li's D statistic.
/// D = (S - eta/a1) / sqrt(u1*S + u2*S^2)
/// where eta = number of mutations on external branches.
pub fn biopython_fu_li_d(sequences : Array[String]) -> Double {
  let n = sequences.length()
  if n < 3 {
    return 0.0
  }
  let s = biopython_segregating_sites(sequences)
  // Simplified: eta ≈ S (external mutations estimated as S)
  let eta = s.to_double()
  let mut a1 = 0.0
  for i = 1; i < n; i = i + 1 {
    a1 = a1 + 1.0 / i.to_double()
  }
  let numerator = s.to_double() - eta / a1
  // Simplified variance
  let v = s.to_double() *
    (n.to_double() - 2.0) /
    ((n.to_double() - 1.0) * (n.to_double() - 1.0))
  if v == 0.0 {
    0.0
  } else {
    numerator / biopython_popgen_sqrt(v)
  }
}

///|
/// Helper: sqrt
fn biopython_popgen_sqrt(x : Double) -> Double {
  @math.pow(x, 0.5)
}
