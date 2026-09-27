///|
/// Vegan — community ecology helpers (R `vegan` package surface).
///
/// Two high-frequency pieces that the existing alpha-diversity and phyloseq
/// modules do not cover: rarefaction (expected species count under
/// without-replacement subsampling) and `decostand` matrix standardization.

///|
/// Sum of a sample's count vector as an Int (species abundances are whole
/// numbers stored as Double).
fn vegan_total(counts : Array[Double]) -> Int {
  let mut s = 0
  for c in counts {
    s = s + c.to_int()
  }
  s
}

///|
/// Expected number of species observed when `k` individuals are drawn without
/// replacement from a sample. This is `vegan::rarefy`: for each species with
/// `n` individuals the probability of being missed is
/// `prod_{j=0}^{k-1} (N-n-j)/(N-j)` (a stable product ratio, no combinatorial
/// overflow), so the expected richness is `sum_i (1 - miss_i)`.
pub fn vegan_rarefy(counts : Array[Double], k : Int) -> Double {
  let n_total = vegan_total(counts)
  if n_total <= 0 {
    return 0.0
  }
  let sample = if k > n_total { n_total } else { k }
  let mut expected = 0.0
  for c in counts {
    let n_i = c.to_int()
    if n_i <= 0 {
      continue
    }
    if sample > n_total - n_i {
      // 必然抽到该物种 → miss 概率 0，贡献 1
      expected = expected + 1.0
    } else {
      let mut miss = 1.0
      for j = 0; j < sample; j = j + 1 {
        miss = miss *
          (n_total - n_i - j).to_double() / (n_total - j).to_double()
      }
      expected = expected + 1.0 - miss
    }
  }
  expected
}

///|
/// Rarefaction curve: expected species richness at increasing depths
/// (`vegan::rarecurve`). Returns `(depth, expected_species)` pairs starting at
/// `step` and ending at full depth (where richness equals observed richness).
pub fn vegan_rarefaction_curve(
  counts : Array[Double],
  step? : Int = 1,
) -> Array[(Int, Double)] {
  let n_total = vegan_total(counts)
  let result : Array[(Int, Double)] = Array(capacity=0)
  if n_total <= 0 {
    return result
  }
  let st = if step < 1 { 1 } else { step }
  let mut k = st
  while k < n_total {
    result.push((k, vegan_rarefy(counts, k)))
    k = k + st
  }
  let observed = vegan_rarefy(counts, n_total)
  result.push((n_total, observed))
  result
}

///|
/// Standardize a community matrix (rows = samples, columns = species) with one
/// of `vegan::decostand`'s methods.
///
/// Supported methods:
/// - `"total"`      — relative abundance (divide by sample total)
/// - `"max"`        — divide by sample maximum
/// - `"hellinger"`  — square-root of relative abundance
/// - `"pa"`         — presence/absence (0/1)
/// - `"normalize"`  — divide by sample Euclidean norm
/// - `"frequency"`  — alias of `"total"`
pub fn vegan_decostand(
  matrix : Array[Array[Double]],
  method : String,
) -> Array[Array[Double]] {
  let out : Array[Array[Double]] = Array(capacity=0)
  for row in matrix {
    let n = row.length()
    let total : Double = vegan_row_sum(row)
    let max_v = vegan_row_max(row)
    let mut norm2 = 0.0
    for v in row {
      norm2 = norm2 + v * v
    }
    let norm = norm2.sqrt()

    let new_row : Array[Double] = Array(capacity=0)
    for j = 0; j < n; j = j + 1 {
      let v = row[j]
      let v_new = if method == "pa" {
        if v > 0.0 { 1.0 } else { 0.0 }
      } else if method == "total" || method == "frequency" {
        if total > 0.0 { v / total } else { 0.0 }
      } else if method == "max" {
        if max_v > 0.0 { v / max_v } else { 0.0 }
      } else if method == "hellinger" {
        if total > 0.0 && v > 0.0 { (v / total).sqrt() } else { 0.0 }
      } else if method == "normalize" {
        if norm > 0.0 { v / norm } else { 0.0 }
      } else {
        v
      }
      new_row.push(v_new)
    }
    out.push(new_row)
  }
  out
}

///|
fn vegan_row_sum(row : Array[Double]) -> Double {
  let mut s = 0.0
  for v in row {
    s = s + v
  }
  s
}

///|
fn vegan_row_max(row : Array[Double]) -> Double {
  let mut m = 0.0
  for v in row {
    if v > m {
      m = v
    }
  }
  m
}
