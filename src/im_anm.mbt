///|
/// ProDy-style imANM (mixed ANM).
///
/// imANM assigns a different force constant (gamma) to each residue,
/// allowing heterogeneous stiffness across the protein.  The effective
/// spring constant for a pair (i, j) is:
///
///   γ_ij = combine(γ_i, γ_j) / d_ij²
///
/// where `combine` is one of five rules (arithmetic mean, geometric
/// mean, product, minimum, harmonic mean).
///
/// Use cases:
/// - Softer springs in loop / coil regions, stiffer in the core.
/// - Incorporate experimental B-factors: γ_i ∝ 1 / B_i.
/// - Model multi-domain proteins with different rigidities.
///
/// Reference: Yang, Bahar (2005) Structure 13:893-904.

///|
/// Combination rule for per-residue gamma values.
pub(all) enum GammaCombine {
  ArithmeticMean // (γ_i + γ_j) / 2
  GeometricMean // √(γ_i × γ_j)
  Product // γ_i × γ_j
  Minimum // min(γ_i, γ_j)
  HarmonicMean // 2 γ_i γ_j / (γ_i + γ_j)
}

///|
/// Result of a mixed ANM calculation.
pub struct ImANM {
  n : Int
  cutoff : Double
  gammas : Array[Double]
  combine : GammaCombine
  hessian : Array[Array[Double]]
  eigenvalues : Array[Double]
  eigenvectors : Array[Array[Double]]
}

///|
/// Combine two per-residue gamma values according to the given rule.
fn combine_gammas(gi : Double, gj : Double, rule : GammaCombine) -> Double {
  match rule {
    ArithmeticMean => (gi + gj) / 2.0
    GeometricMean => {
      let prod = gi * gj
      if prod < 0.0 {
        0.0
      } else {
        prod.sqrt()
      }
    }
    Product => gi * gj
    Minimum => if gi < gj { gi } else { gj }
    HarmonicMean => {
      let s = gi + gj
      if s < 1.0e-15 {
        0.0
      } else {
        2.0 * gi * gj / s
      }
    }
  }
}

///|
/// Build an imANM from a protein chain and per-residue gamma values.
///
/// `gammas` must have one entry per CA atom in the chain.
/// If the lengths don't match, returns `None`.
pub fn ImANM::new(
  chain : Chain,
  gammas : Array[Double],
  cutoff? : Double = 15.0,
  combine? : GammaCombine = ArithmeticMean,
) -> ImANM? {
  let residues = chain.get_residues()
  let ca_coords : Array[Vector3] = Array(capacity=0)
  for res in residues {
    match res.get_ca() {
      Some(ca) => ca_coords.push(ca.coord)
      None => ()
    }
  }
  let n = ca_coords.length()
  if n < 3 {
    return None
  }
  if gammas.length() != n {
    return None
  }
  let dim = 3 * n
  let h : Array[Array[Double]] = Array::makei(dim, { _ => Array::make(dim, 0.0) })
  let mut i = 0
  while i < n {
    let mut j = i + 1
    while j < n {
      let dx = ca_coords[j].x - ca_coords[i].x
      let dy = ca_coords[j].y - ca_coords[i].y
      let dz = ca_coords[j].z - ca_coords[i].z
      let d2 = dx * dx + dy * dy + dz * dz
      let d = d2.sqrt()
      if d <= cutoff && d > 0.001 {
        let g_ij = combine_gammas(gammas[i], gammas[j], combine) / d2
        let hxx = g_ij * dx * dx
        let hxy = g_ij * dx * dy
        let hxz = g_ij * dx * dz
        let hyy = g_ij * dy * dy
        let hyz = g_ij * dy * dz
        let hzz = g_ij * dz * dz
        let i3 = 3 * i
        let j3 = 3 * j
        h[i3][j3] = h[i3][j3] - hxx
        h[i3][j3 + 1] = h[i3][j3 + 1] - hxy
        h[i3][j3 + 2] = h[i3][j3 + 2] - hxz
        h[i3 + 1][j3] = h[i3 + 1][j3] - hxy
        h[i3 + 1][j3 + 1] = h[i3 + 1][j3 + 1] - hyy
        h[i3 + 1][j3 + 2] = h[i3 + 1][j3 + 2] - hyz
        h[i3 + 2][j3] = h[i3 + 2][j3] - hxz
        h[i3 + 2][j3 + 1] = h[i3 + 2][j3 + 1] - hyz
        h[i3 + 2][j3 + 2] = h[i3 + 2][j3 + 2] - hzz
        h[j3][i3] = h[j3][i3] - hxx
        h[j3][i3 + 1] = h[j3][i3 + 1] - hxy
        h[j3][i3 + 2] = h[j3][i3 + 2] - hxz
        h[j3 + 1][i3] = h[j3 + 1][i3] - hxy
        h[j3 + 1][i3 + 1] = h[j3 + 1][i3 + 1] - hyy
        h[j3 + 1][i3 + 2] = h[j3 + 1][i3 + 2] - hyz
        h[j3 + 2][i3] = h[j3 + 2][i3] - hxz
        h[j3 + 2][i3 + 1] = h[j3 + 2][i3 + 1] - hyz
        h[j3 + 2][i3 + 2] = h[j3 + 2][i3 + 2] - hzz
        h[i3][i3] = h[i3][i3] + hxx
        h[i3][i3 + 1] = h[i3][i3 + 1] + hxy
        h[i3][i3 + 2] = h[i3][i3 + 2] + hxz
        h[i3 + 1][i3] = h[i3 + 1][i3] + hxy
        h[i3 + 1][i3 + 1] = h[i3 + 1][i3 + 1] + hyy
        h[i3 + 1][i3 + 2] = h[i3 + 1][i3 + 2] + hyz
        h[i3 + 2][i3] = h[i3 + 2][i3] + hxz
        h[i3 + 2][i3 + 1] = h[i3 + 2][i3 + 1] + hyz
        h[i3 + 2][i3 + 2] = h[i3 + 2][i3 + 2] + hzz
        h[j3][j3] = h[j3][j3] + hxx
        h[j3][j3 + 1] = h[j3][j3 + 1] + hxy
        h[j3][j3 + 2] = h[j3][j3 + 2] + hxz
        h[j3 + 1][j3] = h[j3 + 1][j3] + hxy
        h[j3 + 1][j3 + 1] = h[j3 + 1][j3 + 1] + hyy
        h[j3 + 1][j3 + 2] = h[j3 + 1][j3 + 2] + hyz
        h[j3 + 2][j3] = h[j3 + 2][j3] + hxz
        h[j3 + 2][j3 + 1] = h[j3 + 2][j3 + 1] + hyz
        h[j3 + 2][j3 + 2] = h[j3 + 2][j3 + 2] + hzz
      }
      j = j + 1
    }
    i = i + 1
  }
  let max_iter = if dim > 100 { dim * 20 } else { 200 }
  let (eigenvalues, eigenvectors) = linalg_jacobi(h, dim, max_iter)
  Some(ImANM::{ n, cutoff, gammas, combine, hessian: h, eigenvalues, eigenvectors })
}

///|
/// Get the Hessian matrix.
pub fn ImANM::hessian_matrix(self : ImANM) -> Array[Array[Double]] {
  self.hessian
}

///|
/// Get sorted eigenvalues (ascending, first 6 ~0).
pub fn ImANM::get_eigenvalues(self : ImANM) -> Array[Double] {
  let dim = 3 * self.n
  let vals = Array::make(dim, 0.0)
  for i = 0; i < dim; i = i + 1 {
    vals[i] = self.eigenvalues[i]
  }
  for i = 1; i < dim; i = i + 1 {
    let key = vals[i]
    let mut j = i - 1
    while j >= 0 && vals[j] > key {
      vals[j + 1] = vals[j]
      j = j - 1
    }
    vals[j + 1] = key
  }
  vals
}

///|
/// Get eigenvectors as columns (sorted to match eigenvalues).
pub fn ImANM::get_eigenvectors(self : ImANM) -> Array[Array[Double]] {
  let dim = 3 * self.n
  let pairs : Array[(Double, Array[Double])] = Array(capacity=dim)
  for i = 0; i < dim; i = i + 1 {
    let vec = Array::make(dim, 0.0)
    for j = 0; j < dim; j = j + 1 {
      vec[j] = self.eigenvectors[j][i]
    }
    pairs.push((self.eigenvalues[i], vec))
  }
  for i = 1; i < dim; i = i + 1 {
    let key = pairs[i]
    let mut j = i - 1
    while j >= 0 && pairs[j].0 > key.0 {
      pairs[j + 1] = pairs[j]
      j = j - 1
    }
    pairs[j + 1] = key
  }
  let result : Array[Array[Double]] = Array(capacity=dim)
  for i = 0; i < dim; i = i + 1 {
    result.push(pairs[i].1)
  }
  result
}

///|
/// Predict B-factors from slow modes (excluding 6 zero modes).
pub fn ImANM::predict_bfactors(self : ImANM) -> Array[Double] {
  let dim = 3 * self.n
  let n = self.n
  let bf = Array::make(n, 0.0)
  let pairs : Array[(Double, Array[Double])] = Array(capacity=dim)
  for i = 0; i < dim; i = i + 1 {
    let vec = Array::make(dim, 0.0)
    for j = 0; j < dim; j = j + 1 {
      vec[j] = self.eigenvectors[j][i]
    }
    pairs.push((self.eigenvalues[i], vec))
  }
  for i = 1; i < dim; i = i + 1 {
    let key = pairs[i]
    let mut j = i - 1
    while j >= 0 && pairs[j].0 > key.0 {
      pairs[j + 1] = pairs[j]
      j = j - 1
    }
    pairs[j + 1] = key
  }
  for k = 6; k < dim; k = k + 1 {
    let lambda = pairs[k].0
    if lambda < 1.0e-10 {
      continue
    }
    let vec = pairs[k].1
    for i = 0; i < n; i = i + 1 {
      let i3 = 3 * i
      bf[i] = bf[i] + (vec[i3] * vec[i3] + vec[i3 + 1] * vec[i3 + 1] + vec[i3 + 2] * vec[i3 + 2]) / lambda
    }
  }
  bf
}

///|
/// Get per-residue mean-square fluctuations.
pub fn ImANM::fluctuations(self : ImANM) -> Array[Double] {
  self.predict_bfactors()
}

///|
/// Get the k slowest non-trivial modes.
pub fn ImANM::slow_modes(self : ImANM, k : Int) -> Array[(Int, Double, Array[Double])] {
  let dim = 3 * self.n
  let pairs : Array[(Double, Array[Double])] = Array(capacity=dim)
  for i = 0; i < dim; i = i + 1 {
    let vec = Array::make(dim, 0.0)
    for j = 0; j < dim; j = j + 1 {
      vec[j] = self.eigenvectors[j][i]
    }
    pairs.push((self.eigenvalues[i], vec))
  }
  for i = 1; i < dim; i = i + 1 {
    let key = pairs[i]
    let mut j = i - 1
    while j >= 0 && pairs[j].0 > key.0 {
      pairs[j + 1] = pairs[j]
      j = j - 1
    }
    pairs[j + 1] = key
  }
  let kk = if k < dim - 6 { k } else { dim - 6 }
  let result : Array[(Int, Double, Array[Double])] = Array(capacity=kk)
  for idx = 0; idx < kk; idx = idx + 1 {
    let mode_idx = idx + 6
    result.push((mode_idx, pairs[mode_idx].0, pairs[mode_idx].1))
  }
  result
}

///|
/// Compute cross-correlation matrix from slow modes.
pub fn ImANM::cross_correlation(self : ImANM) -> Array[Array[Double]] {
  let dim = 3 * self.n
  let n = self.n
  let pairs : Array[(Double, Array[Double])] = Array(capacity=dim)
  for i = 0; i < dim; i = i + 1 {
    let vec = Array::make(dim, 0.0)
    for j = 0; j < dim; j = j + 1 {
      vec[j] = self.eigenvectors[j][i]
    }
    pairs.push((self.eigenvalues[i], vec))
  }
  for i = 1; i < dim; i = i + 1 {
    let key = pairs[i]
    let mut j = i - 1
    while j >= 0 && pairs[j].0 > key.0 {
      pairs[j + 1] = pairs[j]
      j = j - 1
    }
    pairs[j + 1] = key
  }
  let corr : Array[Array[Double]] = Array::makei(n, { _ => Array::make(n, 0.0) })
  for k = 6; k < dim; k = k + 1 {
    let lambda = pairs[k].0
    if lambda < 1.0e-10 {
      continue
    }
    let vec = pairs[k].1
    for i = 0; i < n; i = i + 1 {
      let i3 = 3 * i
      for j = 0; j < n; j = j + 1 {
        let j3 = 3 * j
        corr[i][j] = corr[i][j] + (vec[i3] * vec[j3] + vec[i3 + 1] * vec[j3 + 1] + vec[i3 + 2] * vec[j3 + 2]) / lambda
      }
    }
  }
  let mut max_diag = 0.0
  for i = 0; i < n; i = i + 1 {
    if corr[i][i] > max_diag {
      max_diag = corr[i][i]
    }
  }
  if max_diag > 1.0e-15 {
    for i = 0; i < n; i = i + 1 {
      for j = 0; j < n; j = j + 1 {
        corr[i][j] = corr[i][j] / max_diag
      }
    }
  }
  corr
}

///|
/// Number of nodes (residues).
pub fn ImANM::num_nodes(self : ImANM) -> Int {
  self.n
}

///|
/// Get per-residue gamma values.
pub fn ImANM::get_gammas(self : ImANM) -> Array[Double] {
  self.gammas
}

///|
/// Get the cutoff distance.
pub fn ImANM::get_cutoff(self : ImANM) -> Double {
  self.cutoff
}

///|
/// Get the combination rule.
pub fn ImANM::get_combine_rule(self : ImANM) -> GammaCombine {
  self.combine
}

// =========================================================================
// Helper functions for constructing per-residue gamma arrays
// =========================================================================

///|
/// Assign gamma values by residue regions.
///
/// Each region is (start_idx, end_idx, gamma_value) — inclusive range.
/// Residues not covered by any region get `default_gamma`.
pub fn gammas_by_region(
  n : Int,
  regions : Array[(Int, Int, Double)],
  default_gamma? : Double = 1.0,
) -> Array[Double] {
  let gammas = Array::make(n, default_gamma)
  for region in regions {
    let (start, end, g) = region
    let lo = if start >= 0 { start } else { 0 }
    let hi = if end < n { end } else { n - 1 }
    for i = lo; i <= hi; i = i + 1 {
      gammas[i] = g
    }
  }
  gammas
}

///|
/// Assign gamma values from experimental B-factors.
///
/// γ_i = scale / B_i  (inverse B-factor: flexible residues get weak springs)
///
/// B-factors of zero or negative are replaced by the mean of valid entries.
pub fn gammas_from_bfactors(
  bfactors : Array[Double],
  scale? : Double = 1.0,
) -> Array[Double] {
  let n = bfactors.length()
  if n == 0 {
    return []
  }
  // Compute mean of valid B-factors for fallback
  let mut sum = 0.0
  let mut count = 0
  for i = 0; i < n; i = i + 1 {
    if bfactors[i] > 1.0e-10 {
      sum = sum + bfactors[i]
      count = count + 1
    }
  }
  let fallback = if count > 0 { sum / count.to_double() } else { 1.0 }
  let gammas = Array::make(n, 0.0)
  for i = 0; i < n; i = i + 1 {
    let b = if bfactors[i] > 1.0e-10 { bfactors[i] } else { fallback }
    gammas[i] = scale / b
  }
  gammas
}

///|
/// Assign gamma values by secondary structure assignment.
///
/// `ss` is an array of single-character strings: "H" (helix),
/// "E" (sheet), "C" or other (coil/loop).
///
/// Typical values: helix/sheet are rigid (γ=1.0), coil is soft (γ=0.3).
pub fn gammas_by_secondary_structure(
  n : Int,
  ss : Array[String],
  helix_gamma? : Double = 1.0,
  sheet_gamma? : Double = 1.0,
  coil_gamma? : Double = 0.3,
) -> Array[Double] {
  let gammas = Array::make(n, coil_gamma)
  let limit = if ss.length() < n { ss.length() } else { n }
  for i = 0; i < limit; i = i + 1 {
    let s = ss[i]
    if s == "H" || s == "G" || s == "I" {
      // Helix (H=alpha, G=3-10, I=pi)
      gammas[i] = helix_gamma
    } else if s == "E" || s == "B" {
      // Sheet (E=strand, B=bridge)
      gammas[i] = sheet_gamma
    } else {
      // Coil / loop / turn
      gammas[i] = coil_gamma
    }
  }
  gammas
}

///|
/// Create uniform gamma array (all residues same gamma).
///
/// With uniform gamma and ArithmeticMean combine rule, imANM reduces
/// to standard ANM.
pub fn gammas_uniform(n : Int, gamma : Double) -> Array[Double] {
  Array::make(n, gamma)
}

///|
/// Create gamma array with Gaussian profile centered at a residue.
///
/// γ_i = base + amplitude * exp(-(i - center)² / (2σ²))
///
/// Useful for modeling a locally stiffened or softened region.
pub fn gammas_gaussian(
  n : Int,
  center : Int,
  sigma : Double,
  base? : Double = 1.0,
  amplitude? : Double = 1.0,
) -> Array[Double] {
  let gammas = Array::make(n, 0.0)
  let two_sigma2 = 2.0 * sigma * sigma
  for i = 0; i < n; i = i + 1 {
    let diff = (i - center).to_double()
    let g = base + amplitude * @math.exp(-(diff * diff / two_sigma2))
    gammas[i] = g
  }
  gammas
}
