///|
/// DOPE (Discrete Optimized Protein Energy) potential scoring.
///
/// A distance-dependent statistical potential for assessing protein
/// structure quality. For each pair of residues (i, j) with sufficient
/// sequence separation, the Cα-Cα distance is binned and the pair energy
/// is looked up. Lower (more negative) total score = better structure.
///
/// This implementation uses a simplified DOPE-like potential derived from
/// residue physicochemical properties (hydrophobicity, charge, size,
/// aromaticity). The framework supports plugging in real Modeller DOPE
/// tables if available.

///|
/// DOPE scoring result.
pub struct DopeScore {
  total_score : Double
  per_residue : Array[Double]
  num_pairs : Int
  normalized_score : Double
}

///|
/// DOPE potential table.
pub struct DopePotential {
  /// energy[pair_idx][dist_bin] — pair_idx = pair_index(aa1, aa2), 210 pairs
  energy : Array[Array[Double]]
  dist_cutoff : Double
  bin_size : Double
  num_bins : Int
  min_seq_sep : Int
}

// 20 standard amino acids in canonical order
let aa_names = [
  "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE",
  "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL",
]

// Residue properties: (hydrophobicity, charge, aromatic, size_class)
// hydrophobicity: 0=hydrophilic, 1=moderate, 2=hydrophobic
// charge: -1=negative, 0=neutral, 1=positive
// aromatic: 0=no, 1=yes
// size: 0=small, 1=medium, 2=large
let aa_hydro = [1, 0, 0, 0, 1, 0, 0, 0, 0, 2, 2, 0, 1, 2, 1, 0, 0, 2, 1, 2]
let dope_aa_charge = [0, 1, 0, -1, 0, 0, -1, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0]
let aa_arom = [0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 1, 1, 0]
let aa_size = [0, 2, 1, 1, 0, 1, 1, 0, 1, 2, 2, 1, 1, 2, 1, 0, 0, 2, 2, 1]

///|
/// Get amino acid index (0-19) from residue name, or -1 if non-standard.
fn aa_index(resname : String) -> Int {
  let upper = to_uppercase(resname)
  for i = 0; i < 20; i = i + 1 {
    if upper == aa_names[i] {
      return i
    }
  }
  -1
}

///|
/// Pair index for symmetric pair (i, j) with i <= j.
/// For 20 amino acids: pair_index = i * 20 - i*(i-1)/2 + (j - i)
/// Total pairs = 210
fn pair_index(i : Int, j : Int) -> Int {
  let a = if i <= j { i } else { j }
  let b = if i <= j { j } else { i }
  a * 20 - a * (a - 1) / 2 + (b - a)
}

///|
/// Compute pair compatibility score from residue properties.
/// Positive = favorable contact, negative = unfavorable.
fn pair_compatibility(i : Int, j : Int) -> Double {
  if i == j {
    // Same residue type
    if aa_arom[i] == 1 {
      return 1.5 // aromatic stacking
    }
    if aa_hydro[i] == 2 {
      return 2.0 // hydrophobic-hydrophobic
    }
    if dope_aa_charge[i] != 0 {
      return -1.0 // same charge repulsion
    }
    if aa_names[i] == "CYS" {
      return 3.5 // disulfide
    }
    return 0.5
  }
  let mut score = 0.0
  // Charge interaction
  let ci = dope_aa_charge[i]
  let cj = dope_aa_charge[j]
  if ci != 0 && cj != 0 {
    if ci * cj > 0 {
      score = score - 1.5 // same charge
    } else {
      score = score + 2.5 // opposite charge
    }
  }
  // Hydrophobic interaction
  let hi = aa_hydro[i]
  let hj = aa_hydro[j]
  if hi == 2 && hj == 2 {
    score = score + 2.0 // hydrophobic-hydrophobic
  } else if hi == 2 && hj == 0 {
    score = score - 0.5 // hydrophobic-hydrophilic
  } else if hi == 0 && hj == 2 {
    score = score - 0.5
  }
  // Aromatic interaction
  if aa_arom[i] == 1 && aa_arom[j] == 1 {
    score = score + 1.5
  }
  // Cysteine special (disulfide)
  if aa_names[i] == "CYS" && aa_names[j] == "CYS" {
    score = score + 3.0
  }
  // Proline restriction
  if aa_names[i] == "PRO" || aa_names[j] == "PRO" {
    score = score - 0.5
  }
  // Size compatibility
  if aa_size[i] == 2 && aa_size[j] == 2 {
    score = score + 0.3 // large-large packing
  }
  score
}

///|
/// Compute distance-dependent energy for a pair at a given distance.
/// Uses a physics-inspired profile:
///   d < 3.5:  steric clash (+5.0)
///   3.5-5.0:  close contact (pair-dependent)
///   5.0-8.0:  medium contact (moderate pair-dependent)
///   8.0-12.0: long range (weak pair-dependent)
///   d > 12.0: no contact (neutral)
fn dist_energy(d : Double, compat : Double) -> Double {
  if d < 3.5 {
    return 5.0 // steric clash
  }
  if d < 5.0 {
    return -compat * 1.5 // close contact
  }
  if d < 8.0 {
    return -compat * 1.0 // medium contact
  }
  if d < 12.0 {
    return -compat * 0.3 // long range
  }
  0.0 // no contact
}

///|
/// Build the default simplified DOPE potential table.
/// 210 pairs × 30 bins (0-15 Å, 0.5 Å bins).
pub fn DopePotential::default() -> DopePotential {
  let num_bins = 30
  let num_pairs = 210
  let energy : Array[Array[Double]] = Array::makei(num_pairs, { _ =>
    Array::make(num_bins, 0.0)
  })
  // Fill energy table
  for i = 0; i < 20; i = i + 1 {
    for j = i; j < 20; j = j + 1 {
      let pidx = pair_index(i, j)
      let compat = pair_compatibility(i, j)
      for b = 0; b < num_bins; b = b + 1 {
        let dist = b.to_double() * 0.5 + 0.25 // bin center
        energy[pidx][b] = dist_energy(dist, compat)
      }
    }
  }
  DopePotential::{
    energy,
    dist_cutoff: 15.0,
    bin_size: 0.5,
    num_bins,
    min_seq_sep: 5,
  }
}

///|
/// Create a DOPE potential with custom parameters.
pub fn DopePotential::new(
  dist_cutoff? : Double = 15.0,
  bin_size? : Double = 0.5,
  min_seq_sep? : Int = 5,
) -> DopePotential {
  let num_bins = (dist_cutoff / bin_size).to_int()
  let num_pairs = 210
  let energy : Array[Array[Double]] = Array::makei(num_pairs, { _ =>
    Array::make(num_bins, 0.0)
  })
  for i = 0; i < 20; i = i + 1 {
    for j = i; j < 20; j = j + 1 {
      let pidx = pair_index(i, j)
      let compat = pair_compatibility(i, j)
      for b = 0; b < num_bins; b = b + 1 {
        let dist = b.to_double() * bin_size + bin_size / 2.0
        energy[pidx][b] = dist_energy(dist, compat)
      }
    }
  }
  DopePotential::{ energy, dist_cutoff, bin_size, num_bins, min_seq_sep }
}

///|
/// Get the distance bin for a given distance.
fn dist_bin(d : Double, bin_size : Double, num_bins : Int) -> Int {
  let b = (d / bin_size).to_int()
  if b < 0 {
    0
  } else if b >= num_bins {
    num_bins - 1
  } else {
    b
  }
}

///|
/// Score a single residue pair.
pub fn DopePotential::score_pair(
  self : DopePotential,
  resname1 : String,
  resname2 : String,
  dist : Double,
) -> Double {
  let i1 = aa_index(resname1)
  let i2 = aa_index(resname2)
  if i1 < 0 || i2 < 0 {
    return 0.0
  }
  if dist >= self.dist_cutoff {
    return 0.0
  }
  let pidx = pair_index(i1, i2)
  let b = dist_bin(dist, self.bin_size, self.num_bins)
  self.energy[pidx][b]
}

///|
/// Score an entire protein chain.
/// Computes total DOPE energy and per-residue contributions.
pub fn DopePotential::score(self : DopePotential, chain : Chain) -> DopeScore? {
  let residues = chain.get_residues()
  let n = residues.length()
  if n < 2 {
    return None
  }
  // Collect Cα coordinates and residue types
  let coords : Array[Vector3] = Array(capacity=0)
  let aa_idx : Array[Int] = Array(capacity=0)
  for res in residues {
    match res.get_ca() {
      Some(ca) => {
        coords.push(ca.coord)
        aa_idx.push(aa_index(res.resname))
      }
      None => ()
    }
  }
  let m = coords.length()
  if m < 2 {
    return None
  }
  let per_residue : Array[Double] = Array::make(m, 0.0)
  let mut total = 0.0
  let mut num_pairs = 0
  for i = 0; i < m; i = i + 1 {
    for j = i + 1; j < m; j = j + 1 {
      // Sequence separation check
      let seq_sep = j - i
      if seq_sep < self.min_seq_sep {
        continue
      }
      let ai = aa_idx[i]
      let aj = aa_idx[j]
      if ai < 0 || aj < 0 {
        continue
      }
      let dist = coords[i].distance(coords[j])
      if dist >= self.dist_cutoff {
        continue
      }
      let pidx = pair_index(ai, aj)
      let b = dist_bin(dist, self.bin_size, self.num_bins)
      let e = self.energy[pidx][b]
      total = total + e
      per_residue[i] = per_residue[i] + e
      per_residue[j] = per_residue[j] + e
      num_pairs = num_pairs + 1
    }
  }
  let normalized = if num_pairs > 0 {
    total / num_pairs.to_double()
  } else {
    0.0
  }
  Some(DopeScore::{
    total_score: total,
    per_residue,
    num_pairs,
    normalized_score: normalized,
  })
}

///|
/// Score per-residue DOPE profile (for identifying poorly modeled regions).
pub fn DopePotential::score_profile(
  self : DopePotential,
  chain : Chain,
) -> Array[Double]? {
  match self.score(chain) {
    Some(result) => Some(result.per_residue)
    None => None
  }
}

///|
/// Compare DOPE scores of two chains (e.g., model vs reference).
/// Returns (score1, score2, difference). Lower = better.
pub fn DopePotential::compare(
  self : DopePotential,
  chain1 : Chain,
  chain2 : Chain,
) -> (Double, Double, Double)? {
  let s1 = match self.score(chain1) {
    Some(s) => s.total_score
    None => return None
  }
  let s2 = match self.score(chain2) {
    Some(s) => s.total_score
    None => return None
  }
  Some((s1, s2, s1 - s2))
}

///|
/// Rank multiple models by DOPE score (ascending = best first).
pub fn DopePotential::rank_models(
  self : DopePotential,
  chains : Array[Chain],
) -> Array[(Int, Double)] {
  let scored : Array[(Int, Double)] = Array(capacity=0)
  for i = 0; i < chains.length(); i = i + 1 {
    match self.score(chains[i]) {
      Some(s) => scored.push((i, s.total_score))
      None => scored.push((i, 1.0e100))
    }
  }
  // Sort by score ascending (lower = better)
  scored.sort_by(fn(a, b) {
    if a.1 < b.1 {
      -1
    } else if a.1 > b.1 {
      1
    } else {
      0
    }
  })
  scored
}

///|
/// Get the raw energy table value for a specific pair and distance bin.
pub fn DopePotential::get_energy(
  self : DopePotential,
  aa1 : Int,
  aa2 : Int,
  bin : Int,
) -> Double {
  if aa1 < 0 || aa1 >= 20 || aa2 < 0 || aa2 >= 20 {
    return 0.0
  }
  if bin < 0 || bin >= self.num_bins {
    return 0.0
  }
  let pidx = pair_index(aa1, aa2)
  self.energy[pidx][bin]
}

///|
/// Get amino acid index from name (public wrapper).
pub fn dope_aa_index(resname : String) -> Int {
  aa_index(resname)
}

///|
/// Get pair index (public wrapper).
pub fn dope_pair_index(i : Int, j : Int) -> Int {
  pair_index(i, j)
}

///|
/// Get pair compatibility score (public wrapper).
pub fn dope_pair_compatibility(i : Int, j : Int) -> Double {
  pair_compatibility(i, j)
}

// ===== DopeScore accessors =====

///|
pub fn DopeScore::get_total(self : DopeScore) -> Double {
  self.total_score
}

///|
pub fn DopeScore::get_per_residue(self : DopeScore) -> Array[Double] {
  self.per_residue
}

///|
pub fn DopeScore::get_num_pairs(self : DopeScore) -> Int {
  self.num_pairs
}

///|
pub fn DopeScore::get_normalized(self : DopeScore) -> Double {
  self.normalized_score
}

///|
/// Find the worst-scoring residues (highest per-residue energy).
/// Returns indices sorted worst-first.
pub fn DopeScore::worst_residues(self : DopeScore, k : Int) -> Array[Int] {
  let n = self.per_residue.length()
  let indexed : Array[(Int, Double)] = Array(capacity=0)
  for i = 0; i < n; i = i + 1 {
    indexed.push((i, self.per_residue[i]))
  }
  indexed.sort_by(fn(a, b) {
    if a.1 > b.1 {
      -1
    } else if a.1 < b.1 {
      1
    } else {
      0
    }
  })
  let result : Array[Int] = Array(capacity=0)
  let count = if k < indexed.length() { k } else { indexed.length() }
  for i = 0; i < count; i = i + 1 {
    result.push(indexed[i].0)
  }
  result
}

///|
/// Find the best-scoring residues (lowest per-residue energy).
pub fn DopeScore::best_residues(self : DopeScore, k : Int) -> Array[Int] {
  let n = self.per_residue.length()
  let indexed : Array[(Int, Double)] = Array(capacity=0)
  for i = 0; i < n; i = i + 1 {
    indexed.push((i, self.per_residue[i]))
  }
  indexed.sort_by(fn(a, b) {
    if a.1 < b.1 {
      -1
    } else if a.1 > b.1 {
      1
    } else {
      0
    }
  })
  let result : Array[Int] = Array(capacity=0)
  let count = if k < indexed.length() { k } else { indexed.length() }
  for i = 0; i < count; i = i + 1 {
    result.push(indexed[i].0)
  }
  result
}
