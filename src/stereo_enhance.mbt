///|
/// Enhanced stereochemistry analysis:
/// - Chirality detection (CA with 4 substituents)
/// - Cis-peptide detection
/// - Non-standard bond lengths/angles
/// - Peptide planarity (omega angle)

///|
/// Chirality result for a residue.
pub(all) struct ChiralityResult {
  resseq : Int
  resname : String
  has_chirality : Bool
  is_correct : Bool // true = L-amino acid, false = D-amino acid
  volume : Double // signed volume of tetrahedron
}

///|
/// Cis-peptide result.
pub(all) struct CisPeptideResult {
  resseq1 : Int
  resseq2 : Int
  omega : Double // dihedral C(i-1)-N(i)-CA(i)-C(i)
  is_cis : Bool
}

///|
/// Bond validation result.
pub(all) struct BondValidationResult {
  atom1 : String
  atom2 : String
  resseq : Int
  observed : Double
  expected : Double
  deviation : Double
  is_outlier : Bool
}

///|
/// Angle validation result.
pub(all) struct AngleValidationResult {
  atom1 : String
  atom2 : String
  atom3 : String
  resseq : Int
  observed : Double
  expected : Double
  deviation : Double
  is_outlier : Bool
}

///|
/// Peptide planarity result.
pub(all) struct PeptidePlanarityResult {
  resseq : Int
  omega : Double
  deviation_from_180 : Double
  is_planar : Bool
}

///|
/// Comprehensive stereochemistry summary.
pub(all) struct StereoSummary {
  num_residues : Int
  num_chiral : Int
  num_wrong_chirality : Int
  num_cis_peptides : Int
  num_bond_outliers : Int
  num_angle_outliers : Int
  num_nonplanar_peptides : Int
}

///|
/// Detect chirality at CA atom.
/// Computes signed volume of tetrahedron formed by N, CA, C, CB.
/// Positive volume = L-amino acid (correct), negative = D-amino acid.
pub fn detect_chirality(
  resseq : Int,
  resname : String,
  n_pos : Vector3?,
  ca_pos : Vector3?,
  c_pos : Vector3?,
  cb_pos : Vector3?,
) -> ChiralityResult {
  match (n_pos, ca_pos, c_pos, cb_pos) {
    (Some(n), Some(ca), Some(c), Some(cb)) => {
      // Vectors from CA to substituents
      let v1 = Vector3::sub(n, ca)
      let v2 = Vector3::sub(c, ca)
      let v3 = Vector3::sub(cb, ca)
      // Signed volume = (v1 × v2) · v3 / 6
      let cross = Vector3::cross(v1, v2)
      let vol = Vector3::dot(cross, v3) / 6.0
      let has_chirality = vol.abs() > 0.001
      let is_correct = vol > 0.0
      ChiralityResult::{ resseq, resname, has_chirality, is_correct, volume: vol }
    }
    _ =>
      ChiralityResult::{ resseq, resname, has_chirality: false, is_correct: false, volume: 0.0 }
  }
}

///|
/// Count cis-peptides from a list of omega angles.
/// A cis-peptide has |omega| < 30°.
pub fn count_cis_peptides(omega_angles : Array[CisPeptideResult]) -> Int {
  let mut count = 0
  for r in omega_angles {
    if r.is_cis {
      count = count + 1
    }
  }
  count
}

///|
/// Detect cis-peptide from omega dihedral angle.
/// omega = dihedral(C_prev, N, CA, C)
pub fn detect_cis_peptide(
  resseq1 : Int,
  resseq2 : Int,
  c_prev : Vector3?,
  n : Vector3?,
  ca : Vector3?,
  c : Vector3?,
) -> CisPeptideResult {
  match (c_prev, n, ca, c) {
    (Some(cp), Some(n), Some(ca), Some(c)) => {
      let omega = calc_dihedral(cp, n, ca, c)
      let is_cis = omega.abs() < 30.0
      CisPeptideResult::{ resseq1, resseq2, omega, is_cis }
    }
    _ => CisPeptideResult::{ resseq1, resseq2, omega: 0.0, is_cis: false }
  }
}

///|
/// Validate a bond length against expected value.
pub fn validate_bond_length(
  atom1 : String,
  atom2 : String,
  resseq : Int,
  pos1 : Vector3?,
  pos2 : Vector3?,
  expected : Double,
  tolerance : Double,
) -> BondValidationResult {
  match (pos1, pos2) {
    (Some(p1), Some(p2)) => {
      let observed = Vector3::distance_squared(p1, p2).sqrt()
      let deviation = observed - expected
      let is_outlier = deviation.abs() > tolerance
      BondValidationResult::{ atom1, atom2, resseq, observed, expected, deviation, is_outlier }
    }
    _ =>
      BondValidationResult::{ atom1, atom2, resseq, observed: 0.0, expected, deviation: 0.0, is_outlier: false }
  }
}

///|
/// Validate a bond angle against expected value.
pub fn validate_bond_angle(
  atom1 : String,
  atom2 : String,
  atom3 : String,
  resseq : Int,
  pos1 : Vector3?,
  pos2 : Vector3?,
  pos3 : Vector3?,
  expected : Double,
  tolerance : Double,
) -> AngleValidationResult {
  match (pos1, pos2, pos3) {
    (Some(p1), Some(p2), Some(p3)) => {
      let v1 = Vector3::sub(p1, p2)
      let v2 = Vector3::sub(p3, p2)
      let observed = Vector3::angle_degrees(v1, v2)
      let deviation = observed - expected
      let is_outlier = deviation.abs() > tolerance
      AngleValidationResult::{ atom1, atom2, atom3, resseq, observed, expected, deviation, is_outlier }
    }
    _ => AngleValidationResult::{
        atom1, atom2, atom3, resseq, observed: 0.0, expected, deviation: 0.0, is_outlier: false
      }
  }
}

///|
/// Check peptide planarity from omega angle.
/// Trans peptide: omega ≈ ±180°, cis: omega ≈ 0°.
/// Non-planar if deviation from 0 or 180 exceeds threshold.
pub fn peptide_planarity(
  resseq : Int,
  omega : Double,
  threshold : Double,
) -> PeptidePlanarityResult {
  let dev_from_180 = (omega.abs() - 180.0).abs()
  let dev_from_0 = omega.abs()
  let deviation_from_180 = if dev_from_180 < dev_from_0 { dev_from_180 } else { dev_from_0 }
  let is_planar = deviation_from_180 < threshold
  PeptidePlanarityResult::{ resseq, omega, deviation_from_180, is_planar }
}

///|
/// Standard bond lengths for protein backbone (in Angstroms).
pub fn standard_bond_length(atom1 : String, atom2 : String) -> Double {
  let a1 = to_lowercase(atom1)
  let a2 = to_lowercase(atom2)
  let pair = if a1 < a2 { a1 + "-" + a2 } else { a2 + "-" + a1 }
  match pair {
    "c-n" => 1.33 // peptide bond
    "n-ca" => 1.47
    "ca-c" => 1.53
    "ca-cb" => 1.54
    "c-o" => 1.23 // carbonyl
    "c-oxt" => 1.25
    _ => 1.50 // generic
  }
}

///|
/// Standard bond angles for protein backbone (in degrees).
pub fn standard_bond_angle(atom1 : String, atom2 : String, atom3 : String) -> Double {
  let a2 = to_lowercase(atom2) // central atom
  match a2 {
    "n" => 121.0 // C-N-CA
    "ca" => 110.0 // N-CA-C, N-CA-CB, C-CA-CB
    "c" => 121.0 // N-C-O, CA-C-N, CA-C-O
    "cb" => 109.5 // tetrahedral
    _ => 109.5
  }
}

///|
/// Validate bond lengths for a list of results.
pub fn validate_bond_lengths(
  results : Array[BondValidationResult],
) -> Array[BondValidationResult] {
  let outliers : Array[BondValidationResult] = []
  for r in results {
    if r.is_outlier {
      outliers.push(r)
    }
  }
  outliers
}

///|
/// Validate bond angles for a list of results.
pub fn validate_bond_angles(
  results : Array[AngleValidationResult],
) -> Array[AngleValidationResult] {
  let outliers : Array[AngleValidationResult] = []
  for r in results {
    if r.is_outlier {
      outliers.push(r)
    }
  }
  outliers
}

///|
/// Generate comprehensive stereochemistry summary.
pub fn stereo_summary(
  chirality_results : Array[ChiralityResult],
  cis_results : Array[CisPeptideResult],
  bond_results : Array[BondValidationResult],
  angle_results : Array[AngleValidationResult],
  planarity_results : Array[PeptidePlanarityResult],
) -> StereoSummary {
  let num_residues = chirality_results.length()
  let mut num_chiral = 0
  let mut num_wrong = 0
  for r in chirality_results {
    if r.has_chirality {
      num_chiral = num_chiral + 1
      if not(r.is_correct) {
        num_wrong = num_wrong + 1
      }
    }
  }
  let mut num_cis = 0
  for r in cis_results {
    if r.is_cis {
      num_cis = num_cis + 1
    }
  }
  let mut num_bond_outliers = 0
  for r in bond_results {
    if r.is_outlier {
      num_bond_outliers = num_bond_outliers + 1
    }
  }
  let mut num_angle_outliers = 0
  for r in angle_results {
    if r.is_outlier {
      num_angle_outliers = num_angle_outliers + 1
    }
  }
  let mut num_nonplanar = 0
  for r in planarity_results {
    if not(r.is_planar) {
      num_nonplanar = num_nonplanar + 1
    }
  }
  StereoSummary::{
    num_residues,
    num_chiral,
    num_wrong_chirality: num_wrong,
    num_cis_peptides: num_cis,
    num_bond_outliers,
    num_angle_outliers,
    num_nonplanar_peptides: num_nonplanar,
  }
}
