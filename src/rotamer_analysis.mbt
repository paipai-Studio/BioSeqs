///|
/// Rotamer analysis: compute side-chain torsion angles (chi1-chi4)
/// from structure, classify into rotamer bins, detect outliers.

///|
/// Rotamer bin classification.
pub(all) enum RotamerBin {
  Trans // ~180°
  GauchePlus // ~+60°
  GaucheMinus // ~-60°
  Unknown
}

///|
/// A rotamer analysis result for one residue.
pub(all) struct RotamerResult {
  resseq : Int
  resname : String
  chi1 : Double // -180 to 180, NaN if not applicable
  chi2 : Double
  chi3 : Double
  chi4 : Double
  chi1_bin : RotamerBin
  chi2_bin : RotamerBin
  is_outlier : Bool
  deviation : Double // deviation from nearest bin center
}

///|
/// Classify an angle into a rotamer bin.
pub fn classify_rotamer_angle(angle : Double) -> RotamerBin {
  let mut a = angle
  while a > 180.0 {
    a = a - 360.0
  }
  while a < -180.0 {
    a = a + 360.0
  }
  let d_trans = (a - 180.0).abs()
  let d_gp = (a - 60.0).abs()
  let d_gm = (a + 60.0).abs()
  let d_trans2 = (a + 180.0).abs()
  let dt = if d_trans < d_trans2 { d_trans } else { d_trans2 }
  if dt <= d_gp && dt <= d_gm {
    Trans
  } else if d_gp <= d_gm {
    GauchePlus
  } else {
    GaucheMinus
  }
}

///|
/// Deviation of angle from nearest rotamer bin center.
pub fn rotamer_deviation(angle : Double) -> Double {
  let mut a = angle
  while a > 180.0 {
    a = a - 360.0
  }
  while a < -180.0 {
    a = a + 360.0
  }
  let d_trans = (a - 180.0).abs()
  let d_trans2 = (a + 180.0).abs()
  let dt = if d_trans < d_trans2 { d_trans } else { d_trans2 }
  let d_gp = (a - 60.0).abs()
  let d_gm = (a + 60.0).abs()
  let mut min_dev = dt
  if d_gp < min_dev {
    min_dev = d_gp
  }
  if d_gm < min_dev {
    min_dev = d_gm
  }
  min_dev
}

///|
/// Get bin center angle.
pub fn rotamer_bin_center(bin : RotamerBin) -> Double {
  match bin {
    Trans => 180.0
    GauchePlus => 60.0
    GaucheMinus => -60.0
    Unknown => 0.0
  }
}

///|
/// Get bin name as string.
pub fn rotamer_bin_name(bin : RotamerBin) -> String {
  match bin {
    Trans => "trans"
    GauchePlus => "gauche+"
    GaucheMinus => "gauche-"
    Unknown => "unknown"
  }
}

///|
/// Check if a residue type has a chi1 angle.
pub fn has_chi1(resname : String) -> Bool {
  let r = to_lowercase(resname)
  not(r == "gly" || r == "ala")
}

///|
/// Check if a residue type has a chi2 angle.
pub fn has_chi2(resname : String) -> Bool {
  let r = to_lowercase(resname)
  not(r == "gly" || r == "ala" || r == "ser" || r == "cys" || r == "pro" || r == "thr" || r == "val")
}

///|
/// Check if a residue type has a chi3 angle.
pub fn has_chi3(resname : String) -> Bool {
  let r = to_lowercase(resname)
  r == "gln" || r == "glu" || r == "met" || r == "lys" || r == "arg" || r == "ile"
}

///|
/// Check if a residue type has a chi4 angle.
pub fn has_chi4(resname : String) -> Bool {
  let r = to_lowercase(resname)
  r == "lys" || r == "arg"
}

///|
/// Compute chi angles for a residue using Vector3 positions.
/// chi1 = dihedral(N, CA, CB, CG)
/// chi2 = dihedral(CA, CB, CG, CD)
/// chi3 = dihedral(CB, CG, CD, CE)
/// chi4 = dihedral(CG, CD, CE, CZ)
///
/// Returns (chi1, chi2, chi3, chi4) with NaN for missing angles.
pub fn compute_chi_angles(
  n_pos : Vector3?,
  ca_pos : Vector3?,
  cb_pos : Vector3?,
  cg_pos : Vector3?,
  cd_pos : Vector3?,
  ce_pos : Vector3?,
  cz_pos : Vector3?,
) -> (Double, Double, Double, Double) {
  let nan = 0.0 / 0.0
  let chi1 = match (n_pos, ca_pos, cb_pos, cg_pos) {
    (Some(n), Some(ca), Some(cb), Some(cg)) => calc_dihedral(n, ca, cb, cg)
    _ => nan
  }
  let chi2 = match (ca_pos, cb_pos, cg_pos, cd_pos) {
    (Some(ca), Some(cb), Some(cg), Some(cd)) => calc_dihedral(ca, cb, cg, cd)
    _ => nan
  }
  let chi3 = match (cb_pos, cg_pos, cd_pos, ce_pos) {
    (Some(cb), Some(cg), Some(cd), Some(ce)) => calc_dihedral(cb, cg, cd, ce)
    _ => nan
  }
  let chi4 = match (cg_pos, cd_pos, ce_pos, cz_pos) {
    (Some(cg), Some(cd), Some(ce), Some(cz)) => calc_dihedral(cg, cd, ce, cz)
    _ => nan
  }
  (chi1, chi2, chi3, chi4)
}

///|
/// Classify a rotamer and detect outliers.
/// Outlier threshold: deviation > 30° from nearest bin center.
pub fn classify_rotamer(
  resseq : Int,
  resname : String,
  chi1 : Double,
  chi2 : Double,
  chi3 : Double,
  chi4 : Double,
) -> RotamerResult {
  let nan = 0.0 / 0.0
  let chi1_bin = if chi1 == nan || chi1 != chi1 {
    RotamerBin::Unknown
  } else {
    classify_rotamer_angle(chi1)
  }
  let chi2_bin = if chi2 == nan || chi2 != chi2 {
    RotamerBin::Unknown
  } else {
    classify_rotamer_angle(chi2)
  }
  let mut is_outlier = false
  let mut deviation = 0.0
  if chi1 == chi1 && chi1 != nan {
    deviation = rotamer_deviation(chi1)
    if deviation > 30.0 {
      is_outlier = true
    }
  }
  if chi2 == chi2 && chi2 != nan {
    let d2 = rotamer_deviation(chi2)
    if d2 > 30.0 {
      is_outlier = true
    }
    if d2 > deviation {
      deviation = d2
    }
  }
  RotamerResult::{ resseq, resname, chi1, chi2, chi3, chi4, chi1_bin, chi2_bin, is_outlier, deviation }
}

///|
/// Detect rotamer outliers from a list of results.
pub fn rotamer_outliers(
  results : Array[RotamerResult],
) -> Array[RotamerResult] {
  let outliers : Array[RotamerResult] = []
  for r in results {
    if r.is_outlier {
      outliers.push(r)
    }
  }
  outliers
}

///|
/// Summary statistics of rotamer analysis.
pub fn rotamer_summary(
  results : Array[RotamerResult],
) -> (Int, Int, Int, Int, Int) {
  let mut total = 0
  let mut outliers = 0
  let mut trans_count = 0
  let mut gp_count = 0
  let mut gm_count = 0
  let nan = 0.0 / 0.0
  for r in results {
    if r.chi1 == r.chi1 && r.chi1 != nan {
      total = total + 1
      match r.chi1_bin {
        Trans => trans_count = trans_count + 1
        GauchePlus => gp_count = gp_count + 1
        GaucheMinus => gm_count = gm_count + 1
        Unknown => ()
      }
      if r.is_outlier {
        outliers = outliers + 1
      }
    }
  }
  (total, outliers, trans_count, gp_count, gm_count)
}

///|
/// Rotamer library probability (simplified).
pub fn rotamer_probability(
  resname : String,
  chi1_bin : RotamerBin,
  chi2_bin : RotamerBin,
) -> Double {
  let r = to_lowercase(resname)
  let p_chi1 = match chi1_bin {
    Trans => 0.40
    GaucheMinus => 0.35
    GauchePlus => 0.20
    Unknown => 0.05
  }
  let p_chi2 = match chi2_bin {
    Trans => 0.45
    GauchePlus => 0.30
    GaucheMinus => 0.20
    Unknown => 0.05
  }
  if r == "pro" {
    match chi1_bin {
      GaucheMinus => 0.60
      _ => 0.10
    }
  } else {
    p_chi1 * p_chi2
  }
}
