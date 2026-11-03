/// QED (Quantitative Estimate of Drug-likeness).
///
/// RDKit equivalent: `rdkit.Chem.QED`
///
/// QED is a continuous 0-1 score combining 8 physicochemical properties
/// with desirability functions. Unlike Lipinski's binary pass/fail,
/// QED gives a graded drug-likeness score.
///
/// Reference: Bickerton et al., Nature Chem 2012.

///|
/// QED result with component scores
pub struct RdkitQedResult {
  qed : Double
  mw : Double
  logp : Double
  tpsa : Double
  hbd : Int32
  hba : Int32
  rotatable_bonds : Int32
  aromatic_rings : Int32
  alerts : Int32
}

///|
pub fn RdkitQedResult::qed(self : RdkitQedResult) -> Double {
  self.qed
}

///|
pub fn RdkitQedResult::mw(self : RdkitQedResult) -> Double {
  self.mw
}

///|
pub fn RdkitQedResult::logp(self : RdkitQedResult) -> Double {
  self.logp
}

///|
pub fn RdkitQedResult::tpsa(self : RdkitQedResult) -> Double {
  self.tpsa
}

///|
/// Desirability function: Gaussian-like peak.
/// d(x) = exp(-((x - mu)^2) / (2 * sigma^2))
fn qed_gaussian(x : Double, mu : Double, sigma : Double) -> Double {
  let d = x - mu
  @math.exp(-(d * d) / (2.0 * sigma * sigma))
}

///|
/// Desirability function: asymmetric with plateau.
fn qed_asymmetric(
  x : Double,
  mu : Double,
  top : Double,
  bottom : Double,
  left : Double,
) -> Double {
  if x < left {
    qed_gaussian(x, left, top)
  } else {
    qed_gaussian(x, left, bottom)
  }
}

///|
/// Desirability for molecular weight (optimal ~350)
fn qed_d_mw(mw : Double) -> Double {
  qed_asymmetric(mw, 350.0, 150.0, 2.5, 250.0)
}

///|
/// Desirability for LogP (optimal ~2.0)
fn qed_d_logp(logp : Double) -> Double {
  qed_asymmetric(logp, 2.0, 1.0, 0.5, 2.0)
}

///|
/// Desirability for TPSA (optimal ~90)
fn qed_d_tpsa(tpsa : Double) -> Double {
  qed_asymmetric(tpsa, 90.0, 40.0, 0.5, 60.0)
}

///|
/// Desirability for HBD (optimal 0-2)
fn qed_d_hbd(hbd : Int32) -> Double {
  let x = hbd.to_double()
  qed_gaussian(x, 0.0, 2.0)
}

///|
/// Desirability for HBA (optimal ~3)
fn qed_d_hba(hba : Int32) -> Double {
  let x = hba.to_double()
  qed_gaussian(x, 3.0, 3.0)
}

///|
/// Desirability for rotatable bonds (optimal 0-2)
fn qed_d_rot(rot : Int32) -> Double {
  let x = rot.to_double()
  qed_gaussian(x, 0.0, 3.0)
}

///|
/// Desirability for aromatic rings (optimal 0-2)
fn qed_d_arom_rings(rings : Int32) -> Double {
  let x = rings.to_double()
  qed_gaussian(x, 1.0, 2.0)
}

///|
/// Desirability for PAINS alerts (optimal 0)
fn qed_d_alerts(alerts : Int32) -> Double {
  if alerts == 0 {
    1.0
  } else {
    0.0
  }
}

///|
/// Compute QED score for a molecule.
/// Returns 0-1 score (higher = more drug-like).
pub fn rdkit_qed(mol : Molecule) -> RdkitQedResult {
  let mw = mol_molecular_weight(mol)
  let logp = mol_logp(mol)
  let tpsa = mol_tpsa(mol)
  let hbd = mol_hbd_count(mol)
  let hba = mol_hba_count(mol)
  let rot = mol_rotatable_bonds(mol)
  let arom_rings = mol_aromatic_ring_count(mol)
  let alerts = rdkit_pains_count(mol)
  // Compute product of desirability functions
  let d_mw = qed_d_mw(mw)
  let d_logp = qed_d_logp(logp)
  let d_tpsa = qed_d_tpsa(tpsa)
  let d_hbd = qed_d_hbd(hbd)
  let d_hba = qed_d_hba(hba)
  let d_rot = qed_d_rot(rot)
  let d_rings = qed_d_arom_rings(arom_rings)
  let d_alerts = qed_d_alerts(alerts)
  let product = d_mw *
    d_logp *
    d_tpsa *
    d_hbd *
    d_hba *
    d_rot *
    d_rings *
    d_alerts
  // QED = product^(1/8)
  let qed = if product > 0.0 { @math.exp(@math.ln(product) / 8.0) } else { 0.0 }
  {
    qed,
    mw,
    logp,
    tpsa,
    hbd,
    hba,
    rotatable_bonds: rot,
    aromatic_rings: arom_rings,
    alerts,
  }
}

///|
/// Quick QED score (just the number).
pub fn rdkit_qed_score(mol : Molecule) -> Double {
  rdkit_qed(mol).qed
}

///|
/// Is QED score above a threshold (default 0.5)?
pub fn rdkit_qed_passes(mol : Molecule, threshold~ : Double = 0.5) -> Bool {
  rdkit_qed_score(mol) >= threshold
}
