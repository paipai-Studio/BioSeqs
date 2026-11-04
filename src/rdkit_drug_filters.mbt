/// Lipinski Rule of 5 and related drug-likeness filters.
///
/// RDKit equivalent: `rdkit.Chem.Lipinski`, `rdkit.Chem.Filters`
///
/// The Rule of 5 (Ro5) is the most widely used drug-likeness filter:
/// MW ≤ 500, LogP ≤ 5, HBD ≤ 5, HBA ≤ 10.
/// Violations ≤ 1 are generally considered acceptable.

///|
/// Lipinski Rule of 5 result
pub struct RdkitLipinskiResult {
  mw : Double
  logp : Double
  hbd : Int32
  hba : Int32
  violations : Int32
  passes : Bool
}

///|
pub fn RdkitLipinskiResult::mw(self : RdkitLipinskiResult) -> Double {
  self.mw
}

///|
pub fn RdkitLipinskiResult::logp(self : RdkitLipinskiResult) -> Double {
  self.logp
}

///|
pub fn RdkitLipinskiResult::hbd(self : RdkitLipinskiResult) -> Int32 {
  self.hbd
}

///|
pub fn RdkitLipinskiResult::hba(self : RdkitLipinskiResult) -> Int32 {
  self.hba
}

///|
pub fn RdkitLipinskiResult::violations(self : RdkitLipinskiResult) -> Int32 {
  self.violations
}

///|
pub fn RdkitLipinskiResult::passes(self : RdkitLipinskiResult) -> Bool {
  self.passes
}

///|
/// Evaluate Lipinski Rule of 5.
pub fn rdkit_lipinski(mol : Molecule) -> RdkitLipinskiResult {
  let mw = mol_molecular_weight(mol)
  let logp = mol_logp(mol)
  let hbd = mol_hbd_count(mol)
  let hba = mol_hba_count(mol)
  let mut violations = 0
  if mw > 500.0 {
    violations = violations + 1
  }
  if logp > 5.0 {
    violations = violations + 1
  }
  if hbd > 5 {
    violations = violations + 1
  }
  if hba > 10 {
    violations = violations + 1
  }
  { mw, logp, hbd, hba, violations, passes: violations <= 1, }
}

///|
/// Quick boolean check: does molecule pass Lipinski Ro5?
pub fn rdkit_passes_lipinski(mol : Molecule) -> Bool {
  rdkit_lipinski(mol).passes
}

///|
/// Count Lipinski violations (0 = perfect, 4 = all violated).
pub fn rdkit_lipinski_violations(mol : Molecule) -> Int32 {
  rdkit_lipinski(mol).violations
}

// ── Veber filter ──
// Good oral bioavailability: rotatable bonds ≤ 10, TPSA ≤ 140

///|
/// Veber filter: rotatable bonds ≤ 10 and TPSA ≤ 140.
pub fn rdkit_veber_filter(mol : Molecule) -> Bool {
  mol_rotatable_bonds(mol) <= 10 && mol_tpsa(mol) <= 140.0
}

// ── Egan filter ──
// Good oral absorption: LogP ≤ 5.88, TPSA ≤ 131.6

///|
/// Egan filter for oral absorption: LogP ≤ 5.88, TPSA ≤ 131.6.
pub fn rdkit_egan_filter(mol : Molecule) -> Bool {
  mol_logp(mol) <= 5.88 && mol_tpsa(mol) <= 131.6
}

// ── Ghose filter ──
// Drug-like: 160 ≤ MW ≤ 480, -0.4 ≤ LogP ≤ 5.6, 40 ≤ MR ≤ 130, 20 ≤ atoms ≤ 70

///|
/// Ghose filter: 160 ≤ MW ≤ 480, -0.4 ≤ LogP ≤ 5.6, 20 ≤ atoms ≤ 70.
pub fn rdkit_ghose_filter(mol : Molecule) -> Bool {
  let mw = mol_molecular_weight(mol)
  let logp = mol_logp(mol)
  let n = mol_num_atoms(mol)
  mw >= 160.0 &&
  mw <= 480.0 &&
  logp >= -0.4 &&
  logp <= 5.6 &&
  n >= 20 &&
  n <= 70
}

///|
/// Combined drug-likeness report: passes all four filters?
pub struct RdkitDrugLikeness {
  lipinski : Bool
  veber : Bool
  egan : Bool
  ghose : Bool
  all_pass : Bool
}

///|
pub fn RdkitDrugLikeness::lipinski(self : RdkitDrugLikeness) -> Bool {
  self.lipinski
}

///|
pub fn RdkitDrugLikeness::veber(self : RdkitDrugLikeness) -> Bool {
  self.veber
}

///|
pub fn RdkitDrugLikeness::egan(self : RdkitDrugLikeness) -> Bool {
  self.egan
}

///|
pub fn RdkitDrugLikeness::ghose(self : RdkitDrugLikeness) -> Bool {
  self.ghose
}

///|
pub fn RdkitDrugLikeness::all_pass(self : RdkitDrugLikeness) -> Bool {
  self.all_pass
}

///|
/// Evaluate all drug-likeness filters at once.
pub fn rdkit_drug_likeness(mol : Molecule) -> RdkitDrugLikeness {
  let l = rdkit_passes_lipinski(mol)
  let v = rdkit_veber_filter(mol)
  let e = rdkit_egan_filter(mol)
  let g = rdkit_ghose_filter(mol)
  { lipinski: l, veber: v, egan: e, ghose: g, all_pass: l && v && e && g, }
}
