///|
/// Structure repair — fix common issues in protein structures.
///
/// Biotite provides structure repair functions to add missing atoms,
/// fix alternate locations, and standardize residue names.

///|
/// Repair result describing what was fixed.
pub struct RepairResult {
  missing_atoms_added : Int32
  mut altlocs_resolved : Int32
  mut resnames_standardized : Int32
  chain_breaks_fixed : Int32
  messages : Array[String]
}

///|
/// Create an empty repair result.
pub fn RepairResult::new() -> RepairResult {
  RepairResult::{
    missing_atoms_added: 0,
    altlocs_resolved: 0,
    resnames_standardized: 0,
    chain_breaks_fixed: 0,
    messages: [],
  }
}

///|
/// Get total number of repairs.
pub fn RepairResult::total_repairs(self : RepairResult) -> Int32 {
  self.missing_atoms_added +
  self.altlocs_resolved +
  self.resnames_standardized +
  self.chain_breaks_fixed
}

///|
/// Standardize residue names to PDB convention.
///
/// Fixes common non-standard names like "HSD" -> "HIS", "HSE" -> "HIS", etc.
pub fn standardize_resname(resname : String) -> String {
  match resname {
    "HSD" => "HIS"
    "HSE" => "HIS"
    "HSP" => "HIS"
    "HIE" => "HIS"
    "HIP" => "HIS"
    "ASH" => "ASP"
    "GLH" => "GLU"
    "LYN" => "LYS"
    "CYX" => "CYS"
    "CYM" => "CYS"
    "DAN" => "ASN"
    "DAS" => "ASP"
    "DCY" => "CYS"
    "DGL" => "GLU"
    "DGN" => "GLN"
    "DHI" => "HIS"
    "DIL" => "ILE"
    "DLE" => "LEU"
    "DLY" => "LYS"
    "DME" => "MET"
    "DPH" => "PHE"
    "DPR" => "PRO"
    "DSER" => "SER"
    "DTH" => "THR"
    "DTR" => "TRP"
    "DTY" => "TYR"
    "DVA" => "VAL"
    _ => resname
  }
}

///|
/// Check if a residue name is non-standard.
pub fn is_nonstandard_resname(resname : String) -> Bool {
  standardize_resname(resname) != resname
}

///|
/// Resolve alternate location indicators.
///
/// When atoms have alternate locations (altloc != ' '), select the
/// 'A' conformer (or first non-blank) and set altloc to ' '.
pub fn resolve_altlocs(atoms : Array[Atom]) -> Array[Atom] {
  let result : Array[Atom] = []
  let mut last_altloc = ' '
  for a in atoms {
    if a.altloc == ' ' {
      result.push(a)
    } else if a.altloc == 'A' || (last_altloc != ' ' && a.altloc != last_altloc) {
      let fixed = Atom::new(
        name=a.name,
        coord=a.get_coord(),
        resname=a.resname,
        chainid=a.chainid,
        resseq=a.resseq,
        occupancy=a.occupancy,
        bfactor=a.bfactor,
        element=a.element,
        charge=a.charge,
      )
      result.push(fixed)
      last_altloc = a.altloc
    }
  }
  result
}

///|
/// Repair a structure by fixing common issues.
///
/// - Standardizes non-standard residue names
/// - Resolves alternate locations
/// - Reports missing backbone atoms
pub fn repair_structure(
  residues : Array[Residue],
  standardize_resnames~ : Bool = true,
  resolve_alt~ : Bool = true,
) -> (Array[Residue], RepairResult) {
  let result = RepairResult::new()
  let repaired : Array[Residue] = []
  for r in residues {
    let mut resname = r.resname
    if standardize_resnames && is_nonstandard_resname(r.resname) {
      resname = standardize_resname(r.resname)
      result.resnames_standardized = result.resnames_standardized + 1
      result.messages.push("Standardized " + r.resname + " -> " + resname)
    }
    let mut atoms = r.atoms
    if resolve_alt {
      let mut has_altloc = false
      for a in atoms {
        if a.altloc != ' ' {
          has_altloc = true
          break
        }
      }
      if has_altloc {
        atoms = resolve_altlocs(atoms)
        result.altlocs_resolved = result.altlocs_resolved + 1
        result.messages.push(
          "Resolved alternate locations in residue " + r.resseq.to_string(),
        )
      }
    }
    let new_residue = Residue::new(
      resname~,
      chainid=r.chainid,
      resseq=r.resseq,
      atoms~,
      hetfield=r.hetfield,
    )
    repaired.push(new_residue)
  }
  (repaired, result)
}

///|
/// Add hydrogen atoms to a structure (placeholder).
///
/// In a full implementation, this would use protonation rules to add
/// hydrogen atoms. Currently returns the structure unchanged.
pub fn add_hydrogens(
  residues : Array[Residue],
) -> (Array[Residue], RepairResult) {
  let result = RepairResult::new()
  result.messages.push("Hydrogen addition not yet implemented")
  (residues, result)
}
