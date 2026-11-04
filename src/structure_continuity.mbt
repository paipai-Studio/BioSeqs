///|
/// Structure continuity checks — detect gaps and breaks in protein/RNA chains.
///
/// Biotite provides `check_continuity()` and `check_id_continuity()` to
/// validate that structure chains have no missing residues or atoms.

///|
/// Continuity issue type.
pub enum ContinuityIssue {
  MissingResidue(Char, Int32, Int32) // chain, expected, found
  ChainBreak(Char, Char, Int32) // chain1, chain2, resseq
  SpatialGap(Int32, Int32, Double) // resseq1, resseq2, distance
  MissingBackboneAtom(Int32, String) // resseq, atom_name
}

///|
/// Get a human-readable description of a continuity issue.
pub fn ContinuityIssue::to_string(self : ContinuityIssue) -> String {
  match self {
    ContinuityIssue::MissingResidue(chain, expected, found) =>
      "Missing residue in chain " +
      chain.to_string() +
      ": expected " +
      expected.to_string() +
      ", found " +
      found.to_string()
    ContinuityIssue::ChainBreak(chain1, chain2, resseq) =>
      "Chain break: " +
      chain1.to_string() +
      " -> " +
      chain2.to_string() +
      " at resseq " +
      resseq.to_string()
    ContinuityIssue::SpatialGap(resseq1, resseq2, distance) =>
      "Spatial gap between residues " +
      resseq1.to_string() +
      " and " +
      resseq2.to_string() +
      ": " +
      distance.to_string() +
      " A"
    ContinuityIssue::MissingBackboneAtom(resseq, atom_name) =>
      "Missing backbone atom " + atom_name + " in residue " + resseq.to_string()
  }
}

///|
/// Check residue numbering continuity.
///
/// Detects gaps in residue sequence numbers (e.g., 1,2,3,5,6 — missing 4).
pub fn check_residue_continuity(
  chain_id : Char,
  resseqs : Array[Int32],
) -> Array[ContinuityIssue] {
  let issues : Array[ContinuityIssue] = []
  if resseqs.length() < 2 {
    return issues
  }
  for i in 1..<resseqs.length() {
    let expected = resseqs[i - 1] + 1
    if resseqs[i] != expected && resseqs[i] != resseqs[i - 1] {
      issues.push(
        ContinuityIssue::MissingResidue(chain_id, expected, resseqs[i]),
      )
    }
  }
  issues
}

///|
/// Check spatial continuity of C-alpha trace.
///
/// Detects large jumps (> 4.2 Å) between consecutive C-alpha atoms,
/// which typically indicate missing residues or chain breaks.
pub fn check_spatial_continuity(
  ca_coords : Array[Vector3],
  resseqs : Array[Int32],
  max_distance~ : Double = 4.2,
) -> Array[ContinuityIssue] {
  let issues : Array[ContinuityIssue] = []
  let n = ca_coords.length()
  if n < 2 {
    return issues
  }
  for i in 1..<n {
    let dist = ca_coords[i - 1].distance(ca_coords[i])
    if dist > max_distance {
      let r1 = if i - 1 < resseqs.length() { resseqs[i - 1] } else { 0 }
      let r2 = if i < resseqs.length() { resseqs[i] } else { 0 }
      issues.push(ContinuityIssue::SpatialGap(r1, r2, dist))
    }
  }
  issues
}

///|
/// Check for missing backbone atoms.
///
/// For each residue, verifies that N, CA, C, and O atoms are present.
pub fn check_backbone_completeness(
  residues : Array[Residue],
) -> Array[ContinuityIssue] {
  let issues : Array[ContinuityIssue] = []
  let required = ["N", "CA", "C", "O"]
  for r in residues {
    let atom_names : Array[String] = []
    for a in r.get_atoms() {
      atom_names.push(a.name)
    }
    for req in required {
      let mut found = false
      for name in atom_names {
        if name == req {
          found = true
          break
        }
      }
      if !found {
        issues.push(ContinuityIssue::MissingBackboneAtom(r.resseq, req))
      }
    }
  }
  issues
}

///|
/// Comprehensive structure continuity check.
///
/// Combines residue numbering, spatial, and backbone completeness checks.
pub fn check_structure_continuity(
  structure : Structure,
) -> Array[ContinuityIssue] {
  let issues : Array[ContinuityIssue] = []
  let models = structure.get_models()
  if models.length() == 0 {
    return issues
  }
  let chains = models[0].get_chains()
  for chain in chains {
    let residues = chain.get_residues()
    // Collect residue sequence numbers
    let resseqs : Array[Int32] = []
    let ca_coords : Array[Vector3] = []
    for r in residues {
      resseqs.push(r.resseq)
      for a in r.get_atoms() {
        if a.name == "CA" {
          ca_coords.push(a.get_coord())
        }
      }
    }
    // Check residue numbering
    let res_issues = check_residue_continuity(chain.get_id(), resseqs)
    for issue in res_issues {
      issues.push(issue)
    }
    // Check spatial continuity
    let spatial_issues = check_spatial_continuity(ca_coords, resseqs)
    for issue in spatial_issues {
      issues.push(issue)
    }
    // Check backbone completeness
    let backbone_issues = check_backbone_completeness(residues)
    for issue in backbone_issues {
      issues.push(issue)
    }
  }
  issues
}

///|
/// Count issues by type.
pub fn count_continuity_issues(
  issues : Array[ContinuityIssue],
) -> (Int32, Int32, Int32, Int32) {
  let mut missing_res = 0
  let mut chain_break = 0
  let mut spatial_gap = 0
  let mut missing_bb = 0
  for issue in issues {
    match issue {
      ContinuityIssue::MissingResidue(_, _, _) => missing_res = missing_res + 1
      ContinuityIssue::ChainBreak(_, _, _) => chain_break = chain_break + 1
      ContinuityIssue::SpatialGap(_, _, _) => spatial_gap = spatial_gap + 1
      ContinuityIssue::MissingBackboneAtom(_, _) => missing_bb = missing_bb + 1
    }
  }
  (missing_res, chain_break, spatial_gap, missing_bb)
}

///|
fn _continuity_unused_constructors() -> Unit {
  ignore(ContinuityIssue::ChainBreak('A', 'B', 0))
}
