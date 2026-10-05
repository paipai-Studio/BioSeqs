///|
/// Disulfide bond detection and geometric analysis.
///
/// Disulfide bonds are covalent S-S bonds between cysteine residues,
/// typically with SG-SG distance ~2.05 Å. They stabilize protein
/// structure and are characterized by five dihedral angles:
///   chi1  = N-CA-CB-SG
///   chi2  = CA-CB-SG-SG'
///   chiSS = CB-SG-SG'-CB'
///   chi2' = SG-SG'-CB'-CA'
///   chi1' = SG'-CB'-CA'-N'

///|
/// A disulfide bond between two cysteine residues.
pub(all) struct DisulfideBond {
  chainid1 : Char
  resseq1 : Int
  chainid2 : Char
  resseq2 : Int
  sg_sg_distance : Double  // S-S bond length (Å)
  chi1 : Double            // N-CA-CB-SG (degrees)
  chi2 : Double            // CA-CB-SG-SG'
  chi_ss : Double          // CB-SG-SG'-CB'
  chi2p : Double           // SG-SG'-CB'-CA'
  chi1p : Double           // SG'-CB'-CA'-N'
  oxidation_state : String // "oxidized" or "reduced"
}

///|
/// Find all disulfide bonds in a chain.
///
/// Scans all cysteine pairs for SG-SG distance below the cutoff
/// (default 2.5 Å, typical S-S bond is ~2.05 Å).
/// Computes five dihedral angles for each bond.
pub fn find_disulfide_bonds(
  chain : Chain,
  distance_cutoff? : Double = 2.5,
) -> Array[DisulfideBond] {
  let residues = chain.get_residues()
  // Collect cysteine residues with SG atoms
  let cys_list = Array(capacity=0)
  for res in residues {
    if res.resname == "CYS" || res.resname == "Cys" {
      match res.find_atom("SG") {
        Some(_) => cys_list.push(res)
        None => ()
      }
    }
  }
  let bonds = Array(capacity=0)
  let nc = cys_list.length()
  let mut i = 0
  while i < nc {
    let mut j = i + 1
    while j < nc {
      let res1 = cys_list[i]
      let res2 = cys_list[j]
      let sg1 = res1.find_atom("SG").unwrap()
      let sg2 = res2.find_atom("SG").unwrap()
      let d = sg1.get_coord().distance(sg2.get_coord())
      if d <= distance_cutoff {
        let bond = compute_disulfide_bond(res1, res2, d)
        match bond {
          Some(b) => bonds.push(b)
          None => ()
        }
      }
      j = j + 1
    }
    i = i + 1
  }
  bonds
}

///|
/// Compute disulfide bond geometry between two cysteine residues.
fn compute_disulfide_bond(
  res1 : Residue,
  res2 : Residue,
  sg_sg_dist : Double,
) -> DisulfideBond? {
  let sg1 = res1.find_atom("SG")
  let sg2 = res2.find_atom("SG")
  let n1 = res1.find_atom("N")
  let ca1 = res1.find_atom("CA")
  let cb1 = res1.find_atom("CB")
  let n2 = res2.find_atom("N")
  let ca2 = res2.find_atom("CA")
  let cb2 = res2.find_atom("CB")
  // SG atoms required
  if sg1.is_empty() || sg2.is_empty() {
    return None
  }
  let sg1c = sg1.unwrap().get_coord()
  let sg2c = sg2.unwrap().get_coord()
  // Compute dihedral angles (default 0 if atoms missing)
  let chi1 = if n1.is_some() && ca1.is_some() && cb1.is_some() {
    calc_dihedral(
      n1.unwrap().get_coord(),
      ca1.unwrap().get_coord(),
      cb1.unwrap().get_coord(),
      sg1c,
    )
  } else { 0.0 }
  let chi2 = if ca1.is_some() && cb1.is_some() {
    calc_dihedral(
      ca1.unwrap().get_coord(),
      cb1.unwrap().get_coord(),
      sg1c,
      sg2c,
    )
  } else { 0.0 }
  let chi_ss = if cb1.is_some() && cb2.is_some() {
    calc_dihedral(
      cb1.unwrap().get_coord(),
      sg1c,
      sg2c,
      cb2.unwrap().get_coord(),
    )
  } else { 0.0 }
  let chi2p = if ca2.is_some() && cb2.is_some() {
    calc_dihedral(
      sg1c,
      sg2c,
      cb2.unwrap().get_coord(),
      ca2.unwrap().get_coord(),
    )
  } else { 0.0 }
  let chi1p = if n2.is_some() && ca2.is_some() && cb2.is_some() {
    calc_dihedral(
      sg2c,
      cb2.unwrap().get_coord(),
      ca2.unwrap().get_coord(),
      n2.unwrap().get_coord(),
    )
  } else { 0.0 }
  let oxidation = if sg_sg_dist < 3.0 {
    "oxidized"
  } else {
    "reduced"
  }
  Some(DisulfideBond::{
    chainid1: res1.chainid,
    resseq1: res1.resseq,
    chainid2: res2.chainid,
    resseq2: res2.resseq,
    sg_sg_distance: sg_sg_dist,
    chi1,
    chi2,
    chi_ss,
    chi2p,
    chi1p,
    oxidation_state: oxidation,
  })
}

///|
/// Format a disulfide bond as a human-readable string.
pub fn format_disulfide_bond(bond : DisulfideBond) -> String {
  let chi1_str = bond.chi1.to_string()
  let chi2_str = bond.chi2.to_string()
  let chi_ss_str = bond.chi_ss.to_string()
  let chi2p_str = bond.chi2p.to_string()
  let chi1p_str = bond.chi1p.to_string()
  let dist_str = bond.sg_sg_distance.to_string()
  "SS-bond: " +
    bond.chainid1.to_string() + bond.resseq1.to_string() + "-" +
    bond.chainid2.to_string() + bond.resseq2.to_string() +
    "  d=" + dist_str + " Å" +
    "  chi1=" + chi1_str +
    "  chi2=" + chi2_str +
    "  chiSS=" + chi_ss_str +
    "  chi2'=" + chi2p_str +
    "  chi1'=" + chi1p_str +
    "  [" + bond.oxidation_state + "]"
}

///|
/// Summary of all disulfide bonds.
pub fn disulfide_bond_summary(bonds : Array[DisulfideBond]) -> String {
  if bonds.length() == 0 {
    return "No disulfide bonds found."
  }
  let lines = Array(capacity=bonds.length() + 1)
  lines.push("Disulfide bonds (" + bonds.length().to_string() + "):")
  for bond in bonds {
    lines.push("  " + format_disulfide_bond(bond))
  }
  lines.join("\n")
}
