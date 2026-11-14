///|
/// Structure dihedral / torsion angles — backbone (φ/ψ/ω) and
/// side-chain (χ¹..χ⁴) for AtomArray.
///
/// Inspired by Biotite's structure.dihedral module.

// ============================================================
// Generic 4-point dihedral on AtomArray
// ============================================================

///|
/// Compute dihedral angle (degrees) defined by four atom indices.
/// Returns None if any index is out of bounds.
pub fn AtomArray::dihedral(
  self : AtomArray,
  a : Int32,
  b : Int32,
  c : Int32,
  d : Int32,
) -> Double? {
  let n = self.len()
  if a < 0 || b < 0 || c < 0 || d < 0 {
    return None
  }
  if a >= n || b >= n || c >= n || d >= n {
    return None
  }
  let p1 = self.coords[a]
  let p2 = self.coords[b]
  let p3 = self.coords[c]
  let p4 = self.coords[d]
  Some(vector_dihedral_angle_degrees(p1, p2, p3, p4))
}

// ============================================================
// Backbone dihedrals (phi / psi / omega)
// ============================================================

///|
/// Compute backbone φ (phi) for residue at residue index `ri`.
/// φ = C(i-1) — N(i) — CA(i) — C(i)
/// Requires the previous residue to be in the same chain.
pub fn AtomArray::backbone_phi(self : AtomArray, ri : Int32) -> Double? {
  let starts = self.get_residue_starts()
  let n_res = starts.length()
  if ri <= 0 || ri >= n_res {
    return None
  }
  let prev_s = starts[ri - 1]
  let prev_e = starts[ri]
  let cur_s = starts[ri]
  let cur_e = if ri + 1 < starts.length() { starts[ri + 1] } else { self.len() }
  if self.chain_ids[prev_s] != self.chain_ids[cur_s] {
    return None
  }
  let c_prev = find_atom_in_range(self, prev_s, prev_e, "C")
  let n_cur = find_atom_in_range(self, cur_s, cur_e, "N")
  let ca_cur = find_atom_in_range(self, cur_s, cur_e, "CA")
  let c_cur = find_atom_in_range(self, cur_s, cur_e, "C")
  match (c_prev, n_cur, ca_cur, c_cur) {
    (Some(cp), Some(nc), Some(cac), Some(cc)) => self.dihedral(cp, nc, cac, cc)
    _ => None
  }
}

///|
/// Compute backbone ψ (psi) for residue at residue index `ri`.
/// ψ = N(i) — CA(i) — C(i) — N(i+1)
pub fn AtomArray::backbone_psi(self : AtomArray, ri : Int32) -> Double? {
  let starts = self.get_residue_starts()
  let n_res = starts.length()
  if ri < 0 || ri >= n_res - 1 {
    return None
  }
  let cur_s = starts[ri]
  let cur_e = starts[ri + 1]
  let next_s = starts[ri + 1]
  let next_e = if ri + 2 < starts.length() {
    starts[ri + 2]
  } else {
    self.len()
  }
  if self.chain_ids[cur_s] != self.chain_ids[next_s] {
    return None
  }
  let n_cur = find_atom_in_range(self, cur_s, cur_e, "N")
  let ca_cur = find_atom_in_range(self, cur_s, cur_e, "CA")
  let c_cur = find_atom_in_range(self, cur_s, cur_e, "C")
  let n_next = find_atom_in_range(self, next_s, next_e, "N")
  match (n_cur, ca_cur, c_cur, n_next) {
    (Some(nc), Some(cac), Some(cc), Some(nn)) => self.dihedral(nc, cac, cc, nn)
    _ => None
  }
}

///|
/// Compute backbone ω (omega) for the peptide bond between residue ri and ri+1.
/// ω = CA(i) — C(i) — N(i+1) — CA(i+1)
pub fn AtomArray::backbone_omega(self : AtomArray, ri : Int32) -> Double? {
  let starts = self.get_residue_starts()
  let n_res = starts.length()
  if ri < 0 || ri >= n_res - 1 {
    return None
  }
  let cur_s = starts[ri]
  let cur_e = starts[ri + 1]
  let next_s = starts[ri + 1]
  let next_e = if ri + 2 < starts.length() {
    starts[ri + 2]
  } else {
    self.len()
  }
  if self.chain_ids[cur_s] != self.chain_ids[next_s] {
    return None
  }
  let ca_cur = find_atom_in_range(self, cur_s, cur_e, "CA")
  let c_cur = find_atom_in_range(self, cur_s, cur_e, "C")
  let n_next = find_atom_in_range(self, next_s, next_e, "N")
  let ca_next = find_atom_in_range(self, next_s, next_e, "CA")
  match (ca_cur, c_cur, n_next, ca_next) {
    (Some(cac), Some(cc), Some(nn), Some(can)) =>
      self.dihedral(cac, cc, nn, can)
    _ => None
  }
}

///|
/// Find an atom index by name within [start, end).
fn find_atom_in_range(
  arr : AtomArray,
  start : Int32,
  end : Int32,
  name : String,
) -> Int32? {
  let mut i = start
  while i < end && i < arr.len() {
    if arr.names[i] == name {
      return Some(i)
    }
    i = i + 1
  }
  None
}

// ============================================================
// Side-chain χ dihedrals (χ¹..χ⁴)
// ============================================================

///|
/// Residue-specific side-chain χ dihedral atom-name tables.
/// Each χ entry = (atomA, atomB, atomC, atomD).
fn chi_atoms(
  res_name : String,
  chi_idx : Int32,
) -> (String, String, String, String)? {
  let c = chi_idx
  match res_name {
    "GLY" | "ALA" => None
    "VAL" => if c == 1 { Some(("N", "CA", "CB", "CG1")) } else { None }
    "SER" => if c == 1 { Some(("N", "CA", "CB", "OG")) } else { None }
    "THR" => if c == 1 { Some(("N", "CA", "CB", "OG1")) } else { None }
    "CYS" => if c == 1 { Some(("N", "CA", "CB", "SG")) } else { None }
    "PRO" => if c == 1 { Some(("N", "CA", "CB", "CG")) } else { None }
    "PHE" =>
      if c == 1 {
        Some(("N", "CA", "CB", "CG"))
      } else if c == 2 {
        Some(("CA", "CB", "CG", "CD1"))
      } else {
        None
      }
    "TYR" =>
      if c == 1 {
        Some(("N", "CA", "CB", "CG"))
      } else if c == 2 {
        Some(("CA", "CB", "CG", "CD1"))
      } else {
        None
      }
    "TRP" =>
      if c == 1 {
        Some(("N", "CA", "CB", "CG"))
      } else if c == 2 {
        Some(("CA", "CB", "CG", "CD1"))
      } else {
        None
      }
    "HIS" =>
      if c == 1 {
        Some(("N", "CA", "CB", "CG"))
      } else if c == 2 {
        Some(("CA", "CB", "CG", "ND1"))
      } else {
        None
      }
    "ASP" =>
      if c == 1 {
        Some(("N", "CA", "CB", "CG"))
      } else if c == 2 {
        Some(("CA", "CB", "CG", "OD1"))
      } else {
        None
      }
    "ASN" =>
      if c == 1 {
        Some(("N", "CA", "CB", "CG"))
      } else if c == 2 {
        Some(("CA", "CB", "CG", "OD1"))
      } else {
        None
      }
    "LEU" =>
      if c == 1 {
        Some(("N", "CA", "CB", "CG"))
      } else if c == 2 {
        Some(("CA", "CB", "CG", "CD1"))
      } else {
        None
      }
    "ILE" =>
      if c == 1 {
        Some(("N", "CA", "CB", "CG1"))
      } else if c == 2 {
        Some(("CA", "CB", "CG1", "CD1"))
      } else {
        None
      }
    "GLU" =>
      if c == 1 {
        Some(("N", "CA", "CB", "CG"))
      } else if c == 2 {
        Some(("CA", "CB", "CG", "CD"))
      } else if c == 3 {
        Some(("CB", "CG", "CD", "OE1"))
      } else {
        None
      }
    "GLN" =>
      if c == 1 {
        Some(("N", "CA", "CB", "CG"))
      } else if c == 2 {
        Some(("CA", "CB", "CG", "CD"))
      } else if c == 3 {
        Some(("CB", "CG", "CD", "OE1"))
      } else {
        None
      }
    "MET" =>
      if c == 1 {
        Some(("N", "CA", "CB", "CG"))
      } else if c == 2 {
        Some(("CA", "CB", "CG", "SD"))
      } else if c == 3 {
        Some(("CB", "CG", "SD", "CE"))
      } else {
        None
      }
    "ARG" =>
      if c == 1 {
        Some(("N", "CA", "CB", "CG"))
      } else if c == 2 {
        Some(("CA", "CB", "CG", "CD"))
      } else if c == 3 {
        Some(("CB", "CG", "CD", "NE"))
      } else if c == 4 {
        Some(("CG", "CD", "NE", "CZ"))
      } else {
        None
      }
    "LYS" =>
      if c == 1 {
        Some(("N", "CA", "CB", "CG"))
      } else if c == 2 {
        Some(("CA", "CB", "CG", "CD"))
      } else if c == 3 {
        Some(("CB", "CG", "CD", "CE"))
      } else if c == 4 {
        Some(("CG", "CD", "CE", "NZ"))
      } else {
        None
      }
    _ => None
  }
}

///|
/// Compute a side-chain χ dihedral for residue `ri`.
/// `chi_idx` is 1-based (1 = χ¹, 2 = χ², ...).
/// Returns None if the residue has fewer χ angles or any atom is missing.
pub fn AtomArray::dihedral_side_chain(
  self : AtomArray,
  ri : Int32,
  chi_idx : Int32,
) -> Double? {
  let starts = self.get_residue_starts()
  let n_res = starts.length()
  if ri < 0 || ri >= n_res {
    return None
  }
  let s = starts[ri]
  let e = if ri + 1 < starts.length() { starts[ri + 1] } else { self.len() }
  let rn = self.res_names[s]
  match chi_atoms(rn, chi_idx) {
    None => None
    Some((an, bn, cn, dn)) => {
      let a = find_atom_in_range(self, s, e, an)
      let b = find_atom_in_range(self, s, e, bn)
      let c = find_atom_in_range(self, s, e, cn)
      let d = find_atom_in_range(self, s, e, dn)
      match (a, b, c, d) {
        (Some(ai), Some(bi), Some(ci), Some(di)) =>
          self.dihedral(ai, bi, ci, di)
        _ => None
      }
    }
  }
}

///|
/// Compute all side-chain χ dihedrals for residue `ri`.
/// Returns an array of (chi_index, angle_degrees) pairs for each available χ.
pub fn AtomArray::all_side_chain_dihedrals(
  self : AtomArray,
  ri : Int32,
) -> Array[(Int32, Double)] {
  let result : Array[(Int32, Double)] = []
  let starts = self.get_residue_starts()
  let n_res = starts.length()
  if ri < 0 || ri >= n_res {
    return result
  }
  let s = starts[ri]
  let rn = self.res_names[s]
  // Try χ¹..χ⁴ and collect those that work
  let mut ci = 1
  while ci <= 4 {
    match self.dihedral_side_chain(ri, ci) {
      Some(a) => result.push((ci, a))
      None => ()
    }
    ci = ci + 1
  }
  result
}
