/// Biotite `structure.charges` — partial charges for AtomArray structures.
///
/// Implements a Gasteiger electronegativity-equalization scheme over the
/// atom-array bond topology (Gasteiger/Marsili 1980 A/B/C parameters,
/// 23 refinement iterations, matching Biotite's default). Atoms outside
/// the organic subset (H, C, N, O, S, P) receive no charge, mirroring
/// Biotite's NaN convention (None here).

///|
/// Gasteiger A/B/C parameters per element/hybridization.
/// Hybridization is inferred from double/triple/aromatic bonds.
fn biotite_gasteiger_params(
  element : String,
  sp2 : Bool,
) -> (Double, Double, Double) {
  match element {
    "H" => (7.17, -6.14, -0.56)
    "C" => if sp2 { (8.81, -9.98, -0.94) } else { (7.94, -9.22, -0.85) }
    "N" => if sp2 { (15.68, -13.18, -1.93) } else { (11.15, -12.13, -1.45) }
    "O" => if sp2 { (17.68, -16.03, -2.53) } else { (14.18, -14.77, -2.01) }
    "S" => (8.03, -6.33, -0.75)
    "P" => (8.54, -7.86, -0.96)
    _ => (0.0, 0.0, 0.0)
  }
}

///|
/// Whether an element belongs to the organic subset that receives
/// partial charges.
fn biotite_is_organic_element(element : String) -> Bool {
  match element {
    "H" | "C" | "N" | "O" | "S" | "P" => true
    _ => false
  }
}

///|
/// Prepared charge calculation state (template structure).
pub struct BiotitePartialCharges {
  elements : Array[String]
  // adjacency[i] = [(j, bond_type), ...]
  adjacency : Array[Array[(Int32, BondType)]]
  formal : Array[Double]
  organic : Array[Bool]
}

///|
/// Prepare a partial-charge calculation over an atom array.
///
/// `formal_charges` supplies per-atom integer formal charges; defaults
/// to all zero. (The name is suffixed `_gasteiger` to distinguish it
/// from the simpler one-shot `biotite_partial_charges` in
/// biotite_formal_charges.mbt.)
pub fn biotite_partial_charges_gasteiger(
  arr : AtomArray,
  bonds : BondList,
  formal_charges~ : Array[Int32]? = None,
) -> BiotitePartialCharges {
  let n = arr.len()
  let adjacency : Array[Array[(Int32, BondType)]] = Array(capacity=n)
  for _i = 0; _i < n; _i = _i + 1 {
    adjacency.push([])
  }
  for b in bonds.get_bonds() {
    if b.atom1 >= 0 && b.atom1 < n && b.atom2 >= 0 && b.atom2 < n {
      adjacency[b.atom1].push((b.atom2, b.bond_type))
      adjacency[b.atom2].push((b.atom1, b.bond_type))
    }
  }
  let organic : Array[Bool] = Array(capacity=n)
  for i = 0; i < n; i = i + 1 {
    organic.push(biotite_is_organic_element(arr.elements[i]))
  }
  let formal : Array[Double] = Array(capacity=n)
  for i = 0; i < n; i = i + 1 {
    match formal_charges {
      Some(fc) =>
        if i < fc.length() {
          formal.push(fc[i].to_double())
        } else {
          formal.push(0.0)
        }
      None => formal.push(0.0)
    }
  }
  BiotitePartialCharges::{ elements: arr.elements, adjacency, formal, organic, }
}

///|
/// Whether the atom participates in a double/triple/aromatic bond.
fn biotite_atom_is_sp2(
  adjacency : Array[Array[(Int32, BondType)]],
  i : Int32,
) -> Bool {
  for (j, bt) in adjacency[i] {
    let _ = j
    match bt {
      BondType::Double | BondType::Triple | BondType::Aromatic => return true
      _ => ()
    }
  }
  false
}

///|
/// Public bond-type constructor for external packages
/// (BondType is a read-only enum; tests and callers need this helper).
pub fn biotite_bond_double() -> BondType {
  BondType::Double
}

///|
/// Run the Gasteiger iteration and return per-atom charges.
///
/// Organic atoms get a Double; all others get None (Biotite NaN).
pub fn BiotitePartialCharges::compute_charge(
  self : BiotitePartialCharges,
) -> Array[Double?] {
  let n = self.elements.length()
  let q : Array[Double] = Array(capacity=n)
  for i = 0; i < n; i = i + 1 {
    q.push(if self.organic[i] { self.formal[i] } else { 0.0 })
  }
  let max_iter = 23
  for _iter = 0; _iter < max_iter; _iter = _iter + 1 {
    // Electronegativity per atom under the current charge.
    let chi : Array[Double] = Array(capacity=n)
    for i = 0; i < n; i = i + 1 {
      if self.organic[i] {
        let (a, b, c) = biotite_gasteiger_params(
          self.elements[i],
          biotite_atom_is_sp2(self.adjacency, i),
        )
        chi.push(a + b * q[i] + c * q[i] * q[i])
      } else {
        chi.push(0.0)
      }
    }
    // Charge migrates from lower to higher electronegativity.
    let new_q : Array[Double] = q.copy()
    for i = 0; i < n; i = i + 1 {
      if !self.organic[i] {
        continue
      }
      let neighbors = self.adjacency[i]
      let mut delta = 0.0
      for (j, _bt) in neighbors {
        if !self.organic[j] {
          continue
        }
        let (ai, _bi, _ci) = biotite_gasteiger_params(
          self.elements[i],
          biotite_atom_is_sp2(self.adjacency, i),
        )
        let (_aj, _bj, _cj) = biotite_gasteiger_params(
          self.elements[j],
          biotite_atom_is_sp2(self.adjacency, j),
        )
        let denom = ai + ai
        if denom > 0.0 {
          // Equalization step towards the neighbor's electronegativity.
          delta = delta + (chi[j] - chi[i]) / denom
        }
      }
      let deg = if neighbors.length() > 0 { neighbors.length() } else { 1 }
      let new_val = q[i] + delta / deg.to_double()
      // Damping clamp to keep the iteration stable.
      new_q[i] = if new_val > 1.5 {
        1.5
      } else if new_val < -1.5 {
        -1.5
      } else {
        new_val
      }
    }
    for i = 0; i < n; i = i + 1 {
      q[i] = new_q[i]
    }
  }
  // Normalize so that the total equals the sum of formal charges
  // (organic atoms only).
  let mut total = 0.0
  let mut formal_total = 0.0
  let mut n_organic = 0
  for i = 0; i < n; i = i + 1 {
    if self.organic[i] {
      total = total + q[i]
      formal_total = formal_total + self.formal[i]
      n_organic = n_organic + 1
    }
  }
  if n_organic > 0 {
    let correction = (formal_total - total) / n_organic.to_double()
    for i = 0; i < n; i = i + 1 {
      if self.organic[i] {
        q[i] = q[i] + correction
      }
    }
  }
  let out : Array[Double?] = Array(capacity=n)
  for i = 0; i < n; i = i + 1 {
    out.push(if self.organic[i] { Some(q[i]) } else { None })
  }
  out
}

///|
pub fn _biotite_structure_charges_dummy() -> Unit {
  let pc = biotite_partial_charges_gasteiger(
    AtomArray::new(),
    BondList::new(n_atoms=0),
  )
  let _ = pc.compute_charge()
  let _ = pc.elements
  let _ = pc.organic
  let _ = pc.formal
  let _ = pc.adjacency
}
