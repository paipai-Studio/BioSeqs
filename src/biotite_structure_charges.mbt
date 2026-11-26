/// Biotite `structure.charges` — partial charges via the PEOE
/// algorithm of Gasteiger-Marsili.
///
/// Faithful port of the official Rust implementation
/// (`src/rust/structure/charges.rs` in biotite-dev/biotite, backing
/// `biotite.structure.charges.partial_charges`):
/// - per-(element, hybridization) parameters (A, B, C) for the nine
///   parametrized elements H/C/N/O/S/F/Cl/Br/I, chi(q) = A + B*q +
///   C*q^2; carbon sp2 shares B/C with sp3,
/// - hybridization from the highest bond order (triple/quadruple ->
///   sp, double/aromatic -> sp2, single -> sp3); atoms whose bonds
///   are all `Any` fall back to the binding-partner count,
/// - atoms without parameters (unknown element or valence state) get
///   NaN electronegativity: their own charge and that of any bonded
///   partner becomes NaN, exactly like the official implementation,
/// - charge transfer along bonds towards the higher electronegativity
///   atom, divided by the ionization potential of the *donating* atom
///   (a + b + c; hydrogen and deuterium use the constant 20.02 eV),
/// - damping factor halved every iteration starting from 0.5.
/// The official doctest values (fluoromethane, 1 step:
/// C 0.115 / F -0.175 / H 0.020) and the 17 molecules of table 3 of
/// the Gasteiger-Marsili publication are covered by tests.

///|
/// Gasteiger (A, B, C) parameters per element and hybridization.
/// Returns None for elements or valence states without parameters
/// (the caller then uses NaN electronegativity like the official
/// implementation).
pub fn biotite_poeo_params(
  element : String,
  mode : String,
) -> (Double, Double, Double)? {
  match (element, mode) {
    ("H", "*") => Some((7.17, 6.24, -0.56))
    ("C", "sp3") => Some((7.98, 9.18, 1.88))
    ("C", "sp2") => Some((8.79, 9.18, 1.88))
    ("C", "sp") => Some((10.39, 9.45, 0.73))
    ("N", "sp3") => Some((11.54, 10.82, 1.36))
    ("N", "sp2") => Some((12.87, 11.15, 0.85))
    ("N", "sp") => Some((15.68, 11.7, -0.27))
    ("O", "sp3") => Some((14.18, 12.92, 1.39))
    ("O", "sp2") => Some((17.07, 13.79, 0.47))
    ("S", "sp3") => Some((10.14, 9.13, 1.38))
    ("F", "sp3") => Some((14.66, 13.85, 2.31))
    ("Cl", "sp3") => Some((11.00, 9.69, 1.35))
    ("Br", "sp3") => Some((10.08, 8.47, 1.16))
    ("I", "sp3") => Some((9.9, 7.96, 0.96))
    _ => None
  }
}

///|
/// Hybridization mode of an atom based on its bond orders: the
/// highest bond order wins (triple/quadruple -> sp, double/aromatic
/// -> sp2, single -> sp3), `Any` bonds contribute no information.
/// Atoms whose bonds are all `Any` fall back to the binding-partner
/// count (official `hybridization_from_partner_count`); an empty
/// result means "no parameters available".
fn biotite_poeo_mode(
  element : String,
  adjacency : Array[Array[(Int32, BondType)]],
  i : Int32,
) -> String {
  // Hydrogen is always treated as the s-orbital "*".
  if element == "H" {
    return "*"
  }
  let mut rank = 0
  let mut has_info = false
  for (_j, bt) in adjacency[i] {
    match bt {
      BondType::Single => {
        has_info = true
        if rank < 1 {
          rank = 1
        }
      }
      BondType::Double | BondType::Aromatic => {
        has_info = true
        if rank < 2 {
          rank = 2
        }
      }
      BondType::Triple | BondType::Quadruple => {
        has_info = true
        rank = 3
      }
      BondType::Any => ()
    }
  }
  if has_info {
    return match rank {
      1 => "sp3"
      2 => "sp2"
      _ => "sp"
    }
  }
  // All bonds are `Any`: infer the hybridization from the number of
  // binding partners (can be wrong for charged atoms, like upstream).
  let partners = adjacency[i].length()
  match element {
    "H" => if partners == 1 { "sp3" } else { "" }
    "C" =>
      match partners {
        4 => "sp3"
        3 => "sp2"
        2 => "sp"
        _ => ""
      }
    "N" =>
      match partners {
        3 | 4 => "sp3"
        2 => "sp2"
        1 => "sp"
        _ => ""
      }
    "O" =>
      match partners {
        2 => "sp3"
        1 => "sp2"
        _ => ""
      }
    "S" => if partners == 2 { "sp3" } else { "" }
    "F" | "Cl" | "Br" | "I" => if partners == 1 { "sp3" } else { "" }
    _ => ""
  }
}

///|
/// Compute Gasteiger-Marsili (PEOE) partial charges for an AtomArray.
///
/// Mirrors `biotite.structure.charges.partial_charges` (official Rust
/// implementation): `iteration_step_num` defaults to 6 ("described as
/// sufficient" by Gasteiger and Marsili), `formal_charges` seeds the
/// initial charges (defaults to all zero). Structures without bonds
/// return the formal charges unchanged. Atoms without parameters get
/// NaN charges, which propagate to their bonded partners. Total
/// charge is conserved by the symmetric transfer among parametrized
/// atoms, so no normalization is applied.
pub fn biotite_partial_charges_peoe(
  arr : AtomArray,
  bonds : BondList,
  iteration_step_num~ : Int32 = 6,
  formal_charges~ : Array[Int32]? = None,
) -> Array[Double] {
  let n = arr.len()
  let q : Array[Double] = Array(capacity=n)
  for i = 0; i < n; i = i + 1 {
    match formal_charges {
      Some(fc) =>
        if i < fc.length() {
          q.push(fc[i].to_double())
        } else {
          q.push(0.0)
        }
      None => q.push(0.0)
    }
  }
  // Adjacency with bond types for hybridization inference.
  let adjacency : Array[Array[(Int32, BondType)]] = Array(capacity=n)
  for _i = 0; _i < n; _i = _i + 1 {
    adjacency.push([])
  }
  let mut n_bonds = 0
  for b in bonds.get_bonds() {
    if b.atom1 >= 0 && b.atom1 < n && b.atom2 >= 0 && b.atom2 < n {
      adjacency[b.atom1].push((b.atom2, b.bond_type))
      adjacency[b.atom2].push((b.atom1, b.bond_type))
      n_bonds = n_bonds + 1
    }
  }
  if n_bonds == 0 {
    // No bonds along which charge transfer can occur
    // -> partial charges equal the formal charges (e.g. ions).
    return q
  }
  // Per-atom parameters; NaN for unparametrized elements/valence
  // states (their charge becomes NaN and propagates to partners).
  let nan : Double = @double.not_a_number
  let a_p : Array[Double] = Array(capacity=n)
  let b_p : Array[Double] = Array(capacity=n)
  let c_p : Array[Double] = Array(capacity=n)
  let ion_x : Array[Double] = Array(capacity=n)
  for i = 0; i < n; i = i + 1 {
    let mode = biotite_poeo_mode(arr.elements[i], adjacency, i)
    match biotite_poeo_params(arr.elements[i], mode) {
      Some((a, b, c)) => {
        a_p.push(a)
        b_p.push(b)
        c_p.push(c)
        if arr.elements[i] == "H" || arr.elements[i] == "D" {
          // Electronegativity of positively charged hydrogen.
          ion_x.push(20.02)
        } else {
          ion_x.push(a + b + c)
        }
      }
      None => {
        a_p.push(nan)
        b_p.push(nan)
        c_p.push(nan)
        ion_x.push(nan)
      }
    }
  }
  let energ : Array[Double] = Array(capacity=n)
  for _i = 0; _i < n; _i = _i + 1 {
    energ.push(0.0)
  }
  let mut damp = 0.5
  for _it = 0; _it < iteration_step_num; _it = _it + 1 {
    // Electronegativity of every atom under the current charge;
    // NaN parameters propagate automatically through the formula.
    for i = 0; i < n; i = i + 1 {
      energ[i] = a_p[i] + q[i] * (b_p[i] + c_p[i] * q[i])
    }
    // Iterate over the bonds and transfer charges towards the higher
    // electronegativity atom; the divisor is the ionization potential
    // of the donating (lower electronegativity) atom. A NaN
    // electronegativity on either side marks that side's charge NaN.
    for b in bonds.get_bonds() {
      if b.atom1 < 0 || b.atom1 >= n || b.atom2 < 0 || b.atom2 >= n {
        continue
      }
      let en_i = energ[b.atom1]
      let en_j = energ[b.atom2]
      let i_is_nan = en_i != en_i
      let j_is_nan = en_j != en_j
      if i_is_nan {
        q[b.atom1] = nan
      }
      if j_is_nan {
        q[b.atom2] = nan
      }
      if !i_is_nan && !j_is_nan {
        let denom = if en_j > en_i { ion_x[b.atom1] } else { ion_x[b.atom2] }
        let transfer = (en_j - en_i) / denom * damp
        q[b.atom1] = q[b.atom1] + transfer
        q[b.atom2] = q[b.atom2] - transfer
      }
    }
    damp = damp * 0.5
  }
  q
}

///|
/// Public bond-type constructor for external packages
/// (BondType is a read-only enum; tests and callers need this helper).
pub fn biotite_bond_double() -> BondType {
  BondType::Double
}

///|
/// Public bond-type constructor for external packages (single bond).
pub fn biotite_bond_triple() -> BondType {
  BondType::Triple
}

///|
/// Public bond-type constructor for external packages (single bond).
pub fn biotite_bond_single() -> BondType {
  BondType::Single
}

///|
/// Public bond-type constructor for external packages
/// (unspecified bond type; hybridization then falls back to the
/// binding-partner count).
pub fn biotite_bond_any() -> BondType {
  BondType::Any
}

///|
pub fn _biotite_structure_charges_dummy() -> Unit {
  let charges = biotite_partial_charges_peoe(
    AtomArray::new(),
    BondList::new(n_atoms=0),
  )
  let _ = charges
  let _ = biotite_poeo_params("C", "sp3")
}
