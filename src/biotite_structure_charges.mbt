/// Biotite `structure.charges` — partial charges via the PEOE
/// algorithm of Gasteiger-Marsili.
///
/// Faithful port of the RDKit implementation used by Biotite
/// (`biotite.structure.charges.partial_charges`, which delegates to
/// the Rust port of RDKit's Gasteiger charges):
/// - per-(element, hybridization) parameters (A, B, C) from RDKit's
///   GasteigerParams (chi(q) = A + B*q + C*q^2),
/// - ionization potentials ionX = A + B + C (hydrogen: IONXH = 20.02),
/// - charge transfer along bonds towards the higher electronegativity
///   atom, divided by the ionX of the *receiving* atom,
/// - damping factor 0.5 halved every iteration (DAMP / DAMP_SCALE),
/// - unknown elements fall back to the X parameters (0, 0, 0).
/// The reference doctest values (fluoromethane, 1 step:
/// C 0.115 / F -0.175 / H 0.020) are covered by tests.

///|
/// RDKit Gasteiger (A, B, C) parameters per element and hybridization.
/// Returns None for elements without parameters (caller falls back to
/// the X row: 0, 0, 0).
pub fn biotite_poeo_params(
  element : String,
  mode : String,
) -> (Double, Double, Double)? {
  match (element, mode) {
    ("H", "*") => Some((7.17, 6.24, -0.56))
    ("C", "sp3") => Some((7.98, 9.18, 1.88))
    ("C", "sp2") => Some((8.79, 9.32, 1.51))
    ("C", "sp") => Some((10.39, 9.45, 0.73))
    ("N", "sp3") => Some((11.54, 10.82, 1.36))
    ("N", "sp2") => Some((12.87, 11.15, 0.85))
    ("N", "sp") => Some((15.68, 11.7, -0.27))
    ("O", "sp3") => Some((14.18, 12.92, 1.39))
    ("O", "sp2") => Some((17.07, 13.79, 0.47))
    ("F", "sp3") => Some((14.66, 13.85, 2.31))
    ("Cl", "sp3") => Some((11.00, 9.69, 1.35))
    ("Br", "sp3") => Some((10.08, 8.47, 1.16))
    ("I", "sp3") => Some((9.9, 7.96, 0.96))
    ("S", "sp3") => Some((10.14, 9.13, 1.38))
    ("S", "so") => Some((10.14, 9.13, 1.38))
    ("S", "so2") => Some((12.00, 10.81, 1.20))
    ("S", "sp2") => Some((10.88, 9.49, 1.33))
    ("P", "sp3") => Some((8.90, 8.24, 0.96))
    ("P", "sp2") => Some((9.665, 8.530, 0.735))
    ("Si", "sp3") => Some((7.300, 6.567, 0.657))
    ("Si", "sp2") => Some((7.905, 6.748, 0.443))
    ("Si", "sp") => Some((9.065, 7.027, -0.002))
    ("B", "sp3") => Some((5.980, 6.820, 1.605))
    ("B", "sp2") => Some((6.420, 6.807, 1.322))
    ("Be", "sp3") => Some((3.845, 6.755, 3.165))
    ("Be", "sp2") => Some((4.005, 6.725, 3.035))
    ("Mg", "sp2") => Some((3.565, 5.572, 2.197))
    ("Mg", "sp3") => Some((3.300, 5.587, 2.447))
    ("Mg", "sp") => Some((4.040, 5.472, 1.823))
    ("Al", "sp3") => Some((5.375, 4.953, 0.867))
    ("Al", "sp2") => Some((5.795, 5.020, 0.695))
    ("X", "*") => Some((0.0, 0.0, 0.0))
    _ => None
  }
}

///|
/// Hybridization mode of an atom based on its bond orders
/// (RDKit uses the graph hybridization; we infer it from bond types).
/// Sulfur follows RDKit's special case: count the oxygen neighbors.
fn biotite_poeo_mode(
  element : String,
  adjacency : Array[Array[(Int32, BondType)]],
  elements : Array[String],
  i : Int32,
) -> String {
  // RDKit: hydrogen without hybridization info uses the "*" mode.
  if element == "H" {
    return "*"
  }
  let mut double_bond = false
  let mut triple_bond = false
  for (_j, bt) in adjacency[i] {
    match bt {
      BondType::Double | BondType::Aromatic => double_bond = true
      BondType::Triple => triple_bond = true
      _ => ()
    }
  }
  if triple_bond {
    return "sp"
  }
  if double_bond {
    return "sp2"
  }
  if element == "S" {
    // RDKit: sulfur without hybridization info is classified by its
    // oxygen neighbor count (2 -> so2, 1 -> so, else sp3).
    let mut n_o = 0
    for (j, _bt) in adjacency[i] {
      if elements[j] == "O" {
        n_o = n_o + 1
      }
    }
    if n_o == 2 {
      return "so2"
    }
    if n_o == 1 {
      return "so"
    }
  }
  "sp3"
}

///|
/// Resolve the (A, B, C) parameters of an atom; unknown combinations
/// fall back to the X row (0, 0, 0) like RDKit's non-throwing lookup.
fn biotite_poeo_lookup(
  element : String,
  mode : String,
) -> (Double, Double, Double) {
  match biotite_poeo_params(element, mode) {
    Some(p) => p
    None => (0.0, 0.0, 0.0)
  }
}

///|
/// Compute Gasteiger-Marsili (PEOE) partial charges for an AtomArray.
///
/// Mirrors `biotite.structure.charges.partial_charges` (RDKit
/// implementation): `iteration_step_num` defaults to 6 ("described as
/// sufficient" by Gasteiger and Marsili), `formal_charges` seeds the
/// initial charges (defaults to all zero). Structures without bonds
/// return the formal charges unchanged. Total charge is conserved by
/// the symmetric transfer, so no normalization is applied.
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
    // No bonds along which charge transfer can occur.
    return q
  }
  // Per-atom parameters and ionization potentials.
  let a_p : Array[Double] = Array(capacity=n)
  let b_p : Array[Double] = Array(capacity=n)
  let c_p : Array[Double] = Array(capacity=n)
  let ion_x : Array[Double] = Array(capacity=n)
  for i = 0; i < n; i = i + 1 {
    let mode = biotite_poeo_mode(arr.elements[i], adjacency, arr.elements, i)
    let (a, b, c) = biotite_poeo_lookup(arr.elements[i], mode)
    a_p.push(a)
    b_p.push(b)
    c_p.push(c)
    if arr.elements[i] == "H" {
      // RDKit uses the constant IONXH for hydrogen.
      ion_x.push(20.02)
    } else {
      ion_x.push(a + b + c)
    }
  }
  let damp_start = 0.5
  let damp_scale = 0.5
  let energ : Array[Double] = Array(capacity=n)
  for _i = 0; _i < n; _i = _i + 1 {
    energ.push(0.0)
  }
  let mut damp = damp_start
  for _it = 0; _it < iteration_step_num; _it = _it + 1 {
    // Electronegativity of every atom under the current charge.
    for i = 0; i < n; i = i + 1 {
      energ[i] = a_p[i] + q[i] * (b_p[i] + c_p[i] * q[i])
    }
    // Charge migrates towards the higher electronegativity neighbor;
    // the transfer is divided by the ionX of the receiving atom.
    // (energ is frozen for the whole round, so update order is
    // irrelevant; RDKit updates charges in place exactly like this.)
    for i = 0; i < n; i = i + 1 {
      let mut dq = 0.0
      for (j, _bt) in adjacency[i] {
        let dx = energ[j] - energ[i]
        // dx >= 0: this atom receives -> denominator is its own ionX;
        // dx < 0: the neighbor receives -> denominator is the
        // neighbor's ionX.
        let denom = if dx < 0.0 { ion_x[j] } else { ion_x[i] }
        if denom > 0.0 {
          dq = dq + dx / denom
        }
      }
      q[i] = q[i] + damp * dq
    }
    damp = damp * damp_scale
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
pub fn _biotite_structure_charges_dummy() -> Unit {
  let charges = biotite_partial_charges_peoe(
    AtomArray::new(),
    BondList::new(n_atoms=0),
  )
  let _ = charges
  let _ = biotite_poeo_params("C", "sp3")
  let _ = biotite_poeo_lookup("Fe", "sp3")
}
