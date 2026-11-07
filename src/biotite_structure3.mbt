///|
/// biotite_structure3.mbt — Biotite round 3: structure layer
/// 6 features: AtomArrayStack, connect_via_residue_names, connect_via_distances,
/// find_connected, find_rotatable_bonds, infer_bond_types.

// #22 AtomArrayStack — multi-model structure
pub struct Biotite3Atom {
  atom_name : String
  element : String
  res_name : String
  res_id : Int32
  chain_id : String
  coord : (Double, Double, Double)
  b_factor : Double
  occupancy : Double
  hetero : Bool
}

///|
pub fn biotite3_atom(
  name : String,
  element : String,
  res_name : String,
  res_id : Int32,
  chain_id : String,
  coord : (Double, Double, Double),
) -> Biotite3Atom {
  {
    atom_name: name,
    element,
    res_name,
    res_id,
    chain_id,
    coord,
    b_factor: 0.0,
    occupancy: 1.0,
    hetero: false,
  }
}

///|
pub struct Biotite3AtomArrayStack {
  atoms : Array[Biotite3Atom]
  coord_stack : Array[Array[(Double, Double, Double)]]
  n_models : Int32
}

///|
pub fn biotite3_atom_array_stack(
  atoms : Array[Biotite3Atom],
  coord_stack : Array[Array[(Double, Double, Double)]],
) -> Biotite3AtomArrayStack {
  { atoms, coord_stack, n_models: coord_stack.length(), }
}

///|
pub fn biotite3_atom_array_stack_get_model(
  stack : Biotite3AtomArrayStack,
  model_idx : Int32,
) -> Array[Biotite3Atom] {
  if model_idx < 0 || model_idx >= stack.coord_stack.length() {
    return []
  }
  let coords = stack.coord_stack[model_idx]
  let result = Array(capacity=stack.atoms.length())
  for i in 0..<stack.atoms.length() {
    let a = stack.atoms[i]
    let c = if i < coords.length() { coords[i] } else { a.coord }
    result.push({
      atom_name: a.atom_name,
      element: a.element,
      res_name: a.res_name,
      res_id: a.res_id,
      chain_id: a.chain_id,
      coord: c,
      b_factor: a.b_factor,
      occupancy: a.occupancy,
      hetero: a.hetero,
    })
  }
  result
}

///|
pub fn biotite3_atom_array_stack_concatenate(
  stacks : Array[Biotite3AtomArrayStack],
) -> Biotite3AtomArrayStack {
  let all_atoms = Array(capacity=0)
  let all_coords = Array(capacity=0)
  let mut n_models = 0
  for s in stacks {
    for a in s.atoms {
      all_atoms.push(a)
    }
    for c in s.coord_stack {
      all_coords.push(c)
    }
    n_models = n_models + s.n_models
  }
  { atoms: all_atoms, coord_stack: all_coords, n_models, }
}

// #23 connect_via_distances — bond inference from distance thresholds

///|
pub struct Biotite3Bond {
  atom1 : Int32
  atom2 : Int32
  bond_type : Int32
}

///|
pub fn biotite3_connect_via_distances(
  atoms : Array[Biotite3Atom],
) -> Array[Biotite3Bond] {
  let bonds = Array(capacity=0)
  let covalent_radii : Map[String, Double] = Map([
    ("H", 0.31),
    ("C", 0.76),
    ("N", 0.71),
    ("O", 0.66),
    ("F", 0.57),
    ("P", 1.07),
    ("S", 1.05),
    ("CL", 1.02),
    ("BR", 1.20),
    ("I", 1.39),
    ("FE", 1.32),
    ("ZN", 1.22),
    ("MG", 1.30),
    ("CA", 1.74),
    ("MN", 1.35),
    ("CU", 1.32),
    ("NA", 1.54),
    ("K", 2.03),
  ])
  let tolerance = 0.45
  for i in 0..<atoms.length() {
    let _start_j = i + 1
    for j in _start_j..<atoms.length() {
      let r1 = covalent_radii.get(atoms[i].element).unwrap_or(0.7)
      let r2 = covalent_radii.get(atoms[j].element).unwrap_or(0.7)
      let max_dist = r1 + r2 + tolerance
      let (x1, y1, z1) = atoms[i].coord
      let (x2, y2, z2) = atoms[j].coord
      let dx = x1 - x2
      let dy = y1 - y2
      let dz = z1 - z2
      let dist = (dx * dx + dy * dy + dz * dz).sqrt()
      if dist < max_dist && dist > 0.5 {
        bonds.push(Biotite3Bond::{ atom1: i, atom2: j, bond_type: 1, })
      }
    }
  }
  bonds
}

// #24 connect_via_residue_names — bond inference from residue templates

///|
pub fn biotite3_connect_via_residue_names(
  atoms : Array[Biotite3Atom],
) -> Array[Biotite3Bond] {
  let bonds = Array(capacity=0)
  let mut i = 0
  while i < atoms.length() {
    let mut j = i
    while j < atoms.length() &&
          atoms[j].res_id == atoms[i].res_id &&
          atoms[j].chain_id == atoms[i].chain_id {
      j = j + 1
    }
    for a in i..<j {
      let _start_b = a + 1
      for b in _start_b..<j {
        let n1 = atoms[a].atom_name
        let n2 = atoms[b].atom_name
        let is_bond = biotite3_is_standard_bond(n1, n2)
        if is_bond {
          bonds.push(Biotite3Bond::{ atom1: a, atom2: b, bond_type: 1, })
        }
      }
    }
    for a in i..<j {
      if atoms[a].atom_name == "C" && j < atoms.length() {
        for b in j..<atoms.length() {
          if atoms[b].atom_name == "N" &&
            atoms[b].res_id == atoms[a].res_id + 1 &&
            atoms[b].chain_id == atoms[a].chain_id {
            bonds.push(Biotite3Bond::{ atom1: a, atom2: b, bond_type: 1, })
            break
          }
        }
      }
    }
    i = j
  }
  bonds
}

///|
fn biotite3_is_standard_bond(n1 : String, n2 : String) -> Bool {
  let protein_backbone = [
    ("N", "CA"),
    ("CA", "C"),
    ("C", "O"),
    ("CA", "CB"),
    ("CB", "CG"),
    ("CG", "CD"),
    ("CD", "CE"),
    ("CE", "NZ"),
    ("CG", "ND1"),
    ("ND1", "CE1"),
    ("CE1", "NE2"),
    ("NE2", "CD2"),
    ("CD2", "CG"),
    ("CB", "OG"),
    ("CB", "OG1"),
    ("CB", "SG"),
    ("CG", "OD1"),
    ("OD1", "ND2"),
    ("CG", "ND1"),
    ("ND1", "CE1"),
    ("CE1", "NE2"),
    ("NE2", "CD2"),
    ("CD2", "CG"),
    ("CA", "OXT"),
  ]
  for i in 0..<protein_backbone.length() {
    let (a, b) = protein_backbone[i]
    if (n1 == a && n2 == b) || (n1 == b && n2 == a) {
      return true
    }
  }
  false
}

// #25 find_connected — all atoms connected to a root

///|
pub fn biotite3_find_connected(
  bonds : Array[Biotite3Bond],
  n_atoms : Int32,
  root : Int32,
) -> Array[Int32] {
  let visited : Array[Bool] = Array(capacity=n_atoms)
  for _i in 0..<n_atoms {
    visited.push(false)
  }
  let queue = Array(capacity=0)
  queue.push(root)
  visited[root] = true
  let result = Array(capacity=0)
  while queue.length() > 0 {
    let current = queue[0]
    ignore(queue.remove(0))
    result.push(current)
    for i in 0..<bonds.length() {
      let neighbor = if bonds[i].atom1 == current {
        bonds[i].atom2
      } else if bonds[i].atom2 == current {
        bonds[i].atom1
      } else {
        -1
      }
      if neighbor >= 0 && neighbor < n_atoms && !visited[neighbor] {
        visited[neighbor] = true
        queue.push(neighbor)
      }
    }
  }
  result
}

// #26 find_rotatable_bonds — find bonds that can rotate

///|
pub fn biotite3_find_rotatable_bonds(
  bonds : Array[Biotite3Bond],
  atoms : Array[Biotite3Atom],
) -> Array[Biotite3Bond] {
  let result = Array(capacity=0)
  for i in 0..<bonds.length() {
    let a1 = bonds[i].atom1
    let a2 = bonds[i].atom2
    if a1 >= atoms.length() || a2 >= atoms.length() {
      continue
    }
    if atoms[a1].element == "H" || atoms[a2].element == "H" {
      continue
    }
    if bonds[i].bond_type != 1 {
      continue
    }
    let n1 = biotite3_count_bonds(bonds, a1)
    let n2 = biotite3_count_bonds(bonds, a2)
    if n1 >= 2 && n2 >= 2 {
      result.push(bonds[i])
    }
  }
  result
}

///|
fn biotite3_count_bonds(bonds : Array[Biotite3Bond], atom_idx : Int32) -> Int32 {
  let mut count = 0
  for i in 0..<bonds.length() {
    if bonds[i].atom1 == atom_idx || bonds[i].atom2 == atom_idx {
      count = count + 1
    }
  }
  count
}

// #27 infer_bond_types — infer bond order from geometry

///|
pub fn biotite3_infer_bond_types(
  bonds : Array[Biotite3Bond],
  atoms : Array[Biotite3Atom],
) -> Array[Biotite3Bond] {
  let result = Array(capacity=bonds.length())
  for i in 0..<bonds.length() {
    let a1 = bonds[i].atom1
    let a2 = bonds[i].atom2
    if a1 >= atoms.length() || a2 >= atoms.length() {
      result.push(bonds[i])
      continue
    }
    let (x1, y1, z1) = atoms[a1].coord
    let (x2, y2, z2) = atoms[a2].coord
    let dx = x1 - x2
    let dy = y1 - y2
    let dz = z1 - z2
    let dist = (dx * dx + dy * dy + dz * dz).sqrt()
    let e1 = atoms[a1].element
    let e2 = atoms[a2].element
    let bond_type = biotite3_classify_bond(e1, e2, dist)
    result.push({ atom1: a1, atom2: a2, bond_type, })
  }
  result
}

///|
fn biotite3_classify_bond(e1 : String, e2 : String, dist : Double) -> Int32 {
  let pair = if e1 < e2 { e1 + "-" + e2 } else { e2 + "-" + e1 }
  match pair {
    "C-C" =>
      if dist < 1.20 {
        3
      } else if dist < 1.35 {
        2
      } else if dist < 1.60 {
        1
      } else {
        0
      }
    "C-N" =>
      if dist < 1.20 {
        3
      } else if dist < 1.35 {
        2
      } else if dist < 1.55 {
        1
      } else {
        0
      }
    "C-O" =>
      if dist < 1.15 {
        3
      } else if dist < 1.30 {
        2
      } else if dist < 1.50 {
        1
      } else {
        0
      }
    "N-N" =>
      if dist < 1.10 {
        3
      } else if dist < 1.25 {
        2
      } else if dist < 1.40 {
        1
      } else {
        0
      }
    "N-O" => if dist < 1.15 { 2 } else if dist < 1.35 { 1 } else { 0 }
    "O-O" => if dist < 1.25 { 2 } else if dist < 1.50 { 1 } else { 0 }
    "C-S" => if dist < 1.60 { 2 } else if dist < 1.85 { 1 } else { 0 }
    _ => 1
  }
}

///|
pub fn _biotite_structure3_dummy() -> Unit {
  let _ = Biotite3Atom::{
    atom_name: "",
    element: "",
    res_name: "",
    res_id: 0,
    chain_id: "",
    coord: (0.0, 0.0, 0.0),
    b_factor: 0.0,
    occupancy: 0.0,
    hetero: false,
  }
  let _ = Biotite3AtomArrayStack::{ atoms: [], coord_stack: [], n_models: 0, }
  let _ = Biotite3Bond::{ atom1: 0, atom2: 0, bond_type: 0, }
}
