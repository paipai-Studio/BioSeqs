///|
/// BondList — explicit bond connectivity table with bond types.
///
/// Biotite's `BondList` stores bonds as pairs of atom indices with
/// an associated `BondType`. This module provides an equivalent for
/// MoonBit, supporting graph traversal and connectivity queries.

///|
/// Bond type classification.
pub enum BondType {
  Single
  Double
  Triple
  Quadruple
  Aromatic
  Any // unspecified
}

///|
/// Get bond type as a string.
pub fn BondType::to_string(self : BondType) -> String {
  match self {
    BondType::Single => "SINGLE"
    BondType::Double => "DOUBLE"
    BondType::Triple => "TRIPLE"
    BondType::Quadruple => "QUADRUPLE"
    BondType::Aromatic => "AROMATIC"
    BondType::Any => "ANY"
  }
}

///|
/// A single bond between two atoms.
pub struct Bond {
  atom1 : Int32
  atom2 : Int32
  bond_type : BondType
}

///|
/// Construct a Bond.
pub fn Bond::new(
  atom1~ : Int32,
  atom2~ : Int32,
  bond_type~ : BondType = BondType::Single,
) -> Bond {
  Bond::{ atom1, atom2, bond_type, }
}

///|
/// BondList — stores all bonds in a structure.
pub struct BondList {
  bonds : Array[Bond]
  n_atoms : Int32
}

///|
/// Create an empty BondList for a structure with n_atoms.
pub fn BondList::new(n_atoms~ : Int32) -> BondList {
  BondList::{ bonds: [], n_atoms, }
}

///|
/// Add a bond.
pub fn BondList::add(
  self : BondList,
  atom1 : Int32,
  atom2 : Int32,
  bond_type~ : BondType = BondType::Single,
) -> BondList {
  let new_bonds = Array(capacity=self.bonds.length() + 1)
  for b in self.bonds {
    new_bonds.push(b)
  }
  new_bonds.push(Bond::new(atom1~, atom2~, bond_type~))
  BondList::{ bonds: new_bonds, n_atoms: self.n_atoms, }
}

///|
/// Number of bonds.
pub fn BondList::len(self : BondList) -> Int32 {
  self.bonds.length()
}

///|
/// Get all bonds.
pub fn BondList::get_bonds(self : BondList) -> Array[Bond] {
  self.bonds
}

///|
/// Find all atoms bonded to a given atom.
pub fn BondList::get_neighbors(self : BondList, atom : Int32) -> Array[Int32] {
  let neighbors : Array[Int32] = []
  for b in self.bonds {
    if b.atom1 == atom {
      neighbors.push(b.atom2)
    }
    if b.atom2 == atom {
      neighbors.push(b.atom1)
    }
  }
  neighbors
}

///|
/// Check if two atoms are directly bonded.
pub fn BondList::are_bonded(
  self : BondList,
  atom1 : Int32,
  atom2 : Int32,
) -> Bool {
  for b in self.bonds {
    if (b.atom1 == atom1 && b.atom2 == atom2) ||
      (b.atom1 == atom2 && b.atom2 == atom1) {
      return true
    }
  }
  false
}

///|
/// Get the bond type between two atoms, if any.
pub fn BondList::get_bond_type(
  self : BondList,
  atom1 : Int32,
  atom2 : Int32,
) -> BondType? {
  for b in self.bonds {
    if (b.atom1 == atom1 && b.atom2 == atom2) ||
      (b.atom1 == atom2 && b.atom2 == atom1) {
      return Some(b.bond_type)
    }
  }
  None
}

///|
/// Find all atoms connected to a root atom (transitive closure).
///
/// Uses BFS to traverse the bond graph.
pub fn BondList::find_connected(self : BondList, root : Int32) -> Array[Int32] {
  let visited = Array::make(self.n_atoms, false)
  let result : Array[Int32] = []
  let queue : Array[Int32] = [root]
  visited[root] = true
  while queue.length() > 0 {
    let current = queue[0]
    ignore(queue.remove(0))
    result.push(current)
    let neighbors = self.get_neighbors(current)
    for n in neighbors {
      if n >= 0 && n < self.n_atoms && !visited[n] {
        visited[n] = true
        queue.push(n)
      }
    }
  }
  result
}

///|
/// Count the number of connected components (molecules).
pub fn BondList::count_components(self : BondList) -> Int32 {
  let visited = Array::make(self.n_atoms, false)
  let mut count = 0
  for i in 0..<self.n_atoms {
    if !visited[i] {
      count = count + 1
      let connected = self.find_connected(i)
      for c in connected {
        if c >= 0 && c < self.n_atoms {
          visited[c] = true
        }
      }
    }
  }
  count
}

///|
/// Find rotatable bonds (single bonds between non-terminal, non-ring atoms).
///
/// A bond is rotatable if it is a single bond, both atoms have at least
/// 2 neighbors, and it is not part of a ring (simplified: not checked
/// for rings in this version).
pub fn BondList::find_rotatable_bonds(self : BondList) -> Array[Bond] {
  let rotatable : Array[Bond] = []
  for b in self.bonds {
    if b.bond_type is BondType::Single {
      let n1 = self.get_neighbors(b.atom1).length()
      let n2 = self.get_neighbors(b.atom2).length()
      // Both atoms must have at least 2 neighbors (non-terminal)
      if n1 >= 2 && n2 >= 2 {
        rotatable.push(b)
      }
    }
  }
  rotatable
}

///|
/// Generate a bond list from atom coordinates using distance cutoffs.
///
/// Two atoms are bonded if their distance is within [min_dist, max_dist]
/// based on their van der Waals radii.
pub fn bondlist_from_coordinates(
  coords : Array[Vector3],
  elements : Array[String],
) -> BondList {
  let n = coords.length()
  let mut bl = BondList::new(n_atoms=n)
  for i in 0..<n {
    for j in (i + 1)..<n {
      let dist = coords[i].distance(coords[j])
      // Simple distance-based cutoff
      let r1 = simple_vdw_radius(elements[i])
      let r2 = simple_vdw_radius(elements[j])
      let max_dist = r1 + r2 + 0.4
      let min_dist = 0.4
      if dist >= min_dist && dist <= max_dist {
        bl = bl.add(i, j)
      }
    }
  }
  bl
}

///|
/// Simplified van der Waals radius lookup.
fn simple_vdw_radius(element : String) -> Double {
  match element {
    "H" => 1.2
    "C" => 1.7
    "N" => 1.55
    "O" => 1.52
    "S" => 1.8
    "P" => 1.8
    "F" => 1.47
    "CL" => 1.75
    "BR" => 1.85
    "I" => 1.98
    _ => 1.5
  }
}

///|
fn _bond_type_unused_constructors() -> Unit {
  ignore(BondType::Double)
  ignore(BondType::Triple)
  ignore(BondType::Quadruple)
  ignore(BondType::Aromatic)
  ignore(BondType::Any)
}
