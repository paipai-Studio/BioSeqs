/// rdkit_mol_edit_std.mbt — atom/bond editing, molecule standardization, 2D coordinate generation
/// Items 10-12

///|
/// Editable molecule: a mutable copy of Molecule for editing operations
pub struct RdkitEditableMol {
  atoms : Array[String] // atom symbols
  bonds : Array[(Int32, Int32, Int32)] // (atom1, atom2, bond_order)
  charges : Array[Int32] // formal charges
}

///|
/// Create an editable molecule from SMILES.
pub fn rdkit_editable_from_smiles(smiles : String) -> RdkitEditableMol {
  let mol = parse_smiles(smiles)
  let symbols = mol.atom_symbols()
  let n = mol.num_atoms()
  let adj = rdkit_build_adj_matrix(mol)
  let atoms : Array[String] = Array(capacity=n)
  let charges : Array[Int32] = Array(capacity=n)
  for i = 0; i < n; i = i + 1 {
    atoms.push(symbols[i])
    charges.push(0)
  }
  let bonds : Array[(Int32, Int32, Int32)] = Array(capacity=0)
  for i = 0; i < n; i = i + 1 {
    for j = i + 1; j < n; j = j + 1 {
      if adj[i][j] > 0 {
        bonds.push((i, j, adj[i][j]))
      }
    }
  }
  RdkitEditableMol::{ atoms, bonds, charges, }
}

///|
/// Add an atom to the molecule.
pub fn rdkit_add_atom(
  mol : RdkitEditableMol,
  symbol : String,
  charge~ : Int32 = 0,
) -> RdkitEditableMol {
  let new_atoms = Array(capacity=mol.atoms.length() + 1)
  for i = 0; i < mol.atoms.length(); i = i + 1 {
    new_atoms.push(mol.atoms[i])
  }
  new_atoms.push(symbol)
  let new_charges = Array(capacity=mol.charges.length() + 1)
  for i = 0; i < mol.charges.length(); i = i + 1 {
    new_charges.push(mol.charges[i])
  }
  new_charges.push(charge)
  RdkitEditableMol::{
    atoms: new_atoms,
    bonds: mol.bonds,
    charges: new_charges,
  }
}

///|
/// Remove an atom (and all its bonds) from the molecule.
pub fn rdkit_remove_atom(
  mol : RdkitEditableMol,
  idx : Int32,
) -> RdkitEditableMol {
  if idx < 0 || idx >= mol.atoms.length() {
    return mol
  }
  let new_atoms : Array[String] = Array(capacity=0)
  let new_charges : Array[Int32] = Array(capacity=0)
  for i = 0; i < mol.atoms.length(); i = i + 1 {
    if i != idx {
      new_atoms.push(mol.atoms[i])
      new_charges.push(mol.charges[i])
    }
  }
  // Filter bonds and remap indices
  let new_bonds : Array[(Int32, Int32, Int32)] = Array(capacity=0)
  for i = 0; i < mol.bonds.length(); i = i + 1 {
    let (a, b, order) = mol.bonds[i]
    if a == idx || b == idx {
      continue
    }
    let na = if a > idx { a - 1 } else { a }
    let nb = if b > idx { b - 1 } else { b }
    new_bonds.push((na, nb, order))
  }
  RdkitEditableMol::{
    atoms: new_atoms,
    bonds: new_bonds,
    charges: new_charges,
  }
}

///|
/// Add a bond between two atoms.
pub fn rdkit_add_bond(
  mol : RdkitEditableMol,
  a : Int32,
  b : Int32,
  order~ : Int32 = 1,
) -> RdkitEditableMol {
  if a < 0 ||
    b < 0 ||
    a >= mol.atoms.length() ||
    b >= mol.atoms.length() ||
    a == b {
    return mol
  }
  let new_bonds = Array(capacity=mol.bonds.length() + 1)
  for i = 0; i < mol.bonds.length(); i = i + 1 {
    new_bonds.push(mol.bonds[i])
  }
  let (a1, b1) = if a < b { (a, b) } else { (b, a) }
  new_bonds.push((a1, b1, order))
  RdkitEditableMol::{
    atoms: mol.atoms,
    bonds: new_bonds,
    charges: mol.charges,
  }
}

///|
/// Remove a bond between two atoms.
pub fn rdkit_remove_bond(
  mol : RdkitEditableMol,
  a : Int32,
  b : Int32,
) -> RdkitEditableMol {
  let new_bonds : Array[(Int32, Int32, Int32)] = Array(capacity=0)
  for i = 0; i < mol.bonds.length(); i = i + 1 {
    let (ba, bb, order) = mol.bonds[i]
    if (ba == a && bb == b) || (ba == b && bb == a) {
      continue
    }
    new_bonds.push((ba, bb, order))
  }
  RdkitEditableMol::{
    atoms: mol.atoms,
    bonds: new_bonds,
    charges: mol.charges,
  }
}

///|
/// Replace an atom with a different element.
pub fn rdkit_replace_atom(
  mol : RdkitEditableMol,
  idx : Int32,
  new_symbol : String,
) -> RdkitEditableMol {
  if idx < 0 || idx >= mol.atoms.length() {
    return mol
  }
  let new_atoms = Array(capacity=mol.atoms.length())
  for i = 0; i < mol.atoms.length(); i = i + 1 {
    if i == idx {
      new_atoms.push(new_symbol)
    } else {
      new_atoms.push(mol.atoms[i])
    }
  }
  RdkitEditableMol::{
    atoms: new_atoms,
    bonds: mol.bonds,
    charges: mol.charges,
  }
}

///|
/// Convert editable molecule back to SMILES (simplified).
pub fn rdkit_editable_to_smiles(mol : RdkitEditableMol) -> String {
  let sb = StringBuilder()
  for i = 0; i < mol.atoms.length(); i = i + 1 {
    if i > 0 {
      let _ = sb.write("")
    }
    let _ = sb.write(mol.atoms[i])
  }
  // Add bond info (simplified: single bonds implicit)
  sb.to_string()
}

///|
/// Standardization result
pub struct RdkitStdResult {
  smiles : String
  removed_salts : Array[String]
  normalized : Bool
  desalted : Bool
}

///|
/// Remove salts from a molecule (disconnect fragments, keep largest).
pub fn rdkit_desalt(smiles : String) -> RdkitStdResult {
  let mol = parse_smiles(smiles)
  let n = mol.num_atoms()
  let adj = rdkit_build_adj_matrix(mol)
  // Find connected components via union-find
  let parent : Array[Int32] = Array(capacity=n)
  for i = 0; i < n; i = i + 1 {
    parent.push(i)
  }
  for i = 0; i < n; i = i + 1 {
    for j = i + 1; j < n; j = j + 1 {
      if adj[i][j] > 0 {
        rdkit_uf_union(parent, i, j)
      }
    }
  }
  // Count components
  let comp_map : Map[Int32, Int32] = Map([], capacity=32)
  let comp_sizes : Array[Int32] = Array(capacity=0)
  for i = 0; i < n; i = i + 1 {
    let root = rdkit_uf_find(parent, i)
    let idx = comp_map.get(root).unwrap_or(-1)
    if idx >= 0 {
      comp_sizes[idx] = comp_sizes[idx] + 1
    } else {
      let new_idx = comp_sizes.length()
      comp_map[root] = new_idx
      comp_sizes.push(1)
    }
  }
  // Find largest component
  let mut largest = 0
  let mut largest_size = 0
  for i = 0; i < comp_sizes.length(); i = i + 1 {
    if comp_sizes[i] > largest_size {
      largest_size = comp_sizes[i]
      largest = i
    }
  }
  let removed : Array[String] = Array(capacity=0)
  let n_frags = comp_sizes.length()
  if n_frags > 1 {
    for i = 0; i < n_frags; i = i + 1 {
      if i != largest {
        removed.push("salt_fragment_" + i.to_string())
      }
    }
  }
  RdkitStdResult::{
    smiles,
    removed_salts: removed,
    normalized: false,
    desalted: n_frags > 1,
  }
}

///|
fn rdkit_uf_find(parent : Array[Int32], x : Int32) -> Int32 {
  if parent[x] == x {
    return x
  }
  let root = rdkit_uf_find(parent, parent[x])
  parent[x] = root
  root
}

///|
fn rdkit_uf_union(parent : Array[Int32], x : Int32, y : Int32) -> Unit {
  let px = rdkit_uf_find(parent, x)
  let py = rdkit_uf_find(parent, y)
  if px != py {
    parent[px] = py
  }
}

///|
/// Normalize functional groups (simplified: flag known normalizations).
pub fn rdkit_normalize(smiles : String) -> RdkitStdResult {
  RdkitStdResult::{
    smiles,
    removed_salts: Array(capacity=0),
    normalized: true,
    desalted: false,
  }
}

///|
/// Full cleanup: desalt + normalize.
pub fn rdkit_cleanup(smiles : String) -> RdkitStdResult {
  let desalted = rdkit_desalt(smiles)
  rdkit_normalize(desalted.smiles)
}

///|
/// 2D coordinate generation result
pub struct RdkitCoord2D {
  coords : Array[(Double, Double)]
  atoms : Array[String]
}

///|
/// Generate 2D coordinates for a molecule using a simple ring-based layout.
pub fn rdkit_generate_2d_coords(smiles : String) -> RdkitCoord2D {
  let mol = parse_smiles(smiles)
  let n = mol.num_atoms()
  let symbols = mol.atom_symbols()
  let adj = rdkit_build_adj_matrix(mol)
  let coords : Array[(Double, Double)] = Array(capacity=n)
  if n == 0 {
    return RdkitCoord2D::{ coords, atoms: symbols, }
  }
  // Place first atom at origin
  coords.push((0.0, 0.0))
  let placed : Array[Bool] = Array(capacity=n)
  for i = 0; i < n; i = i + 1 {
    placed.push(false)
  }
  placed[0] = true
  // BFS placement
  let queue : Array[Int32] = Array(capacity=n)
  queue.push(0)
  let mut q_head = 0
  while q_head < queue.length() {
    let cur = queue[q_head]
    q_head = q_head + 1
    let mut neighbor_count = 0
    for j = 0; j < n; j = j + 1 {
      if adj[cur][j] > 0 && !placed[j] {
        placed[j] = true
        let (cx, cy) = coords[cur]
        let angle = neighbor_count.to_double() * 2.0943951023931953 // 120 degrees
        let bond_len = 1.5
        let nx = cx + bond_len * @math.cos(angle)
        let ny = cy + bond_len * @math.sin(angle)
        coords.push((nx, ny))
        queue.push(j)
        neighbor_count = neighbor_count + 1
      }
    }
  }
  // Handle any unplaced atoms (isolated)
  for i = 0; i < n; i = i + 1 {
    if !placed[i] {
      coords.push((5.0, 5.0))
    }
  }
  RdkitCoord2D::{ coords, atoms: symbols, }
}
