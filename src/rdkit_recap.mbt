/// rdkit_recap.mbt — RECAP decomposition (RDKit gap fill, item A4)
///
/// Retrosynthetic combinatorial analysis procedure on the bond-order
/// graph: a subset of the RECAP rules (amide, ester, amine, urea, ether,
/// olefin, aromatic C-N) is evaluated on every bond; matched bonds are cut
/// and the resulting fragments are returned as canonical SMILES.

///|
/// A RECAP fragment with the rule that produced it.
pub struct RdkitRecapFragment {
  smiles : String
  rule : String
}

///|
fn recap_adj(mol : Molecule) -> Array[Array[(Int32, Int32)]] {
  let n = mol.num_atoms()
  let neighbors : Array[Array[(Int32, Int32)]] = Array(capacity=0)
  for _k = 0; _k < n; _k = _k + 1 {
    neighbors.push([])
  }
  for k = 0; k < mol.bonds.length(); k = k + 1 {
    let b = mol.bonds[k]
    if b.from < n && b.to < n {
      neighbors[b.from].push((b.to, k))
      neighbors[b.to].push((b.from, k))
    }
  }
  neighbors
}

///|
fn recap_has_double_o(
  mol : Molecule,
  atom : Int32,
  neighbors : Array[Array[(Int32, Int32)]],
) -> Bool {
  for (nb, _) in neighbors[atom] {
    if mol.atoms[nb].symbol == "O" {
      if mol_bond_order(mol, atom, nb) == 2 {
        return true
      }
    }
  }
  false
}

///|
fn recap_n_neighbor_symbol(
  mol : Molecule,
  atom : Int32,
  neighbors : Array[Array[(Int32, Int32)]],
) -> Int32 {
  let mut c = 0
  for (nb, _) in neighbors[atom] {
    if mol.atoms[nb].symbol == "N" {
      c = c + 1
    }
  }
  c
}

///|
/// Classify a bond against the implemented RECAP rules (returns "" when
/// no rule matches).
fn recap_classify_bond(
  mol : Molecule,
  bond_idx : Int32,
  neighbors : Array[Array[(Int32, Int32)]],
) -> String {
  let b = mol.bonds[bond_idx]
  let a = mol.atoms[b.from]
  let c = mol.atoms[b.to]
  let sym_a = a.symbol
  let sym_b = c.symbol
  if b.order == 2 &&
    sym_a == "C" &&
    sym_b == "C" &&
    !mol_is_in_ring(mol, b.from) &&
    !mol_is_in_ring(mol, b.to) {
    return "olefin"
  }
  if b.order == 1 &&
    ((sym_a == "C" && sym_b == "N") || (sym_a == "N" && sym_b == "C")) {
    let c_idx = if sym_a == "C" { b.from } else { b.to }
    let c_has_oco = recap_has_double_o(mol, c_idx, neighbors)
    if c_has_oco && recap_n_neighbor_symbol(mol, c_idx, neighbors) >= 2 {
      return "urea"
    }
    if c_has_oco {
      return "amide"
    }
    if mol_is_in_ring(mol, b.from) && mol_is_in_ring(mol, b.to) {
      return "aromatic_C_N"
    }
    return "amine"
  }
  if b.order == 1 && sym_a == "C" && sym_b == "O" {
    // C-O where the C carries a =O -> ester (or carboxylic acid pattern)
    if recap_has_double_o(mol, b.from, neighbors) {
      return "ester"
    }
  }
  if b.order == 1 && sym_a == "O" && sym_b == "C" {
    if recap_has_double_o(mol, b.to, neighbors) {
      return "ester"
    }
  }
  if b.order == 1 && sym_a == "O" && sym_b == "O" {
    return ""
  }
  if b.order == 1 &&
    sym_a == "O" &&
    sym_b == "C" &&
    neighbors[b.from].length() == 2 {
    // ether: O bonded to two carbons
    let mut both_c = true
    for (nb, _) in neighbors[b.from] {
      if mol.atoms[nb].symbol != "C" {
        both_c = false
      }
    }
    if both_c {
      return "ether"
    }
  }
  if b.order == 1 &&
    sym_a == "C" &&
    sym_b == "O" &&
    neighbors[b.to].length() == 2 {
    let mut both_c = true
    for (nb, _) in neighbors[b.to] {
      if mol.atoms[nb].symbol != "C" {
        both_c = false
      }
    }
    if both_c {
      return "ether"
    }
  }
  ""
}

///|
/// Extract the connected component containing `seed` when the given bond
/// index set is removed. Returns the member atom indices.
fn recap_component(
  mol : Molecule,
  seed : Int32,
  cut_bonds : Set[Int32],
) -> Array[Int32] {
  let n = mol.num_atoms()
  let visited : Array[Bool] = Array::make(n, false)
  let out : Array[Int32] = Array(capacity=0)
  let queue : Array[Int32] = Array(capacity=0)
  queue.push(seed)
  visited[seed] = true
  let mut head = 0
  while head < queue.length() {
    let cur = queue[head]
    head = head + 1
    out.push(cur)
    for b in mol.bonds {
      let other = if b.from == cur {
        b.to
      } else if b.to == cur {
        b.from
      } else {
        -1
      }
      // determine bond index by identity scan
      if other >= 0 && !visited[other] {
        // bond index check
        if recap_bond_in_cut(mol, cur, other, cut_bonds) {
          continue
        }
        visited[other] = true
        queue.push(other)
      }
    }
  }
  out
}

///|
fn recap_bond_in_cut(
  mol : Molecule,
  a : Int32,
  b : Int32,
  cut_bonds : Set[Int32],
) -> Bool {
  for i = 0; i < mol.bonds.length(); i = i + 1 {
    if cut_bonds.contains(i) {
      let bond = mol.bonds[i]
      let matches = (bond.from == a && bond.to == b) ||
        (bond.from == b && bond.to == a)
      if matches {
        return true
      }
    }
  }
  false
}

///|
/// Build a sub-molecule from an atom set (preserving symbols, charges and
/// aromatic flags, keeping bonds whose both ends are in the set).
fn recap_sub_mol(mol : Molecule, atoms : Array[Int32]) -> Molecule {
  let in_set : Array[Bool] = Array::make(mol.num_atoms(), false)
  for a in atoms {
    in_set[a] = true
  }
  let index_map : Map[Int32, Int32] = Map([], capacity=0)
  let new_atoms : Array[MolAtom] = Array(capacity=0)
  for a in atoms {
    index_map[a] = new_atoms.length()
    new_atoms.push(MolAtom::{
      symbol: mol.atoms[a].symbol,
      charge: mol.atoms[a].charge,
      aromatic: mol.atoms[a].aromatic,
    })
  }
  let new_bonds : Array[MolBond] = Array(capacity=0)
  for b in mol.bonds {
    if in_set[b.from] && in_set[b.to] {
      new_bonds.push(MolBond::{
        from: index_map[b.from],
        to: index_map[b.to],
        order: b.order,
      })
    }
  }
  Molecule::{ atoms: new_atoms, bonds: new_bonds, }
}

///|
/// RECAP decomposition: cut every bond matched by the implemented rules
/// and return the unique fragment SMILES with their rule names.
pub fn rdkit_recap_decompose(mol : Molecule) -> Array[RdkitRecapFragment] {
  let out : Array[RdkitRecapFragment] = Array(capacity=0)
  let seen : Set[String] = Set([])
  let neighbors = recap_adj(mol)
  for i = 0; i < mol.bonds.length(); i = i + 1 {
    let rule = recap_classify_bond(mol, i, neighbors)
    if rule.length() == 0 {
      continue
    }
    let cut : Set[Int32] = Set([i])
    let comp_a = recap_component(mol, mol.bonds[i].from, cut)
    let comp_b = recap_component(mol, mol.bonds[i].to, cut)
    let frag_pairs = [
      (recap_sub_mol(mol, comp_a), comp_a.length()),
      (recap_sub_mol(mol, comp_b), comp_b.length()),
    ]
    for pair in frag_pairs {
      let (sub, _) = pair
      let smi = mol_to_smiles(mol_perceive_aromaticity(sub))
      if smi.length() > 0 && !seen.contains(smi) {
        seen.add(smi)
        out.push(RdkitRecapFragment::{ smiles: smi, rule, })
      }
    }
  }
  out
}
