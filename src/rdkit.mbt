///|
/// RDKit — cheminformatics core (SMILES parsing + Morgan/ECFP fingerprints).
///
/// A self-contained molecular graph + a SMILES parser for the common core
/// grammar (elements, aromatic lowercase atoms, bonds, branches, ring
/// closures, charges) plus the Morgan (circular / ECFP) fingerprint.

///|
/// A heavy atom in a molecule graph.
pub struct MolAtom {
  symbol : String
  charge : Int
  aromatic : Bool
}

///|
/// A bond between two atom indices.
pub struct MolBond {
  from : Int
  to : Int
  order : Int // 1 = single, 2 = double, 3 = triple
}

///|
/// A molecular graph.
pub struct Molecule {
  atoms : Array[MolAtom]
  bonds : Array[MolBond]
}

///|
/// Number of atoms.
pub fn Molecule::num_atoms(self : Molecule) -> Int {
  self.atoms.length()
}

///|
/// Number of bonds.
pub fn Molecule::num_bonds(self : Molecule) -> Int {
  self.bonds.length()
}

///|
/// Atom symbols in order.
pub fn Molecule::atom_symbols(self : Molecule) -> Array[String] {
  let out : Array[String] = Array(capacity=0)
  for a in self.atoms {
    out.push(a.symbol)
  }
  out
}

///|
/// Bonds as `(from, to, order)` triples.
pub fn Molecule::bond_triples(self : Molecule) -> Array[(Int, Int, Int)] {
  let out : Array[(Int, Int, Int)] = Array(capacity=0)
  for b in self.bonds {
    out.push((b.from, b.to, b.order))
  }
  out
}

///|
fn mol_is_upper(c : UInt16) -> Bool {
  c >= 65 && c <= 90
}

///|
fn mol_is_lower(c : UInt16) -> Bool {
  c >= 97 && c <= 122
}

///|
fn mol_is_digit(c : UInt16) -> Bool {
  c >= 48 && c <= 57
}

///|
fn mol_digit_value(c : UInt16) -> Int {
  (c - 48).to_int()
}

///|
/// Parse an element symbol starting at `i`; returns (symbol, aromatic, length).
/// Lowercase letters are aromatic atoms (c/n/o/s/p); uppercase starts a normal
/// element (1 or 2 letters), and a bracket `[...]` is handled inline by taking
/// the first element symbol inside.
fn mol_parse_element(s : String, i : Int) -> (String, Bool, Int) {
  let n = s.length()
  let c = s.unsafe_get(i)
  if mol_is_lower(c) {
    // 芳香原子（小写）
    let sym = mol_upper_symbol(c)
    (sym, true, 1)
  } else if mol_is_upper(c) {
    // 两字母元素（Cl/Br/Si...）
    if i + 1 < n && mol_is_lower(s.unsafe_get(i + 1)) {
      let sb = StringBuilder()
      sb.write_char(c.unsafe_to_char())
      sb.write_char(s.unsafe_get(i + 1).unsafe_to_char())
      (sb.to_string(), false, 2)
    } else {
      let sb = StringBuilder()
      sb.write_char(c.unsafe_to_char())
      (sb.to_string(), false, 1)
    }
  } else if c.to_int() == 91 {
    // '[' 括号内：提取第一个大写元素（2 字母）
    let mut k = i + 1
    while k < n && s.unsafe_get(k).to_int() != 93 {
      if mol_is_upper(s.unsafe_get(k)) {
        // 尝试两字母
        if k + 1 < n && mol_is_lower(s.unsafe_get(k + 1)) {
          let sb = StringBuilder()
          sb.write_char(s.unsafe_get(k).unsafe_to_char())
          sb.write_char(s.unsafe_get(k + 1).unsafe_to_char())
          return (sb.to_string(), false, k + 2 - i)
        }
        let sb = StringBuilder()
        sb.write_char(s.unsafe_get(k).unsafe_to_char())
        return (sb.to_string(), false, k + 1 - i)
      }
      k = k + 1
    }
    ("C", false, 1)
  } else {
    ("C", false, 1)
  }
}

///|
fn mol_upper_symbol(c : UInt16) -> String {
  let upper = if c >= 97 && c <= 122 { c - 32 } else { c }
  let sb = StringBuilder()
  sb.write_char(upper.unsafe_to_char())
  sb.to_string()
}

///|
/// Parse a SMILES string into a molecular graph (core grammar).
pub fn parse_smiles(smiles : String) -> Molecule {
  let atoms : Array[MolAtom] = Array(capacity=0)
  let bonds : Array[MolBond] = Array(capacity=0)
  let ring : Map[Int, Int] = Map([], capacity=0)
  let branch_stack : Array[Int] = Array(capacity=0)
  let mut prev = -1
  let mut pending_order = 1

  let n = smiles.length()
  let mut i = 0
  while i < n {
    let c = smiles.unsafe_get(i)
    let ci = c.to_int()
    if mol_is_upper(c) || mol_is_lower(c) || ci == 91 {
      let (symbol, aromatic, len) = mol_parse_element(smiles, i)
      let atom_idx = atoms.length()
      atoms.push(MolAtom::{ symbol, charge: 0, aromatic })
      if prev >= 0 {
        bonds.push(MolBond::{ from: prev, to: atom_idx, order: pending_order })
      }
      prev = atom_idx
      pending_order = 1
      i = i + len
    } else if ci == 40 {
      // (
      branch_stack.push(prev)
      i = i + 1
    } else if ci == 41 {
      // )
      if branch_stack.length() > 0 {
        prev = branch_stack[branch_stack.length() - 1]
        let _ = branch_stack.pop()
      }
      i = i + 1
    } else if ci == 61 {
      // =
      pending_order = 2
      i = i + 1
    } else if ci == 35 {
      // #
      pending_order = 3
      i = i + 1
    } else if ci == 45 {
      // -
      pending_order = 1
      i = i + 1
    } else if mol_is_digit(c) {
      let rn = mol_digit_value(c)
      match ring.get(rn) {
        Some(open_atom) => {
          if prev >= 0 && prev != open_atom {
            bonds.push(MolBond::{ from: prev, to: open_atom, order: pending_order })
          }
          ring.remove(rn)
        }
        None => {
          if prev >= 0 {
            ring.set(rn, prev)
          }
        }
      }
      pending_order = 1
      i = i + 1
    } else {
      i = i + 1
    }
  }
  Molecule::{ atoms, bonds }
}


///|
/// FNV-1a 32-bit hash over a string's code units, reduced to a non-negative
/// 31-bit Int (Int64 intermediate avoids overflow).
fn mol_string_hash(s : String) -> Int {
  let mut h : Int64 = 2166136261L
  let n = s.length()
  let mut i = 0
  while i < n {
    h = h ^ s.unsafe_get(i).to_int().to_int64()
    h = h * 16777619L
    i = i + 1
  }
  h.to_int() & 0x7fffffff
}

///|
fn mol_insertion_sort_strings(arr : Array[String]) -> Array[String] {
  let out = arr.copy()
  let mut i = 1
  while i < out.length() {
    let key = out[i]
    let mut j = i - 1
    while j >= 0 && mol_str_gt(out[j], key) {
      out[j + 1] = out[j]
      j = j - 1
    }
    out[j + 1] = key
    i = i + 1
  }
  out
}

///|
fn mol_str_gt(a : String, b : String) -> Bool {
  // 按字节序比较（短者为 MT 前缀？这里用最短优先的简单 lexicographic）
  let na = a.length()
  let nb = b.length()
  let m = if na < nb { na } else { nb }
  let mut i = 0
  while i < m {
    let ca = a.unsafe_get(i)
    let cb = b.unsafe_get(i)
    if ca > cb {
      return true
    } else if ca < cb {
      return false
    }
    i = i + 1
  }
  na > nb
}

///|
/// Morgan (circular / ECFP) fingerprint.
///
/// Returns the multiset of environment hashes: one hash per atom at every
/// radius from 0 up to `radius`. Identical molecular graphs yield identical
/// multisets; different (sub)structures diverge as the radius grows.
pub fn morgan_fingerprint(mol : Molecule, radius : Int) -> Array[Int] {
  let n = mol.atoms.length()
  if n == 0 {
    return []
  }

  // 邻接表（邻居 + 键序）
  let neighbors : Array[Array[(Int, Int)]] = Array(capacity=0)
  for _k = 0; _k < n; _k = _k + 1 {
    neighbors.push([])
  }
  for b in mol.bonds {
    neighbors[b.from].push((b.to, b.order))
    neighbors[b.to].push((b.from, b.order))
  }

  // 初始 invariant 字符串
  let mut atom_inv : Array[String] = Array(capacity=0)
  for k = 0; k < n; k = k + 1 {
    let a = mol.atoms[k]
    let deg = neighbors[k].length()
    let arom = if a.aromatic { "A" } else { "S" }
    atom_inv.push(
      a.symbol + "|" + deg.to_string() + "|" + a.charge.to_string() + "|" + arom,
    )
  }

  let fingerprint : Array[Int] = Array(capacity=0)
  for k = 0; k < n; k = k + 1 {
    fingerprint.push(mol_string_hash(atom_inv[k]))
  }

  let mut r = 0
  while r < radius {
    let new_inv : Array[String] = Array(capacity=0)
    for k = 0; k < n; k = k + 1 {
      let parts : Array[String] = Array(capacity=0)
      for (nb, order) in neighbors[k] {
        parts.push(order.to_string() + ":" + atom_inv[nb])
      }
      let sorted = mol_insertion_sort_strings(parts)
      let sb = StringBuilder()
      sb.write_string(atom_inv[k])
      sb.write_char('(')
      for p in sorted {
        sb.write_string(p)
        sb.write_char(',')
      }
      sb.write_char(')')
      let sig = sb.to_string()
      new_inv.push(sig)
      fingerprint.push(mol_string_hash(sig))
    }
    atom_inv = new_inv
    r = r + 1
  }

  fingerprint
}
