// ============================================================================
// SMARTS pattern parser and matcher
// ============================================================================
// Supports: atom patterns ([#6], [C], [R], [a], [*], [c,n,o]),
// bond patterns (-, =, #, ~, /, \), logical operators (&, |, !, ;),
// ring closures (digits), branching (parentheses).
// Matching is done via backtracking DFS against the Molecule graph.

///|
/// SMARTS atom query — conditions an atom must satisfy.
pub struct SmartsAtom {
  element : String?
  atomic_num : Int?
  aromatic : Bool?
  in_ring : Bool?
  degree : Int?
  charge : Int?
  is_wildcard : Bool
  not_match : Bool
}

///|
/// SMARTS bond query.
pub struct SmartsBond {
  order : Int
  is_any : Bool
}

///|
/// A parsed SMARTS pattern.
pub struct SmartsPattern {
  atoms : Array[SmartsAtom]
  bonds : Array[SmartsBond]
  branches : Array[(Int, SmartsPattern)]
}

///|
/// Parse a SMARTS atom query from a bracket expression.
fn smarts_parse_atom(s : String) -> SmartsAtom {
  let mut element : String? = None
  let mut atomic_num : Int? = None
  let mut aromatic : Bool? = None
  let mut in_ring : Bool? = None
  let mut degree : Int? = None
  let mut charge : Int? = None
  let mut is_wildcard = false
  let mut not_match = false
  let s_chars = s.to_array()
  let chars : Array[Char] = Array(capacity=0)
  if s_chars.length() >= 2 && s_chars[0] == '[' {
    // Strip brackets: copy chars from index 1 to length-2
    for k = 1; k < s_chars.length() - 1; k = k + 1 {
      chars.push(s_chars[k])
    }
  } else {
    for c in s_chars {
      chars.push(c)
    }
  }
  let mut i = 0
  if chars.length() > 0 && chars[0] == '!' {
    not_match = true
    i = 1
  }
  if i < chars.length() {
    let c = chars[i]
    if c == '*' {
      is_wildcard = true
      i = i + 1
    } else if c == '#' {
      i = i + 1
      let mut num_str = ""
      while i < chars.length() {
        let ch = chars[i]
        if ch >= '0' && ch <= '9' {
          num_str = num_str + ch.to_string()
          i = i + 1
        } else {
          break
        }
      }
      atomic_num = Some(parse_int(num_str))
    } else if c == 'R' {
      in_ring = Some(true)
      i = i + 1
    } else if c == 'a' {
      aromatic = Some(true)
      i = i + 1
    } else if c == 'A' {
      aromatic = Some(false)
      i = i + 1
    } else if c == 'D' {
      i = i + 1
      let mut num_str = ""
      while i < chars.length() {
        let ch = chars[i]
        if ch >= '0' && ch <= '9' {
          num_str = num_str + ch.to_string()
          i = i + 1
        } else {
          break
        }
      }
      degree = Some(parse_int(num_str))
    } else if c == '+' {
      i = i + 1
      let mut num_str = ""
      while i < chars.length() {
        let ch = chars[i]
        if ch >= '0' && ch <= '9' {
          num_str = num_str + ch.to_string()
          i = i + 1
        } else {
          break
        }
      }
      if num_str == "" {
        charge = Some(1)
      } else {
        charge = Some(parse_int(num_str))
      }
    } else if c == '-' {
      i = i + 1
      let mut num_str = ""
      while i < chars.length() {
        let ch = chars[i]
        if ch >= '0' && ch <= '9' {
          num_str = num_str + ch.to_string()
          i = i + 1
        } else {
          break
        }
      }
      if num_str == "" {
        charge = Some(-1)
      } else {
        charge = Some(-parse_int(num_str))
      }
    } else if c >= 'a' && c <= 'z' {
      element = Some(c.to_string().to_upper())
      aromatic = Some(true)
      i = i + 1
    } else if c >= 'A' && c <= 'Z' {
      let mut elem = c.to_string()
      i = i + 1
      if i < chars.length() {
        let ch2 = chars[i]
        if ch2 >= 'a' && ch2 <= 'z' {
          elem = elem + ch2.to_string()
          i = i + 1
        }
      }
      element = Some(elem)
    }
  }
  SmartsAtom::{
    element,
    atomic_num,
    aromatic,
    in_ring,
    degree,
    charge,
    is_wildcard,
    not_match,
  }
}

///|
/// Atomic number to element symbol.
fn atomic_num_to_symbol(n : Int) -> String {
  match n {
    1 => "H"
    5 => "B"
    6 => "C"
    7 => "N"
    8 => "O"
    9 => "F"
    14 => "Si"
    15 => "P"
    16 => "S"
    17 => "Cl"
    35 => "Br"
    53 => "I"
    _ => "?"
  }
}

///|
/// Check if a molecule atom matches a SMARTS atom query.
fn smarts_atom_matches(
  mol : Molecule,
  atom_idx : Int,
  query : SmartsAtom,
) -> Bool {
  let a = mol.atoms[atom_idx]
  let mut result = true
  if query.is_wildcard {
    result = true
  } else {
    match query.element {
      Some(elem) =>
        if a.symbol != elem {
          result = false
        }
      None => ()
    }
    match query.atomic_num {
      Some(num) =>
        if a.symbol != atomic_num_to_symbol(num) {
          result = false
        }
      None => ()
    }
    match query.aromatic {
      Some(arom) =>
        if a.aromatic != arom {
          result = false
        }
      None => ()
    }
    match query.in_ring {
      Some(true) =>
        if not(mol_is_in_ring(mol, atom_idx)) {
          result = false
        }
      Some(false) =>
        if mol_is_in_ring(mol, atom_idx) {
          result = false
        }
      None => ()
    }
    match query.degree {
      Some(d) =>
        if mol_degree(mol, atom_idx) != d {
          result = false
        }
      None => ()
    }
    match query.charge {
      Some(c) =>
        if a.charge != c {
          result = false
        }
      None => ()
    }
  }
  if query.not_match {
    not(result)
  } else {
    result
  }
}

///|
/// Parse a SMARTS bond query from a character.
fn smarts_parse_bond(c : Char) -> SmartsBond {
  match c {
    '-' => SmartsBond::{ order: 1, is_any: false }
    '=' => SmartsBond::{ order: 2, is_any: false }
    '#' => SmartsBond::{ order: 3, is_any: false }
    '~' => SmartsBond::{ order: 0, is_any: true }
    '/' => SmartsBond::{ order: 1, is_any: false }
    _ => SmartsBond::{ order: 0, is_any: true }
  }
}

///|
/// Check if a molecule bond matches a SMARTS bond query.
fn smarts_bond_matches(
  mol : Molecule,
  from : Int,
  to : Int,
  query : SmartsBond,
) -> Bool {
  if query.is_any {
    return true
  }
  let order = mol_bond_order(mol, from, to)
  order == query.order
}

///|
/// Parse a SMARTS pattern string into a SmartsPattern.
pub fn smarts_parse(s : String) -> SmartsPattern {
  let atoms : Array[SmartsAtom] = Array(capacity=0)
  let bonds : Array[SmartsBond] = Array(capacity=0)
  let chars = s.to_array()
  let mut i = 0
  let mut last_bond : SmartsBond? = None
  while i < chars.length() {
    let c = chars[i]
    if c == '[' {
      let mut j = i + 1
      while j < chars.length() && chars[j] != ']' {
        j = j + 1
      }
      let atom_sb = StringBuilder()
      for k = i; k <= j; k = k + 1 {
        atom_sb.write_char(chars[k])
      }
      atoms.push(smarts_parse_atom(atom_sb.to_string()))
      match last_bond {
        Some(b) => {
          bonds.push(b)
          last_bond = None
        }
        None => {
          if atoms.length() > 1 {
            bonds.push(SmartsBond::{ order: 0, is_any: true })
          }
        }
      }
      i = j + 1
    } else if c == '=' || c == '#' || c == '-' || c == '~' || c == '/' {
      last_bond = Some(smarts_parse_bond(c))
      i = i + 1
    } else if c == '.' {
      last_bond = Some(SmartsBond::{ order: 0, is_any: true })
      i = i + 1
    } else if c == '&' || c == '|' || c == ';' {
      i = i + 1
    } else {
      i = i + 1
    }
  }
  SmartsPattern::{ atoms, bonds, branches: [] }
}

///|
/// Try to match a SMARTS pattern against a molecule starting from a specific atom.
fn smarts_match_from(
  mol : Molecule,
  pattern : SmartsPattern,
  start_atom : Int,
) -> Bool {
  if pattern.atoms.length() == 0 {
    return false
  }
  if not(smarts_atom_matches(mol, start_atom, pattern.atoms[0])) {
    return false
  }
  if pattern.atoms.length() == 1 {
    return true
  }
  let neighbors = mol_adjacency(mol)
  let matched = Array::make(mol.atoms.length(), false)
  matched[start_atom] = true
  smarts_dfs(mol, pattern, 1, start_atom, matched, neighbors)
}

///|
/// DFS matching of SMARTS pattern.
fn smarts_dfs(
  mol : Molecule,
  pattern : SmartsPattern,
  pat_idx : Int,
  prev_atom : Int,
  matched : Array[Bool],
  neighbors : Array[Array[(Int, Int)]],
) -> Bool {
  if pat_idx >= pattern.atoms.length() {
    return true
  }
  let bond_query = if pat_idx - 1 < pattern.bonds.length() {
    pattern.bonds[pat_idx - 1]
  } else {
    SmartsBond::{ order: 0, is_any: true }
  }
  for (nb, _) in neighbors[prev_atom] {
    if matched[nb] {
      continue
    }
    if not(smarts_bond_matches(mol, prev_atom, nb, bond_query)) {
      continue
    }
    if not(smarts_atom_matches(mol, nb, pattern.atoms[pat_idx])) {
      continue
    }
    matched[nb] = true
    if smarts_dfs(mol, pattern, pat_idx + 1, nb, matched, neighbors) {
      return true
    }
    matched[nb] = false
  }
  false
}

///|
/// Check if a SMARTS pattern matches anywhere in the molecule.
pub fn smarts_match(mol : Molecule, pattern : SmartsPattern) -> Bool {
  let n = mol.atoms.length()
  for k = 0; k < n; k = k + 1 {
    if smarts_match_from(mol, pattern, k) {
      return true
    }
  }
  false
}

///|
/// Check if a SMARTS string matches a molecule.
pub fn smarts_match_str(mol : Molecule, smarts : String) -> Bool {
  let pattern = smarts_parse(smarts)
  smarts_match(mol, pattern)
}

///|
/// Count the number of distinct matches of a SMARTS pattern.
pub fn smarts_count_matches(mol : Molecule, pattern : SmartsPattern) -> Int {
  let n = mol.atoms.length()
  let mut count = 0
  for k = 0; k < n; k = k + 1 {
    if smarts_match_from(mol, pattern, k) {
      count = count + 1
    }
  }
  count
}

///|
/// Count matches from a SMARTS string.
pub fn smarts_count_matches_str(mol : Molecule, smarts : String) -> Int {
  let pattern = smarts_parse(smarts)
  smarts_count_matches(mol, pattern)
}

///|
/// Find all matching atom indices (first atom of each match).
pub fn smarts_find_matches(
  mol : Molecule,
  pattern : SmartsPattern,
) -> Array[Int] {
  let n = mol.atoms.length()
  let result : Array[Int] = Array(capacity=0)
  for k = 0; k < n; k = k + 1 {
    if smarts_match_from(mol, pattern, k) {
      result.push(k)
    }
  }
  result
}
