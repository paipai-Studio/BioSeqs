// ============================================================================
// MACCS 166-bit fingerprint
// ============================================================================
// MACCS (Molecular ACCess System) 166-bit structural keys.
// Each bit corresponds to a predefined substructure pattern.
// Bits 1-166 (1-based); bit 0 is unused in the original MACCS spec.

///|
/// Check if molecule has at least one atom of the given element.
fn maccs_has_element(mol : Molecule, symbol : String) -> Bool {
  for a in mol.atoms {
    if a.symbol == symbol {
      return true
    }
  }
  false
}

///|
/// Count atoms of a given element.
fn maccs_count_element(mol : Molecule, symbol : String) -> Int {
  let mut count = 0
  for a in mol.atoms {
    if a.symbol == symbol {
      count = count + 1
    }
  }
  count
}

///|
/// Check if molecule has a bond of the given order.
fn maccs_has_bond_order(mol : Molecule, order : Int) -> Bool {
  for b in mol.bonds {
    if b.order == order {
      return true
    }
  }
  false
}

///|
/// Count bonds of a given order.
fn maccs_count_bond_order(mol : Molecule, order : Int) -> Int {
  let mut count = 0
  for b in mol.bonds {
    if b.order == order {
      count = count + 1
    }
  }
  count
}

///|
/// Check if molecule has a ring of the given size.
fn maccs_has_ring_size(mol : Molecule, size : Int) -> Bool {
  let rings = mol_find_rings(mol)
  for ring in rings {
    if ring.length() == size {
      return true
    }
  }
  false
}

///|
/// Check if molecule has an atom with specific degree.
fn maccs_has_degree(mol : Molecule, degree : Int) -> Bool {
  let n = mol.atoms.length()
  for k = 0; k < n; k = k + 1 {
    if mol_degree(mol, k) == degree {
      return true
    }
  }
  false
}

///|
/// Check for a specific atom with specific degree.
fn maccs_has_atom_degree(
  mol : Molecule,
  symbol : String,
  degree : Int,
) -> Bool {
  let n = mol.atoms.length()
  for k = 0; k < n; k = k + 1 {
    if mol.atoms[k].symbol == symbol && mol_degree(mol, k) == degree {
      return true
    }
  }
  false
}

///|
/// Check for aromatic atom.
fn maccs_has_aromatic(mol : Molecule, symbol : String) -> Bool {
  for a in mol.atoms {
    if a.aromatic && a.symbol == symbol {
      return true
    }
  }
  false
}

///|
/// Check for a path of specific length starting from an atom of given element.
fn maccs_has_path(
  mol : Molecule,
  start_symbol : String,
  length : Int,
) -> Bool {
  let n = mol.atoms.length()
  let neighbors = mol_adjacency(mol)
  for k = 0; k < n; k = k + 1 {
    if mol.atoms[k].symbol == start_symbol {
      if maccs_bfs_path(neighbors, k, length) {
        return true
      }
    }
  }
  false
}

///|
/// BFS to check if there's a path of at least `length` bonds from src.
fn maccs_bfs_path(
  neighbors : Array[Array[(Int, Int)]],
  src : Int,
  length : Int,
) -> Bool {
  let n = neighbors.length()
  let dist = Array::make(n, -1)
  dist[src] = 0
  let queue : Array[Int] = [src]
  let mut head = 0
  while head < queue.length() {
    let u = queue[head]
    head = head + 1
    for (v, _) in neighbors[u] {
      if dist[v] == -1 {
        dist[v] = dist[u] + 1
        if dist[v] >= length {
          return true
        }
        queue.push(v)
      }
    }
  }
  false
}

///|
/// Check for X-C-Y pattern (two-hop path with specific elements).
fn maccs_has_xcy(
  mol : Molecule,
  x : String,
  c : String,
  y : String,
) -> Bool {
  let n = mol.atoms.length()
  let neighbors = mol_adjacency(mol)
  for k = 0; k < n; k = k + 1 {
    if mol.atoms[k].symbol == c {
      let mut has_x = false
      let mut has_y = false
      for (nb, _) in neighbors[k] {
        if mol.atoms[nb].symbol == x {
          has_x = true
        }
        if mol.atoms[nb].symbol == y {
          has_y = true
        }
      }
      if has_x && has_y {
        return true
      }
    }
  }
  false
}

///|
/// Count atoms in rings.
fn maccs_count_ring_atoms(mol : Molecule) -> Int {
  let n = mol.atoms.length()
  let mut count = 0
  for k = 0; k < n; k = k + 1 {
    if mol_is_in_ring(mol, k) {
      count = count + 1
    }
  }
  count
}

///|
/// Evaluate a single MACCS key against a molecule.
fn maccs_eval_key(mol : Molecule, bit : Int) -> Bool {
  match bit {
    1 => true
    2 => maccs_has_element(mol, "F")
    3 => maccs_count_element(mol, "Cl") >= 1
    4 => maccs_count_element(mol, "Cl") >= 2
    5 => maccs_count_element(mol, "Cl") >= 3
    6 => maccs_count_element(mol, "Cl") >= 4
    7 => maccs_count_element(mol, "Br") >= 1
    8 => maccs_count_element(mol, "Br") >= 2
    9 => maccs_count_element(mol, "I") >= 1
    10 => maccs_count_element(mol, "I") >= 2
    11 => maccs_count_element(mol, "P") >= 1
    12 => maccs_count_element(mol, "P") >= 2
    13 => maccs_count_element(mol, "S") >= 1
    14 => maccs_count_element(mol, "S") >= 2
    15 => maccs_count_element(mol, "S") >= 3
    16 => maccs_count_element(mol, "S") >= 4
    17 => maccs_has_element(mol, "O")
    18 => maccs_count_element(mol, "O") >= 2
    19 => maccs_count_element(mol, "O") >= 3
    20 => maccs_count_element(mol, "O") >= 4
    21 => maccs_count_element(mol, "O") >= 5
    22 => maccs_count_element(mol, "O") >= 6
    23 => maccs_count_element(mol, "O") >= 7
    24 => maccs_has_element(mol, "N")
    25 => maccs_count_element(mol, "N") >= 2
    26 => maccs_count_element(mol, "N") >= 3
    27 => maccs_count_element(mol, "N") >= 4
    28 => maccs_count_element(mol, "N") >= 5
    29 => maccs_count_element(mol, "N") >= 6
    30 => maccs_count_element(mol, "N") >= 7
    31 => maccs_count_element(mol, "N") >= 8
    32 => maccs_count_element(mol, "C") >= 1
    33 => maccs_count_element(mol, "C") >= 2
    34 => maccs_count_element(mol, "C") >= 3
    35 => maccs_count_element(mol, "C") >= 4
    36 => maccs_count_element(mol, "C") >= 5
    37 => maccs_count_element(mol, "C") >= 6
    38 => maccs_count_element(mol, "C") >= 7
    39 => maccs_count_element(mol, "C") >= 8
    40 => maccs_count_element(mol, "C") >= 9
    41 => maccs_count_element(mol, "C") >= 10
    42 => maccs_count_element(mol, "C") >= 11
    43 => maccs_count_element(mol, "C") >= 12
    44 => maccs_count_element(mol, "C") >= 13
    45 => maccs_count_element(mol, "C") >= 14
    46 => maccs_count_element(mol, "C") >= 15
    47 => maccs_count_element(mol, "C") >= 16
    48 => maccs_count_element(mol, "C") >= 17
    49 => maccs_count_element(mol, "C") >= 18
    50 => maccs_count_element(mol, "C") >= 19
    51 => maccs_count_element(mol, "C") >= 20
    52 => maccs_has_bond_order(mol, 2)
    53 => maccs_count_bond_order(mol, 2) >= 2
    54 => maccs_count_bond_order(mol, 2) >= 3
    55 => maccs_count_bond_order(mol, 2) >= 4
    56 => maccs_has_bond_order(mol, 3)
    57 => maccs_count_bond_order(mol, 3) >= 2
    58 => maccs_has_ring_size(mol, 3)
    59 => maccs_has_ring_size(mol, 4)
    60 => maccs_has_ring_size(mol, 5)
    61 => maccs_has_ring_size(mol, 6)
    62 => mol_ring_count(mol) >= 2
    63 => mol_ring_count(mol) >= 3
    64 => mol_ring_count(mol) >= 4
    65 => mol_aromatic_ring_count(mol) >= 1
    66 => maccs_has_aromatic(mol, "C")
    67 => maccs_has_aromatic(mol, "N")
    68 => maccs_has_aromatic(mol, "O")
    69 => maccs_has_aromatic(mol, "S")
    70 => maccs_has_xcy(mol, "N", "C", "N")
    71 => maccs_has_xcy(mol, "O", "C", "N")
    72 => maccs_has_xcy(mol, "O", "C", "O")
    73 => maccs_has_xcy(mol, "S", "C", "N")
    74 => maccs_has_xcy(mol, "S", "C", "O")
    75 => maccs_has_xcy(mol, "N", "C", "O")
    76 => maccs_has_atom_degree(mol, "C", 4)
    77 => maccs_has_atom_degree(mol, "N", 4)
    78 => maccs_has_atom_degree(mol, "O", 2) && maccs_has_bond_order(mol, 2)
    79 => maccs_has_atom_degree(mol, "N", 3) && maccs_has_bond_order(mol, 2)
    80 => maccs_has_atom_degree(mol, "C", 3) && maccs_has_bond_order(mol, 2)
    81 => maccs_has_degree(mol, 3)
    82 => maccs_count_bond_order(mol, 1) >= 4
    83 => maccs_count_bond_order(mol, 1) >= 6
    84 => maccs_count_bond_order(mol, 1) >= 8
    85 => maccs_count_bond_order(mol, 1) >= 10
    86 => maccs_count_bond_order(mol, 1) >= 12
    87 => maccs_count_bond_order(mol, 1) >= 14
    88 => maccs_count_bond_order(mol, 1) >= 16
    89 => maccs_count_bond_order(mol, 1) >= 18
    90 => maccs_count_bond_order(mol, 1) >= 20
    91 => maccs_count_bond_order(mol, 1) >= 22
    92 => maccs_has_xcy(mol, "C", "C", "N")
    93 => maccs_has_xcy(mol, "C", "C", "O")
    94 => maccs_has_xcy(mol, "C", "N", "C")
    95 => maccs_has_xcy(mol, "C", "O", "C")
    96 => maccs_has_xcy(mol, "N", "N", "C")
    97 => maccs_has_xcy(mol, "N", "O", "C")
    98 => maccs_has_xcy(mol, "O", "O", "C")
    99 => maccs_has_xcy(mol, "C", "C", "C")
    100 => maccs_has_xcy(mol, "C", "C", "S")
    101 => maccs_has_xcy(mol, "C", "S", "C")
    102 => maccs_has_xcy(mol, "C", "S", "S")
    103 => maccs_has_xcy(mol, "N", "C", "S")
    104 => maccs_has_xcy(mol, "O", "C", "S")
    105 => maccs_has_xcy(mol, "S", "S", "C")
    106 => maccs_has_xcy(mol, "C", "C", "F")
    107 => maccs_has_xcy(mol, "C", "C", "Cl")
    108 => maccs_has_xcy(mol, "C", "C", "Br")
    109 => maccs_has_xcy(mol, "C", "C", "I")
    110 => maccs_has_xcy(mol, "C", "N", "N")
    111 => maccs_has_xcy(mol, "C", "O", "N")
    112 => maccs_has_xcy(mol, "C", "O", "O")
    113 => maccs_has_xcy(mol, "N", "C", "C")
    114 => maccs_has_xcy(mol, "O", "C", "C")
    115 => maccs_has_xcy(mol, "C", "C", "P")
    116 => maccs_has_xcy(mol, "C", "P", "C")
    117 => maccs_has_xcy(mol, "N", "C", "P")
    118 => maccs_has_xcy(mol, "O", "C", "P")
    119 => maccs_has_xcy(mol, "P", "C", "P")
    120 => maccs_has_ring_size(mol, 7)
    121 => maccs_has_ring_size(mol, 8)
    122 => mol_ring_count(mol) >= 5
    123 => mol_ring_count(mol) >= 6
    124 => mol_aromatic_ring_count(mol) >= 2
    125 => mol_aromatic_ring_count(mol) >= 3
    126 => mol_saturated_ring_count(mol) >= 1
    127 => mol_saturated_ring_count(mol) >= 2
    128 => maccs_count_ring_atoms(mol) >= 6
    129 => maccs_count_ring_atoms(mol) >= 10
    130 => maccs_count_ring_atoms(mol) >= 15
    131 => maccs_has_path(mol, "C", 3)
    132 => maccs_has_path(mol, "C", 4)
    133 => maccs_has_path(mol, "C", 5)
    134 => maccs_has_path(mol, "C", 6)
    135 => maccs_has_path(mol, "C", 7)
    136 => maccs_has_path(mol, "C", 8)
    137 => maccs_has_path(mol, "C", 9)
    138 => maccs_has_path(mol, "C", 10)
    139 => maccs_has_path(mol, "N", 2)
    140 => maccs_has_path(mol, "O", 2)
    141 => maccs_has_atom_degree(mol, "C", 3)
    142 => maccs_has_atom_degree(mol, "N", 3)
    143 => maccs_has_atom_degree(mol, "O", 1)
    144 => maccs_has_atom_degree(mol, "N", 1)
    145 => maccs_has_atom_degree(mol, "C", 1)
    146 => maccs_count_element(mol, "F") >= 2
    147 => maccs_count_element(mol, "F") >= 3
    148 => maccs_count_element(mol, "F") >= 4
    149 => maccs_has_xcy(mol, "F", "C", "F")
    150 => maccs_has_xcy(mol, "F", "C", "N")
    151 => maccs_has_xcy(mol, "F", "C", "O")
    152 => maccs_has_xcy(mol, "Cl", "C", "Cl")
    153 => maccs_has_xcy(mol, "Br", "C", "Br")
    154 => maccs_has_xcy(mol, "I", "C", "I")
    155 => maccs_count_element(mol, "Cl") >= 5
    156 => maccs_count_element(mol, "Br") >= 3
    157 => maccs_count_element(mol, "I") >= 3
    158 => maccs_has_xcy(mol, "N", "N", "N")
    159 => maccs_has_xcy(mol, "O", "N", "O")
    160 => maccs_has_xcy(mol, "S", "O", "S")
    161 => maccs_has_xcy(mol, "N", "S", "N")
    162 => maccs_has_xcy(mol, "O", "S", "O")
    163 => mol_num_atoms(mol) >= 20
    164 => mol_num_atoms(mol) >= 30
    165 => mol_num_atoms(mol) >= 40
    166 => mol_num_bonds(mol) >= 30
    _ => false
  }
}

///|
/// Compute MACCS 166-bit fingerprint as a bit array (1-based, index 0 unused).
/// Returns array of 167 Bool values (index 0 = false, indices 1-166 = key bits).
pub fn maccs_fingerprint(mol : Molecule) -> Array[Bool] {
  let bits = Array::make(167, false)
  for bit = 1; bit <= 166; bit = bit + 1 {
    bits[bit] = maccs_eval_key(mol, bit)
  }
  bits
}

///|
/// Compute MACCS fingerprint as a 166-bit integer (bits 1-166 packed).
pub fn maccs_fingerprint_int(mol : Molecule) -> Int64 {
  let bits = maccs_fingerprint(mol)
  let mut result : Int64 = 0L
  for bit = 1; bit <= 166; bit = bit + 1 {
    if bits[bit] {
      result = result | (1L << bit)
    }
  }
  result
}

///|
/// MACCS fingerprint as a binary string (166 chars, left=bit1, right=bit166).
pub fn maccs_fingerprint_string(mol : Molecule) -> String {
  let bits = maccs_fingerprint(mol)
  let sb = StringBuilder()
  for bit = 1; bit <= 166; bit = bit + 1 {
    if bits[bit] {
      sb.write_char('1')
    } else {
      sb.write_char('0')
    }
  }
  sb.to_string()
}

///|
/// Count number of set bits in MACCS fingerprint.
pub fn maccs_bit_count(mol : Molecule) -> Int {
  let bits = maccs_fingerprint(mol)
  let mut count = 0
  for bit = 1; bit <= 166; bit = bit + 1 {
    if bits[bit] {
      count = count + 1
    }
  }
  count
}

///|
/// Tanimoto similarity between two MACCS fingerprints.
pub fn maccs_tanimoto(mol1 : Molecule, mol2 : Molecule) -> Double {
  let bits1 = maccs_fingerprint(mol1)
  let bits2 = maccs_fingerprint(mol2)
  let mut intersection = 0
  let mut union = 0
  for bit = 1; bit <= 166; bit = bit + 1 {
    if bits1[bit] && bits2[bit] {
      intersection = intersection + 1
    }
    if bits1[bit] || bits2[bit] {
      union = union + 1
    }
  }
  if union == 0 {
    return 0.0
  }
  intersection.to_double() / union.to_double()
}
