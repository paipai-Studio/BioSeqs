///|
/// Build adjacency matrix from Molecule (helper for rdkit second-round modules).
pub fn rdkit_build_adj_matrix(mol : Molecule) -> Array[Array[Int32]] {
  let n = mol.num_atoms()
  let adj : Array[Array[Int32]] = Array(capacity=n)
  for i = 0; i < n; i = i + 1 {
    let row : Array[Int32] = Array(capacity=n)
    for j = 0; j < n; j = j + 1 {
      row.push(0)
    }
    adj.push(row)
  }
  let adj_list = mol_adjacency(mol)
  for i = 0; i < n; i = i + 1 {
    for j = 0; j < adj_list[i].length(); j = j + 1 {
      let (nb, order) = adj_list[i][j]
      adj[i][nb] = order
    }
  }
  adj
}

/// rdkit_fp_identifiers.mbt — InChI, MACCS, Avalon, topo torsion FP, ECFP variants
/// Items 1-5

///|
/// Generate a simplified InChI string from a SMILES.
/// This is a simplified version — full InChI requires the IUPAC algorithm.
pub fn rdkit_smiles_to_inchi(smiles : String) -> String {
  // Parse the molecule to get formula and connectivity
  let mol = parse_smiles(smiles)
  let _n_atoms = mol.num_atoms()
  // Build formula
  let formula = rdkit_inchi_formula(mol)
  // Build connectivity layer
  let conn = rdkit_inchi_connectivity(mol)
  // Build hydrogen layer
  let h_layer = rdkit_inchi_hydrogen(mol)
  let inchi = "InChI=1S/" + formula + "/c" + conn + "/h" + h_layer
  inchi
}

///|
/// Generate InChIKey (26-char hash) from InChI string.
/// Uses a simplified hash: first 14 chars from formula, rest from connectivity.
pub fn rdkit_inchi_to_inchikey(inchi : String) -> String {
  // Simplified: hash the InChI string
  let mut hash1 = 0
  let mut hash2 = 0
  for i = 0; i < inchi.length(); i = i + 1 {
    let c = inchi[i].to_int()
    hash1 = (hash1 * 31 + c) & 0x3FFFFFFF
    hash2 = (hash2 * 37 + c) & 0xFFFFF
  }
  // Encode as base26 (A-Z)
  let key1 = rdkit_encode_base26(hash1, 14)
  let key2 = rdkit_encode_base26(hash2, 10)
  key1 + "-" + key2 + "-N"
}

///|
fn rdkit_encode_base26(val : Int32, len : Int32) -> String {
  let sb = StringBuilder()
  let mut v = val
  for i = 0; i < len; i = i + 1 {
    let c = v % 26 + 65
    let _ = sb.write_char(c.unsafe_to_char())
    v = v / 26
  }
  sb.to_string()
}

///|
fn rdkit_inchi_formula(mol : Molecule) -> String {
  let symbols = mol.atom_symbols()
  // Count each element
  let elem_counts : Array[(String, Int32)] = Array(capacity=0)
  for i = 0; i < symbols.length(); i = i + 1 {
    let s = symbols[i]
    let mut found = false
    for j = 0; j < elem_counts.length(); j = j + 1 {
      if elem_counts[j].0 == s {
        elem_counts[j] = (s, elem_counts[j].1 + 1)
        found = true
      }
    }
    if !found {
      elem_counts.push((s, 1))
    }
  }
  // Sort: C first, H second, then alphabetical
  let sorted : Array[(String, Int32)] = Array(capacity=0)
  for i = 0; i < elem_counts.length(); i = i + 1 {
    if elem_counts[i].0 == "C" {
      sorted.push(elem_counts[i])
    }
  }
  for i = 0; i < elem_counts.length(); i = i + 1 {
    if elem_counts[i].0 == "H" {
      sorted.push(elem_counts[i])
    }
  }
  for i = 0; i < elem_counts.length(); i = i + 1 {
    if elem_counts[i].0 != "C" && elem_counts[i].0 != "H" {
      sorted.push(elem_counts[i])
    }
  }
  let sb = StringBuilder()
  for i = 0; i < sorted.length(); i = i + 1 {
    let _ = sb.write(sorted[i].0)
    if sorted[i].1 > 1 {
      let _ = sb.write(sorted[i].1.to_string())
    }
  }
  sb.to_string()
}

///|
fn rdkit_inchi_connectivity(mol : Molecule) -> String {
  let n = mol.num_atoms()
  let adj = rdkit_build_adj_matrix(mol)
  let sb = StringBuilder()
  for i = 0; i < n; i = i + 1 {
    for j = i + 1; j < n; j = j + 1 {
      if adj[i][j] > 0 {
        let _ = sb.write((i + 1).to_string())
        let _ = sb.write("-")
        let _ = sb.write((j + 1).to_string())
        if j < n - 1 {
          let _ = sb.write(",")
        }
      }
    }
  }
  sb.to_string()
}

///|
fn rdkit_inchi_hydrogen(mol : Molecule) -> String {
  let n = mol.num_atoms()
  let symbols = mol.atom_symbols()
  let sb = StringBuilder()
  let mut first = true
  for i = 0; i < n; i = i + 1 {
    // Simplified: assume 1 H for C in typical organic molecules
    if symbols[i] != "H" {
      if !first {
        let _ = sb.write(",")
      }
      let _ = sb.write((i + 1).to_string())
      let _ = sb.write("H")
      first = false
    }
  }
  sb.to_string()
}

///|
/// MACCS keys: 166-bit structural fingerprint.
/// Simplified implementation covering key substructure patterns.
pub fn rdkit_maccs_fp(smiles : String) -> Array[Int32] {
  let mol = parse_smiles(smiles)
  let fp = Array(capacity=167)
  for i = 0; i < 167; i = i + 1 {
    fp.push(0)
  }
  // Key 1: molecular weight > 100
  let mw = rdkit_fp_mol_wt(mol)
  if mw > 100.0 {
    fp[1] = 1
  }
  // Key 2: has at least one ring
  if mol_ring_count(mol) > 0 {
    fp[2] = 1
  }
  // Key 3: has aromatic ring
  if mol_aromatic_ring_count(mol) > 0 {
    fp[3] = 1
  }
  // Key 6: has C
  if rdkit_fp_has_element(mol, "C") {
    fp[6] = 1
  }
  // Key 7: has N
  if rdkit_fp_has_element(mol, "N") {
    fp[7] = 1
  }
  // Key 8: has O
  if rdkit_fp_has_element(mol, "O") {
    fp[8] = 1
  }
  // Key 9: has F
  if rdkit_fp_has_element(mol, "F") {
    fp[9] = 1
  }
  // Key 10: has P
  if rdkit_fp_has_element(mol, "P") {
    fp[10] = 1
  }
  // Key 11: has S
  if rdkit_fp_has_element(mol, "S") {
    fp[11] = 1
  }
  // Key 12: has Cl
  if rdkit_fp_has_element(mol, "Cl") {
    fp[12] = 1
  }
  // Key 13: has Br
  if rdkit_fp_has_element(mol, "Br") {
    fp[13] = 1
  }
  // Key 14: has I
  if rdkit_fp_has_element(mol, "I") {
    fp[14] = 1
  }
  // Key 20: >= 2 rings
  if mol_ring_count(mol) >= 2 {
    fp[20] = 1
  }
  // Key 21: >= 3 rings
  if mol_ring_count(mol) >= 3 {
    fp[21] = 1
  }
  // Key 22: >= 4 rings
  if mol_ring_count(mol) >= 4 {
    fp[22] = 1
  }
  // Key 44: >= 6 heavy atoms
  if mol.num_atoms() >= 6 {
    fp[44] = 1
  }
  // Key 45: >= 8 heavy atoms
  if mol.num_atoms() >= 8 {
    fp[45] = 1
  }
  // Key 46: >= 10 heavy atoms
  if mol.num_atoms() >= 10 {
    fp[46] = 1
  }
  // Key 47: >= 12 heavy atoms
  if mol.num_atoms() >= 12 {
    fp[47] = 1
  }
  // Key 48: >= 14 heavy atoms
  if mol.num_atoms() >= 14 {
    fp[48] = 1
  }
  // Key 49: >= 16 heavy atoms
  if mol.num_atoms() >= 16 {
    fp[49] = 1
  }
  // Key 50: >= 18 heavy atoms
  if mol.num_atoms() >= 18 {
    fp[50] = 1
  }
  // Key 51: >= 20 heavy atoms
  if mol.num_atoms() >= 20 {
    fp[51] = 1
  }
  // Key 166: MW > 600
  if mw > 600.0 {
    fp[166] = 1
  }
  fp
}

///|
fn rdkit_fp_has_element(mol : Molecule, elem : String) -> Bool {
  let symbols = mol.atom_symbols()
  for i = 0; i < symbols.length(); i = i + 1 {
    if symbols[i] == elem {
      return true
    }
  }
  false
}

///|
fn rdkit_fp_mol_wt(mol : Molecule) -> Double {
  let symbols = mol.atom_symbols()
  let mut total = 0.0
  for i = 0; i < symbols.length(); i = i + 1 {
    total = total + rdkit_fp_atomic_wt(symbols[i])
  }
  total
}

///|
fn rdkit_fp_atomic_wt(sym : String) -> Double {
  match sym {
    "H" => 1.008
    "C" => 12.011
    "N" => 14.007
    "O" => 15.999
    "F" => 18.998
    "P" => 30.974
    "S" => 32.06
    "Cl" => 35.45
    "Br" => 79.904
    "I" => 126.904
    _ => 0.0
  }
}

///|
/// Avalon fingerprint: 512-bit substructure fingerprint.
/// Simplified: hashes atom environments into 512 bits.
pub fn rdkit_avalon_fp(smiles : String, n_bits~ : Int32 = 512) -> Array[Int32] {
  let mol = parse_smiles(smiles)
  let fp = Array(capacity=n_bits)
  for i = 0; i < n_bits; i = i + 1 {
    fp.push(0)
  }
  // Hash each atom environment
  let n = mol.num_atoms()
  let symbols = mol.atom_symbols()
  let adj = rdkit_build_adj_matrix(mol)
  for i = 0; i < n; i = i + 1 {
    let mut env_hash = 0
    for c = 0; c < symbols[i].length(); c = c + 1 {
      env_hash = env_hash * 31 + symbols[i][c].to_int()
    }
    for j = 0; j < n; j = j + 1 {
      if adj[i][j] > 0 {
        env_hash = env_hash * 17 + adj[i][j]
      }
    }
    let bit = (env_hash & 0x7FFFFFFF) % n_bits
    fp[bit] = 1
  }
  // Set structure-based bits
  if mol_ring_count(mol) > 0 {
    fp[0] = 1
  }
  if mol_aromatic_ring_count(mol) > 0 {
    fp[1] = 1
  }
  fp
}

///|
/// Topological torsion fingerprint: 4-atom path fingerprint.
pub fn rdkit_topo_torsion_fp(
  smiles : String,
  n_bits~ : Int32 = 2048,
) -> Array[Int32] {
  let mol = parse_smiles(smiles)
  let fp = Array(capacity=n_bits)
  for i = 0; i < n_bits; i = i + 1 {
    fp.push(0)
  }
  let n = mol.num_atoms()
  let symbols = mol.atom_symbols()
  let adj = rdkit_build_adj_matrix(mol)
  // Enumerate all 4-atom paths
  for a = 0; a < n; a = a + 1 {
    for b = 0; b < n; b = b + 1 {
      if a == b || adj[a][b] == 0 {
        continue
      }
      for c = 0; c < n; c = c + 1 {
        if c == a || c == b || adj[b][c] == 0 {
          continue
        }
        for d = c + 1; d < n; d = d + 1 {
          if d == a || d == b || adj[c][d] == 0 {
            continue
          }
          // Hash the torsion (a,b,c,d)
          let mut h = rdkit_fp_hash_sym(symbols[a].to_string()) * 7
          h = h + rdkit_fp_hash_sym(symbols[b].to_string()) * 11
          h = h + rdkit_fp_hash_sym(symbols[c].to_string()) * 13
          h = h + rdkit_fp_hash_sym(symbols[d].to_string()) * 17
          let bit = (h & 0x7FFFFFFF) % n_bits
          fp[bit] = 1
        }
      }
    }
  }
  fp
}

///|
fn rdkit_fp_hash_sym(s : String) -> Int32 {
  let mut h = 0
  for i = 0; i < s.length(); i = i + 1 {
    h = h * 31 + s[i].to_int()
  }
  h
}

///|
/// ECFP fingerprint with configurable radius and bit count.
/// ECFP4 = radius 2, ECFP6 = radius 3.
pub fn rdkit_ecfp_fp(
  smiles : String,
  radius~ : Int32 = 2, // 2=ECFP4, 3=ECFP6
  n_bits~ : Int32 = 2048,
) -> Array[Int32] {
  let mol = parse_smiles(smiles)
  let fp = Array(capacity=n_bits)
  for i = 0; i < n_bits; i = i + 1 {
    fp.push(0)
  }
  let n = mol.num_atoms()
  let symbols = mol.atom_symbols()
  let adj = rdkit_build_adj_matrix(mol)
  // For each atom, compute its environment hash at each radius
  let env_hashes : Array[Int32] = Array(capacity=n)
  for i = 0; i < n; i = i + 1 {
    env_hashes.push(rdkit_fp_hash_sym(symbols[i].to_string()))
  }
  // Iteratively expand
  for r = 0; r < radius; r = r + 1 {
    let new_hashes : Array[Int32] = Array(capacity=n)
    for i = 0; i < n; i = i + 1 {
      let mut h = env_hashes[i] * 31
      for j = 0; j < n; j = j + 1 {
        if adj[i][j] > 0 {
          h = h + env_hashes[j] * adj[i][j]
        }
      }
      new_hashes.push(h & 0x7FFFFFFF)
    }
    for i = 0; i < n; i = i + 1 {
      env_hashes[i] = new_hashes[i]
    }
  }
  // Fold into bit vector
  for i = 0; i < n; i = i + 1 {
    let bit = env_hashes[i] % n_bits
    fp[bit] = 1
  }
  fp
}

///|
/// Compute Tanimoto similarity between two bit vectors.
pub fn rdkit_fp_tanimoto(fp1 : Array[Int32], fp2 : Array[Int32]) -> Double {
  let n = if fp1.length() < fp2.length() { fp1.length() } else { fp2.length() }
  let mut intersection = 0
  let mut union = 0
  for i = 0; i < n; i = i + 1 {
    if fp1[i] > 0 && fp2[i] > 0 {
      intersection = intersection + 1
      union = union + 1
    } else if fp1[i] > 0 || fp2[i] > 0 {
      union = union + 1
    }
  }
  if union == 0 {
    return 0.0
  }
  intersection.to_double() / union.to_double()
}
