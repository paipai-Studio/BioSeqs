/// rdkit_substruct_feature.mbt — substructure matching + feature/pharmacophore fingerprint
/// Items 6-7

///|
/// Check if a molecule contains a substructure (given as SMILES).
pub fn rdkit_has_substruct(mol_smiles : String, query_smiles : String) -> Bool {
  let mol = parse_smiles(mol_smiles)
  let query = parse_smiles(query_smiles)
  let mol_n = mol.num_atoms()
  let query_n = query.num_atoms()
  if query_n > mol_n {
    return false
  }
  if query_n == 0 {
    return true
  }
  let mol_adj = rdkit_build_adj_matrix(mol)
  let query_adj = rdkit_build_adj_matrix(query)
  let mol_sym = mol.atom_symbols()
  let query_sym = query.atom_symbols()
  // Try all possible mappings
  let mapping : Array[Int32] = Array(capacity=query_n)
  for i = 0; i < query_n; i = i + 1 {
    mapping.push(-1)
  }
  rdkit_substruct_backtrack(
    0, query_n, mol_n, mapping, mol_sym, query_sym, mol_adj, query_adj,
  )
}

///|
fn rdkit_substruct_backtrack(
  depth : Int32,
  query_n : Int32,
  mol_n : Int32,
  mapping : Array[Int32],
  mol_sym : Array[String],
  query_sym : Array[String],
  mol_adj : Array[Array[Int32]],
  query_adj : Array[Array[Int32]],
) -> Bool {
  if depth == query_n {
    return true
  }
  for cand = 0; cand < mol_n; cand = cand + 1 {
    // Check if already used
    let mut used = false
    for i = 0; i < depth; i = i + 1 {
      if mapping[i] == cand {
        used = true
      }
    }
    if used {
      continue
    }
    // Check atom type match
    if mol_sym[cand] != query_sym[depth] {
      continue
    }
    // Check bond consistency
    let mut ok = true
    for i = 0; i < depth; i = i + 1 {
      if query_adj[i][depth] > 0 {
        if mol_adj[mapping[i]][cand] == 0 {
          ok = false
        }
      }
    }
    if !ok {
      continue
    }
    mapping[depth] = cand
    if rdkit_substruct_backtrack(
        depth + 1,
        query_n,
        mol_n,
        mapping,
        mol_sym,
        query_sym,
        mol_adj,
        query_adj,
      ) {
      return true
    }
    mapping[depth] = -1
  }
  false
}

///|
/// Get all substructure matches (returns list of atom index mappings).
pub fn rdkit_get_substruct_matches(
  mol_smiles : String,
  query_smiles : String,
  max_matches~ : Int32 = 100,
) -> Array[Array[Int32]] {
  let mol = parse_smiles(mol_smiles)
  let query = parse_smiles(query_smiles)
  let mol_n = mol.num_atoms()
  let query_n = query.num_atoms()
  let results : Array[Array[Int32]] = Array(capacity=0)
  if query_n > mol_n {
    return results
  }
  let mol_adj = rdkit_build_adj_matrix(mol)
  let query_adj = rdkit_build_adj_matrix(query)
  let mol_sym = mol.atom_symbols()
  let query_sym = query.atom_symbols()
  let mapping : Array[Int32] = Array(capacity=query_n)
  for i = 0; i < query_n; i = i + 1 {
    mapping.push(-1)
  }
  rdkit_substruct_collect(
    0, query_n, mol_n, mapping, mol_sym, query_sym, mol_adj, query_adj, results,
    max_matches,
  )
  results
}

///|
fn rdkit_substruct_collect(
  depth : Int32,
  query_n : Int32,
  mol_n : Int32,
  mapping : Array[Int32],
  mol_sym : Array[String],
  query_sym : Array[String],
  mol_adj : Array[Array[Int32]],
  query_adj : Array[Array[Int32]],
  results : Array[Array[Int32]],
  max_matches : Int32,
) -> Unit {
  if results.length() >= max_matches {
    return
  }
  if depth == query_n {
    let copy : Array[Int32] = Array(capacity=query_n)
    for i = 0; i < query_n; i = i + 1 {
      copy.push(mapping[i])
    }
    results.push(copy)
    return
  }
  for cand = 0; cand < mol_n; cand = cand + 1 {
    let mut used = false
    for i = 0; i < depth; i = i + 1 {
      if mapping[i] == cand {
        used = true
      }
    }
    if used {
      continue
    }
    if mol_sym[cand] != query_sym[depth] {
      continue
    }
    let mut ok = true
    for i = 0; i < depth; i = i + 1 {
      if query_adj[i][depth] > 0 {
        if mol_adj[mapping[i]][cand] == 0 {
          ok = false
        }
      }
    }
    if !ok {
      continue
    }
    mapping[depth] = cand
    rdkit_substruct_collect(
      depth + 1,
      query_n,
      mol_n,
      mapping,
      mol_sym,
      query_sym,
      mol_adj,
      query_adj,
      results,
      max_matches,
    )
    mapping[depth] = -1
  }
}

///|
/// Filter a list of molecules by substructure query.
pub fn rdkit_filter_by_substruct(
  mol_smiles_list : Array[String],
  query_smiles : String,
) -> Array[String] {
  let result : Array[String] = Array(capacity=0)
  for i = 0; i < mol_smiles_list.length(); i = i + 1 {
    if rdkit_has_substruct(mol_smiles_list[i], query_smiles) {
      result.push(mol_smiles_list[i])
    }
  }
  result
}

///|
// Feature types for pharmacophore fingerprint
pub enum FcfpFeature {
  FcfpHBD // hydrogen bond donor
  FcfpHBA // hydrogen bond acceptor
  FcfpHydrophobic
  FcfpAromatic
  FcfpPos // positive charge
  FcfpNeg // negative charge
}

///|
/// Identify pharmacophore features in a molecule.
pub fn rdkit_identify_features(smiles : String) -> Array[(Int32, FcfpFeature)] {
  let mol = parse_smiles(smiles)
  let symbols = mol.atom_symbols()
  let n = mol.num_atoms()
  let adj = rdkit_build_adj_matrix(mol)
  let features : Array[(Int32, FcfpFeature)] = Array(capacity=0)
  for i = 0; i < n; i = i + 1 {
    let s = symbols[i]
    // HBD: N or O with implicit H
    if s == "N" || s == "O" {
      features.push((i, FcfpHBD))
    }
    // HBA: N or O
    if s == "N" || s == "O" {
      features.push((i, FcfpHBA))
    }
    // Hydrophobic: C
    if s == "C" {
      features.push((i, FcfpHydrophobic))
    }
    // Aromatic: C in aromatic ring (simplified)
    if s == "C" {
      if mol_is_in_ring(mol, i) && mol_aromatic_ring_count(mol) > 0 {
        features.push((i, FcfpAromatic))
      }
    }
    // Positive: N with explicit positive charge or amine
    if s == "N" {
      features.push((i, FcfpPos))
    }
    // Negative: O in COOH or similar
    if s == "O" {
      // Check if bonded to another O (carboxyl)
      let mut has_o_neighbor = false
      for j = 0; j < n; j = j + 1 {
        if adj[i][j] > 0 && symbols[j] == "O" {
          has_o_neighbor = true
        }
      }
      if has_o_neighbor {
        features.push((i, FcfpNeg))
      }
    }
  }
  features
}

///|
/// Generate feature fingerprint (FCFP): bit vector based on pharmacophore features.
pub fn rdkit_feature_fp(
  smiles : String,
  n_bits~ : Int32 = 2048,
) -> Array[Int32] {
  let features = rdkit_identify_features(smiles)
  let fp = Array(capacity=n_bits)
  for i = 0; i < n_bits; i = i + 1 {
    fp.push(0)
  }
  // Hash each feature into bits
  for i = 0; i < features.length(); i = i + 1 {
    let atom_idx = features[i].0
    let feat_hash = rdkit_feat_hash(features[i].1)
    let bit = ((feat_hash * 31 + atom_idx) & 0x7FFFFFFF) % n_bits
    fp[bit] = 1
  }
  // Also hash pairs of features (2-point pharmacophores)
  for i = 0; i < features.length(); i = i + 1 {
    for j = i + 1; j < features.length(); j = j + 1 {
      let h1 = rdkit_feat_hash(features[i].1)
      let h2 = rdkit_feat_hash(features[j].1)
      let dist = features[i].0 * 7 + features[j].0 * 11
      let bit = ((h1 * 31 + h2 * 37 + dist) & 0x7FFFFFFF) % n_bits
      fp[bit] = 1
    }
  }
  fp
}

///|
fn rdkit_feat_hash(f : FcfpFeature) -> Int32 {
  match f {
    FcfpHBD => 1
    FcfpHBA => 2
    FcfpHydrophobic => 3
    FcfpAromatic => 5
    FcfpPos => 7
    FcfpNeg => 11
  }
}

///|
/// Count features by type. Returns array of (feature_name, count).
pub fn rdkit_count_features(smiles : String) -> Array[(String, Int32)] {
  let features = rdkit_identify_features(smiles)
  let names : Array[String] = [
    "HBD", "HBA", "Hydrophobic", "Aromatic", "Pos", "Neg",
  ]
  let counts : Array[Int32] = [0, 0, 0, 0, 0, 0]
  for i = 0; i < features.length(); i = i + 1 {
    let idx = rdkit_feat_hash(features[i].1)
    if idx >= 1 && idx <= 6 {
      counts[idx - 1] = counts[idx - 1] + 1
    }
  }
  let result : Array[(String, Int32)] = Array(capacity=6)
  for i = 0; i < 6; i = i + 1 {
    result.push((names[i], counts[i]))
  }
  result
}

///|
pub fn rdkit_substruct_feature_dummy() -> Unit {
  ignore(FcfpHBD)
  ignore(FcfpHBA)
  ignore(FcfpHydrophobic)
  ignore(FcfpAromatic)
  ignore(FcfpPos)
  ignore(FcfpNeg)
}
