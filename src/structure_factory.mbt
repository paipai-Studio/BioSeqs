///|
/// Factory and IO convenience methods for AtomArray.
///
/// Bridge between the high-level Structure/Chain/Residue model and
/// the flat AtomArray representation.

// ============================================================
// P5-1: AtomArray factory — Structure → AtomArray
// ============================================================

///|
/// Build an AtomArray from a pdb::Structure (first model).
pub fn AtomArray::from_structure(
  structure : Structure,
  model_index~ : Int32 = 0,
) -> AtomArray {
  let models = structure.get_models()
  if models.length() == 0 {
    return AtomArray::new()
  }
  let mi = if model_index < models.length() { model_index } else { 0 }
  let model = models[mi]
  let all_atoms = model.get_atoms()
  AtomArray::from_atoms(all_atoms)
}

///|
/// Build an AtomArray from PDB text.
pub fn AtomArray::from_pdb_text(text : String) -> AtomArray {
  let structure = parse_pdb(text)
  AtomArray::from_structure(structure)
}

///|
/// Build an AtomArray from a pdb::Chain.
pub fn AtomArray::from_chain(chain : Chain) -> AtomArray {
  let atoms : Array[Atom] = []
  for residue in chain.get_residues() {
    for atom in residue.atoms {
      atoms.push(atom)
    }
  }
  AtomArray::from_atoms(atoms)
}

///|
/// Build an AtomArray from Array[Residue].
pub fn AtomArray::from_residues(residues : Array[Residue]) -> AtomArray {
  let atoms : Array[Atom] = []
  for residue in residues {
    for atom in residue.atoms {
      atoms.push(atom)
    }
  }
  AtomArray::from_atoms(atoms)
}

// ============================================================
// P5-2: AtomArray → Chain (flat → hierarchical, for calc_sasa etc.)
// ============================================================

///|
/// Convert AtomArray → Array[Chain] grouped by chain_id + res_id.
pub fn AtomArray::to_chains(self : AtomArray) -> Array[Chain] {
  let starts = self.get_residue_starts()
  let n_res = starts.length()
  if n_res == 0 {
    return []
  }
  let chain_starts = self.get_chain_starts()
  let result : Array[Chain] = Array(capacity=chain_starts.length())
  let mut ci = 0
  while ci < chain_starts.length() {
    let cs = chain_starts[ci]
    let ce = if ci + 1 < chain_starts.length() {
      chain_starts[ci + 1]
    } else {
      n_res
    }
    let cid = self.chain_ids[cs]
    let residues : Array[Residue] = Array(capacity=ce - cs)
    let mut ri = cs
    while ri < ce {
      let rs = starts[ri]
      let re = if ri + 1 < n_res { starts[ri + 1] } else { self.len() }
      let rn = self.res_names[rs]
      let rid = self.res_ids[rs]
      let atoms : Array[Atom] = Array(capacity=re - rs)
      let mut ai = rs
      while ai < re {
        atoms.push(
          Atom::new(
            name=self.names[ai],
            coord=self.coords[ai],
            resname=self.res_names[ai],
            chainid=self.chain_ids[ai],
            resseq=self.res_ids[ai],
            occupancy=self.occupancies[ai],
            bfactor=self.b_factors[ai],
            element=self.elements[ai],
          ),
        )
        ai = ai + 1
      }
      residues.push(Residue::new(resname=rn, chainid=cid, resseq=rid, atoms~))
      ri = ri + 1
    }
    result.push(Chain::new(id=cid, residues~))
    ci = ci + 1
  }
  result
}

///|
/// Convert AtomArray → pdb::Structure (single model).
pub fn AtomArray::to_structure(self : AtomArray) -> Structure {
  let chains = self.to_chains()
  let model = Model::new(id=1, chains~)
  Structure::new(id="", models=[model])
}

// ============================================================
// P5-3: SASA for AtomArray (wrap existing pdb_analysis.calc_sasa)
// ============================================================

///|
/// Compute per-atom SASA using Shrake-Rupley from pdb_analysis.
/// Returns Array[Double], same length as self.
pub fn AtomArray::calc_sasa(
  self : AtomArray,
  probe_radius~ : Double = 1.4,
  n_points~ : Int32 = 100,
) -> Array[Double] {
  let n = self.len()
  let result = Array::make(n, 0.0)
  let chains = self.to_chains()
  let mut offset = 0
  for chain in chains {
    let chain_id = chain.get_id()
    let per_atom = calc_sasa(chain, probe_radius~, n_points~)
    let mut pi = 0
    while pi < per_atom.length() {
      let (key, sasa_val) = per_atom[pi]
      // key = "chainid_resseq_atomname"
      let mut best_idx = -1
      let mut k = offset
      while k < n {
        if self.chain_ids[k] != chain_id {
          break
        }
        let expected_key = self.chain_ids[k].to_string() +
          "_" +
          self.res_ids[k].to_string() +
          "_" +
          self.names[k]
        if expected_key == key {
          best_idx = k
          break
        }
        k = k + 1
      }
      if best_idx >= 0 {
        result[best_idx] = sasa_val
      }
      pi = pi + 1
    }
    // Advance offset past this chain
    let mut cnt = 0
    let mut k2 = offset
    while k2 < n && self.chain_ids[k2] == chain_id {
      cnt = cnt + 1
      k2 = k2 + 1
    }
    offset = offset + cnt
  }
  result
}

///|
/// Compute total SASA (sum of per-atom).
pub fn AtomArray::calc_total_sasa(
  self : AtomArray,
  probe_radius~ : Double = 1.4,
  n_points~ : Int32 = 100,
) -> Double {
  let per_atom = self.calc_sasa(probe_radius~, n_points~)
  let mut sum = 0.0
  for v in per_atom {
    sum = sum + v
  }
  sum
}

///|
/// Compute per-residue SASA.
pub fn AtomArray::calc_residue_sasa(
  self : AtomArray,
  probe_radius~ : Double = 1.4,
  n_points~ : Int32 = 100,
) -> Array[Double] {
  let per_atom = self.calc_sasa(probe_radius~, n_points~)
  let starts = self.get_residue_starts()
  let n_res = starts.length()
  let result : Array[Double] = Array(capacity=n_res)
  let mut i = 0
  while i < n_res {
    let s = starts[i]
    let e = if i + 1 < n_res { starts[i + 1] } else { self.len() }
    let mut sum = 0.0
    let mut j = s
    while j < e {
      sum = sum + per_atom[j]
      j = j + 1
    }
    result.push(sum)
    i = i + 1
  }
  result
}

// ============================================================
// P5-4: Convenience utilities
// ============================================================

///|
/// Count unique residues (same as get_residue_starts length).
pub fn AtomArray::count_unique_residues(self : AtomArray) -> Int32 {
  self.get_residue_starts().length()
}

///|
/// Get unique chain IDs in appearance order.
pub fn AtomArray::get_unique_chain_ids(self : AtomArray) -> Array[Char] {
  let result : Array[Char] = []
  let seen : Array[Char] = []
  for id in self.chain_ids {
    let mut found = false
    for s in seen {
      if s == id {
        found = true
        break
      }
    }
    if !found {
      seen.push(id)
      result.push(id)
    }
  }
  result
}

///|
/// Per-residue average B-factor.
pub fn AtomArray::per_residue_bfactor(self : AtomArray) -> Array[Double] {
  let starts = self.get_residue_starts()
  let n_res = starts.length()
  let result : Array[Double] = Array(capacity=n_res)
  let mut i = 0
  while i < n_res {
    let s = starts[i]
    let e = if i + 1 < n_res { starts[i + 1] } else { self.len() }
    let mut sum = 0.0
    let mut cnt = 0
    let mut j = s
    while j < e {
      sum = sum + self.b_factors[j]
      cnt = cnt + 1
      j = j + 1
    }
    result.push(if cnt > 0 { sum / cnt.to_double() } else { 0.0 })
    i = i + 1
  }
  result
}

///|
/// Per-residue average occupancy.
pub fn AtomArray::per_residue_occupancy(self : AtomArray) -> Array[Double] {
  let starts = self.get_residue_starts()
  let n_res = starts.length()
  let result : Array[Double] = Array(capacity=n_res)
  let mut i = 0
  while i < n_res {
    let s = starts[i]
    let e = if i + 1 < n_res { starts[i + 1] } else { self.len() }
    let mut sum = 0.0
    let mut cnt = 0
    let mut j = s
    while j < e {
      sum = sum + self.occupancies[j]
      cnt = cnt + 1
      j = j + 1
    }
    result.push(if cnt > 0 { sum / cnt.to_double() } else { 0.0 })
    i = i + 1
  }
  result
}

///|
/// Per-chain average B-factor.
pub fn AtomArray::per_chain_bfactor(self : AtomArray) -> Array[Double] {
  let chain_starts = self.get_chain_starts()
  let n_chains = chain_starts.length()
  let result : Array[Double] = Array(capacity=n_chains)
  let mut i = 0
  while i < n_chains {
    let s = chain_starts[i]
    let e = if i + 1 < n_chains { chain_starts[i + 1] } else { self.len() }
    let mut sum = 0.0
    let mut cnt = 0
    let mut j = s
    while j < e {
      sum = sum + self.b_factors[j]
      cnt = cnt + 1
      j = j + 1
    }
    result.push(if cnt > 0 { sum / cnt.to_double() } else { 0.0 })
    i = i + 1
  }
  result
}

///|
/// Per-residue geometric center.
pub fn AtomArray::per_residue_center(self : AtomArray) -> Array[Vector3] {
  let starts = self.get_residue_starts()
  let n_res = starts.length()
  let result : Array[Vector3] = Array(capacity=n_res)
  let mut i = 0
  while i < n_res {
    let s = starts[i]
    let e = if i + 1 < n_res { starts[i + 1] } else { self.len() }
    let mut sx = 0.0
    let mut sy = 0.0
    let mut sz = 0.0
    let mut cnt = 0
    let mut j = s
    while j < e {
      sx = sx + self.coords[j].x
      sy = sy + self.coords[j].y
      sz = sz + self.coords[j].z
      cnt = cnt + 1
      j = j + 1
    }
    let nf = cnt.to_double()
    result.push(Vector3::new(sx / nf, sy / nf, sz / nf))
    i = i + 1
  }
  result
}
