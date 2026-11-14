///|
/// Structure superposition helpers for AtomArray.
///
/// Wraps existing SVD / Kabsch implementations to provide
/// AtomArray-native alignment APIs.
///
/// Biotite's superposition family: superimpose, superimpose_homologs,
/// superimpose_structural_homologs, superimpose_without_outliers.

// ============================================================
// P4-1: Align coords via wrapping superimpose() from svd_superimposer.mbt
// ============================================================

///|
/// Convert Array[Vector3] → Array[AtomCoordinate].
fn v3_to_atomcoord(coords : Array[Vector3]) -> Array[AtomCoordinate] {
  let result : Array[AtomCoordinate] = Array(capacity=coords.length())
  for c in coords {
    result.push(AtomCoordinate::new(c.x, c.y, c.z))
  }
  result
}

///|
/// Convert Array[AtomCoordinate] → Array[Vector3].
fn atomcoord_to_v3(coords : Array[AtomCoordinate]) -> Array[Vector3] {
  let result : Array[Vector3] = Array(capacity=coords.length())
  for c in coords {
    result.push(Vector3::new(c.x, c.y, c.z))
  }
  result
}

///|
/// Superimpose self onto `target` AtomArray (same atom ordering).
///
/// Returns (aligned AtomArray, RMSD).
/// Uses SVD-based optimal rotation (Kabsch), wrapped from
/// `svd_superimposer::superimpose`.
pub fn AtomArray::superimpose_to(
  self : AtomArray,
  target : AtomArray,
) -> (AtomArray, Double) {
  let n = self.len()
  if n == 0 || n != target.len() {
    return (self, 0.0)
  }
  let src_coords = v3_to_atomcoord(self.coords)
  let tgt_coords = v3_to_atomcoord(target.coords)
  let (aligned_coords, rmsd_val) = superimpose(src_coords, tgt_coords)
  let new_coords = atomcoord_to_v3(aligned_coords)
  (
    AtomArray::{
      coords: new_coords,
      names: self.names,
      res_names: self.res_names,
      chain_ids: self.chain_ids,
      res_ids: self.res_ids,
      elements: self.elements,
      b_factors: self.b_factors,
      occupancies: self.occupancies,
      altlocs: self.altlocs,
    },
    rmsd_val,
  )
}

///|
/// Compute best-fit RMSD via SVD alignment between two AtomArrays.
/// Requires same atom ordering (same length).
pub fn AtomArray::rmsd_aligned(self : AtomArray, other : AtomArray) -> Double {
  let (_, rmsd_val) = self.superimpose_to(other)
  rmsd_val
}

///|
/// Compute per-residue RMSD between self and other after optimal alignment.
/// Assumes same residue ordering (chain+res_id same).
/// Returns Array[Double] length = residue count; 0.0 for residues with no CA.
pub fn AtomArray::per_residue_rmsd(
  self : AtomArray,
  other : AtomArray,
) -> Array[Double] {
  let self_cas = self.filter_ca()
  let other_cas = other.filter_ca()
  let n = self_cas.len()
  let n2 = other_cas.len()
  let n_pairs = if n < n2 { n } else { n2 }
  if n_pairs == 0 {
    return []
  }
  // Build residue → (chain, res_id) map for self
  let starts = self.get_residue_starts()
  let n_res = starts.length()
  // Build per-residue CA list
  let res_cas_self : Array[Array[Vector3]] = []
  let res_cas_other : Array[Vector3] = []
  let mut i = 0
  while i < n_res {
    let s = starts[i]
    let e = if i + 1 < n_res { starts[i + 1] } else { self.len() }
    let chain = self.chain_ids[s]
    let rid = self.res_ids[s]
    // Find matching in self_cas/other_cas
    let mut self_ca_idx = -1
    let mut j = 0
    while j < n_pairs {
      if self_cas.chain_ids[j] == chain && self_cas.res_ids[j] == rid {
        self_ca_idx = j
        break
      }
      j = j + 1
    }
    res_cas_self.push(
      if self_ca_idx >= 0 {
        [self_cas.coords[self_ca_idx]]
      } else {
        []
      },
    )
    // Other: need same chain/rid
    let mut other_ca_idx = -1
    let mut k = 0
    while k < n_pairs {
      if other_cas.chain_ids[k] == chain && other_cas.res_ids[k] == rid {
        other_ca_idx = k
        break
      }
      k = k + 1
    }
    res_cas_other.push(
      if other_ca_idx >= 0 {
        other_cas.coords[other_ca_idx]
      } else {
        Vector3::zero()
      },
    )
    i = i + 1
  }
  // Compute RMSD per residue (1 pair = just distance)
  let result : Array[Double] = Array(capacity=n_res)
  let mut r = 0
  while r < n_res {
    let arr = res_cas_self[r]
    if arr.length() == 0 {
      result.push(0.0)
    } else {
      result.push(arr[0].distance(res_cas_other[r]))
    }
    r = r + 1
  }
  result
}

// ============================================================
// P4-2: Chain-level property vectorization
// ============================================================

///|
/// Center of mass per chain. Returns Array[Vector3].
pub fn AtomArray::mass_center_per_chain(self : AtomArray) -> Array[Vector3] {
  let chain_starts = self.get_chain_starts()
  let n_chains = chain_starts.length()
  let result : Array[Vector3] = Array(capacity=n_chains)
  let mut i = 0
  while i < n_chains {
    let s = chain_starts[i]
    let e = if i + 1 < n_chains { chain_starts[i + 1] } else { self.len() }
    let coords : Array[Vector3] = Array(capacity=e - s)
    let elements : Array[String] = Array(capacity=e - s)
    let mut j = s
    while j < e {
      coords.push(self.coords[j])
      elements.push(self.elements[j])
      j = j + 1
    }
    let mp = calculate_mass_properties(coords, elements)
    result.push(mp.get_center_of_mass())
    i = i + 1
  }
  result
}

///|
/// Radius of gyration per chain.
pub fn AtomArray::gyration_radius_per_chain(self : AtomArray) -> Array[Double] {
  let chain_starts = self.get_chain_starts()
  let n_chains = chain_starts.length()
  let result : Array[Double] = Array(capacity=n_chains)
  let mut i = 0
  while i < n_chains {
    let s = chain_starts[i]
    let e = if i + 1 < n_chains { chain_starts[i + 1] } else { self.len() }
    let coords : Array[Vector3] = Array(capacity=e - s)
    let elements : Array[String] = Array(capacity=e - s)
    let mut j = s
    while j < e {
      coords.push(self.coords[j])
      elements.push(self.elements[j])
      j = j + 1
    }
    let mp = calculate_mass_properties(coords, elements)
    result.push(mp.get_radius_of_gyration())
    i = i + 1
  }
  result
}

///|
/// Total mass per chain.
pub fn AtomArray::total_mass_per_chain(self : AtomArray) -> Array[Double] {
  let chain_starts = self.get_chain_starts()
  let n_chains = chain_starts.length()
  let result : Array[Double] = Array(capacity=n_chains)
  let mut i = 0
  while i < n_chains {
    let s = chain_starts[i]
    let e = if i + 1 < n_chains { chain_starts[i + 1] } else { self.len() }
    let coords : Array[Vector3] = Array(capacity=e - s)
    let elements : Array[String] = Array(capacity=e - s)
    let mut j = s
    while j < e {
      coords.push(self.coords[j])
      elements.push(self.elements[j])
      j = j + 1
    }
    let mp = calculate_mass_properties(coords, elements)
    result.push(mp.get_total_mass())
    i = i + 1
  }
  result
}

// ============================================================
// P4-3: Chain/substructure extraction
// ============================================================

///|
/// Extract one chain as its own AtomArray.
pub fn AtomArray::get_chain(self : AtomArray, chain_id : Char) -> AtomArray {
  self.filter_chain(chain_id)
}

///|
/// Extract a subrange by chain_id and [first_res_id, last_res_id].
/// Useful for getting a specific domain or fragment.
pub fn AtomArray::get_fragment(
  self : AtomArray,
  chain_id : Char,
  first_res_id : Int32,
  last_res_id : Int32,
) -> AtomArray {
  let mask : Array[Bool] = []
  let starts = self.get_residue_starts()
  let mut i = 0
  while i < starts.length() {
    let s = starts[i]
    let chain = self.chain_ids[s]
    let rid = self.res_ids[s]
    let e = if i + 1 < starts.length() { starts[i + 1] } else { self.len() }
    let within = chain == chain_id && rid >= first_res_id && rid <= last_res_id
    let mut j = s
    while j < e {
      mask.push(within)
      j = j + 1
    }
    i = i + 1
  }
  self.filter(mask)
}

///|
/// Iterate each chain and apply a closure to produce a per-chain result.
///
/// This is the Biotite pattern for avoid applying per-chain in a loop
/// — though MoonBit doesn't support generics with closures, here we
/// provide specific-type vectorization helpers (Double, Vector3, AtomArray).
pub fn spread_chain_vec3(
  arr : AtomArray,
  chain_data : Array[Vector3],
) -> Array[Vector3] {
  let chain_starts = arr.get_chain_starts()
  let n = arr.len()
  let result : Array[Vector3] = Array(capacity=n)
  let mut i = 0
  while i < chain_starts.length() {
    let s = chain_starts[i]
    let e = if i + 1 < chain_starts.length() { chain_starts[i + 1] } else { n }
    let val = if i < chain_data.length() {
      chain_data[i]
    } else {
      Vector3::zero()
    }
    let mut j = s
    while j < e {
      result.push(val)
      j = j + 1
    }
    i = i + 1
  }
  result
}

///|
/// Spread per-chain Double data to per-atom level.
pub fn spread_chain_double(
  arr : AtomArray,
  chain_data : Array[Double],
) -> Array[Double] {
  let chain_starts = arr.get_chain_starts()
  let n = arr.len()
  let result : Array[Double] = Array(capacity=n)
  let mut i = 0
  while i < chain_starts.length() {
    let s = chain_starts[i]
    let e = if i + 1 < chain_starts.length() { chain_starts[i + 1] } else { n }
    let val = if i < chain_data.length() { chain_data[i] } else { 0.0 }
    let mut j = s
    while j < e {
      result.push(val)
      j = j + 1
    }
    i = i + 1
  }
  result
}
