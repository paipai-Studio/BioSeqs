///|
/// structure_similarity.mbt — Bio.PDB.Dice / Contact / ResidueDepth
///
/// Structure-level similarity & quality metrics:
///   - Dice / Tanimoto / Jaccard coefficient (binding site similarity)
///   - Contact order (topological complexity)
///   - Residue depth (per-residue distance to nearest cross-chain CA)

// ============================================================================
// Dice / Tanimoto / Jaccard
// ============================================================================

///|
/// Collect all (i, j) residue pairs (i < j) whose CAs are within `cutoff`.
fn collect_contacts(arr : AtomArray, cutoff : Double) -> Array[(Int32, Int32)] {
  let res_starts = arr.get_residue_starts()
  let res_count = res_starts.length()
  if res_count < 2 {
    return Array(capacity=0)
  }
  let ca_coords : Array[(Bool, Vector3)] = Array(capacity=res_count)
  let mut i = 0
  while i < res_count {
    let start = res_starts[i]
    let end = if i + 1 < res_count { res_starts[i + 1] } else { arr.len() }
    let mut found = false
    let mut coord = Vector3::new(0.0, 0.0, 0.0)
    let mut j = start
    while j < end && !found {
      let name = arr.names[j]
      if name == "CA" {
        coord = arr.coords[j]
        found = true
      }
      j = j + 1
    }
    ca_coords.push((found, coord))
    i = i + 1
  }
  let r2 = cutoff * cutoff
  let pairs : Array[(Int32, Int32)] = Array(capacity=res_count * 8)
  let mut a = 0
  while a < res_count {
    let (ok_a, ca_a) = ca_coords[a]
    if !ok_a {
      a = a + 1
      continue
    }
    let mut b = a + 1
    while b < res_count {
      let (ok_b, ca_b) = ca_coords[b]
      if ok_b {
        let dx = ca_a.x - ca_b.x
        let dy = ca_a.y - ca_b.y
        let dz = ca_a.z - ca_b.z
        if dx * dx + dy * dy + dz * dz <= r2 {
          pairs.push((a, b))
        }
      }
      b = b + 1
    }
    a = a + 1
  }
  pairs
}

///|
/// Count identical pairs shared between two contact pair arrays.
fn contact_pair_intersection(
  a : Array[(Int32, Int32)],
  b : Array[(Int32, Int32)],
) -> Int32 {
  let mut count = 0
  let mut i = 0
  let mut j = 0
  let la = a.length()
  let lb = b.length()
  while i < la && j < lb {
    let (a1, a2) = a[i]
    let (b1, b2) = b[j]
    if a1 == b1 && a2 == b2 {
      count = count + 1
      i = i + 1
      j = j + 1
    } else if a1 < b1 || (a1 == b1 && a2 < b2) {
      i = i + 1
    } else {
      j = j + 1
    }
  }
  count
}

///|
/// Dice coefficient between two structures based on residue-residue
/// contacts within `cutoff` (default 5.0 Å).
///
///   Dice = 2 * |A ∩ B| / (|A| + |B|)
///
/// Ranges [0.0, 1.0]. 1.0 means identical contact patterns.
pub fn dice_similarity(
  arr_a : AtomArray,
  arr_b : AtomArray,
  cutoff~ : Double = 5.0,
) -> Double {
  let pairs_a = collect_contacts(arr_a, cutoff)
  let pairs_b = collect_contacts(arr_b, cutoff)
  let na = pairs_a.length()
  let nb = pairs_b.length()
  if na == 0 && nb == 0 {
    return 1.0
  }
  if na == 0 || nb == 0 {
    return 0.0
  }
  let inter = contact_pair_intersection(pairs_a, pairs_b)
  2.0 * inter.to_double() / (na.to_double() + nb.to_double())
}

///|
/// Tanimoto (Jaccard) coefficient between two structures based on
/// residue-residue contacts within `cutoff`.
///
///   Tanimoto = |A ∩ B| / |A ∪ B|
pub fn tanimoto_similarity(
  arr_a : AtomArray,
  arr_b : AtomArray,
  cutoff~ : Double = 5.0,
) -> Double {
  let pairs_a = collect_contacts(arr_a, cutoff)
  let pairs_b = collect_contacts(arr_b, cutoff)
  let na = pairs_a.length()
  let nb = pairs_b.length()
  if na == 0 && nb == 0 {
    return 1.0
  }
  if na == 0 || nb == 0 {
    return 0.0
  }
  let inter = contact_pair_intersection(pairs_a, pairs_b)
  let uni = na + nb - inter
  inter.to_double() / uni.to_double()
}

///|
/// Jaccard coefficient — alias for Tanimoto.
pub fn jaccard_similarity(
  arr_a : AtomArray,
  arr_b : AtomArray,
  cutoff~ : Double = 5.0,
) -> Double {
  tanimoto_similarity(arr_a, arr_b, cutoff~)
}

// ============================================================================
// Contact Order
// ============================================================================

///|
/// Topological contact order of a protein.
///
///   CO = (Σ |i − j| over all intra-chain contacts within `cutoff`)
///        / N_contacts
///
/// A measure of structural complexity (lower CO = more compact / helical).
pub fn AtomArray::contact_order(
  self : AtomArray,
  cutoff~ : Double = 5.0,
) -> Double {
  let res_starts = self.get_residue_starts()
  let res_count = res_starts.length()
  if res_count < 2 {
    return 0.0
  }
  let chain_ids = self.chain_ids
  let mut co_sum = 0.0
  let mut n_contacts = 0
  let mut i = 0
  while i < res_count {
    let start_i = res_starts[i]
    let end_i = if i + 1 < res_count { res_starts[i + 1] } else { self.len() }
    let mut ok_i = false
    let mut ca_i = Vector3::new(0.0, 0.0, 0.0)
    let mut j = start_i
    while j < end_i && !ok_i {
      if self.names[j] == "CA" {
        ca_i = self.coords[j]
        ok_i = true
      }
      j = j + 1
    }
    if !ok_i {
      i = i + 1
      continue
    }
    let chain_i = chain_ids[start_i]
    let r2 = cutoff * cutoff
    let mut a = i + 2
    while a < res_count {
      let start_a = res_starts[a]
      if chain_ids[start_a] != chain_i {
        a = a + 1
        continue
      }
      let end_a = if a + 1 < res_count { res_starts[a + 1] } else { self.len() }
      let mut ok_a = false
      let mut ca_a = Vector3::new(0.0, 0.0, 0.0)
      let mut b = start_a
      while b < end_a && !ok_a {
        if self.names[b] == "CA" {
          ca_a = self.coords[b]
          ok_a = true
        }
        b = b + 1
      }
      if ok_a {
        let dx = ca_i.x - ca_a.x
        let dy = ca_i.y - ca_a.y
        let dz = ca_i.z - ca_a.z
        if dx * dx + dy * dy + dz * dz <= r2 {
          co_sum = co_sum + (a - i).to_double()
          n_contacts = n_contacts + 1
        }
      }
      a = a + 1
    }
    i = i + 1
  }
  if n_contacts == 0 {
    0.0
  } else {
    co_sum / n_contacts.to_double()
  }
}

// ============================================================================
// Residue Depth
// ============================================================================

///|
/// Per-residue depth: distance (Å) from each residue's Cα to the
/// nearest Cα of a residue in a *different chain*.
///
/// Residues at chain interfaces have low depth; residues isolated in
/// their own chain have high depth. Missing residues / single-chain
/// setups return −1.0 for convenience (filter with `d >= 0.0`).
pub fn AtomArray::residue_depth(self : AtomArray) -> Array[Double] {
  let res_starts = self.get_residue_starts()
  let res_count = res_starts.length()
  let result : Array[Double] = Array(capacity=res_count)
  if res_count == 0 {
    return result
  }
  let chain_ids = self.chain_ids
  let per_res_chain : Array[Char] = Array(capacity=res_count)
  let per_res_ca : Array[(Bool, Vector3)] = Array(capacity=res_count)
  let mut i = 0
  while i < res_count {
    let start = res_starts[i]
    let end = if i + 1 < res_count { res_starts[i + 1] } else { self.len() }
    let chain = chain_ids[start]
    let mut ok = false
    let mut ca = Vector3::new(0.0, 0.0, 0.0)
    let mut j = start
    while j < end && !ok {
      if self.names[j] == "CA" {
        ca = self.coords[j]
        ok = true
      }
      j = j + 1
    }
    per_res_chain.push(chain)
    per_res_ca.push((ok, ca))
    i = i + 1
  }
  let mut a = 0
  while a < res_count {
    let (ok_a, ca_a) = per_res_ca[a]
    if !ok_a {
      result.push(-1.0)
      a = a + 1
      continue
    }
    let chain_a = per_res_chain[a]
    let mut best = 1.0e30
    let mut c = 0
    while c < res_count {
      let chain_c = per_res_chain[c]
      if chain_c != chain_a {
        let (ok_c, ca_c) = per_res_ca[c]
        if ok_c {
          let dx = ca_a.x - ca_c.x
          let dy = ca_a.y - ca_c.y
          let dz = ca_a.z - ca_c.z
          let d2 = dx * dx + dy * dy + dz * dz
          if d2 < best {
            best = d2
          }
        }
      }
      c = c + 1
    }
    if best >= 1.0e30 {
      result.push(-1.0)
    } else {
      result.push(best.sqrt())
    }
    a = a + 1
  }
  result
}

///|
/// Mean residue depth over residues with a valid cross-chain depth.
pub fn AtomArray::mean_residue_depth(self : AtomArray) -> Double {
  let depths = self.residue_depth()
  let mut sum = 0.0
  let mut n = 0
  for d in depths {
    if d >= 0.0 {
      sum = sum + d
      n = n + 1
    }
  }
  if n == 0 {
    0.0
  } else {
    sum / n.to_double()
  }
}
