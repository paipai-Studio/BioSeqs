///|
/// Structure quality assessment — lddt (Local Distance Difference Test)
/// per-residue and global score.
///
/// lddt is the standard CASP metric for structure prediction quality,
/// and is the basis of AlphaFold's pLDDT confidence score.
///
/// Reference: Mariani et al. (2013) "lDDT: a local superposition-free
/// score to assess the quality of protein structure predictions."

///|
/// lddt result container: per-residue scores + global average.
pub struct LddtResult {
  /// Per-residue lddt score (same length as residue count, 0.0 for skipped)
  per_residue : Array[Double]
  /// Global lddt average over all residues with at least 1 pair
  global : Double
  /// Number of residues evaluated
  n_residues : Int32
}

///|
/// Get per-residue lddt scores.
pub fn LddtResult::get_per_residue(self : LddtResult) -> Array[Double] {
  self.per_residue
}

///|
/// Get global lddt score.
pub fn LddtResult::get_global(self : LddtResult) -> Double {
  self.global
}

///|
/// Get number of residues.
pub fn LddtResult::get_n_residues(self : LddtResult) -> Int32 {
  self.n_residues
}

// --- internal helpers ---

///|
/// Find CA atom index for each residue. Returns array[Int32], -1 if missing.
fn collect_ca_indices(arr : AtomArray) -> Array[Int32] {
  let starts = arr.get_residue_starts()
  let n_res = starts.length()
  let cas : Array[Int32] = Array(capacity=n_res)
  let mut i = 0
  while i < n_res {
    let s = starts[i]
    let e = if i + 1 < n_res { starts[i + 1] } else { arr.len() }
    let mut found = -1
    let mut j = s
    while j < e {
      if arr.names[j] == "CA" {
        found = j
        break
      }
      j = j + 1
    }
    cas.push(found)
    i = i + 1
  }
  cas
}

///|
/// Compute lddt between reference and model structures.
/// Both AtomArrays are assumed to have the same residue ordering
/// (same chain_id / res_id sequence).
///
/// Uses the standard 15Å cutoff on the reference structure and
/// the four thresholds (0.5, 1.0, 2.0, 4.0 Å).
pub fn lddt(reference : AtomArray, model : AtomArray) -> LddtResult {
  let ref_starts = reference.get_residue_starts()
  let mod_starts = model.get_residue_starts()
  let n_ref = ref_starts.length()
  let n_mod = mod_starts.length()
  let n = if n_ref < n_mod { n_ref } else { n_mod }
  let cas_ref = collect_ca_indices(reference)
  let cas_mod = collect_ca_indices(model)
  let per_residue : Array[Double] = Array(capacity=n)
  // Precompute reference CA coords
  let ref_coords : Array[Vector3] = Array(capacity=n)
  let mod_coords : Array[Vector3] = Array(capacity=n)
  let mut k = 0
  while k < n {
    if cas_ref[k] >= 0 && cas_mod[k] >= 0 {
      ref_coords.push(reference.coords[cas_ref[k]])
      mod_coords.push(model.coords[cas_mod[k]])
    } else {
      ref_coords.push(Vector3::zero())
      mod_coords.push(Vector3::zero())
    }
    k = k + 1
  }
  // lddt thresholds
  let t0 = 0.5
  let t1 = 1.0
  let t2 = 2.0
  let t3 = 4.0
  let mut sum_global = 0.0
  let mut n_counted = 0
  let mut i = 0
  while i < n {
    let cas_i = cas_ref[i]
    // Skip if this residue or its model counterpart has no CA
    if cas_i < 0 || cas_mod[i] < 0 {
      per_residue.push(0.0)
      i = i + 1
      continue
    }
    let chain_i = reference.chain_ids[ref_starts[i]]
    let ri = ref_coords[i]
    let mi = mod_coords[i]
    // Accumulators for 4 thresholds
    let mut c0 = 0 // count within 0.5Å
    let mut c1 = 0 // within 1.0Å
    let mut c2 = 0 // within 2.0Å
    let mut c3 = 0 // within 4.0Å
    let mut total_pairs = 0
    let mut j = 0
    while j < n {
      if j == i {
        j = j + 1
        continue
      }
      // Same chain check
      let chain_j = reference.chain_ids[ref_starts[j]]
      if chain_i != chain_j {
        j = j + 1
        continue
      }
      let cas_j_ref = cas_ref[j]
      let cas_j_mod = cas_mod[j]
      if cas_j_ref < 0 || cas_j_mod < 0 {
        j = j + 1
        continue
      }
      let rj = ref_coords[j]
      let mj = mod_coords[j]
      // Reference distance — must be within 15Å
      let d_ref = ri.distance(rj)
      if d_ref > 15.0 {
        j = j + 1
        continue
      }
      let d_model = mi.distance(mj)
      let delta = (d_model - d_ref).abs()
      total_pairs = total_pairs + 1
      if delta < t0 {
        c0 = c0 + 1
      }
      if delta < t1 {
        c1 = c1 + 1
      }
      if delta < t2 {
        c2 = c2 + 1
      }
      if delta < t3 {
        c3 = c3 + 1
      }
      j = j + 1
    }
    if total_pairs == 0 {
      per_residue.push(0.0)
      i = i + 1
      continue
    }
    let tp = total_pairs.to_double()
    let score = (
        c0.to_double() / tp +
        c1.to_double() / tp +
        c2.to_double() / tp +
        c3.to_double() / tp
      ) /
      4.0
    per_residue.push(score)
    sum_global = sum_global + score
    n_counted = n_counted + 1
    i = i + 1
  }
  let global = if n_counted > 0 {
    sum_global / n_counted.to_double()
  } else {
    0.0
  }
  LddtResult::{ per_residue, global, n_residues: n, }
}

///|
/// Compute lddt of this AtomArray against a reference AtomArray.
pub fn AtomArray::lddt_against(
  self : AtomArray,
  reference : AtomArray,
) -> LddtResult {
  lddt(reference, self)
}
