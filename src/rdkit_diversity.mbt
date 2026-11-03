/// Diversity picking: MaxMin algorithm.
///
/// RDKit equivalent: `rdkit.SimDivFilters.MaxMinPicker`
///
/// The MaxMin algorithm selects a diverse subset from a compound library:
/// at each step, pick the molecule that maximizes the minimum distance
/// to all already-selected molecules. This ensures maximum spread.

///|
/// Pick a diverse subset using the MaxMin algorithm.
///
/// fingerprints: array of fingerprint arrays (one per molecule)
/// n_to_pick: number of molecules to select
/// Returns indices of selected molecules.
pub fn rdkit_maxmin_pick(
  fingerprints : Array[Array[Int32]],
  n_to_pick : Int32,
) -> Array[Int32] {
  let n = fingerprints.length()
  if n == 0 || n_to_pick <= 0 {
    return []
  }
  if n_to_pick >= n {
    let all : Array[Int32] = Array(capacity=0)
    for i = 0; i < n; i = i + 1 {
      all.push(i)
    }
    return all
  }
  // Tanimoto distance = 1 - Tanimoto similarity
  let selected : Array[Int32] = Array(capacity=0)
  // min_dist[i] = minimum distance from molecule i to any selected molecule
  let min_dist : Array[Double] = Array::make(n, 10000000000.0)
  // Start with molecule 0 (or the one with highest diversity potential)
  selected.push(0)
  for i = 0; i < n; i = i + 1 {
    min_dist[i] = rdkit_tanimoto_distance(fingerprints[0], fingerprints[i])
  }
  // Greedy MaxMin selection
  while selected.length() < n_to_pick {
    let mut best_idx = -1
    let mut best_dist = -1.0
    for i = 0; i < n; i = i + 1 {
      // Skip already selected
      let mut is_selected = false
      for s in selected {
        if s == i {
          is_selected = true
          break
        }
      }
      if !is_selected && min_dist[i] > best_dist {
        best_dist = min_dist[i]
        best_idx = i
      }
    }
    if best_idx < 0 {
      break
    }
    selected.push(best_idx)
    // Update min distances
    for i = 0; i < n; i = i + 1 {
      let d = rdkit_tanimoto_distance(fingerprints[best_idx], fingerprints[i])
      if d < min_dist[i] {
        min_dist[i] = d
      }
    }
  }
  selected
}

///|
/// Tanimoto distance between two integer fingerprints.
/// distance = 1 - similarity
fn rdkit_tanimoto_distance(fp1 : Array[Int32], fp2 : Array[Int32]) -> Double {
  1.0 - rdkit_tanimoto_sim(fp1, fp2)
}

///|
/// Tanimoto similarity between two integer fingerprints.
/// T = |A ∩ B| / |A ∪ B| = |A ∩ B| / (|A| + |B| - |A ∩ B|)
fn rdkit_tanimoto_sim(fp1 : Array[Int32], fp2 : Array[Int32]) -> Double {
  let len = if fp1.length() < fp2.length() {
    fp1.length()
  } else {
    fp2.length()
  }
  let mut intersection = 0
  let mut union = 0
  for i = 0; i < len; i = i + 1 {
    if fp1[i] != 0 || fp2[i] != 0 {
      union = union + 1
    }
    if fp1[i] != 0 && fp2[i] != 0 {
      intersection = intersection + 1
    }
  }
  if union == 0 {
    return 0.0
  }
  intersection.to_double() / union.to_double()
}

///|
/// Sphere exclusion picking: faster alternative to MaxMin.
/// Selects molecules ensuring minimum distance threshold between picks.
pub fn rdkit_sphere_exclusion_pick(
  fingerprints : Array[Array[Int32]],
  min_distance~ : Double = 0.3,
) -> Array[Int32] {
  let n = fingerprints.length()
  if n == 0 {
    return []
  }
  let selected : Array[Int32] = Array(capacity=0)
  let excluded : Array[Bool] = Array::make(n, false)
  // Process in order, skip excluded
  for i = 0; i < n; i = i + 1 {
    if excluded[i] {
      continue
    }
    selected.push(i)
    // Exclude all molecules within min_distance
    for j = i + 1; j < n; j = j + 1 {
      if !excluded[j] {
        let d = rdkit_tanimoto_distance(fingerprints[i], fingerprints[j])
        if d < min_distance {
          excluded[j] = true
        }
      }
    }
  }
  selected
}

///|
/// Compute diversity of a selected subset (average pairwise distance).
pub fn rdkit_diversity_score(
  fingerprints : Array[Array[Int32]],
  selected_indices : Array[Int32],
) -> Double {
  let m = selected_indices.length()
  if m <= 1 {
    return 0.0
  }
  let mut total = 0.0
  let mut count = 0.0
  for i = 0; i < m; i = i + 1 {
    for j = i + 1; j < m; j = j + 1 {
      let d = rdkit_tanimoto_distance(
        fingerprints[selected_indices[i]],
        fingerprints[selected_indices[j]],
      )
      total = total + d
      count = count + 1.0
    }
  }
  if count > 0.0 {
    total / count
  } else {
    0.0
  }
}
