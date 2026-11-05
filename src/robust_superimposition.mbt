///|
/// Outlier-robust superimposition — iteratively removes conformational
/// outliers (flexible regions) to find the optimal structural alignment
/// of the rigid core.
///
/// Biotite provides `superimpose_without_outliers()` and
/// `superimpose_homologs()` for this purpose. This module implements
/// an iterative approach:
/// 1. Superimpose all atoms
/// 2. Compute per-atom deviations
/// 3. Remove atoms with deviation > threshold
/// 4. Repeat until convergence
///
/// Reference: Biotite `superimpose_without_outliers` documentation.

///|
/// Result of outlier-robust superimposition.
pub struct RobustSuperimpositionResult {
  transformed : Array[AtomCoordinate]
  rmsd : Double
  kept_indices : Array[Int32]
  removed_indices : Array[Int32]
  n_iterations : Int32
}

///|
/// Get the transformed coordinates.
pub fn RobustSuperimpositionResult::get_transformed(
  self : RobustSuperimpositionResult,
) -> Array[AtomCoordinate] {
  self.transformed
}

///|
/// Get the final RMSD over the kept atoms.
pub fn RobustSuperimpositionResult::get_rmsd(
  self : RobustSuperimpositionResult,
) -> Double {
  self.rmsd
}

///|
/// Get indices of atoms kept in the final alignment.
pub fn RobustSuperimpositionResult::get_kept_indices(
  self : RobustSuperimpositionResult,
) -> Array[Int32] {
  self.kept_indices
}

///|
/// Get indices of atoms removed as outliers.
pub fn RobustSuperimpositionResult::get_removed_indices(
  self : RobustSuperimpositionResult,
) -> Array[Int32] {
  self.removed_indices
}

///|
/// Get the number of iterations performed.
pub fn RobustSuperimpositionResult::get_n_iterations(
  self : RobustSuperimpositionResult,
) -> Int32 {
  self.n_iterations
}

///|
/// Compute per-atom squared deviation after superimposition.
fn per_atom_deviations(
  transformed : Array[AtomCoordinate],
  target : Array[AtomCoordinate],
) -> Array[Double] {
  let deviations : Array[Double] = []
  for i in 0..<transformed.length() {
    let d = transformed[i].distance(target[i])
    deviations.push(d * d)
  }
  deviations
}

///|
/// Compute mean of an array.
fn mean(arr : Array[Double]) -> Double {
  if arr.length() == 0 {
    return 0.0
  }
  let mut sum = 0.0
  for v in arr {
    sum = sum + v
  }
  sum / arr.length().to_double()
}

///|
/// Compute standard deviation of an array.
fn stddev(arr : Array[Double]) -> Double {
  if arr.length() <= 1 {
    return 0.0
  }
  let m = mean(arr)
  let mut sum_sq = 0.0
  for v in arr {
    let diff = v - m
    sum_sq = sum_sq + diff * diff
  }
  (sum_sq / (arr.length() - 1).to_double()).sqrt()
}

///|
/// Superimpose without outliers — iteratively removes atoms whose
/// deviation exceeds `threshold_sigma` standard deviations from the
/// mean, re-superimposing after each removal.
///
/// Parameters:
/// - `threshold_sigma`: outlier cutoff in units of sigma (default 2.5)
/// - `max_iterations`: maximum iterative removal rounds (default 10)
/// - `min_atoms`: minimum atoms to keep (default 3)
pub fn superimpose_without_outliers(
  atoms1 : Array[AtomCoordinate],
  atoms2 : Array[AtomCoordinate],
  threshold_sigma~ : Double = 2.5,
  max_iterations~ : Int32 = 10,
  min_atoms~ : Int32 = 3,
) -> RobustSuperimpositionResult {
  let n = atoms1.length()
  if n < min_atoms {
    let (transformed, rmsd) = superimpose(atoms1, atoms2)
    return RobustSuperimpositionResult::{
      transformed,
      rmsd,
      kept_indices: {
        let a : Array[Int32] = []
        for i in 0..<n {
          a.push(i)
        }
        a
      },
      removed_indices: [],
      n_iterations: 0,
    }
  }
  // Start with all atoms
  let mut kept : Array[Int32] = []
  for i in 0..<n {
    kept.push(i)
  }
  let removed : Array[Int32] = []
  let mut iteration = 0

  while iteration < max_iterations && kept.length() > min_atoms {
    // Extract subsets
    let subset1 : Array[AtomCoordinate] = []
    let subset2 : Array[AtomCoordinate] = []
    for idx in kept {
      subset1.push(atoms1[idx])
      subset2.push(atoms2[idx])
    }
    // Superimpose
    let (transformed, _) = superimpose(subset1, subset2)

    // Compute per-atom deviations
    let deviations = per_atom_deviations(transformed, subset2)
    let m = mean(deviations)
    let sd = stddev(deviations)
    if sd < 0.0000000001 {
      // All deviations are the same — converged
      break
    }
    // Identify outliers
    let threshold = m + threshold_sigma * sd
    let new_kept : Array[Int32] = []
    for i in 0..<kept.length() {
      if deviations[i] <= threshold {
        new_kept.push(kept[i])
      } else {
        removed.push(kept[i])
      }
    }
    if new_kept.length() == kept.length() {
      // No outliers found — converged
      break
    }
    if new_kept.length() < min_atoms {
      // Would remove too many atoms — stop
      break
    }
    kept = new_kept
    iteration = iteration + 1
  }
  // Final superimposition with kept atoms
  let subset1 : Array[AtomCoordinate] = []
  let subset2 : Array[AtomCoordinate] = []
  for idx in kept {
    subset1.push(atoms1[idx])
    subset2.push(atoms2[idx])
  }
  let (final_transformed, final_rmsd) = superimpose(subset1, subset2)
  RobustSuperimpositionResult::{
    transformed: final_transformed,
    rmsd: final_rmsd,
    kept_indices: kept,
    removed_indices: removed,
    n_iterations: iteration,
  }
}

///|
/// Superimpose homologs — handles sequences of different lengths by
/// first finding the correspondence (matching residues by index), then
/// applying outlier-robust superimposition.
///
/// This is a simplified version that assumes atoms1 and atoms2 are
/// already corresponding (same length). For sequences of different
/// lengths, the caller should first align and extract matching positions.
pub fn superimpose_homologs(
  atoms1 : Array[AtomCoordinate],
  atoms2 : Array[AtomCoordinate],
  threshold_sigma~ : Double = 3.0,
  max_iterations~ : Int32 = 10,
) -> RobustSuperimpositionResult {
  // For same-length sequences, delegate to outlier-robust
  superimpose_without_outliers(
    atoms1,
    atoms2,
    threshold_sigma~,
    max_iterations~,
  )
}
