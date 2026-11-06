/// biotite_permutation.mbt — Permutation operations for MSA order optimization
/// Inspired by biotite.sequence.align.Permutation

///|
/// A permutation represents a reordering of sequence indices.
/// Used to optimize the order of sequences in progressive MSA.
pub struct BiotitePermutation {
  order : Array[Int32] // order[i] = original index of the i-th sequence
}

///|
/// Create identity permutation.
pub fn BiotitePermutation::identity(n : Int32) -> BiotitePermutation {
  let order : Array[Int32] = Array::make(n, 0)
  for i = 0; i < n; i = i + 1 {
    order[i] = i
  }
  { order, }
}

///|
/// Create permutation from an explicit order array.
pub fn BiotitePermutation::from_order(
  order : Array[Int32],
) -> BiotitePermutation {
  { order, }
}

///|
/// Apply permutation to an array.
pub fn[T] BiotitePermutation::apply(
  self : BiotitePermutation,
  arr : Array[T],
) -> Array[T] {
  let result : Array[T] = Array(capacity=0)
  for i = 0; i < self.order.length(); i = i + 1 {
    if self.order[i] < arr.length() {
      result.push(arr[self.order[i]])
    }
  }
  result
}

///|
/// Get the inverse permutation.
pub fn BiotitePermutation::inverse(
  self : BiotitePermutation,
) -> BiotitePermutation {
  let n = self.order.length()
  let inv : Array[Int32] = Array::make(n, 0)
  for i = 0; i < n; i = i + 1 {
    if self.order[i] < n {
      inv[self.order[i]] = i
    }
  }
  { order: inv, }
}

///|
/// Compose two permutations: (self ∘ other).
pub fn BiotitePermutation::compose(
  self : BiotitePermutation,
  other : BiotitePermutation,
) -> BiotitePermutation {
  let n = self.order.length()
  let result : Array[Int32] = Array::make(n, 0)
  for i = 0; i < n; i = i + 1 {
    if other.order[i] < n && self.order[other.order[i]] < n {
      result[i] = self.order[other.order[i]]
    }
  }
  { order: result, }
}

///|
/// Optimize sequence order for progressive alignment using a guide tree.
/// Sequences closer in the guide tree are aligned first.
/// Uses a simple heuristic: sort by average pairwise distance.
pub fn BiotitePermutation::from_distance_matrix(
  dist_matrix : Array[Array[Double]],
) -> BiotitePermutation {
  let n = dist_matrix.length()
  if n == 0 {
    return { order: [], }
  }
  // Calculate average distance for each sequence
  let avg_dist : Array[(Int32, Double)] = Array::make(n, (0, 0.0))
  for i = 0; i < n; i = i + 1 {
    let mut total = 0.0
    let mut count = 0
    for j = 0; j < n; j = j + 1 {
      if i != j && j < dist_matrix[i].length() {
        total = total + dist_matrix[i][j]
        count = count + 1
      }
    }
    avg_dist[i] = if count > 0 {
      (i, total / count.to_double())
    } else {
      (i, 0.0)
    }
  }
  // Sort by average distance (ascending)
  avg_dist.sort_by(fn(a, b) {
    if a.1 < b.1 {
      -1
    } else if a.1 > b.1 {
      1
    } else {
      0
    }
  })
  let order : Array[Int32] = Array::make(n, 0)
  for i = 0; i < n; i = i + 1 {
    order[i] = avg_dist[i].0
  }
  { order, }
}

///|
/// Convert permutation to a series of transpositions (swaps).
/// Returns array of (i, j) pairs meaning "swap positions i and j".
pub fn BiotitePermutation::to_transpositions(
  self : BiotitePermutation,
) -> Array[(Int32, Int32)] {
  let n = self.order.length()
  let current : Array[Int32] = Array::make(n, 0)
  for i = 0; i < n; i = i + 1 {
    current[i] = i
  }
  let swaps : Array[(Int32, Int32)] = Array(capacity=0)
  for target = 0; target < n; target = target + 1 {
    // Find where the target element currently is
    let mut pos = target
    for i = target; i < n; i = i + 1 {
      if current[i] == self.order[target] {
        pos = i
        break
      }
    }
    if pos != target {
      // Swap
      let tmp = current[target]
      current[target] = current[pos]
      current[pos] = tmp
      swaps.push((target, pos))
    }
  }
  swaps
}

///|
/// Check if permutation is the identity.
pub fn BiotitePermutation::is_identity(self : BiotitePermutation) -> Bool {
  for i = 0; i < self.order.length(); i = i + 1 {
    if self.order[i] != i {
      return false
    }
  }
  true
}
