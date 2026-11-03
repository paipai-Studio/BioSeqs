/// Find k nearest intervals.
///
/// Bioframe equivalent: `closest()` with k parameter
/// (generalization of bedtools closest -k)

///|
/// Result of k-closest query
pub struct BfClosestKResult {
  query : BedInterval
  /// Up to k nearest intervals, sorted by distance (ascending)
  hits : Array[(BedInterval, Int32)]
}

///|
pub fn BfClosestKResult::query(self : BfClosestKResult) -> BedInterval {
  self.query
}

///|
pub fn BfClosestKResult::hits(
  self : BfClosestKResult,
) -> Array[(BedInterval, Int32)] {
  self.hits
}

///|
/// Compute genomic distance between two intervals on the same chromosome.
/// Returns 0 if overlapping, otherwise the gap size.
fn bf_interval_distance(a : BedInterval, b : BedInterval) -> Int32 {
  if a.chrom != b.chrom {
    return 2147483647 // max Int32 as "infinity"
  }
  if a.start < b.end && b.start < a.end {
    return 0
  }
  if a.end <= b.start {
    return b.start - a.end
  }
  a.start - b.end
}

///|
/// For each interval in `query`, find the k nearest intervals in `subject`.
///
/// k: number of nearest neighbors to return (default 1)
/// Returns one result per query interval, each containing up to k (interval, distance) pairs.
pub fn bf_closest_k(
  query : Array[BedInterval],
  subject : Array[BedInterval],
  k~ : Int32 = 1,
) -> Array[BfClosestKResult] {
  let results : Array[BfClosestKResult] = Array(capacity=0)
  for qi in query {
    // Collect all (interval, distance) pairs
    let candidates : Array[(BedInterval, Int32)] = Array(capacity=0)
    for si in subject {
      let d = bf_interval_distance(qi, si)
      if d < 2147483647 {
        candidates.push((si, d))
      }
    }
    // Sort by distance (selection sort for simplicity, k is usually small)
    let sorted : Array[(BedInterval, Int32)] = Array(capacity=0)
    let n = candidates.length()
    let used : Array[Bool] = Array::make(n, false)
    while k > 0 && sorted.length() < n {
      let mut best_idx = -1
      let mut best_dist : Int32 = 2147483647
      for j = 0; j < n; j = j + 1 {
        if !used[j] && candidates[j].1 < best_dist {
          best_dist = candidates[j].1
          best_idx = j
        }
      }
      if best_idx < 0 {
        break
      }
      used[best_idx] = true
      sorted.push(candidates[best_idx])
      if sorted.length() >= k {
        break
      }
    }
    results.push({ query: qi, hits: sorted, })
  }
  results
}

///|
/// For each query interval, return only the distances to k nearest (no interval objects).
pub fn bf_closest_k_distances(
  query : Array[BedInterval],
  subject : Array[BedInterval],
  k~ : Int32 = 1,
) -> Array[Array[Int32]] {
  let closest = bf_closest_k(query, subject, k~)
  let result : Array[Array[Int32]] = Array(capacity=0)
  for r in closest {
    let dists : Array[Int32] = Array(capacity=0)
    for h in r.hits {
      dists.push(h.1)
    }
    result.push(dists)
  }
  result
}
