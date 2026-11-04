/// BEDPE pairs set operations.
///
/// Bioframe equivalent: `pairs_closest`, `pairs_overlap`, `pairs_subtract`,
/// `pairs_cluster`
///
/// These operate on BedpePair (3D genomic interval pairs) for Hi-C loop analysis.

///|
/// Check if two BEDPE pairs overlap: both anchors must overlap.
pub fn bf_pairs_overlap(a : BedpePair, b : BedpePair) -> Bool {
  let anchor1_overlap = a.chrom1 == b.chrom1 &&
    a.start1 < b.end1 &&
    b.start1 < a.end1
  let anchor2_overlap = a.chrom2 == b.chrom2 &&
    a.start2 < b.end2 &&
    b.start2 < a.end2
  // Also check reverse orientation (anchor1 matches anchor2)
  let rev_overlap = a.chrom1 == b.chrom2 &&
    a.start1 < b.end2 &&
    b.start2 < a.end1 &&
    a.chrom2 == b.chrom1 &&
    a.start2 < b.end1 &&
    b.start1 < a.end2
  (anchor1_overlap && anchor2_overlap) || rev_overlap
}

///|
/// Find overlapping pairs between two sets.
/// Returns array of (index_a, index_b) pairs.
pub fn bf_pairs_overlap_pairs(
  a : Array[BedpePair],
  b : Array[BedpePair],
) -> Array[(Int32, Int32)] {
  let result : Array[(Int32, Int32)] = Array(capacity=0)
  for i = 0; i < a.length(); i = i + 1 {
    for j = 0; j < b.length(); j = j + 1 {
      if bf_pairs_overlap(a[i], b[j]) {
        result.push((i, j))
      }
    }
  }
  result
}

///|
/// Subtract pairs in b from a: keep only pairs in a that don't overlap any in b.
pub fn bf_pairs_subtract(
  a : Array[BedpePair],
  b : Array[BedpePair],
) -> Array[BedpePair] {
  let result : Array[BedpePair] = Array(capacity=0)
  for ai in a {
    let mut found = false
    for bi in b {
      if bf_pairs_overlap(ai, bi) {
        found = true
        break
      }
    }
    if !found {
      result.push(ai)
    }
  }
  result
}

///|
/// Cluster overlapping pairs: assign cluster IDs to pairs where
/// any two overlapping pairs are in the same cluster.
pub fn bf_pairs_cluster(pairs : Array[BedpePair]) -> Array[(BedpePair, Int32)] {
  let n = pairs.length()
  // Union-find
  let parent : Array[Int32] = Array::make(n, 0)
  for i = 0; i < n; i = i + 1 {
    parent[i] = i
  }
  fn find(parent : Array[Int32], x : Int32) -> Int32 {
    if parent[x] == x {
      x
    } else {
      let root = find(parent, parent[x])
      parent[x] = root
      root
    }
  }
  for i = 0; i < n; i = i + 1 {
    for j = i + 1; j < n; j = j + 1 {
      if bf_pairs_overlap(pairs[i], pairs[j]) {
        let ri = find(parent, i)
        let rj = find(parent, j)
        if ri != rj {
          parent[ri] = rj
        }
      }
    }
  }
  let result : Array[(BedpePair, Int32)] = Array(capacity=0)
  for i = 0; i < n; i = i + 1 {
    result.push((pairs[i], find(parent, i)))
  }
  result
}

///|
/// Find closest pairs in b for each pair in a.
/// Distance is the sum of distances for both anchors.
pub struct BfPairsClosestResult {
  query : BedpePair
  closest : BedpePair?
  distance : Int32
}

///|
pub fn BfPairsClosestResult::query(self : BfPairsClosestResult) -> BedpePair {
  self.query
}

///|
pub fn BfPairsClosestResult::closest(self : BfPairsClosestResult) -> BedpePair? {
  self.closest
}

///|
pub fn BfPairsClosestResult::distance(self : BfPairsClosestResult) -> Int32 {
  self.distance
}

///|
fn bf_pair_anchor_distance(
  chrom1 : String,
  start1 : Int32,
  end1 : Int32,
  chrom2 : String,
  start2 : Int32,
  end2 : Int32,
) -> Int32 {
  if chrom1 != chrom2 {
    return 2147483647
  }
  if start1 < end2 && start2 < end1 {
    return 0
  }
  if end1 <= start2 {
    return start2 - end1
  }
  start1 - end2
}

///|
pub fn bf_pairs_closest(
  a : Array[BedpePair],
  b : Array[BedpePair],
) -> Array[BfPairsClosestResult] {
  let results : Array[BfPairsClosestResult] = Array(capacity=0)
  for ai in a {
    let mut best : BedpePair? = None
    let mut best_dist : Int32 = 2147483647
    for bi in b {
      let d1 = bf_pair_anchor_distance(
        ai.chrom1,
        ai.start1,
        ai.end1,
        bi.chrom1,
        bi.start1,
        bi.end1,
      )
      let d2 = bf_pair_anchor_distance(
        ai.chrom2,
        ai.start2,
        ai.end2,
        bi.chrom2,
        bi.start2,
        bi.end2,
      )
      let total = if d1 == 2147483647 || d2 == 2147483647 {
        2147483647
      } else {
        d1 + d2
      }
      if total < best_dist {
        best_dist = total
        best = Some(bi)
      }
    }
    results.push({ query: ai, closest: best, distance: best_dist, })
  }
  results
}
