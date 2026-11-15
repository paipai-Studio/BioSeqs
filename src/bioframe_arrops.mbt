/// bioframe_arrops.mbt — Pure-array interval arithmetic (no BedInterval wrapper).
///
/// Bioframe equivalent: `bioframe.arrops`.
///
/// Motivation: BedInterval is ergonomic but allocates a struct per interval.
/// When the working set is thousands of rows, a plain Int array of starts and
/// ends is both faster and smaller in memory.  The function signatures here
/// are intentionally mirror those of bioframe.arrops so that existing callers
/// of the BedInterval-layer API can drop in without a new vocabulary.
///
/// All functions take parallel `starts` / `ends` arrays and return the
/// same representation — no allocation of structs, no chromosome-aware
/// grouping (caller is expected to pre-split by chromosome if needed).

// ================================================================
// 1. merge_intervals  (pure-array)
// ================================================================

///|
/// Merge overlapping or book-ended intervals within the same chromosome.
///
/// Inputs must be **sorted by start** (this function does not sort —
/// use the BedInterval-layer `bf_merge` if your inputs aren't sorted yet).
///
/// Half-open: [s1, e1) and [s2, e2) merge iff s2 <= e1.
///
/// Returns (merged_starts, merged_ends).  Empty input → two empty arrays.
pub fn bf_arrops_merge_intervals(
  starts : Array[Int],
  ends : Array[Int],
) -> (Array[Int], Array[Int]) {
  let n = starts.length()
  if n == 0 {
    return (Array(capacity=0), Array(capacity=0))
  }
  let out_s : Array[Int] = Array(capacity=n)
  let out_e : Array[Int] = Array(capacity=n)
  let mut cur_s = starts[0]
  let mut cur_e = ends[0]
  let mut i = 1
  while i < n {
    let s = starts[i]
    let e = ends[i]
    if s <= cur_e {
      // overlap or book-ended — extend
      if e > cur_e { cur_e = e }
    } else {
      out_s.push(cur_s)
      out_e.push(cur_e)
      cur_s = s
      cur_e = e
    }
    i = i + 1
  }
  out_s.push(cur_s)
  out_e.push(cur_e)
  (out_s, out_e)
}

// ================================================================
// 2. complement_intervals
// ================================================================

///|
/// Complement intervals within a known chromosome extent.
///
/// Given a set of intervals sorted by start, return the gaps between them,
/// clipped to `[0, chrom_end)`.
///
///   starts = [100, 300]
///   ends   = [200, 400]
///   chrom_end = 500
///   → complement: [0, 100), [200, 300), [400, 500)
pub fn bf_arrops_complement_intervals(
  starts : Array[Int],
  ends : Array[Int],
  chrom_end : Int,
) -> (Array[Int], Array[Int]) {
  let n = starts.length()
  let out_s : Array[Int] = Array(capacity=n + 1)
  let out_e : Array[Int] = Array(capacity=n + 1)
  let mut prev_end = 0
  let mut i = 0
  while i < n {
    let s = starts[i]
    let e = ends[i]
    if prev_end < s && prev_end < chrom_end {
      let cs = prev_end
      let ce = if s < chrom_end { s } else { chrom_end }
      if cs < ce {
        out_s.push(cs)
        out_e.push(ce)
      }
    }
    if e > prev_end { prev_end = e }
    i = i + 1
  }
  // final tail gap
  if prev_end < chrom_end {
    out_s.push(prev_end)
    out_e.push(chrom_end)
  }
  (out_s, out_e)
}

// ================================================================
// 3. overlap_intervals  (inner join)
// ================================================================

///|
/// Pair every interval in A with every interval in B that it overlaps.
///
/// Inner join semantics: only pairs where both intervals have ≥1bp overlap
/// are returned.  A in [10, 20) overlaps B in [15, 25) → matched.
/// A in [10, 15) does NOT overlap B in [15, 25) (half-open convention).
///
/// Returns two parallel index arrays `idx_a` and `idx_b` where each
/// pair (idx_a[k], idx_b[k]) is an overlapping (A[i], B[j]) pair.
pub fn bf_arrops_overlap_intervals(
  starts_a : Array[Int],
  ends_a : Array[Int],
  starts_b : Array[Int],
  ends_b : Array[Int],
) -> (Array[Int], Array[Int]) {
  let out_a : Array[Int] = Array(capacity=0)
  let out_b : Array[Int] = Array(capacity=0)
  let na = starts_a.length()
  let nb = starts_b.length()
  let mut i = 0
  while i < na {
    let sa = starts_a[i]
    let ea = ends_a[i]
    let mut j = 0
    while j < nb {
      let sb = starts_b[j]
      let eb = ends_b[j]
      if sa < eb && sb < ea {
        out_a.push(i)
        out_b.push(j)
      }
      j = j + 1
    }
    i = i + 1
  }
  (out_a, out_b)
}

// ================================================================
// 4. closest_intervals
// ================================================================

///|
/// Find the single closest B interval for each A interval.
///
/// Distance is measured between the A interval midpoint and the B interval
/// midpoint.  If two B intervals tie on distance, the earlier (lower index)
/// one wins.  Intervals in A that don't overlap any B interval are still
/// matched to their closest neighbour (unlike overlap_intervals which only
/// returns actual overlaps).
///
/// Returns two parallel arrays:
///   - `closest_idx_b`: for each A interval, index of its closest B interval
///   - `distance`:       bp distance between the two midpoints
pub fn bf_arrops_closest_intervals(
  starts_a : Array[Int],
  ends_a : Array[Int],
  starts_b : Array[Int],
  ends_b : Array[Int],
) -> (Array[Int], Array[Int]) {
  let na = starts_a.length()
  let nb = starts_b.length()
  let out_idx : Array[Int] = Array(capacity=na)
  let out_dist : Array[Int] = Array(capacity=na)
  let mut i = 0
  while i < na {
    let sa = starts_a[i]
    let ea = ends_a[i]
    let mut best_j : Int = 0
    let mut best_dist : Int = bf_big_int()
    let mut j = 0
    while j < nb {
      let sb = starts_b[j]
      let eb = ends_b[j]
      // midpoint distance — symmetric
      let mid_a = (sa + ea) / 2
      let mid_b = (sb + eb) / 2
      let d = if mid_a > mid_b { mid_a - mid_b } else { mid_b - mid_a }
      if d < best_dist {
        best_dist = d
        best_j = j
      }
      j = j + 1
    }
    out_idx.push(best_j)
    out_dist.push(best_dist)
    i = i + 1
  }
  (out_idx, out_dist)
}

/// A sentinel "very big number" for distance initialization.
fn bf_big_int() -> Int { 1000000000 }

// ================================================================
// 5. arange_multi  (per-chromosome range expansion)
// ================================================================

///|
/// For each chromosome length, produce the equivalent of `arange(chrom_len)`
/// concatenated into a single flat offset array.
///
/// This is the pure-array equivalent of Bioframe's `arrops.arange_multi`
/// — converts per-chromosome (contig-length) info into an indexing layout
/// for coordinate-to-global-offset tables.
///
/// Example:
///   lengths = [100, 200, 50]
///   → [0..99, 0..199, 0..49]   ← 400 elements total
pub fn bf_arrops_arange_multi(lengths : Array[Int]) -> Array[Int] {
  let mut total : Int = 0
  let mut i = 0
  while i < lengths.length() {
    total = total + lengths[i]
    i = i + 1
  }
  let out : Array[Int] = Array(capacity=total)
  i = 0
  while i < lengths.length() {
    let len = lengths[i]
    let mut j : Int = 0
    while j < len {
      out.push(j)
      j = j + 1
    }
    i = i + 1
  }
  out
}

// ================================================================
// 6. interweave  (merge two arrays into interval pairs)
// ================================================================

///|
/// Interleave two equal-length arrays into an interval pair array.
///
///   starts = [10, 30, 50]
///   ends   = [20, 40, 60]
///   → [10, 20, 30, 40, 50, 60]
///
/// Useful when you have parallel start/end arrays from different sources
/// that need to be flattened.  The input lengths must match.
pub fn bf_arrops_interweave(starts : Array[Int], ends : Array[Int]) -> Array[Int] {
  let n = starts.length()
  let m = ends.length()
  if n != m {
    // Best-effort: use the shorter length
    let n2 = if n < m { n } else { m }
    let out : Array[Int] = Array(capacity=(n2 * 2))
    let mut i = 0
    while i < n2 {
      out.push(starts[i])
      out.push(ends[i])
      i = i + 1
    }
    return out
  }
  let out : Array[Int] = Array(capacity=(n * 2))
  let mut i = 0
  while i < n {
    out.push(starts[i])
    out.push(ends[i])
    i = i + 1
  }
  out
}

// ================================================================
// 7. sum_slices  (batch sum of slices through an array)
// ================================================================

///|
/// For each interval [start, end) in `slices`, compute `sum(arr[start..end])`.
///
/// Inputs:
///   - arr:   source array of values (can be Int or Double — here we stick to Int)
///   - slice_starts / slice_ends: parallel start/end arrays of half-open intervals
///
/// Returns an array of the same length as `slice_starts` with each element
/// being the sum of the corresponding slice.
///
///   arr         = [10, 20, 30, 40, 50]
///   slice_starts = [0, 2]
///   slice_ends   = [2, 5]
///   → [30, 120]   i.e. (10+20) and (30+40+50)
pub fn bf_arrops_sum_slices(
  arr : Array[Int],
  slice_starts : Array[Int],
  slice_ends : Array[Int],
) -> Array[Int] {
  // Build prefix sums once — O(N) precompute, then each slice is O(1)
  let n = arr.length()
  let prefix : Array[Int] = Array(capacity=n + 1)
  prefix.push(0)
  let mut psum : Int = 0
  let mut i = 0
  while i < n {
    psum = psum + arr[i]
    prefix.push(psum)
    i = i + 1
  }
  // Answer each query in O(1)
  let nslices = slice_starts.length()
  let out : Array[Int] = Array(capacity=nslices)
  i = 0
  while i < nslices {
    let s = slice_starts[i]
    let e = slice_ends[i]
    let mut ss = s
    let mut ee = e
    if ss < 0 { ss = 0 }
    if ee > n { ee = n }
    if ss >= ee {
      out.push(0)
    } else {
      out.push(prefix[ee] - prefix[ss])
    }
    i = i + 1
  }
  out
}
