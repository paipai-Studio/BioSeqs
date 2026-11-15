///|
/// Tests for bioframe_arrops.mbt — pure-array interval arithmetic.

// ---- merge_intervals ------------------------------------------------------

test "arrops_merge: empty input → empty" {
  let (s, e) = @src.bf_arrops_merge_intervals([], [])
  assert_eq(s.length(), 0)
  assert_eq(e.length(), 0)
}

test "arrops_merge: no overlap → unchanged" {
  let starts = [10, 50, 100]
  let ends = [20, 60, 110]
  let (s, e) = @src.bf_arrops_merge_intervals(starts, ends)
  assert_eq(s.length(), 3)
  assert_eq(s[0], 10); assert_eq(e[0], 20)
  assert_eq(s[1], 50); assert_eq(e[1], 60)
  assert_eq(s[2], 100); assert_eq(e[2], 110)
}

test "arrops_merge: fully overlapping → single interval" {
  let starts = [0, 20, 50]
  let ends = [100, 80, 70]
  let (s, e) = @src.bf_arrops_merge_intervals(starts, ends)
  assert_eq(s.length(), 1)
  assert_eq(s[0], 0); assert_eq(e[0], 100)
}

test "arrops_merge: book-ended intervals merge" {
  // [0, 100) and [100, 200): s2=100 <= e1=100 → merge
  let (s, e) = @src.bf_arrops_merge_intervals([0, 100], [100, 200])
  assert_eq(s.length(), 1)
  assert_eq(s[0], 0); assert_eq(e[0], 200)
}

test "arrops_merge: multiple partial overlaps chain" {
  // [0,50), [30,80), [70,120) → [0,120)
  let (s, e) = @src.bf_arrops_merge_intervals([0, 30, 70], [50, 80, 120])
  assert_eq(s.length(), 1)
  assert_eq(s[0], 0); assert_eq(e[0], 120)
}

test "arrops_merge: some merge some don't" {
  // [0,50) (from [0,30)+[20,50)), [100,170) (from [100,130)+[120,170)) → 2 merged
  let starts = [0, 20, 100, 120]
  let ends = [30, 50, 130, 170]
  let (s, e) = @src.bf_arrops_merge_intervals(starts, ends)
  assert_eq(s.length(), 2)
  assert_eq(s[0], 0); assert_eq(e[0], 50)
  assert_eq(s[1], 100); assert_eq(e[1], 170)
}

// ---- complement_intervals -------------------------------------------------

test "arrops_complement: empty input → full range" {
  let (s, e) = @src.bf_arrops_complement_intervals([], [], 500)
  assert_eq(s.length(), 1)
  assert_eq(s[0], 0); assert_eq(e[0], 500)
}

test "arrops_complement: single interval in middle" {
  let (s, e) = @src.bf_arrops_complement_intervals([200], [300], 500)
  assert_eq(s.length(), 2)
  assert_eq(s[0], 0); assert_eq(e[0], 200)
  assert_eq(s[1], 300); assert_eq(e[1], 500)
}

test "arrops_complement: multiple intervals → gaps" {
  let starts = [100, 300]
  let ends = [200, 400]
  let (s, e) = @src.bf_arrops_complement_intervals(starts, ends, 500)
  assert_eq(s.length(), 3)
  assert_eq(s[0], 0); assert_eq(e[0], 100)
  assert_eq(s[1], 200); assert_eq(e[1], 300)
  assert_eq(s[2], 400); assert_eq(e[2], 500)
}

test "arrops_complement: intervals cover entire range → empty" {
  let (s, e) = @src.bf_arrops_complement_intervals([0, 500], [500, 1000], 1000)
  assert_eq(s.length(), 0)
}

test "arrops_complement: clips to chrom_end" {
  let (s, e) = @src.bf_arrops_complement_intervals([0], [600], 500)
  assert_eq(s.length(), 0)
}

// ---- overlap_intervals ----------------------------------------------------

test "arrops_overlap: empty inputs → empty" {
  let (ia, ib) = @src.bf_arrops_overlap_intervals([], [], [], [])
  assert_eq(ia.length(), 0)
}

test "arrops_overlap: no overlap → empty" {
  let (ia, ib) = @src.bf_arrops_overlap_intervals(
    [0, 100], [50, 150],
    [200, 300], [250, 350],
  )
  assert_eq(ia.length(), 0)
}

test "arrops_overlap: half-open boundary does NOT count" {
  // A=[10,15), B=[15,25): 15 is exclusive → no overlap
  let (ia, ib) = @src.bf_arrops_overlap_intervals([10], [15], [15], [25])
  assert_eq(ia.length(), 0)
}

test "arrops_overlap: inner join with multiple matches" {
  let (ia, ib) = @src.bf_arrops_overlap_intervals(
    [10, 50, 100], [60, 80, 150],
    [30, 70], [55, 90],
  )
  // A[0]=[10,60) overlaps B[0]=[30,55): yes
  // A[0]=[10,60) overlaps B[1]=[70,90)? 60<=70: no
  // A[1]=[50,80) overlaps B[0]=[30,55): yes
  // A[1]=[50,80) overlaps B[1]=[70,90): yes
  // A[2]=[100,150) overlaps nothing
  assert_eq(ia.length(), 3)
  assert_eq(ia[0], 0); assert_eq(ib[0], 0)
  assert_eq(ia[1], 1); assert_eq(ib[1], 0)
  assert_eq(ia[2], 1); assert_eq(ib[2], 1)
}

// ---- closest_intervals ----------------------------------------------------

test "arrops_closest: single pair" {
  let (idx, dist) = @src.bf_arrops_closest_intervals([10], [20], [100], [200])
  // mid_a=15, mid_b=150 → dist=135
  assert_eq(idx.length(), 1)
  assert_eq(dist.length(), 1)
  assert_eq(idx[0], 0)
  assert_eq(dist[0], 135)
}

test "arrops_closest: ties → earlier wins" {
  let (idx, dist) = @src.bf_arrops_closest_intervals(
    [50], [150],          // mid=100
    [0, 200], [100, 300], // mids=50, 250
  )
  assert_eq(idx[0], 0)  // both distance 50 → first wins
  assert_eq(dist[0], 50)
}

test "arrops_closest: overlapping still computes midpoint dist" {
  let (idx, dist) = @src.bf_arrops_closest_intervals([0], [100], [25], [75])
  assert_eq(dist[0], 0)
}

// ---- arange_multi ---------------------------------------------------------

test "arange_multi: lengths → flat offset array" {
  let r = @src.bf_arrops_arange_multi([3, 2])
  assert_eq(r.length(), 5)
  assert_eq(r[0], 0); assert_eq(r[1], 1); assert_eq(r[2], 2)
  assert_eq(r[3], 0); assert_eq(r[4], 1)
}

test "arange_multi: empty → empty" {
  let r = @src.bf_arrops_arange_multi([])
  assert_eq(r.length(), 0)
}

test "arange_multi: zero length skipped" {
  let r = @src.bf_arrops_arange_multi([0, 3, 0, 2])
  assert_eq(r.length(), 5)
  assert_eq(r[0], 0); assert_eq(r[1], 1); assert_eq(r[2], 2)
  assert_eq(r[3], 0); assert_eq(r[4], 1)
}

// ---- interweave -----------------------------------------------------------

test "arrops_interweave: equal length" {
  let r = @src.bf_arrops_interweave([10, 30, 50], [20, 40, 60])
  assert_eq(r.length(), 6)
  assert_eq(r[0], 10); assert_eq(r[1], 20)
  assert_eq(r[2], 30); assert_eq(r[3], 40)
  assert_eq(r[4], 50); assert_eq(r[5], 60)
}

test "arrops_interweave: unequal → truncated to shorter" {
  let r = @src.bf_arrops_interweave([1, 2, 3, 4], [10, 20])
  assert_eq(r.length(), 4)
  assert_eq(r[0], 1); assert_eq(r[1], 10)
  assert_eq(r[2], 2); assert_eq(r[3], 20)
}

// ---- sum_slices -----------------------------------------------------------

test "arrops_sum_slices: basic prefix sum" {
  let arr = [10, 20, 30, 40, 50]
  let r = @src.bf_arrops_sum_slices(arr, [0, 2], [2, 5])
  assert_eq(r.length(), 2)
  assert_eq(r[0], 30)   // 10 + 20
  assert_eq(r[1], 120)  // 30 + 40 + 50
}

test "arrops_sum_slices: OOB corrected to [0, n]" {
  let arr = [1, 2, 3, 4, 5]
  let r = @src.bf_arrops_sum_slices(arr, [-1], [10])
  assert_eq(r.length(), 1)
  assert_eq(r[0], 15)
}

test "arrops_sum_slices: empty slices → empty" {
  let r = @src.bf_arrops_sum_slices([1, 2, 3], [], [])
  assert_eq(r.length(), 0)
}

test "arrops_sum_slices: ss >= ee → 0" {
  let arr = [1, 2, 3, 4, 5]
  let r = @src.bf_arrops_sum_slices(arr, [3], [2])
  assert_eq(r.length(), 1)
  assert_eq(r[0], 0)
}
