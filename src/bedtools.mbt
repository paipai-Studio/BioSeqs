///|
/// bedtools-style genome interval arithmetic: an `IntervalTree` data structure
/// plus `merge` (book-ended union) and `intersect` (≥1bp overlap) operations.
///
/// Semantics follow bedtools defaults (0-based half-open coordinates):
/// `merge` combines overlapping or "book-ended" intervals (0 bp gap); `intersect`
/// reports an A interval if it shares ≥1 bp with any B interval.

///|
pub struct BedInterval {
  chrom : String
  start : Int
  end : Int
}

///|
pub fn BedInterval::new(chrom : String, start : Int, end : Int) -> BedInterval {
  BedInterval::{ chrom, start, end }
}

///|
enum ItreeNode {
  ItreeNil
  ItreeNode(Int, ItreeNode, ItreeNode, Array[Int], Array[Int])
}

///|
pub struct IntervalTree {
  root : ItreeNode
  starts : Array[Int]
  ends : Array[Int]
}

///|
fn itree_sort_by_start(starts : Array[Int], idxs : Array[Int]) -> Array[Int] {
  let out = idxs.copy()
  let mut i = 1
  while i < out.length() {
    let key = out[i]
    let mut j = i - 1
    while j >= 0 && starts[out[j]] > starts[key] {
      out[j + 1] = out[j]
      j = j - 1
    }
    out[j + 1] = key
    i = i + 1
  }
  out
}

///|
fn itree_sort_by_end_desc(ends : Array[Int], idxs : Array[Int]) -> Array[Int] {
  let out = idxs.copy()
  let mut i = 1
  while i < out.length() {
    let key = out[i]
    let mut j = i - 1
    while j >= 0 && ends[out[j]] < ends[key] {
      out[j + 1] = out[j]
      j = j - 1
    }
    out[j + 1] = key
    i = i + 1
  }
  out
}

///|
fn itree_build(
  starts : Array[Int],
  ends : Array[Int],
  idxs : Array[Int],
) -> ItreeNode {
  let n = idxs.length()
  if n == 0 {
    return ItreeNil
  }
  let sorted = itree_sort_by_start(starts, idxs)
  let center = starts[sorted[n / 2]]
  let left_idxs : Array[Int] = Array::new()
  let right_idxs : Array[Int] = Array::new()
  let span_idxs : Array[Int] = Array::new()
  for idx in sorted {
    let s = starts[idx]
    let e = ends[idx]
    if e <= center {
      left_idxs.push(idx)
    } else if s > center {
      right_idxs.push(idx)
    } else {
      span_idxs.push(idx) // s <= center < e
    }
  }
  let by_start = itree_sort_by_start(starts, span_idxs)
  let by_end = itree_sort_by_end_desc(ends, span_idxs)
  ItreeNode(
    center,
    itree_build(starts, ends, left_idxs),
    itree_build(starts, ends, right_idxs),
    by_start,
    by_end,
  )
}

///|
/// Build a centered interval tree from parallel `starts`/`ends` arrays
/// (half-open `[start, end)` intervals).
pub fn interval_tree_new(starts : Array[Int], ends : Array[Int]) -> IntervalTree {
  let idxs : Array[Int] = Array::new()
  for i = 0; i < starts.length(); i = i + 1 {
    idxs.push(i)
  }
  IntervalTree::{ root: itree_build(starts, ends, idxs), starts, ends }
}

///|
fn itree_query(
  node : ItreeNode,
  starts : Array[Int],
  ends : Array[Int],
  a : Int,
  b : Int,
  out : Array[Int],
) -> Unit {
  match node {
    ItreeNil => ()
    ItreeNode(center, left, right, by_start, by_end) => {
      if b <= center {
        let mut k = 0
        while k < by_start.length() && starts[by_start[k]] < b {
          out.push(by_start[k])
          k = k + 1
        }
        itree_query(left, starts, ends, a, b, out)
      } else if a >= center {
        let mut k = 0
        while k < by_end.length() && ends[by_end[k]] > a {
          out.push(by_end[k])
          k = k + 1
        }
        itree_query(right, starts, ends, a, b, out)
      } else {
        for idx in by_start {
          out.push(idx)
        }
        itree_query(left, starts, ends, a, b, out)
        itree_query(right, starts, ends, a, b, out)
      }
    }
  }
}

///|
/// Indices of intervals overlapping the half-open interval `[a, b)`.
pub fn IntervalTree::query(self : IntervalTree, a : Int, b : Int) -> Array[Int] {
  let out : Array[Int] = Array::new()
  itree_query(self.root, self.starts, self.ends, a, b, out)
  out
}

///|
fn bedinterval_less(a : BedInterval, b : BedInterval) -> Bool {
  if a.chrom != b.chrom {
    a.chrom < b.chrom
  } else if a.start != b.start {
    a.start < b.start
  } else {
    a.end < b.end
  }
}

///|
/// Sort intervals by (chrom, start, end).
fn sort_bed_intervals(intervals : Array[BedInterval]) -> Array[BedInterval] {
  let out = intervals.copy()
  let mut i = 1
  while i < out.length() {
    let key = out[i]
    let mut j = i - 1
    while j >= 0 && bedinterval_less(key, out[j]) {
      out[j + 1] = out[j]
      j = j - 1
    }
    out[j + 1] = key
    i = i + 1
  }
  out
}

///|
/// Merge overlapping or book-ended intervals (per chromosome). Matching
/// `bedtools merge -d 0`: adjacent intervals with 0 bp gap are merged.
pub fn bedtools_merge(intervals : Array[BedInterval]) -> Array[BedInterval] {
  let sorted = sort_bed_intervals(intervals)
  let out : Array[BedInterval] = Array::new()
  let mut cur : BedInterval? = None
  for iv in sorted {
    match cur {
      None => cur = Some(iv)
      Some(c) => {
        if iv.chrom == c.chrom && iv.start <= c.end {
          let ne = if iv.end > c.end { iv.end } else { c.end }
          cur = Some(BedInterval::{ chrom: c.chrom, start: c.start, end: ne })
        } else {
          out.push(c)
          cur = Some(iv)
        }
      }
    }
  }
  match cur {
    Some(c) => out.push(c)
    None => ()
  }
  out
}

///|
/// Report A intervals that overlap any B interval by ≥1 bp (same chromosome).
/// Matching `bedtools intersect` default (`-wa`-style full-A output).
pub fn bedtools_intersect(
  a : Array[BedInterval],
  b : Array[BedInterval],
) -> Array[BedInterval] {
  // group B indices by chrom
  let chroms : Array[String] = Array::new()
  let b_per_chrom : Array[Array[Int]] = Array::new()
  for i = 0; i < b.length(); i = i + 1 {
    let c = b[i].chrom
    let mut found = -1
    for k = 0; k < chroms.length(); k = k + 1 {
      if chroms[k] == c {
        found = k
        break
      }
    }
    if found < 0 {
      chroms.push(c)
      b_per_chrom.push(Array::new())
      b_per_chrom[b_per_chrom.length() - 1].push(i)
    } else {
      b_per_chrom[found].push(i)
    }
  }
  // build one interval tree per chrom
  let tree_chroms : Array[String] = Array::new()
  let trees : Array[IntervalTree] = Array::new()
  for k = 0; k < chroms.length(); k = k + 1 {
    let idxs = b_per_chrom[k]
    let starts : Array[Int] = Array::new()
    let ends : Array[Int] = Array::new()
    for i in idxs {
      starts.push(b[i].start)
      ends.push(b[i].end)
    }
    tree_chroms.push(chroms[k])
    trees.push(interval_tree_new(starts, ends))
  }
  let out : Array[BedInterval] = Array::new()
  for iv in a {
    let mut hits = false
    for k = 0; k < tree_chroms.length(); k = k + 1 {
      if tree_chroms[k] == iv.chrom {
        let q = trees[k].query(iv.start, iv.end)
        if q.length() > 0 {
          hits = true
        }
        break
      }
    }
    if hits {
      out.push(iv)
    }
  }
  out
}
///|
pub struct ClosestResult {
  a : BedInterval
  b : BedInterval
  distance : Int
}

///|
fn bed_distance(a : BedInterval, b : BedInterval) -> Int {
  // Signed distance matching bedtools closest: overlap 0, downstream positive,
  // upstream negative, with a +1 offset (adjacent features => distance 1).
  if a.start < b.end && b.start < a.end {
    0
  } else if b.start >= a.end {
    b.start - a.end + 1
  } else {
    -(a.start - b.end) - 1
  }
}

///|
/// For each A, the closest B on the same chromosome (minimum genomic distance,
/// 0 when overlapping). Ties are broken by the smaller B start (matching
/// `bedtools closest -t first`).
pub fn bedtools_closest(
  a : Array[BedInterval],
  b : Array[BedInterval],
) -> Array[ClosestResult] {
  let out : Array[ClosestResult] = Array::new()
  for ai in a {
    let overlaps : Array[BedInterval] = Array::new()
    let up_recs : Array[BedInterval] = Array::new()
    let down_recs : Array[BedInterval] = Array::new()
    let mut up_dist = 0
    let mut down_dist = 0
    let mut has_up = false
    let mut has_down = false
    for bi in b {
      if ai.chrom != bi.chrom {
        continue
      }
      if ai.start < bi.end && bi.start < ai.end {
        overlaps.push(bi)
      } else if bi.start >= ai.end {
        let d = bi.start - ai.end + 1
        if not(has_down) || d < down_dist {
          down_dist = d
          down_recs.clear()
          down_recs.push(bi)
          has_down = true
        } else if d == down_dist {
          down_recs.push(bi)
        }
      } else if ai.start >= bi.end {
        let d = ai.start - bi.end + 1
        if not(has_up) || d < up_dist {
          up_dist = d
          up_recs.clear()
          up_recs.push(bi)
          has_up = true
        } else if d == up_dist {
          up_recs.push(bi)
        }
      }
    }
    if overlaps.length() > 0 {
      for bb in overlaps {
        out.push(ClosestResult::{ a: ai, b: bb, distance: 0 })
      }
    } else {
      // ALL_TIES: report the closer side(s); a tie reports both
      let use_up = has_up && (not(has_down) || up_dist <= down_dist)
      let use_down = has_down && (not(has_up) || down_dist <= up_dist)
      if use_up {
        for bb in up_recs {
          out.push(ClosestResult::{ a: ai, b: bb, distance: 0 - up_dist })
        }
      }
      if use_down {
        for bb in down_recs {
          out.push(ClosestResult::{ a: ai, b: bb, distance: down_dist })
        }
      }
    }
  }
  out
}

///|
/// Report A features that have a B feature within `w` bp of A's start (on both
/// sides), matching `bedtools window -w w`.
pub fn bedtools_window(
  a : Array[BedInterval],
  b : Array[BedInterval],
  w : Int,
) -> Array[BedInterval] {
  let out : Array[BedInterval] = Array::new()
  for ai in a {
    // window extends w bp upstream of A.start and w bp downstream of A.end
    let lo = ai.start - w
    let hi = ai.end + w
    let mut hit = false
    for bi in b {
      if ai.chrom == bi.chrom && bi.start < hi && bi.end > lo {
        hit = true
        break
      }
    }
    if hit {
      out.push(ai)
    }
  }
  out
}

///|
pub struct CoverageResult {
  a : BedInterval
  count : Int
  bases_covered : Int
  length : Int
  fraction : Double
}

///|
/// For each A, count the overlapping B features and report bases covered,
/// length, and covered fraction (matching `bedtools coverage` default).
pub fn bedtools_coverage(
  a : Array[BedInterval],
  b : Array[BedInterval],
) -> Array[CoverageResult] {
  let out : Array[CoverageResult] = Array::new()
  for ai in a {
    let len = ai.end - ai.start
    let depth : Array[Int] = Array::make(len, 0)
    let mut count = 0
    for bi in b {
      if ai.chrom == bi.chrom && ai.start < bi.end && bi.start < ai.end {
        let ms = if ai.start > bi.start { ai.start } else { bi.start }
        let me = if ai.end < bi.end { ai.end } else { bi.end }
        let mut p = ms
        while p < me {
          depth[p - ai.start] = depth[p - ai.start] + 1
          p = p + 1
        }
        count = count + 1
      }
    }
    let mut covered = 0
    for i = 0; i < len; i = i + 1 {
      if depth[i] >= 1 {
        covered = covered + 1
      }
    }
    let fraction = covered.to_double() / len.to_double()
    out.push(CoverageResult::{ a: ai, count, bases_covered: covered, length: len, fraction })
  }
  out
}

///|
pub struct GenomecovSegment {
  chrom : String
  start : Int
  end : Int
  depth : Int
}

///|
fn sort_events(events : Array[(Int, Int)]) -> Array[(Int, Int)] {
  let out = events.copy()
  let mut i = 1
  while i < out.length() {
    let key = out[i]
    let mut j = i - 1
    while j >= 0 && out[j].0 > key.0 {
      out[j + 1] = out[j]
      j = j - 1
    }
    out[j + 1] = key
    i = i + 1
  }
  out
}

///|
/// Sweep one chromosome's intervals into bedgraph segments (depth >= 1),
/// matching `bedtools genomecov -bg`.
fn genomecov_bg_one_chrom(ivals : Array[(Int, Int)]) -> Array[(Int, Int, Int)] {
  let events : Array[(Int, Int)] = Array::new()
  for (s, e) in ivals {
    events.push((s, 1))
    events.push((e, -1))
  }
  let sorted = sort_events(events)
  let out : Array[(Int, Int, Int)] = Array::new()
  let mut cur_depth = 0
  let mut i = 0
  while i < sorted.length() {
    let pos = sorted[i].0
    let mut delta = 0
    while i < sorted.length() && sorted[i].0 == pos {
      delta = delta + sorted[i].1
      i = i + 1
    }
    cur_depth = cur_depth + delta
    if i < sorted.length() && cur_depth >= 1 {
      out.push((pos, sorted[i].0, cur_depth))
    }
  }
  out
}

///|
/// Report contiguous depth segments (depth >= 1) via differential-event sweep.
pub fn bedtools_genomecov_bg(intervals : Array[BedInterval]) -> Array[GenomecovSegment] {
  let chroms : Array[String] = Array::new()
  let per_chrom : Array[Array[BedInterval]] = Array::new()
  for iv in intervals {
    let mut found = -1
    for k = 0; k < chroms.length(); k = k + 1 {
      if chroms[k] == iv.chrom {
        found = k
        break
      }
    }
    if found < 0 {
      chroms.push(iv.chrom)
      per_chrom.push(Array::new())
      per_chrom[per_chrom.length() - 1].push(iv)
    } else {
      per_chrom[found].push(iv)
    }
  }
  let out : Array[GenomecovSegment] = Array::new()
  for k = 0; k < chroms.length(); k = k + 1 {
    let ivs : Array[(Int, Int)] = Array::new()
    for iv in per_chrom[k] {
      ivs.push((iv.start, iv.end))
    }
    for (s, e, d) in genomecov_bg_one_chrom(ivs) {
      out.push(GenomecovSegment::{ chrom: chroms[k], start: s, end: e, depth: d })
    }
  }
  out
}

///|
fn clamp_int(v : Int, lo : Int, hi : Int) -> Int {
  if v < lo {
    lo
  } else if v > hi {
    hi
  } else {
    v
  }
}

///|
fn chrom_size_lookup(chrom_sizes : Array[(String, Int)], chrom : String) -> Int {
  for (c, s) in chrom_sizes {
    if c == chrom {
      return s
    }
  }
  -1
}

///|
/// Remove the portions of each A that overlap any B (matching `bedtools subtract`).
pub fn bedtools_subtract(
  a : Array[BedInterval],
  b : Array[BedInterval],
) -> Array[BedInterval] {
  let out : Array[BedInterval] = Array::new()
  for ai in a {
    let len = ai.end - ai.start
    let uncovered = Array::make(len, true)
    for bi in b {
      if ai.chrom == bi.chrom && ai.start < bi.end && bi.start < ai.end {
        let ms = if ai.start > bi.start { ai.start } else { bi.start }
        let me = if ai.end < bi.end { ai.end } else { bi.end }
        let mut p = ms
        while p < me {
          uncovered[p - ai.start] = false
          p = p + 1
        }
      }
    }
    let mut i = 0
    while i < len {
      if uncovered[i] {
        let mut j = i
        while j < len && uncovered[j] {
          j = j + 1
        }
        out.push(BedInterval::{ chrom: ai.chrom, start: ai.start + i, end: ai.start + j })
        i = j
      } else {
        i = i + 1
      }
    }
  }
  out
}

///|
/// Report the genome gaps not covered by A (chrom lengths from `chrom_sizes`),
/// matching `bedtools complement`.
pub fn bedtools_complement(
  a : Array[BedInterval],
  chrom_sizes : Array[(String, Int)],
) -> Array[BedInterval] {
  let merged = bedtools_merge(a)
  let out : Array[BedInterval] = Array::new()
  for (chrom, size) in chrom_sizes {
    let mut cur = 0
    for iv in merged {
      if iv.chrom == chrom {
        if iv.start > cur {
          out.push(BedInterval::{ chrom, start: cur, end: iv.start })
        }
        if iv.end > cur {
          cur = iv.end
        }
      }
    }
    if cur < size {
      out.push(BedInterval::{ chrom, start: cur, end: size })
    }
  }
  out
}

///|
/// Flanking regions of each A: `[start-L, start)` and `[end, end+R)`, clamped.
pub fn bedtools_flank(
  a : Array[BedInterval],
  chrom_sizes : Array[(String, Int)],
  left : Int,
  right : Int,
) -> Array[BedInterval] {
  let out : Array[BedInterval] = Array::new()
  for ai in a {
    let size = chrom_size_lookup(chrom_sizes, ai.chrom)
    if size < 0 {
      continue
    }
    if left > 0 && ai.start != 0 {
      let s = if ai.start - left > 0 { ai.start - left } else { 0 }
      out.push(BedInterval::{ chrom: ai.chrom, start: s, end: ai.start })
    }
    if right > 0 && ai.end < size {
      let e = if ai.end + right < size { ai.end + right } else { size }
      out.push(BedInterval::{ chrom: ai.chrom, start: ai.end, end: e })
    }
  }
  out
}

///|
/// Shift each A by signed `s` bp, clamped to the chromosome bounds.
pub fn bedtools_shift(
  a : Array[BedInterval],
  chrom_sizes : Array[(String, Int)],
  s : Int,
) -> Array[BedInterval] {
  let out : Array[BedInterval] = Array::new()
  for ai in a {
    let size = chrom_size_lookup(chrom_sizes, ai.chrom)
    if size < 0 {
      continue
    }
    let ns = clamp_int(ai.start + s, 0, size - 1)
    let ne = clamp_int(ai.end + s, 1, size)
    out.push(BedInterval::{ chrom: ai.chrom, start: ns, end: ne })
  }
  out
}

///|
/// Extend each A by `left` bp upstream and `right` bp downstream, clamped.
pub fn bedtools_slop(
  a : Array[BedInterval],
  chrom_sizes : Array[(String, Int)],
  left : Int,
  right : Int,
) -> Array[BedInterval] {
  let out : Array[BedInterval] = Array::new()
  for ai in a {
    let size = chrom_size_lookup(chrom_sizes, ai.chrom)
    if size < 0 {
      continue
    }
    let ns = clamp_int(ai.start - left, 0, size)
    let ne = clamp_int(ai.end + right, 0, size)
    out.push(BedInterval::{ chrom: ai.chrom, start: ns, end: ne })
  }
  out
}

///|
fn collect_overlap_vals(
  ai : BedInterval,
  b : Array[BedInterval],
  b_values : Array[Double],
) -> Array[Double] {
  let vals : Array[Double] = Array::new()
  for k = 0; k < b.length(); k = k + 1 {
    if ai.chrom == b[k].chrom && ai.start < b[k].end && b[k].start < ai.end {
      vals.push(b_values[k])
    }
  }
  vals
}

///|
/// Number of B features overlapping each A (matching `bedtools map -c 5 -o count`).
pub fn bedtools_map_count(a : Array[BedInterval], b : Array[BedInterval]) -> Array[Int] {
  let out : Array[Int] = Array::new()
  for ai in a {
    let mut c = 0
    for bi in b {
      if ai.chrom == bi.chrom && ai.start < bi.end && bi.start < ai.end {
        c = c + 1
      }
    }
    out.push(c)
  }
  out
}

///|
pub fn bedtools_map_sum(
  a : Array[BedInterval],
  b : Array[BedInterval],
  b_values : Array[Double],
) -> Array[Double] {
  let out : Array[Double] = Array::new()
  for ai in a {
    let vals = collect_overlap_vals(ai, b, b_values)
    let mut s = 0.0
    for v in vals {
      s = s + v
    }
    out.push(s)
  }
  out
}

///|
pub fn bedtools_map_mean(
  a : Array[BedInterval],
  b : Array[BedInterval],
  b_values : Array[Double],
) -> Array[Double] {
  let out : Array[Double] = Array::new()
  for ai in a {
    let vals = collect_overlap_vals(ai, b, b_values)
    if vals.length() == 0 {
      out.push(0.0)
    } else {
      let mut s = 0.0
      for v in vals {
        s = s + v
      }
      out.push(s / vals.length().to_double())
    }
  }
  out
}

///|
pub fn bedtools_map_min(
  a : Array[BedInterval],
  b : Array[BedInterval],
  b_values : Array[Double],
) -> Array[Double] {
  let out : Array[Double] = Array::new()
  for ai in a {
    let vals = collect_overlap_vals(ai, b, b_values)
    if vals.length() == 0 {
      out.push(0.0)
    } else {
      let mut m = vals[0]
      for v in vals {
        if v < m {
          m = v
        }
      }
      out.push(m)
    }
  }
  out
}

///|
pub fn bedtools_map_max(
  a : Array[BedInterval],
  b : Array[BedInterval],
  b_values : Array[Double],
) -> Array[Double] {
  let out : Array[Double] = Array::new()
  for ai in a {
    let vals = collect_overlap_vals(ai, b, b_values)
    if vals.length() == 0 {
      out.push(0.0)
    } else {
      let mut m = vals[0]
      for v in vals {
        if v > m {
          m = v
        }
      }
      out.push(m)
    }
  }
  out
}


///|
/// Report A if a B overlaps at least `frac` of A's length (matching `bedtools
/// intersect -f frac`).
pub fn bedtools_intersect_frac(
  a : Array[BedInterval],
  b : Array[BedInterval],
  frac : Double,
) -> Array[BedInterval] {
  let out : Array[BedInterval] = Array::new()
  for ai in a {
    let mut hit = false
    for bi in b {
      if ai.chrom == bi.chrom {
        let ov = (if ai.end < bi.end { ai.end } else { bi.end }) -
          (if ai.start > bi.start { ai.start } else { bi.start })
        if ov > 0 && ov.to_double() / (ai.end - ai.start).to_double() >= frac {
          hit = true
          break
        }
      }
    }
    if hit {
      out.push(ai)
    }
  }
  out
}

///|
/// Report A that overlaps a same-strand B (matching `bedtools intersect -s`).
/// Strand values are "+", "-", or "." (unknown).
pub fn bedtools_intersect_stranded(
  a : Array[BedInterval],
  a_strand : Array[String],
  b : Array[BedInterval],
  b_strand : Array[String],
) -> Array[BedInterval] {
  let out : Array[BedInterval] = Array::new()
  for i = 0; i < a.length(); i = i + 1 {
    let ai = a[i]
    let mut hit = false
    for k = 0; k < b.length(); k = k + 1 {
      let bi = b[k]
      if ai.chrom == bi.chrom &&
        ai.start < bi.end &&
        bi.start < ai.end &&
        a_strand[i] == b_strand[k] &&
        a_strand[i] != "." {
        hit = true
        break
      }
    }
    if hit {
      out.push(ai)
    }
  }
  out
}
