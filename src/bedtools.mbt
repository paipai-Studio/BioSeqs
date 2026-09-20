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