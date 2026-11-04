/// bioframe_set_ops.mbt — DataFrame-style complement, subtract, sort, sanitize
/// Items 2, 3, 5, 6: bf_complement, bf_subtract, bf_sort, bf_sanitize

///|
/// Compute complement of intervals within chromsizes.
/// Returns intervals covering the regions NOT in the input.
pub fn bf_complement(
  intervals : Array[BedInterval],
  chromsizes : Array[BfChromSize],
) -> Array[BedInterval] {
  // Group intervals by chrom
  let results : Array[BedInterval] = Array(capacity=0)
  for ci = 0; ci < chromsizes.length(); ci = ci + 1 {
    let cs = chromsizes[ci]
    // Collect intervals on this chrom
    let on_chrom : Array[BedInterval] = Array(capacity=0)
    for i = 0; i < intervals.length(); i = i + 1 {
      if intervals[i].chrom == cs.name {
        on_chrom.push(intervals[i])
      }
    }
    // Sort by start
    on_chrom.sort_by(fn(a, b) {
      if a.start < b.start {
        -1
      } else if a.start > b.start {
        1
      } else {
        0
      }
    })
    // Walk and emit gaps
    let mut pos : Int32 = 0
    for i = 0; i < on_chrom.length(); i = i + 1 {
      let iv = on_chrom[i]
      // Clamp to chromsize
      let s = if iv.start < 0 { 0 } else { iv.start }
      let e = if iv.end > cs.length { cs.length } else { iv.end }
      if s > pos {
        results.push(BedInterval::new(cs.name, pos, s))
      }
      if e > pos {
        pos = e
      }
    }
    // Final gap to end of chrom
    if pos < cs.length {
      results.push(BedInterval::new(cs.name, pos, cs.length))
    }
  }
  results
}

///|
/// Subtract intervals B from intervals A.
/// Returns regions in A that do not overlap any interval in B.
pub fn bf_subtract(
  a : Array[BedInterval],
  b : Array[BedInterval],
) -> Array[BedInterval] {
  let results : Array[BedInterval] = Array(capacity=0)
  // Group B by chrom for efficiency
  let b_chroms : Array[String] = Array(capacity=0)
  let b_by_chrom : Array[Array[(Int32, Int32)]] = Array(capacity=0)
  for i = 0; i < b.length(); i = i + 1 {
    let bv = b[i]
    let mut found = false
    for j = 0; j < b_chroms.length(); j = j + 1 {
      if b_chroms[j] == bv.chrom {
        b_by_chrom[j].push((bv.start, bv.end))
        found = true
      }
    }
    if !found {
      b_chroms.push(bv.chrom)
      let arr : Array[(Int32, Int32)] = Array(capacity=0)
      arr.push((bv.start, bv.end))
      b_by_chrom.push(arr)
    }
  }
  // Sort each chrom's B intervals
  for j = 0; j < b_by_chrom.length(); j = j + 1 {
    b_by_chrom[j].sort_by(fn(a, b) {
      if a.0 < b.0 {
        -1
      } else if a.0 > b.0 {
        1
      } else {
        0
      }
    })
  }
  // For each A interval, subtract overlapping B intervals
  for i = 0; i < a.length(); i = i + 1 {
    let av = a[i]
    // Find B intervals on same chrom
    let mut b_idx = -1
    for j = 0; j < b_chroms.length(); j = j + 1 {
      if b_chroms[j] == av.chrom {
        b_idx = j
      }
    }
    if b_idx < 0 {
      // No B on this chrom — keep entire A
      results.push(BedInterval::new(av.chrom, av.start, av.end))
      continue
    }
    // Walk through A interval, cutting out B overlaps
    let b_list = b_by_chrom[b_idx]
    let segments : Array[(Int32, Int32)] = Array(capacity=0)
    segments.push((av.start, av.end))
    for j = 0; j < b_list.length(); j = j + 1 {
      let (b_start, b_end) = b_list[j]
      let new_segs : Array[(Int32, Int32)] = Array(capacity=0)
      for k = 0; k < segments.length(); k = k + 1 {
        let (s, e) = segments[k]
        if b_end <= s || b_start >= e {
          // No overlap
          new_segs.push((s, e))
        } else {
          // Overlap — split
          if b_start > s {
            new_segs.push((s, b_start))
          }
          if b_end < e {
            new_segs.push((b_end, e))
          }
        }
      }
      // Replace segments
      segments.clear()
      for k = 0; k < new_segs.length(); k = k + 1 {
        segments.push(new_segs[k])
      }
    }
    for k = 0; k < segments.length(); k = k + 1 {
      let (s, e) = segments[k]
      if e > s {
        results.push(BedInterval::new(av.chrom, s, e))
      }
    }
  }
  results
}

///|
/// Sort intervals by (chrom, start, end) in genomic order.
pub fn bf_sort(intervals : Array[BedInterval]) -> Array[BedInterval] {
  let arr = Array(capacity=intervals.length())
  for i = 0; i < intervals.length(); i = i + 1 {
    arr.push(intervals[i])
  }
  arr.sort_by(fn(a, b) {
    if a.chrom < b.chrom {
      -1
    } else if a.chrom > b.chrom {
      1
    } else if a.start < b.start {
      -1
    } else if a.start > b.start {
      1
    } else if a.end < b.end {
      -1
    } else if a.end > b.end {
      1
    } else {
      0
    }
  })
  arr
}

///|
/// Sanitize intervals: sort, merge overlapping, clamp to chromsizes.
/// Returns a clean, non-overlapping, sorted set of intervals.
pub fn bf_sanitize(
  intervals : Array[BedInterval],
  chromsizes : Array[BfChromSize],
) -> Array[BedInterval] {
  // Build chromsize lookup
  let cs_map_chroms : Array[String] = Array(capacity=0)
  let cs_map_lengths : Array[Int32] = Array(capacity=0)
  for i = 0; i < chromsizes.length(); i = i + 1 {
    cs_map_chroms.push(chromsizes[i].name)
    cs_map_lengths.push(chromsizes[i].length)
  }
  // Clamp each interval to chromsize
  let clamped : Array[BedInterval] = Array(capacity=0)
  for i = 0; i < intervals.length(); i = i + 1 {
    let iv = intervals[i]
    let mut max_len : Int32 = iv.end
    for j = 0; j < cs_map_chroms.length(); j = j + 1 {
      if cs_map_chroms[j] == iv.chrom {
        max_len = cs_map_lengths[j]
      }
    }
    let s = if iv.start < 0 { 0 } else { iv.start }
    let e = if iv.end > max_len { max_len } else { iv.end }
    if e > s {
      clamped.push(BedInterval::new(iv.chrom, s, e))
    }
  }
  // Sort
  let sorted = bf_sort(clamped)
  // Merge overlapping (gap=0)
  let merged = bf_merge(sorted)
  // Convert back to BedInterval
  let results : Array[BedInterval] = Array(capacity=0)
  for i = 0; i < merged.length(); i = i + 1 {
    let m = merged[i]
    results.push(BedInterval::new(m.chrom, m.start, m.end))
  }
  results
}
