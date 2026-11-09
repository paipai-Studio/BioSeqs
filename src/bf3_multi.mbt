// Bioframe third-round: multi-set operations (items 30-32)
// multi_intersect, multi_union, multi_coverage

///|
pub struct Bf3MultiIntersectResult {
  chrom : String
  start : Int32
  end_ : Int32
  count : Int32 // how many sets contain this interval
  set_ids : Array[Int32] // which sets contain this interval
}

///|
pub fn bf3_multi_intersect(
  sets : Array[Array[(String, Int32, Int32)]],
  min_count : Int32,
) -> Array[Bf3MultiIntersectResult] {
  // Find regions present in at least min_count sets
  // Collect all boundary points
  let all_events : Array[(String, Int32, Int32, Int32)] = [] // (chrom, pos, set_id, delta)
  for s in 0..<sets.length() {
    let set_intervals = sets[s]
    for i in 0..<set_intervals.length() {
      let (chrom, start, end_) = set_intervals[i]
      all_events.push((chrom, start, s, 1))
      all_events.push((chrom, end_, s, -1))
    }
  }
  // Sort by chrom, then position
  all_events.sort_by(fn(a, b) {
    let (ca, pa, _, _) = a
    let (cb, pb, _, _) = b
    if ca < cb {
      -1
    } else if ca > cb {
      1
    } else if pa < pb {
      -1
    } else if pa > pb {
      1
    } else {
      0
    }
  })
  // Sweep line
  let result : Array[Bf3MultiIntersectResult] = []
  let active_sets : Array[Int32] = [] // set_ids currently active
  let mut prev_pos = -1
  let mut prev_chrom = ""
  for e in 0..<all_events.length() {
    let (chrom, pos, set_id, delta) = all_events[e]
    // Emit region from prev_pos to pos if count >= min_count
    if prev_pos >= 0 && prev_pos < pos && prev_chrom == chrom {
      let count = active_sets.length()
      if count >= min_count {
        result.push({
          chrom: prev_chrom,
          start: prev_pos,
          end_: pos,
          count,
          set_ids: active_sets.copy(),
        })
      }
    }
    // Update active sets
    if delta > 0 {
      active_sets.push(set_id)
    } else {
      // Remove set_id
      let mut idx = -1
      for k in 0..<active_sets.length() {
        if active_sets[k] == set_id {
          idx = k
          break
        }
      }
      if idx >= 0 {
        let _ = active_sets.remove(idx)
      }
    }
    prev_pos = pos
    prev_chrom = chrom
  }
  result
}

///|
pub struct Bf3MultiUnionResult {
  chrom : String
  start : Int32
  end_ : Int32
  source_sets : Array[Int32]
}

///|
pub fn bf3_multi_union(
  sets : Array[Array[(String, Int32, Int32)]],
) -> Array[Bf3MultiUnionResult] {
  // Union of all sets, with source tracking
  let all_events : Array[(String, Int32, Int32, Int32)] = [] // (chrom, pos, set_id, delta)
  for s in 0..<sets.length() {
    let set_intervals = sets[s]
    for i in 0..<set_intervals.length() {
      let (chrom, start, end_) = set_intervals[i]
      all_events.push((chrom, start, s, 1))
      all_events.push((chrom, end_, s, -1))
    }
  }
  all_events.sort_by(fn(a, b) {
    let (ca, pa, _, _) = a
    let (cb, pb, _, _) = b
    if ca < cb {
      -1
    } else if ca > cb {
      1
    } else if pa < pb {
      -1
    } else if pa > pb {
      1
    } else {
      0
    }
  })
  let result : Array[Bf3MultiUnionResult] = []
  let active : Array[Int32] = []
  let mut prev_pos = -1
  let mut prev_chrom = ""
  for e in 0..<all_events.length() {
    let (chrom, pos, set_id, delta) = all_events[e]
    if prev_pos >= 0 &&
      prev_pos < pos &&
      prev_chrom == chrom &&
      not(active.is_empty()) {
      result.push({
        chrom: prev_chrom,
        start: prev_pos,
        end_: pos,
        source_sets: active.copy(),
      })
    }
    if delta > 0 {
      active.push(set_id)
    } else {
      let mut idx = -1
      for k in 0..<active.length() {
        if active[k] == set_id {
          idx = k
          break
        }
      }
      if idx >= 0 {
        let _ = active.remove(idx)
      }
    }
    prev_pos = pos
    prev_chrom = chrom
  }
  result
}

///|
pub struct Bf3MultiCoverageResult {
  chrom : String
  start : Int32
  end_ : Int32
  coverage : Array[Int32] // coverage[i] = number of intervals from set i covering this region
}

///|
pub fn bf3_multi_coverage(
  sets : Array[Array[(String, Int32, Int32)]],
) -> Array[Bf3MultiCoverageResult] {
  // Layered coverage: for each region, how many intervals from each set cover it
  let num_sets = sets.length()
  // Collect all boundary points
  let all_points : Array[(String, Int32)] = []
  for s in 0..<num_sets {
    let set_intervals = sets[s]
    for i in 0..<set_intervals.length() {
      let (chrom, start, end_) = set_intervals[i]
      all_points.push((chrom, start))
      all_points.push((chrom, end_))
    }
  }
  // Sort and deduplicate
  all_points.sort_by(fn(a, b) {
    let (ca, pa) = a
    let (cb, pb) = b
    if ca < cb {
      -1
    } else if ca > cb {
      1
    } else if pa < pb {
      -1
    } else if pa > pb {
      1
    } else {
      0
    }
  })
  let result : Array[Bf3MultiCoverageResult] = []
  let num_pts = all_points.length() - 1
  for i in 0..<num_pts {
    let (chrom, start) = all_points[i]
    let (next_chrom, next_pos) = all_points[i + 1]
    if chrom != next_chrom || start == next_pos {
      continue
    }
    let end_ = next_pos
    // Count coverage from each set
    let coverage : Array[Int32] = []
    for s in 0..<num_sets {
      let mut cnt = 0
      let set_intervals = sets[s]
      for j in 0..<set_intervals.length() {
        let (c, s_pos, e_pos) = set_intervals[j]
        if c == chrom && s_pos <= start && e_pos >= end_ {
          cnt = cnt + 1
        }
      }
      coverage.push(cnt)
    }
    result.push({ chrom, start, end_, coverage, })
  }
  result
}
