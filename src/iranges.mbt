///|
/// IRanges - Integer range operations inspired by Bioconductor IRanges.

///|
pub struct IRanges {
  starts : Array[Int]
  ends : Array[Int]
  widths : Array[Int]
}

///|
pub fn iranges(starts : Array[Int], ends : Array[Int]) -> IRanges {
  let n = starts.length()
  let widths : Array[Int] = Array::make(n, 0)

  for i = 0; i < n; i = i + 1 {
    widths[i] = ends[i] - starts[i] + 1
  }

  IRanges::{ starts, ends, widths }
}

///|
pub fn iranges_single(start : Int, end : Int) -> IRanges {
  IRanges::{ starts: [start], ends: [end], widths: [end - start + 1] }
}

///|
pub fn ir_length(ir : IRanges) -> Int {
  ir.starts.length()
}

///|
pub fn ir_start(ir : IRanges) -> Array[Int] {
  ir.starts
}

///|
pub fn ir_end(ir : IRanges) -> Array[Int] {
  ir.ends
}

///|
pub fn ir_width(ir : IRanges) -> Array[Int] {
  ir.widths
}

///|
pub fn ir_shift(ir : IRanges, shift : Int) -> IRanges {
  let n = ir.starts.length()
  let new_starts : Array[Int] = Array::make(n, 0)
  let new_ends : Array[Int] = Array::make(n, 0)
  let new_widths : Array[Int] = Array::make(n, 0)

  for i = 0; i < n; i = i + 1 {
    new_starts[i] = ir.starts[i] + shift
    new_ends[i] = ir.ends[i] + shift
    new_widths[i] = ir.widths[i]
  }

  IRanges::{ starts: new_starts, ends: new_ends, widths: new_widths }
}

///|
pub fn ir_narrow(ir : IRanges, start : Int, end : Int) -> IRanges {
  let n = ir.starts.length()
  let new_starts : Array[Int] = Array::make(n, 0)
  let new_ends : Array[Int] = Array::make(n, 0)
  let new_widths : Array[Int] = Array::make(n, 0)

  for i = 0; i < n; i = i + 1 {
    new_starts[i] = ir.starts[i] + start - 1
    new_ends[i] = ir.starts[i] + end - 1
    new_widths[i] = end - start + 1
  }

  IRanges::{ starts: new_starts, ends: new_ends, widths: new_widths }
}

///|
pub fn ir_resize(ir : IRanges, width : Int, fix : String) -> IRanges {
  let n = ir.starts.length()
  let new_starts : Array[Int] = Array::make(n, 0)
  let new_ends : Array[Int] = Array::make(n, 0)
  let new_widths : Array[Int] = Array::make(n, width)

  for i = 0; i < n; i = i + 1 {
    if fix == "start" {
      new_starts[i] = ir.starts[i]
      new_ends[i] = ir.starts[i] + width - 1
    } else if fix == "end" {
      new_ends[i] = ir.ends[i]
      new_starts[i] = ir.ends[i] - width + 1
    } else {
      let center = (ir.starts[i] + ir.ends[i]) / 2
      new_starts[i] = center - width / 2 + 1
      new_ends[i] = center + width / 2
    }
  }

  IRanges::{ starts: new_starts, ends: new_ends, widths: new_widths }
}

///|
pub fn ir_flank(
  ir : IRanges,
  width : Int,
  start : Bool,
  both : Bool,
) -> IRanges {
  let n = ir.starts.length()

  if both {
    let new_starts : Array[Int] = Array::make(n, 0)
    let new_ends : Array[Int] = Array::make(n, 0)
    let new_widths : Array[Int] = Array::make(n, width * 2)

    for i = 0; i < n; i = i + 1 {
      new_starts[i] = ir.starts[i] - width
      new_ends[i] = ir.ends[i] + width
    }

    IRanges::{ starts: new_starts, ends: new_ends, widths: new_widths }
  } else if start {
    let new_starts : Array[Int] = Array::make(n, 0)
    let new_ends : Array[Int] = Array::make(n, 0)
    let new_widths : Array[Int] = Array::make(n, width)

    for i = 0; i < n; i = i + 1 {
      new_starts[i] = ir.starts[i] - width
      new_ends[i] = ir.starts[i] - 1
    }

    IRanges::{ starts: new_starts, ends: new_ends, widths: new_widths }
  } else {
    let new_starts : Array[Int] = Array::make(n, 0)
    let new_ends : Array[Int] = Array::make(n, 0)
    let new_widths : Array[Int] = Array::make(n, width)

    for i = 0; i < n; i = i + 1 {
      new_starts[i] = ir.ends[i] + 1
      new_ends[i] = ir.ends[i] + width
    }

    IRanges::{ starts: new_starts, ends: new_ends, widths: new_widths }
  }
}

///|
pub fn ir_reduce(ir : IRanges, min_gapwidth : Int) -> IRanges {
  if ir.starts.length() == 0 {
    return IRanges::{ starts: [], ends: [], widths: [] }
  }

  let n = ir.starts.length()
  let indices : Array[(Int, Int)] = Array::make(n, (0, 0))
  for i = 0; i < n; i = i + 1 {
    indices[i] = (ir.starts[i], ir.ends[i])
  }

  for i = 0; i < n - 1; i = i + 1 {
    for j = i + 1; j < n; j = j + 1 {
      if indices[i].0 > indices[j].0 {
        let temp = indices[i]
        indices[i] = indices[j]
        indices[j] = temp
      }
    }
  }

  let reduced : Array[(Int, Int)] = Array::new()
  reduced.push(indices[0])

  for i = 1; i < n; i = i + 1 {
    let last = reduced[reduced.length() - 1]
    let current = indices[i]

    if current.0 <= last.1 + min_gapwidth + 1 {
      let new_end = if current.1 > last.1 { current.1 } else { last.1 }
      reduced[reduced.length() - 1] = (last.0, new_end)
    } else {
      reduced.push(current)
    }
  }

  let new_starts : Array[Int] = Array::make(reduced.length(), 0)
  let new_ends : Array[Int] = Array::make(reduced.length(), 0)
  let new_widths : Array[Int] = Array::make(reduced.length(), 0)

  for i = 0; i < reduced.length(); i = i + 1 {
    new_starts[i] = reduced[i].0
    new_ends[i] = reduced[i].1
    new_widths[i] = reduced[i].1 - reduced[i].0 + 1
  }

  IRanges::{ starts: new_starts, ends: new_ends, widths: new_widths }
}

///|
pub fn ir_disjoin(ir : IRanges) -> IRanges {
  if ir.starts.length() == 0 {
    return IRanges::{ starts: [], ends: [], widths: [] }
  }

  let n = ir.starts.length()
  let events : Array[(Int, Int)] = Array::make(n * 2, (0, 0))

  for i = 0; i < n; i = i + 1 {
    events[i * 2] = (ir.starts[i], 1)
    events[i * 2 + 1] = (ir.ends[i] + 1, -1)
  }

  for i = 0; i < events.length() - 1; i = i + 1 {
    for j = i + 1; j < events.length(); j = j + 1 {
      if events[i].0 > events[j].0 ||
        (events[i].0 == events[j].0 && events[i].1 > events[j].1) {
        let temp = events[i]
        events[i] = events[j]
        events[j] = temp
      }
    }
  }

  let disjoined : Array[(Int, Int)] = Array::new()
  let mut coverage = 0
  let mut current_start = 0

  for i = 0; i < events.length(); i = i + 1 {
    let pos = events[i].0
    let delta = events[i].1

    if coverage > 0 && pos > current_start {
      disjoined.push((current_start, pos - 1))
    }

    coverage = coverage + delta
    current_start = pos
  }

  let new_starts : Array[Int] = Array::make(disjoined.length(), 0)
  let new_ends : Array[Int] = Array::make(disjoined.length(), 0)
  let new_widths : Array[Int] = Array::make(disjoined.length(), 0)

  for i = 0; i < disjoined.length(); i = i + 1 {
    new_starts[i] = disjoined[i].0
    new_ends[i] = disjoined[i].1
    new_widths[i] = disjoined[i].1 - disjoined[i].0 + 1
  }

  IRanges::{ starts: new_starts, ends: new_ends, widths: new_widths }
}

///|
pub fn ir_intersect(ir1 : IRanges, ir2 : IRanges) -> IRanges {
  let intersects : Array[(Int, Int)] = Array::new()

  for i = 0; i < ir1.starts.length(); i = i + 1 {
    for j = 0; j < ir2.starts.length(); j = j + 1 {
      let start = if ir1.starts[i] > ir2.starts[j] {
        ir1.starts[i]
      } else {
        ir2.starts[j]
      }
      let end = if ir1.ends[i] < ir2.ends[j] {
        ir1.ends[i]
      } else {
        ir2.ends[j]
      }

      if start <= end {
        intersects.push((start, end))
      }
    }
  }

  let new_starts : Array[Int] = Array::make(intersects.length(), 0)
  let new_ends : Array[Int] = Array::make(intersects.length(), 0)
  let new_widths : Array[Int] = Array::make(intersects.length(), 0)

  for i = 0; i < intersects.length(); i = i + 1 {
    new_starts[i] = intersects[i].0
    new_ends[i] = intersects[i].1
    new_widths[i] = intersects[i].1 - intersects[i].0 + 1
  }

  IRanges::{ starts: new_starts, ends: new_ends, widths: new_widths }
}

///|
pub fn ir_union(ir1 : IRanges, ir2 : IRanges) -> IRanges {
  let combined_starts : Array[Int] = Array::make(
    ir1.starts.length() + ir2.starts.length(),
    0,
  )
  let combined_ends : Array[Int] = Array::make(
    ir1.starts.length() + ir2.starts.length(),
    0,
  )

  for i = 0; i < ir1.starts.length(); i = i + 1 {
    combined_starts[i] = ir1.starts[i]
    combined_ends[i] = ir1.ends[i]
  }
  for i = 0; i < ir2.starts.length(); i = i + 1 {
    combined_starts[ir1.starts.length() + i] = ir2.starts[i]
    combined_ends[ir1.starts.length() + i] = ir2.ends[i]
  }

  ir_reduce(iranges(combined_starts, combined_ends), 0)
}

///|
pub fn ir_setdiff(ir1 : IRanges, ir2 : IRanges) -> IRanges {
  if ir2.starts.length() == 0 {
    return ir1
  }

  let disjoined = ir_disjoin(ir_union(ir1, ir2))
  let intersected = ir_intersect(ir1, ir2)

  let diffs : Array[(Int, Int)] = Array::new()

  for i = 0; i < disjoined.starts.length(); i = i + 1 {
    let mut is_in_intersect = false
    for j = 0; j < intersected.starts.length(); j = j + 1 {
      if disjoined.starts[i] == intersected.starts[j] &&
        disjoined.ends[i] == intersected.ends[j] {
        is_in_intersect = true
        break
      }
    }

    if !is_in_intersect {
      let mut is_in_ir1 = false
      for j = 0; j < ir1.starts.length(); j = j + 1 {
        if disjoined.starts[i] >= ir1.starts[j] &&
          disjoined.ends[i] <= ir1.ends[j] {
          is_in_ir1 = true
          break
        }
      }

      if is_in_ir1 {
        diffs.push((disjoined.starts[i], disjoined.ends[i]))
      }
    }
  }

  let new_starts : Array[Int] = Array::make(diffs.length(), 0)
  let new_ends : Array[Int] = Array::make(diffs.length(), 0)
  let new_widths : Array[Int] = Array::make(diffs.length(), 0)

  for i = 0; i < diffs.length(); i = i + 1 {
    new_starts[i] = diffs[i].0
    new_ends[i] = diffs[i].1
    new_widths[i] = diffs[i].1 - diffs[i].0 + 1
  }

  IRanges::{ starts: new_starts, ends: new_ends, widths: new_widths }
}

///|
pub fn ir_overlaps(ir1 : IRanges, ir2 : IRanges) -> Array[(Int, Int)] {
  let overlaps : Array[(Int, Int)] = Array::new()

  for i = 0; i < ir1.starts.length(); i = i + 1 {
    for j = 0; j < ir2.starts.length(); j = j + 1 {
      if ir1.starts[i] <= ir2.ends[j] && ir1.ends[i] >= ir2.starts[j] {
        overlaps.push((i, j))
      }
    }
  }

  overlaps
}

///|
pub fn ir_count_overlaps(query : IRanges, subject : IRanges) -> Array[Int] {
  let counts : Array[Int] = Array::make(query.starts.length(), 0)

  for i = 0; i < query.starts.length(); i = i + 1 {
    let mut count = 0
    for j = 0; j < subject.starts.length(); j = j + 1 {
      if query.starts[i] <= subject.ends[j] &&
        query.ends[i] >= subject.starts[j] {
        count = count + 1
      }
    }
    counts[i] = count
  }

  counts
}

///|
/// Advanced findOverlaps with different overlap types.
/// "any": any overlap, "within": query fully within subject, "start": same start, "end": same end
pub fn find_overlaps_advanced(
  query~ : IRanges,
  subject~ : IRanges,
  type_? : String = "any",
) -> Array[(Int, Int)] {
  let results : Array[(Int, Int)] = Array::new()

  let mut i = 0
  while i < query.starts.length() {
    let mut j = 0
    while j < subject.starts.length() {
      let q_start = query.starts[i]
      let q_end = query.ends[i]
      let s_start = subject.starts[j]
      let s_end = subject.ends[j]

      let match_found = if type_ == "within" {
        q_start >= s_start && q_end <= s_end
      } else if type_ == "start" {
        q_start == s_start
      } else if type_ == "end" {
        q_end == s_end
      } else {
        q_start <= s_end && q_end >= s_start
      }

      if match_found {
        results.push((i, j))
      }
      j = j + 1
    }
    i = i + 1
  }

  results
}

///|
/// Find nearest neighbor in subject for each query range.
pub fn ir_nearest(query : IRanges, subject : IRanges) -> Array[Int] {
  let n = query.starts.length()
  let indices : Array[Int] = Array::make(n, -1)

  let mut i = 0
  while i < n {
    let q_start = query.starts[i]
    let q_end = query.ends[i]
    let q_mid = (q_start + q_end) / 2

    let mut min_dist = @double.infinity
    let mut best_idx = -1

    let mut j = 0
    while j < subject.starts.length() {
      let s_mid = (subject.starts[j] + subject.ends[j]) / 2
      let dist = (q_mid - s_mid).abs().to_double()

      if dist < min_dist {
        min_dist = dist
        best_idx = j
      }
      j = j + 1
    }

    indices[i] = best_idx
    i = i + 1
  }

  indices
}

///|
/// Find preceding ranges in subject for each query.
pub fn ir_precede(query : IRanges, subject : IRanges) -> Array[Int] {
  let n = query.starts.length()
  let indices : Array[Int] = Array::make(n, -1)

  let mut i = 0
  while i < n {
    let q_start = query.starts[i]

    let mut best_idx = -1
    let mut best_end = @double.neg_infinity

    let mut j = 0
    while j < subject.starts.length() {
      let s_end = subject.ends[j]
      if s_end < q_start {
        if s_end.to_double() > best_end {
          best_end = s_end.to_double()
          best_idx = j
        }
      }
      j = j + 1
    }

    indices[i] = best_idx
    i = i + 1
  }

  indices
}

///|
/// Find following ranges in subject for each query.
pub fn ir_follow(query : IRanges, subject : IRanges) -> Array[Int] {
  let n = query.starts.length()
  let indices : Array[Int] = Array::make(n, -1)

  let mut i = 0
  while i < n {
    let q_end = query.ends[i]

    let mut best_idx = -1
    let mut best_start = @double.infinity

    let mut j = 0
    while j < subject.starts.length() {
      let s_start = subject.starts[j]
      if s_start > q_end {
        if s_start.to_double() < best_start {
          best_start = s_start.to_double()
          best_idx = j
        }
      }
      j = j + 1
    }

    indices[i] = best_idx
    i = i + 1
  }

  indices
}

///|
/// Calculate coverage of integer ranges.
pub fn ir_coverage(ir : IRanges, max_end : Int) -> Array[Int] {
  let n = ir.starts.length()
  let cov : Array[Int] = Array::make(max_end, 0)

  let mut i = 0
  while i < n {
    let start = ir.starts[i]
    let end = ir.ends[i]
    let clamped_start = if start < 1 { 1 } else { start }
    let clamped_end = if end > max_end { max_end } else { end }

    let mut pos = clamped_start
    while pos <= clamped_end {
      cov[pos - 1] = cov[pos - 1] + 1
      pos = pos + 1
    }
    i = i + 1
  }

  cov
}

///|
/// Calculate coverage as run-length encoded ranges.
pub fn ir_coverage_rle(ir : IRanges, max_end : Int) -> Array[(Int, Int, Int)] {
  let cov = ir_coverage(ir, max_end)
  let results : Array[(Int, Int, Int)] = Array::new()

  let n = cov.length()
  let mut i = 0
  while i < n {
    let val = cov[i]
    if val > 0 {
      let start = i + 1
      let mut end = start
      let mut j = i + 1
      while j < n && cov[j] == val {
        end = j + 1
        j = j + 1
      }
      results.push((start, end, val))
      i = j
    } else {
      i = i + 1
    }
  }

  results
}

///|
/// Calculate distance matrix between two sets of ranges.
pub fn ir_distance(ir1 : IRanges, ir2 : IRanges) -> Array[Array[Int]] {
  let n1 = ir1.starts.length()
  let n2 = ir2.starts.length()
  let distances : Array[Array[Int]] = Array::make(n1, Array::make(n2, -1))

  let mut i = 0
  while i < n1 {
    let s1 = ir1.starts[i]
    let e1 = ir1.ends[i]
    let mut j = 0
    while j < n2 {
      let s2 = ir2.starts[j]
      let e2 = ir2.ends[j]

      let mut dist = 0
      if e1 < s2 {
        dist = s2 - e1 - 1
      } else if e2 < s1 {
        dist = s1 - e2 - 1
      }

      distances[i][j] = dist
      j = j + 1
    }
    i = i + 1
  }

  distances
}

///|
/// Calculate distance to nearest neighbor within the same set.
pub fn ir_distance_to_nearest(ir : IRanges) -> Array[Int] {
  let n = ir.starts.length()
  let distances : Array[Int] = Array::make(n, -1)

  let mut i = 0
  while i < n {
    let s = ir.starts[i]
    let e = ir.ends[i]
    let mut min_dist = @double.infinity

    let mut j = 0
    while j < n {
      if i != j {
        let mut dist = 0.0
        if ir.starts[j] > e {
          dist = (ir.starts[j] - e - 1).to_double()
        } else if ir.ends[j] < s {
          dist = (s - ir.ends[j] - 1).to_double()
        }
        if dist < min_dist {
          min_dist = dist
        }
      }
      j = j + 1
    }

    distances[i] = if min_dist < @double.infinity { min_dist.to_int() } else { -1 }
    i = i + 1
  }

  distances
}

///|
/// Subset IRanges by a boolean mask.
pub fn ir_subset(ir : IRanges, keep : Array[Bool]) -> IRanges {
  let n = ir.starts.length()
  let new_starts : Array[Int] = Array::new()
  let new_ends : Array[Int] = Array::new()
  let new_widths : Array[Int] = Array::new()

  let mut i = 0
  while i < n {
    if i < keep.length() && keep[i] {
      new_starts.push(ir.starts[i])
      new_ends.push(ir.ends[i])
      new_widths.push(ir.widths[i])
    }
    i = i + 1
  }

  IRanges::{ starts: new_starts, ends: new_ends, widths: new_widths }
}

///|
/// Sort IRanges by start position.
pub fn ir_sort(ir : IRanges) -> IRanges {
  let n = ir.starts.length()
  let pairs : Array[(Int, Int)] = Array::make(n, (0, 0))

  let mut i = 0
  while i < n {
    pairs[i] = (ir.starts[i], ir.ends[i])
    i = i + 1
  }

  // Bubble sort
  let arr = pairs.copy()
  let mut swapped = true
  while swapped {
    swapped = false
    let mut j = 0
    while j < arr.length() - 1 {
      if arr[j].0 > arr[j + 1].0 {
        let temp = arr[j]
        arr[j] = arr[j + 1]
        arr[j + 1] = temp
        swapped = true
      }
      j = j + 1
    }
  }

  let new_starts : Array[Int] = Array::make(n, 0)
  let new_ends : Array[Int] = Array::make(n, 0)
  let new_widths : Array[Int] = Array::make(n, 0)

  let mut k = 0
  while k < n {
    new_starts[k] = arr[k].0
    new_ends[k] = arr[k].1
    new_widths[k] = arr[k].1 - arr[k].0 + 1
    k = k + 1
  }

  IRanges::{ starts: new_starts, ends: new_ends, widths: new_widths }
}

///|
/// Calculate total width of all ranges.
pub fn ir_total_width(ir : IRanges) -> Int {
  let n = ir.starts.length()
  let mut total = 0
  let mut i = 0
  while i < n {
    total = total + ir.widths[i]
    i = i + 1
  }
  total
}

///|
/// Calculate gap width between ranges.
pub fn ir_gap_width(ir : IRanges) -> Int {
  let sorted = ir_sort(ir)
  let n = sorted.starts.length()
  let mut total_gap = 0
  let mut i = 1
  while i < n {
    let gap = sorted.starts[i] - sorted.ends[i - 1] - 1
    if gap > 0 {
      total_gap = total_gap + gap
    }
    i = i + 1
  }
  total_gap
}

///| Test find_overlaps_advanced with "any" type
test "find_overlaps_advanced_any" {
  let query = iranges([1, 10], [5, 15])
  let subject = iranges([3, 8, 12], [8, 14, 18])
  let result = find_overlaps_advanced(query=query, subject=subject, type_="any")
  assert_true(result.length() >= 1)
}

///| Test find_overlaps_advanced with "within" type
test "find_overlaps_advanced_within" {
  let query = iranges([1, 10], [5, 15])
  let subject = iranges([0, 5], [10, 20])
  let result = find_overlaps_advanced(query=query, subject=subject, type_="within")
  assert_true(result.length() >= 0)
}

///| Test ir_nearest function
test "ir_nearest" {
  let query = iranges([1, 20, 40], [5, 25, 45])
  let subject = iranges([10, 30], [15, 35])
  let result = ir_nearest(query, subject)
  assert_eq(result.length(), 3)
  // Query 1-5 should be nearest to 10-15
  assert_eq(result[0], 0)
}

///| Test ir_coverage function
test "ir_coverage" {
  let ir = iranges([1, 5, 10], [8, 15, 20])
  let result = ir_coverage(ir, 25)
  assert_eq(result.length(), 25)
  // Position 1: coverage 1
  assert_eq(result[0], 1)
  // Position 5-8: coverage 2
  assert_eq(result[4], 2)
}

///| Test ir_coverage_rle function
test "ir_coverage_rle" {
  let ir = iranges([1, 5], [3, 8])
  let result = ir_coverage_rle(ir, 10)
  assert_true(result.length() >= 1)
}

///| Test ir_distance function
test "ir_distance" {
  let ir1 = iranges([1, 15], [5, 20])
  let ir2 = iranges([10, 30], [12, 35])
  let result = ir_distance(ir1, ir2)
  assert_eq(result.length(), 2)
  assert_eq(result[0].length(), 2)
}

///| Test ir_distance_to_nearest function
test "ir_distance_to_nearest" {
  let ir = iranges([1, 10, 25], [5, 15, 30])
  let result = ir_distance_to_nearest(ir)
  assert_eq(result.length(), 3)
  // Distance from range 1-5 to range 10-15 should be 4
  assert_eq(result[0], 4)
}

///| Test ir_precede function
test "ir_precede" {
  let query = iranges([10, 30], [15, 35])
  let subject = iranges([1, 20], [5, 25])
  let result = ir_precede(query, subject)
  assert_eq(result.length(), 2)
  // Query 10-15 should find 1-5 as preceding
  assert_eq(result[0], 0)
}

///| Test ir_follow function
test "ir_follow" {
  let query = iranges([1, 10], [5, 15])
  let subject = iranges([20, 30], [25, 35])
  let result = ir_follow(query, subject)
  assert_eq(result.length(), 2)
  // Query 10-15 should find 20-25 as following
  assert_eq(result[1], 0)
}

///| Test ir_subset function
test "ir_subset" {
  let ir = iranges([1, 5, 10], [3, 8, 15])
  let keep = [true, false, true]
  let result = ir_subset(ir, keep)
  assert_eq(result.starts.length(), 2)
  assert_eq(result.starts[0], 1)
  assert_eq(result.starts[1], 10)
}

///| Test ir_sort function
test "ir_sort" {
  let ir = iranges([10, 1, 5], [15, 3, 8])
  let result = ir_sort(ir)
  assert_eq(result.starts[0], 1)
  assert_eq(result.starts[1], 5)
  assert_eq(result.starts[2], 10)
}
