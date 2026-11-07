// Bioframe third-round: visualization helpers (items 6-9)
// pile, staircase, snap_to_grid, midpoints

///|
pub struct Bf3PileRow {
  chrom : String
  start : Int32
  end_ : Int32
  row : Int32 // vertical stacking level (0-based)
  name : String
}

///|
pub fn bf3_pile(
  intervals : Array[(String, Int32, Int32, String)],
) -> Array[Bf3PileRow] {
  // Stack intervals vertically: each row contains non-overlapping intervals
  // Greedy algorithm: assign each interval to the lowest available row
  let sorted = intervals.copy()
  sorted.sort_by(fn(a, b) {
    let (_, sa, _, _) = a
    let (_, sb, _, _) = b
    if sa < sb {
      -1
    } else if sa > sb {
      1
    } else {
      0
    }
  })
  // Track the end position of the last interval in each row per chrom
  let row_ends : Array[(String, Int32)] = [] // (chrom, end) per row
  let result : Array[Bf3PileRow] = []
  for i in 0..<sorted.length() {
    let (chrom, start, end_, name) = sorted[i]
    let mut assigned_row = -1
    for r in 0..<row_ends.length() {
      let (r_chrom, r_end) = row_ends[r]
      if r_chrom == chrom && r_end <= start {
        assigned_row = r
        row_ends[r] = (chrom, end_)
        break
      }
    }
    if assigned_row < 0 {
      assigned_row = row_ends.length()
      row_ends.push((chrom, end_))
    }
    result.push({ chrom, start, end_, row: assigned_row, name, })
  }
  result
}

///|
pub struct Bf3StairStep {
  pos : Int32
  value : Double // cumulative coverage at this position
  step_type : String // "up" or "down"
}

///|
pub fn bf3_staircase(
  intervals : Array[(Int32, Int32)],
  genome_size : Int32,
) -> Array[Bf3StairStep] {
  // Generate staircase (cumulative coverage) from intervals
  // Each interval start = +1 step, each end = -1 step
  let events : Array[(Int32, Int32)] = [] // (pos, delta)
  for i in 0..<intervals.length() {
    let (start, end_) = intervals[i]
    events.push((start, 1))
    events.push((end_, -1))
  }
  events.sort_by(fn(a, b) {
    let (pa, _) = a
    let (pb, _) = b
    if pa < pb {
      -1
    } else if pa > pb {
      1
    } else {
      0
    }
  })
  let steps : Array[Bf3StairStep] = []
  let mut cumulative = 0.0
  for i in 0..<events.length() {
    let (pos, delta) = events[i]
    cumulative = cumulative + delta.to_double()
    let step_type = if delta > 0 { "up" } else { "down" }
    steps.push({ pos, value: cumulative, step_type, })
  }
  // Ensure coverage at genome end is 0
  if steps.length() > 0 {
    let last = steps[steps.length() - 1]
    if last.pos < genome_size {
      steps.push({ pos: genome_size, value: 0.0, step_type: "down", })
    }
  }
  steps
}

///|
pub fn bf3_snap_to_grid(
  intervals : Array[(String, Int32, Int32)],
  grid_size : Int32,
) -> Array[(String, Int32, Int32)] {
  // Snap interval boundaries to grid (bin edges)
  let result : Array[(String, Int32, Int32)] = []
  for i in 0..<intervals.length() {
    let (chrom, start, end_) = intervals[i]
    let snapped_start = start / grid_size * grid_size
    let snapped_end = (end_ + grid_size - 1) / grid_size * grid_size
    result.push((chrom, snapped_start, snapped_end))
  }
  result
}

///|
pub fn bf3_midpoints(
  intervals : Array[(String, Int32, Int32)],
) -> Array[(String, Int32)] {
  // Calculate midpoint of each interval
  let result : Array[(String, Int32)] = []
  for i in 0..<intervals.length() {
    let (chrom, start, end_) = intervals[i]
    result.push((chrom, (start + end_) / 2))
  }
  result
}

///|
pub fn bf3_midpoints_with_offset(
  intervals : Array[(String, Int32, Int32)],
  offset : Double,
) -> Array[(String, Int32)] {
  // Calculate midpoint with fractional offset (0.0 = start, 0.5 = middle, 1.0 = end)
  let result : Array[(String, Int32)] = []
  for i in 0..<intervals.length() {
    let (chrom, start, end_) = intervals[i]
    let mid = start.to_double() + (end_ - start).to_double() * offset
    result.push((chrom, mid.floor().to_int()))
  }
  result
}
