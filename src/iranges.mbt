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
