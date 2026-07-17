///|
/// GenomicRanges - Genomic range operations inspired by R/Bioconductor GenomicRanges.

///|
pub enum Strand {
  Plus
  Minus
  Star
} derive(Eq)

///|
pub fn strand_plus() -> Strand {
  Strand::Plus
}

///|
pub fn strand_minus() -> Strand {
  Strand::Minus
}

///|
pub fn strand_star() -> Strand {
  Strand::Star
}

///|
pub struct GRanges {
  seqnames : Array[String]
  starts : Array[Int]
  ends : Array[Int]
  strands : Array[Strand]
  widths : Array[Int]
}

///|
pub fn granges(
  seqnames : Array[String],
  ranges : Array[(Int, Int)],
  strands : Array[Strand],
) -> GRanges {
  let n = seqnames.length()
  let starts : Array[Int] = Array::make(n, 0)
  let ends : Array[Int] = Array::make(n, 0)
  let widths : Array[Int] = Array::make(n, 0)

  let mut i = 0
  while i < n {
    starts[i] = ranges[i].0
    ends[i] = ranges[i].1
    widths[i] = ranges[i].1 - ranges[i].0 + 1
    i = i + 1
  }

  GRanges::{ seqnames, starts, ends, strands, widths }
}

///|
pub fn granges_single(
  seqname : String,
  start : Int,
  end : Int,
  strand : Strand,
) -> GRanges {
  GRanges::{
    seqnames: [seqname],
    starts: [start],
    ends: [end],
    strands: [strand],
    widths: [end - start + 1],
  }
}

///|
pub fn granges_length(gr : GRanges) -> Int {
  gr.seqnames.length()
}

///|
pub fn granges_shift(gr : GRanges, shift : Int) -> GRanges {
  let n = gr.seqnames.length()
  let new_starts : Array[Int] = Array::make(n, 0)
  let new_ends : Array[Int] = Array::make(n, 0)

  let mut i = 0
  while i < n {
    new_starts[i] = gr.starts[i] + shift
    new_ends[i] = gr.ends[i] + shift
    i = i + 1
  }

  GRanges::{
    seqnames: gr.seqnames.clone(),
    starts: new_starts,
    ends: new_ends,
    strands: gr.strands.clone(),
    widths: gr.widths.clone(),
  }
}

///|
pub fn granges_narrow(gr : GRanges, start : Int, end : Int) -> GRanges {
  let n = gr.seqnames.length()
  let new_starts : Array[Int] = Array::make(n, 0)
  let new_ends : Array[Int] = Array::make(n, 0)
  let new_widths : Array[Int] = Array::make(n, 0)

  let mut i = 0
  while i < n {
    let r_start = gr.starts[i]
    let r_end = gr.ends[i]
    let new_start = if start < 1 { r_start } else { r_start + start - 1 }
    let new_end = if end < 1 { r_end } else { r_start + end - 1 }
    let clamped_start = if new_start < r_start { r_start } else { new_start }
    let clamped_end = if new_end > r_end { r_end } else { new_end }
    new_starts[i] = clamped_start
    new_ends[i] = clamped_end
    new_widths[i] = clamped_end - clamped_start + 1
    i = i + 1
  }

  GRanges::{
    seqnames: gr.seqnames.clone(),
    starts: new_starts,
    ends: new_ends,
    strands: gr.strands.clone(),
    widths: new_widths,
  }
}

///|
pub fn granges_resize(gr : GRanges, width : Int, fix : String) -> GRanges {
  let n = gr.seqnames.length()
  let new_starts : Array[Int] = Array::make(n, 0)
  let new_ends : Array[Int] = Array::make(n, 0)
  let new_widths : Array[Int] = Array::make(n, 0)

  let mut i = 0
  while i < n {
    let start = gr.starts[i]
    let end = gr.ends[i]
    let old_width = end - start + 1

    if width <= 0 {
      new_starts[i] = start
      new_ends[i] = start - 1
      new_widths[i] = 0
    } else if width == old_width {
      new_starts[i] = start
      new_ends[i] = end
      new_widths[i] = width
    } else if fix == "start" {
      new_starts[i] = start
      new_ends[i] = start + width - 1
      new_widths[i] = width
    } else if fix == "end" {
      new_starts[i] = end - width + 1
      new_ends[i] = end
      new_widths[i] = width
    } else {
      let diff = width - old_width
      let left = diff / 2
      let right = diff - left
      new_starts[i] = start - left
      new_ends[i] = end + right
      new_widths[i] = width
    }
    i = i + 1
  }

  GRanges::{
    seqnames: gr.seqnames.clone(),
    starts: new_starts,
    ends: new_ends,
    strands: gr.strands.clone(),
    widths: new_widths,
  }
}

///|
pub fn granges_flank(
  gr : GRanges,
  width : Int,
  start : Bool,
  both : Bool,
) -> GRanges {
  let n = gr.seqnames.length()
  let new_starts : Array[Int] = Array::make(n, 0)
  let new_ends : Array[Int] = Array::make(n, 0)
  let new_widths : Array[Int] = Array::make(n, 0)

  let mut i = 0
  while i < n {
    let r_start = gr.starts[i]
    let r_end = gr.ends[i]

    if both {
      new_starts[i] = r_start - width
      new_ends[i] = r_end + width
      new_widths[i] = r_end - r_start + 1 + 2 * width
    } else if start {
      new_starts[i] = r_start - width
      new_ends[i] = r_start - 1
      new_widths[i] = width
    } else {
      new_starts[i] = r_end + 1
      new_ends[i] = r_end + width
      new_widths[i] = width
    }
    i = i + 1
  }

  GRanges::{
    seqnames: gr.seqnames.clone(),
    starts: new_starts,
    ends: new_ends,
    strands: gr.strands.clone(),
    widths: new_widths,
  }
}

///|
fn sort_ranges_by_start(ranges : Array[(Int, Int)]) -> Array[(Int, Int)] {
  let n = ranges.length()
  let sorted = ranges.clone()

  let mut i = 0
  while i < n - 1 {
    let mut j = 0
    while j < n - i - 1 {
      if sorted[j].0 > sorted[j + 1].0 {
        let temp = sorted[j]
        sorted[j] = sorted[j + 1]
        sorted[j + 1] = temp
      }
      j = j + 1
    }
    i = i + 1
  }

  sorted
}

///|
fn sort_ints(arr : Array[Int]) -> Array[Int] {
  let n = arr.length()
  let sorted = arr.clone()

  let mut i = 0
  while i < n - 1 {
    let mut j = 0
    while j < n - i - 1 {
      if sorted[j] > sorted[j + 1] {
        let temp = sorted[j]
        sorted[j] = sorted[j + 1]
        sorted[j + 1] = temp
      }
      j = j + 1
    }
    i = i + 1
  }

  sorted
}

///|
fn merge_overlapping_ranges(ranges : Array[(Int, Int)]) -> Array[(Int, Int)] {
  if ranges.length() == 0 {
    return []
  }

  let sorted = sort_ranges_by_start(ranges)
  let merged : Array[(Int, Int)] = [sorted[0]]

  let mut i = 1
  while i < sorted.length() {
    let last = merged[merged.length() - 1]
    let current = sorted[i]

    if current.0 <= last.1 + 1 {
      let new_end = if current.1 > last.1 { current.1 } else { last.1 }
      merged[merged.length() - 1] = (last.0, new_end)
    } else {
      merged.push(current)
    }
    i = i + 1
  }

  merged
}

///|
pub fn granges_reduce(gr : GRanges, _min_gapwidth : Int) -> GRanges {
  let results : GRanges = GRanges::{
    seqnames: [],
    starts: [],
    ends: [],
    strands: [],
    widths: [],
  }

  let seen = Map([], capacity=16)
  let mut i = 0
  while i < gr.seqnames.length() {
    let seqname = gr.seqnames[i]
    if !seen.contains(seqname) {
      seen[seqname] = ()

      let ranges : Array[(Int, Int)] = []
      let mut j = 0
      while j < gr.seqnames.length() {
        if gr.seqnames[j] == seqname {
          ranges.push((gr.starts[j], gr.ends[j]))
        }
        j = j + 1
      }

      let merged = merge_overlapping_ranges(ranges)
      for r in merged {
        results.seqnames.push(seqname)
        results.starts.push(r.0)
        results.ends.push(r.1)
        results.strands.push(Strand::Star)
        results.widths.push(r.1 - r.0 + 1)
      }
    }
    i = i + 1
  }

  results
}

///|
pub fn granges_disjoin(gr : GRanges) -> GRanges {
  let results : GRanges = GRanges::{
    seqnames: [],
    starts: [],
    ends: [],
    strands: [],
    widths: [],
  }

  let seen = Map([], capacity=16)
  let mut i = 0
  while i < gr.seqnames.length() {
    let seqname = gr.seqnames[i]
    if !seen.contains(seqname) {
      seen[seqname] = ()

      let breakpoints : Array[Int] = []
      let mut j = 0
      while j < gr.seqnames.length() {
        if gr.seqnames[j] == seqname {
          breakpoints.push(gr.starts[j])
          breakpoints.push(gr.ends[j] + 1)
        }
        j = j + 1
      }

      let sorted = sort_ints(breakpoints)
      let unique : Array[Int] = []
      let mut k = 0
      while k < sorted.length() {
        if k == 0 || sorted[k] != sorted[k - 1] {
          unique.push(sorted[k])
        }
        k = k + 1
      }

      let mut l = 0
      while l < unique.length() - 1 {
        let start = unique[l]
        let end = unique[l + 1] - 1
        if start <= end {
          results.seqnames.push(seqname)
          results.starts.push(start)
          results.ends.push(end)
          results.strands.push(Strand::Star)
          results.widths.push(end - start + 1)
        }
        l = l + 1
      }
    }
    i = i + 1
  }

  results
}

///|
pub fn granges_union(gr1 : GRanges, gr2 : GRanges) -> GRanges {
  let combined : GRanges = GRanges::{
    seqnames: [],
    starts: [],
    ends: [],
    strands: [],
    widths: [],
  }

  let mut i = 0
  while i < gr1.seqnames.length() {
    combined.seqnames.push(gr1.seqnames[i])
    combined.starts.push(gr1.starts[i])
    combined.ends.push(gr1.ends[i])
    combined.strands.push(gr1.strands[i])
    combined.widths.push(gr1.widths[i])
    i = i + 1
  }

  let mut j = 0
  while j < gr2.seqnames.length() {
    combined.seqnames.push(gr2.seqnames[j])
    combined.starts.push(gr2.starts[j])
    combined.ends.push(gr2.ends[j])
    combined.strands.push(gr2.strands[j])
    combined.widths.push(gr2.widths[j])
    j = j + 1
  }

  granges_reduce(combined, 0)
}

///|
pub fn granges_intersect(gr1 : GRanges, gr2 : GRanges) -> GRanges {
  let results : GRanges = GRanges::{
    seqnames: [],
    starts: [],
    ends: [],
    strands: [],
    widths: [],
  }

  let mut i = 0
  while i < gr1.seqnames.length() {
    let seqname1 = gr1.seqnames[i]
    let start1 = gr1.starts[i]
    let end1 = gr1.ends[i]

    let mut j = 0
    while j < gr2.seqnames.length() {
      let seqname2 = gr2.seqnames[j]
      let start2 = gr2.starts[j]
      let end2 = gr2.ends[j]

      if seqname1 == seqname2 {
        let overlap_start = if start1 > start2 { start1 } else { start2 }
        let overlap_end = if end1 < end2 { end1 } else { end2 }

        if overlap_start <= overlap_end {
          results.seqnames.push(seqname1)
          results.starts.push(overlap_start)
          results.ends.push(overlap_end)
          results.strands.push(Strand::Star)
          results.widths.push(overlap_end - overlap_start + 1)
        }
      }
      j = j + 1
    }
    i = i + 1
  }

  results
}

///|
pub fn granges_setdiff(gr1 : GRanges, gr2 : GRanges) -> GRanges {
  let results : GRanges = GRanges::{
    seqnames: [],
    starts: [],
    ends: [],
    strands: [],
    widths: [],
  }

  let mut i = 0
  while i < gr1.seqnames.length() {
    let seqname = gr1.seqnames[i]
    let start = gr1.starts[i]
    let end = gr1.ends[i]

    let mut remaining : Array[(Int, Int)] = [(start, end)]

    let mut j = 0
    while j < gr2.seqnames.length() {
      if gr2.seqnames[j] == seqname {
        let exclude_start = gr2.starts[j]
        let exclude_end = gr2.ends[j]

        let new_remaining : Array[(Int, Int)] = []
        for r in remaining {
          if r.1 < exclude_start || r.0 > exclude_end {
            new_remaining.push(r)
          } else {
            if r.0 < exclude_start {
              new_remaining.push((r.0, exclude_start - 1))
            }
            if r.1 > exclude_end {
              new_remaining.push((exclude_end + 1, r.1))
            }
          }
        }
        remaining = new_remaining
      }
      j = j + 1
    }

    for r in remaining {
      results.seqnames.push(seqname)
      results.starts.push(r.0)
      results.ends.push(r.1)
      results.strands.push(Strand::Star)
      results.widths.push(r.1 - r.0 + 1)
    }
    i = i + 1
  }

  results
}

///|
pub fn count_overlaps(query : GRanges, subject : GRanges) -> Array[Int] {
  let n = granges_length(query)
  let counts : Array[Int] = Array::make(n, 0)

  let mut i = 0
  while i < n {
    let seqname_q = query.seqnames[i]
    let start_q = query.starts[i]
    let end_q = query.ends[i]

    let mut count = 0
    let mut j = 0
    while j < granges_length(subject) {
      if subject.seqnames[j] == seqname_q {
        let start_s = subject.starts[j]
        let end_s = subject.ends[j]
        if start_q <= end_s && end_q >= start_s {
          count = count + 1
        }
      }
      j = j + 1
    }
    counts[i] = count
    i = i + 1
  }

  counts
}

///|
pub fn find_overlaps(query : GRanges, subject : GRanges) -> Array[(Int, Int)] {
  let results : Array[(Int, Int)] = []

  let mut i = 0
  while i < granges_length(query) {
    let seqname_q = query.seqnames[i]
    let start_q = query.starts[i]
    let end_q = query.ends[i]

    let mut j = 0
    while j < granges_length(subject) {
      if subject.seqnames[j] == seqname_q {
        let start_s = subject.starts[j]
        let end_s = subject.ends[j]
        if start_q <= end_s && end_q >= start_s {
          results.push((i, j))
        }
      }
      j = j + 1
    }
    i = i + 1
  }

  results
}

///|
pub fn granges_distance(gr1 : GRanges, gr2 : GRanges) -> Array[Double] {
  let n = granges_length(gr1)
  let distances : Array[Double] = Array::make(n, 0.0)

  let mut i = 0
  while i < n {
    let seqname1 = gr1.seqnames[i]
    let start1 = gr1.starts[i]
    let end1 = gr1.ends[i]

    let mut min_dist = @double.infinity
    let mut found = false

    let mut j = 0
    while j < granges_length(gr2) {
      if gr2.seqnames[j] == seqname1 {
        found = true
        let start2 = gr2.starts[j]
        let end2 = gr2.ends[j]

        if end1 < start2 {
          let dist = (start2 - end1 - 1).to_double()
          if dist < min_dist {
            min_dist = dist
          }
        } else if end2 < start1 {
          let dist = (start1 - end2 - 1).to_double()
          if dist < min_dist {
            min_dist = dist
          }
        } else {
          min_dist = 0.0
        }
      }
      j = j + 1
    }

    distances[i] = if found { min_dist } else { @double.not_a_number }
    i = i + 1
  }

  distances
}

///|
pub fn nearest(gr1 : GRanges, gr2 : GRanges) -> Array[Int] {
  let n = granges_length(gr1)
  let indices : Array[Int] = Array::make(n, -1)

  let mut i = 0
  while i < n {
    let seqname1 = gr1.seqnames[i]
    let mid1 = (gr1.starts[i] + gr1.ends[i]) / 2

    let mut min_dist = @double.infinity
    let mut min_idx = -1

    let mut j = 0
    while j < granges_length(gr2) {
      if gr2.seqnames[j] == seqname1 {
        let mid2 = (gr2.starts[j] + gr2.ends[j]) / 2
        let dist = (mid1 - mid2).abs().to_double()

        if dist < min_dist {
          min_dist = dist
          min_idx = j
        }
      }
      j = j + 1
    }

    indices[i] = min_idx
    i = i + 1
  }

  indices
}

///|
pub fn strand_to_char(s : Strand) -> String {
  match s {
    Strand::Plus => "+"
    Strand::Minus => "-"
    Strand::Star => "*"
  }
}
