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

///|
/// Find ranges in gr2 that precede ranges in gr1 (gr2 end < gr1 start).
pub fn precede(gr1 : GRanges, gr2 : GRanges) -> Array[Int] {
  let n = granges_length(gr1)
  let indices : Array[Int] = Array::make(n, -1)

  let mut i = 0
  while i < n {
    let seqname1 = gr1.seqnames[i]
    let start1 = gr1.starts[i]

    let mut best_idx = -1
    let mut best_end = @double.neg_infinity

    let mut j = 0
    while j < granges_length(gr2) {
      if gr2.seqnames[j] == seqname1 {
        let end2 = gr2.ends[j]
        if end2 < start1 {
          if end2.to_double() > best_end {
            best_end = end2.to_double()
            best_idx = j
          }
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
/// Find ranges in gr2 that follow ranges in gr1 (gr2 start > gr1 end).
pub fn follow(gr1 : GRanges, gr2 : GRanges) -> Array[Int] {
  let n = granges_length(gr1)
  let indices : Array[Int] = Array::make(n, -1)

  let mut i = 0
  while i < n {
    let seqname1 = gr1.seqnames[i]
    let end1 = gr1.ends[i]

    let mut best_idx = -1
    let mut best_start = @double.infinity

    let mut j = 0
    while j < granges_length(gr2) {
      if gr2.seqnames[j] == seqname1 {
        let start2 = gr2.starts[j]
        if start2 > end1 {
          if start2.to_double() < best_start {
            best_start = start2.to_double()
            best_idx = j
          }
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
/// Calculate coverage of genomic ranges.
/// Returns a map of seqname to coverage array.
pub fn coverage_granges(gr : GRanges, seqlengths : Map[String, Int]) -> Map[String, Array[Int]] {
  let result : Map[String, Array[Int]] = Map([], capacity=10)

  let seqnames = gr.seqnames
  let n = granges_length(gr)
  let keys = seqlengths.keys().collect()

  let mut seq_idx = 0
  while seq_idx < keys.length() {
    let seqname = keys[seq_idx]
    let seqlen = seqlengths[seqname]
    let cov = Array::make(seqlen, 0)

    let mut i = 0
    while i < n {
      if seqnames[i] == seqname {
        let start = gr.starts[i]
        let end = gr.ends[i]
        let clamped_start = if start < 1 { 1 } else { start }
        let clamped_end = if end > seqlen { seqlen } else { end }
        let mut pos = clamped_start
        while pos <= clamped_end {
          cov[pos - 1] = cov[pos - 1] + 1
          pos = pos + 1
        }
      }
      i = i + 1
    }

    result[seqname] = cov
    seq_idx = seq_idx + 1
  }

  result
}

///|
/// Calculate coverage as a GRanges object (run-length encoding).
pub fn coverage_as_granges(gr : GRanges, seqlengths : Map[String, Int]) -> GRanges {
  let result_seqnames : Array[String] = Array::new()
  let result_starts : Array[Int] = Array::new()
  let result_ends : Array[Int] = Array::new()
  let result_strands : Array[Strand] = Array::new()
  let result_widths : Array[Int] = Array::new()

  let cov_map = coverage_granges(gr, seqlengths)

  let seqnames = seqlengths.keys().collect()
  let mut seq_idx = 0
  while seq_idx < seqnames.length() {
    let seqname = seqnames[seq_idx]
    match cov_map.get(seqname) {
      Some(cov) => {
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
            result_seqnames.push(seqname)
            result_starts.push(start)
            result_ends.push(end)
            result_strands.push(Strand::Star)
            result_widths.push(end - start + 1)
            i = j
          } else {
            i = i + 1
          }
        }
      }
      None => ()
    }
    seq_idx = seq_idx + 1
  }

  GRanges::{
    seqnames: result_seqnames,
    starts: result_starts,
    ends: result_ends,
    strands: result_strands,
    widths: result_widths,
  }
}

///|
/// Calculate the distance to nearest neighbor (sorted by position).
pub fn distance_to_nearest(gr : GRanges) -> Array[Int] {
  let n = granges_length(gr)
  let distances : Array[Int] = Array::make(n, -1)

  let mut i = 0
  while i < n {
    let seqname = gr.seqnames[i]
    let start = gr.starts[i]
    let end = gr.ends[i]

    let mut min_dist = @double.infinity

    let mut j = 0
    while j < n {
      if i != j && gr.seqnames[j] == seqname {
        let mut dist = 0.0
        if gr.starts[j] > end {
          dist = (gr.starts[j] - end - 1).to_double()
        } else if gr.ends[j] < start {
          dist = (start - gr.ends[j] - 1).to_double()
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
/// Subset GRanges by a logical condition (as boolean indices).
pub fn granges_subset(gr : GRanges, keep : Array[Bool]) -> GRanges {
  let n = gr.seqnames.length()
  let new_seqnames : Array[String] = Array::new()
  let new_starts : Array[Int] = Array::new()
  let new_ends : Array[Int] = Array::new()
  let new_strands : Array[Strand] = Array::new()
  let new_widths : Array[Int] = Array::new()

  let mut i = 0
  while i < n {
    if i < keep.length() && keep[i] {
      new_seqnames.push(gr.seqnames[i])
      new_starts.push(gr.starts[i])
      new_ends.push(gr.ends[i])
      new_strands.push(gr.strands[i])
      new_widths.push(gr.widths[i])
    }
    i = i + 1
  }

  GRanges::{
    seqnames: new_seqnames,
    starts: new_starts,
    ends: new_ends,
    strands: new_strands,
    widths: new_widths,
  }
}

///|
/// Sort GRanges by seqname then start position.
pub fn granges_sort(gr : GRanges) -> GRanges {
  let n = gr.seqnames.length()
  let indices = Array::make(n, 0)
  let mut i = 0
  while i < n {
    indices[i] = i
    i = i + 1
  }

  let mut sorted_indices : Array[Int] = Array::new()
  // Simple bubble sort for demonstration
  let arr = indices.copy()
  let mut swapped = true
  while swapped {
    swapped = false
    let mut j = 0
    while j < arr.length() - 1 {
      let idx1 = arr[j]
      let idx2 = arr[j + 1]
      let cmp = compare_granges_pos(gr, idx1, idx2)
      if cmp > 0 {
        let temp = arr[j]
        arr[j] = arr[j + 1]
        arr[j + 1] = temp
        swapped = true
      }
      j = j + 1
    }
  }

  sorted_indices = arr

  let new_seqnames : Array[String] = Array::new()
  let new_starts : Array[Int] = Array::new()
  let new_ends : Array[Int] = Array::new()
  let new_strands : Array[Strand] = Array::new()
  let new_widths : Array[Int] = Array::new()

  let mut k = 0
  while k < sorted_indices.length() {
    let idx = sorted_indices[k]
    new_seqnames.push(gr.seqnames[idx])
    new_starts.push(gr.starts[idx])
    new_ends.push(gr.ends[idx])
    new_strands.push(gr.strands[idx])
    new_widths.push(gr.widths[idx])
    k = k + 1
  }

  GRanges::{
    seqnames: new_seqnames,
    starts: new_starts,
    ends: new_ends,
    strands: new_strands,
    widths: new_widths,
  }
}

fn compare_granges_pos(gr : GRanges, idx1 : Int, idx2 : Int) -> Int {
  if gr.seqnames[idx1] < gr.seqnames[idx2] {
    -1
  } else if gr.seqnames[idx1] > gr.seqnames[idx2] {
    1
  } else if gr.starts[idx1] < gr.starts[idx2] {
    -1
  } else if gr.starts[idx1] > gr.starts[idx2] {
    1
  } else {
    0
  }
}

///| Test precede function
test "precede" {
  let gr1 = granges(["chr1", "chr1"], [(10, 20), (30, 40)], [strand_plus(), strand_plus()])
  let gr2 = granges(["chr1"], [(5, 8)], [strand_plus()])
  let result = precede(gr1, gr2)
  assert_eq(result.length(), 2)
  // First range in gr1: 10-20, should find 5-8 as preceding
  assert_eq(result[0], 0)
}

///| Test follow function
test "follow" {
  let gr1 = granges(["chr1", "chr1"], [(10, 20), (30, 40)], [strand_plus(), strand_plus()])
  let gr2 = granges(["chr1"], [(50, 60)], [strand_plus()])
  let result = follow(gr1, gr2)
  assert_eq(result.length(), 2)
  // Second range in gr1: 30-40, should find 50-60 as following
  assert_eq(result[1], 0)
}

///| Test coverage_granges function
test "coverage_granges" {
  let gr = granges(["chr1", "chr1"], [(1, 5), (3, 8)], [strand_plus(), strand_plus()])
  let seqlengths : Map[String, Int] = Map([], capacity=1)
  seqlengths["chr1"] = 10
  let result = coverage_granges(gr, seqlengths)
  assert_true(result.contains("chr1"))
  let cov = result["chr1"]
  assert_eq(cov.length(), 10)
  // Position 1-2: coverage 1
  assert_eq(cov[0], 1)
  // Position 3-5: coverage 2
  assert_eq(cov[2], 2)
  // Position 6-8: coverage 1
  assert_eq(cov[5], 1)
}

///| Test distance_to_nearest function
test "distance_to_nearest" {
  let gr = granges(["chr1", "chr1", "chr1"], [(1, 5), (15, 20), (30, 35)], [
    strand_plus(),
    strand_plus(),
    strand_plus(),
  ])
  let result = distance_to_nearest(gr)
  assert_eq(result.length(), 3)
  // Distance from range 1 to range 2 should be 9 (15-5-1)
  assert_eq(result[0], 9)
}

///| Test coverage_as_granges function
test "coverage_as_granges" {
  let gr = granges(["chr1", "chr1"], [(1, 5), (3, 8)], [strand_plus(), strand_plus()])
  let seqlengths : Map[String, Int] = Map([], capacity=1)
  seqlengths["chr1"] = 10
  let result = coverage_as_granges(gr, seqlengths)
  assert_true(granges_length(result) >= 2)
}
