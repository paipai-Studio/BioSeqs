///|
/// GenomicRanges Gaps - Advanced GRanges operations.
/// Inspired by Bioconductor GenomicRanges package.
///
/// Provides local types and functions for:
/// - gaps: Find gaps (uncovered regions) between ranges
/// - reduce: Merge overlapping or adjacent ranges
/// - disjoin: Split ranges into non-overlapping segments
/// - setdiff: Set difference between two GRanges
/// - intersect: Intersection of two GRanges
/// - union: Union of two GRanges
/// - compute_coverage: Compute per-base coverage vector
/// - trim: Trim ranges to chromosome bounds
/// - promoters: Compute promoter regions

///|
/// Local GRanges-compatible structure using String strands.
pub struct GRangesLocal {
  seqnames : Array[String]
  starts : Array[Int]
  ends : Array[Int]
  strands : Array[String]
}

///|
/// Result of grg_gaps: gaps between ranges with lengths.
pub struct GapResult {
  seqnames : Array[String]
  starts : Array[Int]
  ends : Array[Int]
  gap_lengths : Array[Int]
}

///|
/// Result of grg_disjoin: disjoined segments with segment ids.
pub struct DisjoinResult {
  seqnames : Array[String]
  starts : Array[Int]
  ends : Array[Int]
  segment_ids : Array[Int]
}

///|
/// Result of grg_reduce: merged ranges with merge count.
pub struct GrgReduceResult {
  seqnames : Array[String]
  starts : Array[Int]
  ends : Array[Int]
  merged_count : Array[Int]
}

///|
/// Construct a GRangesLocal from sequence names, starts, ends, and optional strands.
/// When strands are not provided, all ranges default to "+".
pub fn grg_granges(
  seqnames : Array[String],
  starts : Array[Int],
  ends : Array[Int],
  strands? : Array[String],
) -> GRangesLocal {
  let n = seqnames.length()
  let final_strands : Array[String] = Array::make(n, "+")
  match strands {
    Some(s) => {
      if s.length() == n {
        let mut i = 0
        while i < n {
          final_strands[i] = s[i]
          i = i + 1
        }
      }
    }
    None => ()
  }
  GRangesLocal::{ seqnames, starts, ends, strands: final_strands }
}

///|
/// Create an empty GRangesLocal.
pub fn grg_empty() -> GRangesLocal {
  GRangesLocal::{ seqnames: [], starts: [], ends: [], strands: [] }
}

///|
/// Return the number of ranges in a GRangesLocal.
pub fn grg_length(gr : GRangesLocal) -> Int {
  gr.seqnames.length()
}

///|
/// Helper: sort an array of (start, end) pairs by start position.
fn grg_sort_ranges(ranges : Array[(Int, Int)]) -> Array[(Int, Int)] {
  let n = ranges.length()
  let sorted = ranges.copy()
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
/// Helper: sort an array of integers (bubble sort).
fn grg_sort_ints(arr : Array[Int]) -> Array[Int] {
  let n = arr.length()
  let sorted = arr.copy()
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
/// Find gaps between ranges within each sequence.
///
/// Parameters:
/// - gr: Input GRangesLocal
/// - seqlengths: Map from seqname to chromosome length
///
/// Returns a GapResult with seqnames, starts, ends and gap_lengths.
/// Gaps are the regions not covered by any range in gr, bounded by 1
/// and the seqlength for that chromosome.
pub fn grg_gaps(gr : GRangesLocal, seqlengths : Map[String, Int]) -> GapResult {
  let result_seqnames : Array[String] = []
  let result_starts : Array[Int] = []
  let result_ends : Array[Int] = []
  let result_lengths : Array[Int] = []

  if gr.seqnames.length() == 0 {
    return GapResult::{
      seqnames: result_seqnames,
      starts: result_starts,
      ends: result_ends,
      gap_lengths: result_lengths,
    }
  }

  let seen = Map([], capacity=16)
  let mut i = 0
  while i < gr.seqnames.length() {
    let seqname = gr.seqnames[i]
    if !seen.contains(seqname) {
      seen[seqname] = ()

      let seqlen : Int = if seqlengths.contains(seqname) { seqlengths[seqname] } else { 0 }

      if seqlen > 0 {
        let ranges : Array[(Int, Int)] = []
        let mut j = 0
        while j < gr.seqnames.length() {
          if gr.seqnames[j] == seqname {
            ranges.push((gr.starts[j], gr.ends[j]))
          }
          j = j + 1
        }

        let sorted = grg_sort_ranges(ranges)
        let merged : Array[(Int, Int)] = []
        let mut k = 0
        while k < sorted.length() {
          let current = sorted[k]
          if merged.length() == 0 {
            merged.push(current)
          } else {
            let last = merged[merged.length() - 1]
            if current.0 <= last.1 + 1 {
              let new_end = if current.1 > last.1 { current.1 } else { last.1 }
              merged[merged.length() - 1] = (last.0, new_end)
            } else {
              merged.push(current)
            }
          }
          k = k + 1
        }

        let mut m = 0
        while m < merged.length() {
          let gap_start = if m == 0 { 1 } else { merged[m - 1].1 + 1 }
          let gap_end = merged[m].0 - 1
          if gap_start <= gap_end && gap_start >= 1 {
            result_seqnames.push(seqname)
            result_starts.push(gap_start)
            result_ends.push(gap_end)
            result_lengths.push(gap_end - gap_start + 1)
          }
          m = m + 1
        }

        let last_end = if merged.length() == 0 { 0 } else { merged[merged.length() - 1].1 }
        let final_gap_start = last_end + 1
        let final_gap_end = seqlen
        if final_gap_start <= final_gap_end {
          result_seqnames.push(seqname)
          result_starts.push(final_gap_start)
          result_ends.push(final_gap_end)
          result_lengths.push(final_gap_end - final_gap_start + 1)
        }
      }
    }
    i = i + 1
  }

  GapResult::{
    seqnames: result_seqnames,
    starts: result_starts,
    ends: result_ends,
    gap_lengths: result_lengths,
  }
}

///|
/// Merge overlapping or adjacent ranges within each sequence.
///
/// Parameters:
/// - gr: Input GRangesLocal
/// - min_gap: Minimum gap between ranges to merge (default 0).
///   Ranges closer than min_gap+1 bases apart will be merged.
///
/// Returns a GrgReduceResult where each merged region has a count
/// of how many original ranges were merged into it.
pub fn grg_reduce(gr : GRangesLocal, min_gap? : Int) -> GrgReduceResult {
  let result_seqnames : Array[String] = []
  let result_starts : Array[Int] = []
  let result_ends : Array[Int] = []
  let result_counts : Array[Int] = []

  if gr.seqnames.length() == 0 {
    return GrgReduceResult::{
      seqnames: result_seqnames,
      starts: result_starts,
      ends: result_ends,
      merged_count: result_counts,
    }
  }

  let gap = match min_gap {
    Some(v) => v
    None => 0
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

      let sorted = grg_sort_ranges(ranges)
      let merged_starts : Array[Int] = []
      let merged_ends : Array[Int] = []
      let merged_counts : Array[Int] = []

      let mut k = 0
      while k < sorted.length() {
        let current = sorted[k]
        let cnt : Int = 1
        if merged_starts.length() == 0 {
          merged_starts.push(current.0)
          merged_ends.push(current.1)
          merged_counts.push(cnt)
        } else {
          let last_end = merged_ends[merged_ends.length() - 1]
          if current.0 <= last_end + gap + 1 {
            let new_end = if current.1 > last_end { current.1 } else { last_end }
            merged_ends[merged_ends.length() - 1] = new_end
            merged_counts[merged_counts.length() - 1] =
              merged_counts[merged_counts.length() - 1] + cnt
          } else {
            merged_starts.push(current.0)
            merged_ends.push(current.1)
            merged_counts.push(cnt)
          }
        }
        k = k + 1
      }

      let mut m = 0
      while m < merged_starts.length() {
        result_seqnames.push(seqname)
        result_starts.push(merged_starts[m])
        result_ends.push(merged_ends[m])
        result_counts.push(merged_counts[m])
        m = m + 1
      }
    }
    i = i + 1
  }

  GrgReduceResult::{
    seqnames: result_seqnames,
    starts: result_starts,
    ends: result_ends,
    merged_count: result_counts,
  }
}

///|
/// Split ranges into non-overlapping segments (disjoin operation).
///
/// Parameters:
/// - gr: Input GRangesLocal
///
/// Returns a DisjoinResult where each segment is a maximal
/// non-overlapping piece. Segments are numbered sequentially
/// (segment_ids) within each sequence.
pub fn grg_disjoin(gr : GRangesLocal) -> DisjoinResult {
  let result_seqnames : Array[String] = []
  let result_starts : Array[Int] = []
  let result_ends : Array[Int] = []
  let result_ids : Array[Int] = []

  if gr.seqnames.length() == 0 {
    return DisjoinResult::{
      seqnames: result_seqnames,
      starts: result_starts,
      ends: result_ends,
      segment_ids: result_ids,
    }
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

      let sorted = grg_sort_ints(breakpoints)
      let unique : Array[Int] = []
      let mut k = 0
      while k < sorted.length() {
        if k == 0 || sorted[k] != sorted[k - 1] {
          unique.push(sorted[k])
        }
        k = k + 1
      }

      let mut seg_id = 0
      let mut l = 0
      while l < unique.length() - 1 {
        let start = unique[l]
        let end = unique[l + 1] - 1
        if start <= end {
          result_seqnames.push(seqname)
          result_starts.push(start)
          result_ends.push(end)
          result_ids.push(seg_id)
          seg_id = seg_id + 1
        }
        l = l + 1
      }
    }
    i = i + 1
  }

  DisjoinResult::{
    seqnames: result_seqnames,
    starts: result_starts,
    ends: result_ends,
    segment_ids: result_ids,
  }
}

///|
/// Compute set difference (gr1 - gr2).
///
/// For each range in gr1, subtract any overlapping parts from gr2,
/// producing zero or more sub-ranges that belong exclusively to gr1.
pub fn grg_setdiff(gr1 : GRangesLocal, gr2 : GRangesLocal) -> GRangesLocal {
  let seqnames : Array[String] = []
  let starts : Array[Int] = []
  let ends : Array[Int] = []
  let strands : Array[String] = []

  if gr1.seqnames.length() == 0 {
    return GRangesLocal::{ seqnames, starts, ends, strands }
  }

  let mut i = 0
  while i < gr1.seqnames.length() {
    let seqname = gr1.seqnames[i]
    let start = gr1.starts[i]
    let end = gr1.ends[i]
    let strand = gr1.strands[i]

    let mut remaining : Array[(Int, Int)] = [(start, end)]

    let mut j = 0
    while j < gr2.seqnames.length() {
      if gr2.seqnames[j] == seqname {
        let exclude_start = gr2.starts[j]
        let exclude_end = gr2.ends[j]

        let new_remaining : Array[(Int, Int)] = []
        let mut k = 0
        while k < remaining.length() {
          let r = remaining[k]
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
          k = k + 1
        }
        remaining = new_remaining
      }
      j = j + 1
    }

    let mut m = 0
    while m < remaining.length() {
      let r = remaining[m]
      seqnames.push(seqname)
      starts.push(r.0)
      ends.push(r.1)
      strands.push(strand)
      m = m + 1
    }
    i = i + 1
  }

  GRangesLocal::{ seqnames, starts, ends, strands }
}

///|
/// Compute the intersection of two GRangesLocal objects.
///
/// For each pair of ranges (one from gr1, one from gr2) that overlap
/// on the same sequence, emit the overlapping sub-range.
pub fn grg_intersect(gr1 : GRangesLocal, gr2 : GRangesLocal) -> GRangesLocal {
  let seqnames : Array[String] = []
  let starts : Array[Int] = []
  let ends : Array[Int] = []
  let strands : Array[String] = []

  if gr1.seqnames.length() == 0 || gr2.seqnames.length() == 0 {
    return GRangesLocal::{ seqnames, starts, ends, strands }
  }

  let mut i = 0
  while i < gr1.seqnames.length() {
    let seqname1 = gr1.seqnames[i]
    let start1 = gr1.starts[i]
    let end1 = gr1.ends[i]
    let strand1 = gr1.strands[i]

    let mut j = 0
    while j < gr2.seqnames.length() {
      let seqname2 = gr2.seqnames[j]
      if seqname1 == seqname2 {
        let start2 = gr2.starts[j]
        let end2 = gr2.ends[j]

        let overlap_start = if start1 > start2 { start1 } else { start2 }
        let overlap_end = if end1 < end2 { end1 } else { end2 }

        if overlap_start <= overlap_end {
          seqnames.push(seqname1)
          starts.push(overlap_start)
          ends.push(overlap_end)
          strands.push(strand1)
        }
      }
      j = j + 1
    }
    i = i + 1
  }

  GRangesLocal::{ seqnames, starts, ends, strands }
}

///|
/// Compute the union of two GRangesLocal objects.
///
/// Concatenates ranges from both inputs and reduces the result
/// so that the output contains no overlapping ranges.
pub fn grg_union(gr1 : GRangesLocal, gr2 : GRangesLocal) -> GRangesLocal {
  let combined_seqnames : Array[String] = []
  let combined_starts : Array[Int] = []
  let combined_ends : Array[Int] = []
  let combined_strands : Array[String] = []

  let mut i = 0
  while i < gr1.seqnames.length() {
    combined_seqnames.push(gr1.seqnames[i])
    combined_starts.push(gr1.starts[i])
    combined_ends.push(gr1.ends[i])
    combined_strands.push(gr1.strands[i])
    i = i + 1
  }

  let mut j = 0
  while j < gr2.seqnames.length() {
    combined_seqnames.push(gr2.seqnames[j])
    combined_starts.push(gr2.starts[j])
    combined_ends.push(gr2.ends[j])
    combined_strands.push(gr2.strands[j])
    j = j + 1
  }

  let combined = GRangesLocal::{
    seqnames: combined_seqnames,
    starts: combined_starts,
    ends: combined_ends,
    strands: combined_strands,
  }

  let reduced = grg_reduce(combined)
  GRangesLocal::{
    seqnames: reduced.seqnames,
    starts: reduced.starts,
    ends: reduced.ends,
    strands: combined_strands,
  }
}

///|
/// Compute a per-base coverage vector for each sequence.
///
/// Parameters:
/// - gr: Input GRangesLocal
/// - seqlengths: Map from seqname to chromosome length
///
/// Returns a Map from seqname to Array[Int] of length seqlength,
/// where index i contains the number of ranges covering position i+1.
pub fn grg_compute_coverage(gr : GRangesLocal, seqlengths : Map[String, Int]) -> Map[String, Array[Int]] {
  let result : Map[String, Array[Int]] = Map([], capacity=10)

  let keys = seqlengths.keys().collect()
  let mut seq_idx = 0
  while seq_idx < keys.length() {
    let seqname = keys[seq_idx]
    let seqlen = seqlengths[seqname]
    let cov : Array[Int] = Array::make(seqlen, 0)

    let mut i = 0
    while i < gr.seqnames.length() {
      if gr.seqnames[i] == seqname {
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
/// Trim ranges to chromosome bounds.
///
/// Parameters:
/// - gr: Input GRangesLocal
/// - seqlengths: Map from seqname to chromosome length
///
/// Returns a GRangesLocal where ranges that extend beyond the
/// chromosome boundaries are clipped to [1, seqlength]. Ranges
/// that fall entirely outside are removed.
pub fn grg_trim(gr : GRangesLocal, seqlengths : Map[String, Int]) -> GRangesLocal {
  let seqnames : Array[String] = []
  let starts : Array[Int] = []
  let ends : Array[Int] = []
  let strands : Array[String] = []

  let mut i = 0
  while i < gr.seqnames.length() {
    let seqname = gr.seqnames[i]
    if seqlengths.contains(seqname) {
      let seqlen = seqlengths[seqname]
      let start = gr.starts[i]
      let end = gr.ends[i]
      let strand = gr.strands[i]

      let clamped_start = if start < 1 { 1 } else { start }
      let clamped_end = if end > seqlen { seqlen } else { end }

      if clamped_start <= clamped_end {
        seqnames.push(seqname)
        starts.push(clamped_start)
        ends.push(clamped_end)
        strands.push(strand)
      }
    }
    i = i + 1
  }

  GRangesLocal::{ seqnames, starts, ends, strands }
}

///|
/// Compute promoter regions for each range.
///
/// Parameters:
/// - gr: Input GRangesLocal
/// - upstream: Bases upstream of the start position (default 2000)
/// - downstream: Bases downstream of the start position (default 200)
///
/// For "+" strand, the promoter is [start - upstream, start + downstream - 1].
/// For "-" strand, the promoter is [end - downstream + 1, end + upstream].
/// For "*" strand, both strands' promoters are combined.
pub fn grg_promoters(
  gr : GRangesLocal,
  upstream? : Int,
  downstream? : Int,
) -> GRangesLocal {
  let seqnames : Array[String] = []
  let starts : Array[Int] = []
  let ends : Array[Int] = []
  let strands : Array[String] = []

  let up = match upstream {
    Some(v) => v
    None => 2000
  }
  let down = match downstream {
    Some(v) => v
    None => 200
  }

  let mut i = 0
  while i < gr.seqnames.length() {
    let seqname = gr.seqnames[i]
    let start = gr.starts[i]
    let end = gr.ends[i]
    let strand = gr.strands[i]

    if strand == "+" {
      let p_start = start - up
      let p_end = start + down - 1
      seqnames.push(seqname)
      starts.push(p_start)
      ends.push(p_end)
      strands.push(strand)
    } else if strand == "-" {
      let p_start = end - down + 1
      let p_end = end + up
      seqnames.push(seqname)
      starts.push(p_start)
      ends.push(p_end)
      strands.push(strand)
    } else {
      let p_start_plus = start - up
      let p_end_plus = start + down - 1
      seqnames.push(seqname)
      starts.push(p_start_plus)
      ends.push(p_end_plus)
      strands.push("+")

      let p_start_minus = end - down + 1
      let p_end_minus = end + up
      seqnames.push(seqname)
      starts.push(p_start_minus)
      ends.push(p_end_minus)
      strands.push("-")
    }
    i = i + 1
  }

  GRangesLocal::{ seqnames, starts, ends, strands }
}

///|
test "grg_granges_constructor" {
  let gr = grg_granges(["chr1", "chr1"], [1, 10], [10, 20])
  assert_eq(gr.seqnames.length(), 2)
  assert_eq(gr.starts[0], 1)
  assert_eq(gr.ends[1], 20)
}

///|
test "grg_granges_default_strands" {
  let gr = grg_granges(["chr1"], [1], [10])
  assert_eq(gr.strands[0], "+")
}

///|
test "grg_length_empty" {
  let gr = grg_empty()
  assert_eq(grg_length(gr), 0)
}

///|
test "grg_gaps_no_gaps_single_range" {
  let gr = grg_granges(["chr1"], [1], [100])
  let sl : Map[String, Int] = Map([], capacity=1)
  sl["chr1"] = 100
  let gaps = grg_gaps(gr, sl)
  assert_eq(gaps.seqnames.length(), 0)
}

///|
test "grg_gaps_finds_gaps_between_ranges" {
  let gr = grg_granges(
    ["chr1", "chr1"],
    [1, 30],
    [10, 50],
  )
  let sl : Map[String, Int] = Map([], capacity=1)
  sl["chr1"] = 100
  let gaps = grg_gaps(gr, sl)
  assert_eq(gaps.seqnames.length(), 2)
  assert_eq(gaps.seqnames[0], "chr1")
  assert_eq(gaps.starts[0], 11)
  assert_eq(gaps.ends[0], 29)
  assert_eq(gaps.gap_lengths[0], 19)
  assert_eq(gaps.starts[1], 51)
  assert_eq(gaps.ends[1], 100)
  assert_eq(gaps.gap_lengths[1], 50)
}

///|
test "grg_gaps_empty_input" {
  let gr = grg_empty()
  let sl : Map[String, Int] = Map([], capacity=1)
  sl["chr1"] = 100
  let gaps = grg_gaps(gr, sl)
  assert_eq(gaps.seqnames.length(), 0)
}

///|
test "grg_gaps_overlapping_ranges_merged" {
  let gr = grg_granges(
    ["chr1", "chr1", "chr1"],
    [1, 5, 20],
    [10, 15, 30],
  )
  let sl : Map[String, Int] = Map([], capacity=1)
  sl["chr1"] = 50
  let gaps = grg_gaps(gr, sl)
  assert_eq(gaps.seqnames.length(), 2)
  assert_eq(gaps.starts[0], 16)
  assert_eq(gaps.ends[0], 19)
  assert_eq(gaps.starts[1], 31)
  assert_eq(gaps.ends[1], 50)
}

///|
test "grg_reduce_merges_overlapping" {
  let gr = grg_granges(
    ["chr1", "chr1"],
    [1, 5],
    [10, 15],
  )
  let result = grg_reduce(gr)
  assert_eq(result.seqnames.length(), 1)
  assert_eq(result.starts[0], 1)
  assert_eq(result.ends[0], 15)
  assert_eq(result.merged_count[0], 2)
}

///|
test "grg_reduce_merges_adjacent" {
  let gr = grg_granges(
    ["chr1", "chr1"],
    [1, 11],
    [10, 20],
  )
  let result = grg_reduce(gr)
  assert_eq(result.seqnames.length(), 1)
  assert_eq(result.starts[0], 1)
  assert_eq(result.ends[0], 20)
}

///|
test "grg_reduce_no_overlap" {
  let gr = grg_granges(
    ["chr1", "chr1"],
    [1, 30],
    [10, 50],
  )
  let result = grg_reduce(gr)
  assert_eq(result.seqnames.length(), 2)
  assert_eq(result.starts[0], 1)
  assert_eq(result.ends[0], 10)
  assert_eq(result.starts[1], 30)
  assert_eq(result.ends[1], 50)
}

///|
test "grg_reduce_empty_input" {
  let gr = grg_empty()
  let result = grg_reduce(gr)
  assert_eq(result.seqnames.length(), 0)
}

///|
test "grg_disjoin_splits_overlapping" {
  let gr = grg_granges(
    ["chr1", "chr1"],
    [1, 5],
    [10, 15],
  )
  let result = grg_disjoin(gr)
  assert_eq(result.seqnames.length(), 3)
  assert_eq(result.starts[0], 1)
  assert_eq(result.ends[0], 4)
  assert_eq(result.starts[1], 5)
  assert_eq(result.ends[1], 10)
  assert_eq(result.starts[2], 11)
  assert_eq(result.ends[2], 15)
}

///|
test "grg_disjoin_no_overlap" {
  let gr = grg_granges(
    ["chr1", "chr1"],
    [1, 20],
    [10, 30],
  )
  let result = grg_disjoin(gr)
  assert_eq(result.seqnames.length(), 3)
  assert_eq(result.starts[0], 1)
  assert_eq(result.ends[0], 10)
  assert_eq(result.starts[1], 11)
  assert_eq(result.ends[1], 19)
  assert_eq(result.starts[2], 20)
  assert_eq(result.ends[2], 30)
}

///|
test "grg_disjoin_empty_input" {
  let gr = grg_empty()
  let result = grg_disjoin(gr)
  assert_eq(result.seqnames.length(), 0)
}

///|
test "grg_setdiff_removes_overlap" {
  let gr1 = grg_granges(["chr1"], [1], [30])
  let gr2 = grg_granges(["chr1"], [10], [20])
  let result = grg_setdiff(gr1, gr2)
  assert_eq(result.seqnames.length(), 2)
  assert_eq(result.starts[0], 1)
  assert_eq(result.ends[0], 9)
  assert_eq(result.starts[1], 21)
  assert_eq(result.ends[1], 30)
}

///|
test "grg_setdiff_no_overlap" {
  let gr1 = grg_granges(["chr1"], [1], [10])
  let gr2 = grg_granges(["chr1"], [20], [30])
  let result = grg_setdiff(gr1, gr2)
  assert_eq(result.seqnames.length(), 1)
  assert_eq(result.starts[0], 1)
  assert_eq(result.ends[0], 10)
}

///|
test "grg_setdiff_full_overlap" {
  let gr1 = grg_granges(["chr1"], [1], [10])
  let gr2 = grg_granges(["chr1"], [1], [20])
  let result = grg_setdiff(gr1, gr2)
  assert_eq(result.seqnames.length(), 0)
}

///|
test "grg_intersect_finds_overlap" {
  let gr1 = grg_granges(["chr1"], [1], [30])
  let gr2 = grg_granges(["chr1"], [10], [20])
  let result = grg_intersect(gr1, gr2)
  assert_eq(result.seqnames.length(), 1)
  assert_eq(result.starts[0], 10)
  assert_eq(result.ends[0], 20)
}

///|
test "grg_intersect_no_overlap" {
  let gr1 = grg_granges(["chr1"], [1], [10])
  let gr2 = grg_granges(["chr1"], [20], [30])
  let result = grg_intersect(gr1, gr2)
  assert_eq(result.seqnames.length(), 0)
}

///|
test "grg_union_merges_ranges" {
  let gr1 = grg_granges(["chr1"], [1], [10])
  let gr2 = grg_granges(["chr1"], [5], [20])
  let result = grg_union(gr1, gr2)
  assert_eq(result.seqnames.length(), 1)
  assert_eq(result.starts[0], 1)
  assert_eq(result.ends[0], 20)
}

///|
test "grg_union_disjoint_ranges" {
  let gr1 = grg_granges(["chr1"], [1], [10])
  let gr2 = grg_granges(["chr1"], [20], [30])
  let result = grg_union(gr1, gr2)
  assert_eq(result.seqnames.length(), 2)
}

///|
test "grg_compute_coverage_basic" {
  let gr = grg_granges(
    ["chr1", "chr1"],
    [1, 3],
    [5, 8],
  )
  let sl : Map[String, Int] = Map([], capacity=1)
  sl["chr1"] = 10
  let cov = grg_compute_coverage(gr, sl)
  let arr = cov["chr1"]
  assert_eq(arr.length(), 10)
  assert_eq(arr[0], 1)
  assert_eq(arr[2], 2)
  assert_eq(arr[5], 1)
  assert_eq(arr[9], 0)
}

///|
test "grg_compute_coverage_no_ranges" {
  let gr = grg_empty()
  let sl : Map[String, Int] = Map([], capacity=1)
  sl["chr1"] = 5
  let cov = grg_compute_coverage(gr, sl)
  let arr = cov["chr1"]
  assert_eq(arr.length(), 5)
  assert_eq(arr[0], 0)
  assert_eq(arr[4], 0)
}

///|
test "grg_trim_clips_beyond_chromosome" {
  let gr = grg_granges(
    ["chr1", "chr1"],
    [-5, 8],
    [10, 30],
  )
  let sl : Map[String, Int] = Map([], capacity=1)
  sl["chr1"] = 20
  let result = grg_trim(gr, sl)
  assert_eq(result.seqnames.length(), 2)
  assert_eq(result.starts[0], 1)
  assert_eq(result.ends[0], 10)
  assert_eq(result.starts[1], 8)
  assert_eq(result.ends[1], 20)
}

///|
test "grg_trim_removes_completely_outside" {
  let gr = grg_granges(
    ["chr1", "chr1"],
    [1, 25],
    [10, 30],
  )
  let sl : Map[String, Int] = Map([], capacity=1)
  sl["chr1"] = 20
  let result = grg_trim(gr, sl)
  assert_eq(result.seqnames.length(), 1)
  assert_eq(result.starts[0], 1)
  assert_eq(result.ends[0], 10)
}

///|
test "grg_promoters_plus_strand" {
  let gr = grg_granges(
    ["chr1"],
    [100],
    [200],
    strands=["+"],
  )
  let result = grg_promoters(gr, upstream=50, downstream=20)
  assert_eq(result.seqnames.length(), 1)
  assert_eq(result.starts[0], 50)
  assert_eq(result.ends[0], 119)
}

///|
test "grg_promoters_minus_strand" {
  let gr = grg_granges(
    ["chr1"],
    [100],
    [200],
    strands=["-"],
  )
  let result = grg_promoters(gr, upstream=50, downstream=20)
  assert_eq(result.seqnames.length(), 1)
  assert_eq(result.starts[0], 181)
  assert_eq(result.ends[0], 250)
}

///|
test "grg_promoters_star_strand_two_outputs" {
  let gr = grg_granges(
    ["chr1"],
    [100],
    [200],
    strands=["*"],
  )
  let result = grg_promoters(gr, upstream=50, downstream=20)
  assert_eq(result.seqnames.length(), 2)
  assert_eq(result.strands[0], "+")
  assert_eq(result.strands[1], "-")
}

///|
test "grg_setdiff_different_chromosomes" {
  let gr1 = grg_granges(["chr1", "chr2"], [1, 1], [10, 10])
  let gr2 = grg_granges(["chr1"], [5], [15])
  let result = grg_setdiff(gr1, gr2)
  assert_eq(result.seqnames.length(), 2)
  assert_eq(result.seqnames[0], "chr1")
  assert_eq(result.starts[0], 1)
  assert_eq(result.ends[0], 4)
  assert_eq(result.seqnames[1], "chr2")
  assert_eq(result.starts[1], 1)
  assert_eq(result.ends[1], 10)
}

///|
test "grg_gaps_multiple_chromosomes" {
  let gr = grg_granges(
    ["chr1", "chr2"],
    [1, 1],
    [10, 10],
  )
  let sl : Map[String, Int] = Map([], capacity=2)
  sl["chr1"] = 20
  sl["chr2"] = 30
  let gaps = grg_gaps(gr, sl)
  assert_eq(gaps.seqnames.length(), 2)
  assert_eq(gaps.seqnames[0], "chr1")
  assert_eq(gaps.starts[0], 11)
  assert_eq(gaps.ends[0], 20)
  assert_eq(gaps.seqnames[1], "chr2")
  assert_eq(gaps.starts[1], 11)
  assert_eq(gaps.ends[1], 30)
}
