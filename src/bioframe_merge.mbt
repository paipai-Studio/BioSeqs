/// bioframe_merge.mbt — DataFrame-style merge and close (gap-tolerant merge)
/// Items 1, 4: bf_merge, bf_close

///|
/// Result of a merge operation: merged interval + cluster_id of source intervals.
pub struct BfMergeResult {
  chrom : String
  start : Int32
  end : Int32
  cluster_id : Int32
  n_merged : Int32 // how many source intervals were merged into this one
}

///|
/// Merge overlapping or book-ended intervals.
/// Returns merged intervals with cluster IDs.
/// If gap > 0, intervals within `gap` bp are also merged (like bedtools merge -d).
pub fn bf_merge(
  intervals : Array[BedInterval],
  gap~ : Int32 = 0,
) -> Array[BfMergeResult] {
  if intervals.length() == 0 {
    return Array(capacity=0)
  }
  // Sort by (chrom, start, end)
  let sorted = bf_merge_sort(intervals)
  let results : Array[BfMergeResult] = Array(capacity=0)
  let mut cur_chrom = sorted[0].chrom
  let mut cur_start = sorted[0].start
  let mut cur_end = sorted[0].end
  let mut cluster_id = 0
  let mut n_merged = 1
  for i = 1; i < sorted.length(); i = i + 1 {
    let iv = sorted[i]
    if iv.chrom == cur_chrom && iv.start <= cur_end + gap {
      // Overlapping or within gap — extend
      if iv.end > cur_end {
        cur_end = iv.end
      }
      n_merged = n_merged + 1
    } else {
      // Flush current
      results.push(BfMergeResult::{
        chrom: cur_chrom,
        start: cur_start,
        end: cur_end,
        cluster_id,
        n_merged,
      })
      cluster_id = cluster_id + 1
      cur_chrom = iv.chrom
      cur_start = iv.start
      cur_end = iv.end
      n_merged = 1
    }
  }
  // Flush last
  results.push(BfMergeResult::{
    chrom: cur_chrom,
    start: cur_start,
    end: cur_end,
    cluster_id,
    n_merged,
  })
  results
}

///|
/// Close: merge intervals within a specified gap distance.
/// This is a convenience wrapper around bf_merge with gap.
/// Unlike bf_merge(gap=0) which only merges overlapping/book-ended,
/// bf_close(gap=d) merges intervals that are up to d bp apart.
pub fn bf_close(
  intervals : Array[BedInterval],
  gap~ : Int32 = 0,
) -> Array[BfMergeResult] {
  bf_merge(intervals, gap~)
}

///|
/// Internal sort by (chrom, start, end).
fn bf_merge_sort(intervals : Array[BedInterval]) -> Array[BedInterval] {
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
