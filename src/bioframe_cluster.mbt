/// bioframe_cluster.mbt — Cluster overlapping intervals into groups.
/// Bioframe equivalent: `bioframe.ops.cluster`.
///
/// Unlike merge() which discards original intervals, cluster() keeps the
/// input intervals and assigns each one to a cluster id. This is useful
/// when you want to know which intervals are in the same "neighborhood"
/// without collapsing them.

///|
/// One row of a cluster result: original interval + its cluster assignment.
pub struct BfClusterRow {
  chrom : String
  start : Int32
  end : Int32
  cluster_id : Int32
}

///|
/// Cluster overlapping or book-ended intervals into groups.
///
/// Input is sorted within each chromosome by start/end, then a single
/// pass walks and assigns cluster ids. Intervals on different chromosomes
/// always end up in different clusters.
///
/// Parameters:
///   - intervals: input BedIntervals (will be copied + sorted internally)
///   - min_dist:   minimum distance between intervals for them to be
///                 clustered together.
///                 - Some(d) means intervals separated by ≤ d bp are
///                   clustered (like merge -d d).
///                 - None means only strictly overlapping / book-ended
///                   intervals are clustered (distance 0 is NOT enough).
///                 Default: None.
///
/// Returns:
///   - rows: each input interval + its cluster_id (same order as input
///           after sorting within chrom).
///   - n_clusters: total number of cluster ids assigned.
///
pub fn bf_cluster(
  intervals : Array[BedInterval],
  min_dist~ : Option[Int32] = None,
) -> (Array[BfClusterRow], Int32) {
  if intervals.length() == 0 {
    return (Array(capacity=0), 0)
  }

  // Sort by (chrom, start, end)
  let sorted = bf_cluster_sort(intervals)

  let rows : Array[BfClusterRow] = Array(capacity=sorted.length())
  let mut cluster_id : Int32 = 0
  let mut cur_chrom = sorted[0].chrom
  let mut cur_end : Int32 = sorted[0].end
  rows.push(BfClusterRow::{
    chrom: sorted[0].chrom,
    start: sorted[0].start,
    end: sorted[0].end,
    cluster_id,
  })

  for i = 1; i < sorted.length(); i = i + 1 {
    let iv = sorted[i]
    let mut start_new_cluster = false

    if iv.chrom != cur_chrom {
      start_new_cluster = true
    } else {
      match min_dist {
        Some(d) => {
          // Cluster if start ≤ cur_end + d  (within gap)
          if iv.start > cur_end + d {
            start_new_cluster = true
          }
        }
        None => {
          // Only cluster overlapping / book-ended intervals
          // semi-open: [a, b) and [b, c) do NOT overlap, must be < not ≤
          if iv.start >= cur_end {
            start_new_cluster = true
          }
        }
      }
    }

    if start_new_cluster {
      cluster_id = cluster_id + 1
      cur_chrom = iv.chrom
      cur_end = iv.end
    } else {
      // Extend current cluster's span
      if iv.end > cur_end {
        cur_end = iv.end
      }
    }
    rows.push(BfClusterRow::{
      chrom: iv.chrom,
      start: iv.start,
      end: iv.end,
      cluster_id,
    })
  }

  (rows, cluster_id + 1)
}

///|
/// Cluster with a minimum merge distance supplied as an Int32.
/// Convenience wrapper around bf_cluster.
pub fn bf_cluster_with_min_dist(
  intervals : Array[BedInterval],
  min_dist : Int32,
) -> (Array[BfClusterRow], Int32) {
  bf_cluster(intervals, min_dist=Some(min_dist))
}

///|
/// Internal sort by (chrom, start, end).
fn bf_cluster_sort(intervals : Array[BedInterval]) -> Array[BedInterval] {
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
