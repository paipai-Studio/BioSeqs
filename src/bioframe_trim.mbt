/// Trim and squash intervals to view bounds.
///
/// Bioframe equivalent: `trim_to_view`, `squash_to_view`
///
/// trim: clip intervals that extend beyond view boundaries
/// squash: remove intervals outside views, clip those that straddle

///|
/// Trim intervals to fit within view regions.
/// Intervals partially outside a view are clipped to the view boundary.
/// Intervals fully outside all views are dropped.
pub fn bf_trim_to_view(
  intervals : Array[BedInterval],
  view : BfViewFrame,
) -> Array[BedInterval] {
  let result : Array[BedInterval] = Array(capacity=0)
  for iv in intervals {
    for vr in view.regions {
      if iv.chrom == vr.chrom {
        // Check if there's any overlap at all
        if iv.start < vr.end && vr.start < iv.end {
          let new_start = if iv.start < vr.start { vr.start } else { iv.start }
          let new_end = if iv.end > vr.end { vr.end } else { iv.end }
          if new_end > new_start {
            result.push(BedInterval::{
              chrom: iv.chrom,
              start: new_start,
              end: new_end,
            })
          }
        }
      }
    }
  }
  result
}

///|
/// Squash intervals to view: same as trim but also shifts coordinates
/// to be relative to the view start (local coordinates).
pub fn bf_squash_to_view(
  intervals : Array[BedInterval],
  view : BfViewFrame,
) -> Array[BedInterval] {
  let result : Array[BedInterval] = Array(capacity=0)
  for iv in intervals {
    for vr in view.regions {
      if iv.chrom == vr.chrom {
        if iv.start < vr.end && vr.start < iv.end {
          let new_start = if iv.start < vr.start {
            0
          } else {
            iv.start - vr.start
          }
          let new_end = if iv.end > vr.end {
            vr.end - vr.start
          } else {
            iv.end - vr.start
          }
          if new_end > new_start {
            result.push(BedInterval::{
              chrom: vr.name,
              start: new_start,
              end: new_end,
            })
          }
        }
      }
    }
  }
  result
}

///|
/// Trim intervals to chromsizes bounds.
/// Intervals extending beyond chromosome length are clipped.
pub fn bf_trim_to_chromsizes(
  intervals : Array[BedInterval],
  chromsizes : Array[(String, Int32)],
) -> Array[BedInterval] {
  let bound_map : Map[String, Int32] = Map([])
  for cs in chromsizes {
    bound_map[cs.0] = cs.1
  }
  let result : Array[BedInterval] = Array(capacity=0)
  for iv in intervals {
    match bound_map.get(iv.chrom) {
      Some(length) => {
        let new_start = if iv.start < 0 { 0 } else { iv.start }
        let new_end = if iv.end > length { length } else { iv.end }
        if new_end > new_start {
          result.push(BedInterval::{
            chrom: iv.chrom,
            start: new_start,
            end: new_end,
          })
        }
      }
      None => ()
    }
  }
  result
}

///|
/// Check if an interval is within view bounds.
pub fn bf_is_within_view(iv : BedInterval, view : BfViewFrame) -> Bool {
  for vr in view.regions {
    if iv.chrom == vr.chrom && vr.start <= iv.start && iv.end <= vr.end {
      return true
    }
  }
  false
}
