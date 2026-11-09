/// Viewframe operations.
///
/// Bioframe equivalent: `make_viewframe`, `assign_view`, `select_view`
///
/// A "view" is a named collection of genomic regions that defines the
/// scope of operations. Views are used to restrict interval operations
/// to specific chromosomes or regions (e.g., only chr1, or only gene bodies).

///|
/// A view region: named genomic region
pub struct BfViewRegion {
  name : String
  chrom : String
  start : Int32
  end : Int32
}

///|
pub fn BfViewRegion::name(self : BfViewRegion) -> String {
  self.name
}

///|
pub fn BfViewRegion::chrom(self : BfViewRegion) -> String {
  self.chrom
}

///|
pub fn BfViewRegion::start(self : BfViewRegion) -> Int32 {
  self.start
}

///|
pub fn BfViewRegion::end(self : BfViewRegion) -> Int32 {
  self.end
}

///|
/// A viewframe: collection of view regions
pub struct BfViewFrame {
  regions : Array[BfViewRegion]
}

///|
pub fn BfViewFrame::regions(self : BfViewFrame) -> Array[BfViewRegion] {
  self.regions
}

///|
pub fn BfViewFrame::length(self : BfViewFrame) -> Int32 {
  self.regions.length()
}

///|
/// Make a viewframe from chromsizes.
/// Each chromosome becomes a view region spanning the entire chromosome.
pub fn bf_make_viewframe(chromsizes : Array[(String, Int32)]) -> BfViewFrame {
  let regions : Array[BfViewRegion] = Array(capacity=0)
  for cs in chromsizes {
    regions.push({ name: cs.0, chrom: cs.0, start: 0, end: cs.1, })
  }
  { regions, }
}

///|
/// Make a viewframe from named regions.
pub fn bf_make_viewframe_regions(
  regions : Array[(String, String, Int32, Int32)],
) -> BfViewFrame {
  let out : Array[BfViewRegion] = Array(capacity=0)
  for r in regions {
    out.push({ name: r.0, chrom: r.1, start: r.2, end: r.3, })
  }
  { regions: out, }
}

///|
/// Assign each interval to the view region that contains it.
/// Returns array of (interval_index, view_name) pairs.
/// Intervals not contained in any view are assigned "NA".
pub fn bf_assign_view(
  intervals : Array[BedInterval],
  view : BfViewFrame,
) -> Array[(Int32, String)] {
  let result : Array[(Int32, String)] = Array(capacity=0)
  for i = 0; i < intervals.length(); i = i + 1 {
    let iv = intervals[i]
    let mut assigned = "NA"
    for vr in view.regions {
      if iv.chrom == vr.chrom && vr.start <= iv.start && iv.end <= vr.end {
        assigned = vr.name
        break
      }
    }
    result.push((i, assigned))
  }
  result
}

///|
/// Select only intervals that fall within any view region.
pub fn bf_select_view(
  intervals : Array[BedInterval],
  view : BfViewFrame,
) -> Array[BedInterval] {
  let result : Array[BedInterval] = Array(capacity=0)
  for iv in intervals {
    for vr in view.regions {
      if iv.chrom == vr.chrom && vr.start <= iv.start && iv.end <= vr.end {
        result.push(iv)
        break
      }
    }
  }
  result
}

///|
/// Select intervals by view name.
pub fn bf_select_by_view_name(
  intervals : Array[BedInterval],
  view : BfViewFrame,
  name : String,
) -> Array[BedInterval] {
  let result : Array[BedInterval] = Array(capacity=0)
  for vr in view.regions {
    if vr.name == name {
      for iv in intervals {
        if iv.chrom == vr.chrom && vr.start <= iv.start && iv.end <= vr.end {
          result.push(iv)
        }
      }
    }
  }
  result
}

///|
/// Get the chromsizes represented by a viewframe.
pub fn bf_viewframe_chromsizes(view : BfViewFrame) -> Array[(String, Int32)] {
  let result : Array[(String, Int32)] = Array(capacity=0)
  let seen : Map[String, Int32] = Map([])
  for vr in view.regions {
    match seen.get(vr.chrom) {
      Some(prev) => if vr.end > prev { seen[vr.chrom] = vr.end }
      None => seen[vr.chrom] = vr.end
    }
  }
  let keys = seen.keys().to_array()
  for k in keys {
    result.push((k, seen[k]))
  }
  result
}
