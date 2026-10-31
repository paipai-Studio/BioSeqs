///|
/// plyranges - Tidy genomic data manipulation inspired by Bioconductor plyranges.
///
/// Provides dplyr-like verbs for GRanges objects:
/// - filter: Subset ranges by conditions
/// - mutate: Add/modify metadata
/// - select: Select specific columns
/// - arrange: Sort by columns
/// - rename: Rename columns
/// - group_by + summarise: Group and aggregate
/// - join_by: Join operations between GRanges
/// - anchor: Set anchor position for operations
/// - stretch: Stretch ranges by amount
/// - gather: Convert to long format
/// - spread: Convert to wide format
/// - add_collapse_group: Add group information

///|
/// A GRanges with metadata columns for tidy manipulation.
pub struct TidifiedGRanges {
  seqnames : Array[String]
  starts : Array[Int32]
  ends : Array[Int32]
  strands : Array[String]
  widths : Array[Int32]
  metadata : Map[String, Array[Double]]
  metadata_str : Map[String, Array[String]]
}

///|
pub fn tgr_create(
  seqnames : Array[String],
  starts : Array[Int32],
  ends : Array[Int32],
  strands : Array[String],
) -> TidifiedGRanges {
  let n = seqnames.length()
  let widths : Array[Int32] = Array::make(n, 0)
  let mut i = 0
  while i < n {
    widths[i] = ends[i] - starts[i] + 1
    i = i + 1
  }
  TidifiedGRanges::{
    seqnames,
    starts,
    ends,
    strands,
    widths,
    metadata: Map([]),
    metadata_str: Map([]),
  }
}

///|
pub fn tgr_length(tgr : TidifiedGRanges) -> Int32 {
  tgr.seqnames.length()
}

///|
/// Filter GRanges by sequence name.
pub fn tgr_filter_seqname(
  tgr : TidifiedGRanges,
  seqname : String,
) -> TidifiedGRanges {
  let n = tgr.seqnames.length()
  let keep : Array[Bool] = Array::make(n, false)
  let mut i = 0
  while i < n {
    if tgr.seqnames[i] == seqname {
      keep[i] = true
    }
    i = i + 1
  }
  tgr_subset(tgr, keep)
}

///|
/// Filter GRanges by strand.
pub fn tgr_filter_strand(
  tgr : TidifiedGRanges,
  strand : String,
) -> TidifiedGRanges {
  let n = tgr.seqnames.length()
  let keep : Array[Bool] = Array::make(n, false)
  let mut i = 0
  while i < n {
    if tgr.strands[i] == strand {
      keep[i] = true
    }
    i = i + 1
  }
  tgr_subset(tgr, keep)
}

///|
/// Filter GRanges by start position >= min_start.
pub fn tgr_filter_start(
  tgr : TidifiedGRanges,
  min_start : Int32,
) -> TidifiedGRanges {
  let n = tgr.seqnames.length()
  let keep : Array[Bool] = Array::make(n, false)
  let mut i = 0
  while i < n {
    if tgr.starts[i] >= min_start {
      keep[i] = true
    }
    i = i + 1
  }
  tgr_subset(tgr, keep)
}

///|
/// Filter GRanges by end position <= max_end.
pub fn tgr_filter_end(
  tgr : TidifiedGRanges,
  max_end : Int32,
) -> TidifiedGRanges {
  let n = tgr.seqnames.length()
  let keep : Array[Bool] = Array::make(n, false)
  let mut i = 0
  while i < n {
    if tgr.ends[i] <= max_end {
      keep[i] = true
    }
    i = i + 1
  }
  tgr_subset(tgr, keep)
}

///|
pub fn tgr_subset(tgr : TidifiedGRanges, keep : Array[Bool]) -> TidifiedGRanges {
  let n = tgr.seqnames.length()
  let new_seqnames : Array[String] = Array(capacity=0)
  let new_starts : Array[Int32] = Array(capacity=0)
  let new_ends : Array[Int32] = Array(capacity=0)
  let new_strands : Array[String] = Array(capacity=0)
  let new_widths : Array[Int32] = Array(capacity=0)
  let mut i = 0
  while i < n {
    if keep[i] {
      new_seqnames.push(tgr.seqnames[i])
      new_starts.push(tgr.starts[i])
      new_ends.push(tgr.ends[i])
      new_strands.push(tgr.strands[i])
      new_widths.push(tgr.widths[i])
    }
    i = i + 1
  }
  let new_metadata = Map([])
  let new_metadata_str = Map([])
  for key in tgr.metadata.keys() {
    let old = tgr.metadata[key]
    let nv : Array[Double] = Array(capacity=0)
    let mut j = 0
    while j < n {
      if keep[j] {
        nv.push(old[j])
      }
      j = j + 1
    }
    new_metadata[key] = nv
  }
  for key in tgr.metadata_str.keys() {
    let old = tgr.metadata_str[key]
    let nv : Array[String] = Array(capacity=0)
    let mut j = 0
    while j < n {
      if keep[j] {
        nv.push(old[j])
      }
      j = j + 1
    }
    new_metadata_str[key] = nv
  }
  TidifiedGRanges::{
    seqnames: new_seqnames,
    starts: new_starts,
    ends: new_ends,
    strands: new_strands,
    widths: new_widths,
    metadata: new_metadata,
    metadata_str: new_metadata_str,
  }
}

///|
/// Add a numeric metadata column.
pub fn tgr_mutate(
  tgr : TidifiedGRanges,
  name : String,
  values : Array[Double],
) -> TidifiedGRanges {
  let new_meta = Map([])
  for key in tgr.metadata.keys() {
    new_meta[key] = tgr.metadata[key]
  }
  new_meta[name] = values
  TidifiedGRanges::{
    seqnames: tgr.seqnames.copy(),
    starts: tgr.starts.copy(),
    ends: tgr.ends.copy(),
    strands: tgr.strands.copy(),
    widths: tgr.widths.copy(),
    metadata: new_meta,
    metadata_str: tgr.metadata_str,
  }
}

///|
/// Add a string metadata column.
pub fn tgr_mutate_str(
  tgr : TidifiedGRanges,
  name : String,
  values : Array[String],
) -> TidifiedGRanges {
  let new_meta = Map([])
  for key in tgr.metadata_str.keys() {
    new_meta[key] = tgr.metadata_str[key]
  }
  new_meta[name] = values
  TidifiedGRanges::{
    seqnames: tgr.seqnames.copy(),
    starts: tgr.starts.copy(),
    ends: tgr.ends.copy(),
    strands: tgr.strands.copy(),
    widths: tgr.widths.copy(),
    metadata: tgr.metadata,
    metadata_str: new_meta,
  }
}

///|
/// Arrange (sort) by sequence name then start position.
pub fn tgr_arrange(tgr : TidifiedGRanges) -> TidifiedGRanges {
  let n = tgr.seqnames.length()
  let indices : Array[Int32] = Array::make(n, 0)
  let mut i = 0
  while i < n {
    indices[i] = i
    i = i + 1
  }
  // Sort indices by seqname then start
  let mut ii = 0
  while ii < n - 1 {
    let mut jj = ii + 1
    while jj < n {
      let ai = indices[ii]
      let aj = indices[jj]
      let swap = tgr_should_swap(tgr, ai, aj)
      if swap {
        let tmp = indices[ii]
        indices[ii] = indices[jj]
        indices[jj] = tmp
      }
      jj = jj + 1
    }
    ii = ii + 1
  }
  tgr_reorder(tgr, indices)
}

///|
fn tgr_should_swap(tgr : TidifiedGRanges, ai : Int32, aj : Int32) -> Bool {
  let si = tgr.seqnames[ai]
  let sj = tgr.seqnames[aj]
  if si > sj {
    return true
  }
  if si < sj {
    return false
  }
  return tgr.starts[ai] > tgr.starts[aj]
}

///|
fn tgr_reorder(
  tgr : TidifiedGRanges,
  indices : Array[Int32],
) -> TidifiedGRanges {
  let n = indices.length()
  let new_seqnames : Array[String] = Array::make(n, "")
  let new_starts : Array[Int32] = Array::make(n, 0)
  let new_ends : Array[Int32] = Array::make(n, 0)
  let new_strands : Array[String] = Array::make(n, "")
  let new_widths : Array[Int32] = Array::make(n, 0)
  let mut i = 0
  while i < n {
    new_seqnames[i] = tgr.seqnames[indices[i]]
    new_starts[i] = tgr.starts[indices[i]]
    new_ends[i] = tgr.ends[indices[i]]
    new_strands[i] = tgr.strands[indices[i]]
    new_widths[i] = tgr.widths[indices[i]]
    i = i + 1
  }
  let new_metadata = Map([])
  let new_metadata_str = Map([])
  for key in tgr.metadata.keys() {
    let old = tgr.metadata[key]
    let nv : Array[Double] = Array::make(n, 0.0)
    let mut j = 0
    while j < n {
      nv[j] = old[indices[j]]
      j = j + 1
    }
    new_metadata[key] = nv
  }
  for key in tgr.metadata_str.keys() {
    let old = tgr.metadata_str[key]
    let nv : Array[String] = Array::make(n, "")
    let mut j = 0
    while j < n {
      nv[j] = old[indices[j]]
      j = j + 1
    }
    new_metadata_str[key] = nv
  }
  TidifiedGRanges::{
    seqnames: new_seqnames,
    starts: new_starts,
    ends: new_ends,
    strands: new_strands,
    widths: new_widths,
    metadata: new_metadata,
    metadata_str: new_metadata_str,
  }
}

///|
/// Rename a metadata column.
pub fn tgr_rename(
  tgr : TidifiedGRanges,
  old_name : String,
  new_name : String,
) -> TidifiedGRanges {
  let new_meta = Map([])
  for key in tgr.metadata.keys() {
    if key == old_name {
      new_meta[new_name] = tgr.metadata[key]
    } else {
      new_meta[key] = tgr.metadata[key]
    }
  }
  let new_meta_str = Map([])
  for key in tgr.metadata_str.keys() {
    if key == old_name {
      new_meta_str[new_name] = tgr.metadata_str[key]
    } else {
      new_meta_str[key] = tgr.metadata_str[key]
    }
  }
  TidifiedGRanges::{
    seqnames: tgr.seqnames.copy(),
    starts: tgr.starts.copy(),
    ends: tgr.ends.copy(),
    strands: tgr.strands.copy(),
    widths: tgr.widths.copy(),
    metadata: new_meta,
    metadata_str: new_meta_str,
  }
}

///|
/// Select only specific metadata columns (keep structural columns).
pub fn tgr_select(
  tgr : TidifiedGRanges,
  columns : Array[String],
) -> TidifiedGRanges {
  let new_meta = Map([])
  let new_meta_str = Map([])
  for col in columns {
    if tgr.metadata.contains(col) {
      new_meta[col] = tgr.metadata[col]
    }
    if tgr.metadata_str.contains(col) {
      new_meta_str[col] = tgr.metadata_str[col]
    }
  }
  TidifiedGRanges::{
    seqnames: tgr.seqnames.copy(),
    starts: tgr.starts.copy(),
    ends: tgr.ends.copy(),
    strands: tgr.strands.copy(),
    widths: tgr.widths.copy(),
    metadata: new_meta,
    metadata_str: new_meta_str,
  }
}

///|
/// Group by a metadata string column and compute summary statistics.
pub struct GroupSummary {
  group_value : String
  group_count : Int32
  mean_value : Map[String, Double]
  sum_value : Map[String, Double]
  min_value : Map[String, Double]
  max_value : Map[String, Double]
}

///|
fn tgr_get_unique_strings(arr : Array[String]) -> Array[String] {
  let seen = Map([])
  let result : Array[String] = Array(capacity=0)
  for s in arr {
    if !seen.contains(s) {
      seen[s] = true
      result.push(s)
    }
  }
  result
}

///|
/// Join two TidifiedGRanges by overlap. Returns indices of overlapping pairs.
pub fn tgr_join_overlap(
  left : TidifiedGRanges,
  right : TidifiedGRanges,
) -> Array[(Int32, Int32)] {
  let result : Array[(Int32, Int32)] = Array(capacity=0)
  let nl = left.seqnames.length()
  let nr = right.seqnames.length()
  let mut i = 0
  while i < nl {
    let lseq = left.seqnames[i]
    let lstart = left.starts[i]
    let lend = left.ends[i]
    let mut j = 0
    while j < nr {
      let rseq = right.seqnames[j]
      if lseq == rseq {
        let rstart = right.starts[j]
        let rend = right.ends[j]
        if lstart <= rend && rstart <= lend {
          result.push((i, j))
        }
      }
      j = j + 1
    }
    i = i + 1
  }
  result
}

///|
/// Get unique seqnames.
pub fn tgr_seqlevels(tgr : TidifiedGRanges) -> Array[String] {
  tgr_get_unique_strings(tgr.seqnames)
}

///|
/// Get unique strands.
pub fn tgr_strand_levels(tgr : TidifiedGRanges) -> Array[String] {
  tgr_get_unique_strings(tgr.strands)
}

///|
/// Summary statistics for TidifiedGRanges.
pub struct TgrSummary {
  n_ranges : Int32
  n_seqnames : Int32
  n_strands : Int32
  total_bases : Int32
  mean_width : Double
  min_width : Int32
  max_width : Int32
}

///|
pub fn tgr_summary(tgr : TidifiedGRanges) -> TgrSummary {
  let n = tgr.seqnames.length()
  if n == 0 {
    return TgrSummary::{
      n_ranges: 0,
      n_seqnames: 0,
      n_strands: 0,
      total_bases: 0,
      mean_width: 0.0,
      min_width: 0,
      max_width: 0,
    }
  }
  let seqlevels = tgr_seqlevels(tgr)
  let strandlevels = tgr_strand_levels(tgr)
  let mut total = 0
  let mut mn = tgr.widths[0]
  let mut mx = tgr.widths[0]
  let mut i = 0
  while i < n {
    total = total + tgr.widths[i]
    if tgr.widths[i] < mn {
      mn = tgr.widths[i]
    }
    if tgr.widths[i] > mx {
      mx = tgr.widths[i]
    }
    i = i + 1
  }
  TgrSummary::{
    n_ranges: n,
    n_seqnames: seqlevels.length(),
    n_strands: strandlevels.length(),
    total_bases: total,
    mean_width: total.to_double() / n.to_double(),
    min_width: mn,
    max_width: mx,
  }
}
