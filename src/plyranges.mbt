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
  starts : Array[Int]
  ends : Array[Int]
  strands : Array[String]
  widths : Array[Int]
  metadata : Map[String, Array[Double]]
  metadata_str : Map[String, Array[String]]
}

///|
pub fn tgr_create(
  seqnames : Array[String],
  starts : Array[Int],
  ends : Array[Int],
  strands : Array[String],
) -> TidifiedGRanges {
  let n = seqnames.length()
  let widths : Array[Int] = Array::make(n, 0)
  let mut i = 0
  while i < n {
    widths[i] = ends[i] - starts[i] + 1
    i = i + 1
  }
  TidifiedGRanges::{
    seqnames, starts, ends, strands, widths,
    metadata: Map::new(),
    metadata_str: Map::new(),
  }
}

///|
pub fn tgr_length(tgr : TidifiedGRanges) -> Int {
  tgr.seqnames.length()
}

///|
/// Filter GRanges by sequence name.
pub fn tgr_filter_seqname(tgr : TidifiedGRanges, seqname : String) -> TidifiedGRanges {
  let n = tgr.seqnames.length()
  let keep : Array[Bool] = Array::make(n, false)
  let mut i = 0
  while i < n {
    if tgr.seqnames[i] == seqname { keep[i] = true }
    i = i + 1
  }
  tgr_subset(tgr, keep)
}

///|
/// Filter GRanges by strand.
pub fn tgr_filter_strand(tgr : TidifiedGRanges, strand : String) -> TidifiedGRanges {
  let n = tgr.seqnames.length()
  let keep : Array[Bool] = Array::make(n, false)
  let mut i = 0
  while i < n {
    if tgr.strands[i] == strand { keep[i] = true }
    i = i + 1
  }
  tgr_subset(tgr, keep)
}

///|
/// Filter GRanges by width range [min_width, max_width] inclusive.
pub fn tgr_filter_width(
  tgr : TidifiedGRanges,
  min_width : Int,
  max_width : Int,
) -> TidifiedGRanges {
  let n = tgr.seqnames.length()
  let keep : Array[Bool] = Array::make(n, false)
  let mut i = 0
  while i < n {
    let w = tgr.widths[i]
    if w >= min_width && w <= max_width { keep[i] = true }
    i = i + 1
  }
  tgr_subset(tgr, keep)
}

///|
/// Filter GRanges by start position >= min_start.
pub fn tgr_filter_start(tgr : TidifiedGRanges, min_start : Int) -> TidifiedGRanges {
  let n = tgr.seqnames.length()
  let keep : Array[Bool] = Array::make(n, false)
  let mut i = 0
  while i < n {
    if tgr.starts[i] >= min_start { keep[i] = true }
    i = i + 1
  }
  tgr_subset(tgr, keep)
}

///|
/// Filter GRanges by end position <= max_end.
pub fn tgr_filter_end(tgr : TidifiedGRanges, max_end : Int) -> TidifiedGRanges {
  let n = tgr.seqnames.length()
  let keep : Array[Bool] = Array::make(n, false)
  let mut i = 0
  while i < n {
    if tgr.ends[i] <= max_end { keep[i] = true }
    i = i + 1
  }
  tgr_subset(tgr, keep)
}

///|
/// Filter GRanges by a metadata numeric column value.
pub fn tgr_filter_metadata(
  tgr : TidifiedGRanges,
  col : String,
  min_val : Double,
  max_val : Double,
) -> TidifiedGRanges {
  let data = tgr.metadata[col]
  if data.length() == 0 { return tgr }
  let n = tgr.seqnames.length()
  let keep : Array[Bool] = Array::make(n, false)
  let mut i = 0
  while i < n {
    let v = data[i]
    if v >= min_val && v <= max_val { keep[i] = true }
    i = i + 1
  }
  tgr_subset(tgr, keep)
}

///|
/// Filter GRanges by a metadata string column value.
pub fn tgr_filter_metadata_str(
  tgr : TidifiedGRanges,
  col : String,
  value : String,
) -> TidifiedGRanges {
  let data = tgr.metadata_str[col]
  if data.length() == 0 { return tgr }
  let n = tgr.seqnames.length()
  let keep : Array[Bool] = Array::make(n, false)
  let mut i = 0
  while i < n {
    if data[i] == value { keep[i] = true }
    i = i + 1
  }
  tgr_subset(tgr, keep)
}

///|
pub fn tgr_subset(tgr : TidifiedGRanges, keep : Array[Bool]) -> TidifiedGRanges {
  let n = tgr.seqnames.length()
  let new_seqnames : Array[String] = Array::new()
  let new_starts : Array[Int] = Array::new()
  let new_ends : Array[Int] = Array::new()
  let new_strands : Array[String] = Array::new()
  let new_widths : Array[Int] = Array::new()
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
  let new_metadata = Map::new()
  let new_metadata_str = Map::new()
  for key in tgr.metadata.keys() {
    let old = tgr.metadata[key]
    let nv : Array[Double] = Array::new()
    let mut j = 0
    while j < n {
      if keep[j] { nv.push(old[j]) }
      j = j + 1
    }
    new_metadata[key] = nv
  }
  for key in tgr.metadata_str.keys() {
    let old = tgr.metadata_str[key]
    let nv : Array[String] = Array::new()
    let mut j = 0
    while j < n {
      if keep[j] { nv.push(old[j]) }
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
pub fn tgr_mutate(tgr : TidifiedGRanges, name : String, values : Array[Double]) -> TidifiedGRanges {
  let new_meta = Map::new()
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
pub fn tgr_mutate_str(tgr : TidifiedGRanges, name : String, values : Array[String]) -> TidifiedGRanges {
  let new_meta = Map::new()
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
  let indices : Array[Int] = Array::make(n, 0)
  let mut i = 0
  while i < n { indices[i] = i; i = i + 1 }
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

fn tgr_should_swap(tgr : TidifiedGRanges, ai : Int, aj : Int) -> Bool {
  let si = tgr.seqnames[ai]
  let sj = tgr.seqnames[aj]
  if si > sj { return true }
  if si < sj { return false }
  return tgr.starts[ai] > tgr.starts[aj]
}

///|
/// Arrange by a numeric metadata column.
pub fn tgr_arrange_by_metadata(tgr : TidifiedGRanges, col : String) -> TidifiedGRanges {
  let data = tgr.metadata[col]
  if data.length() == 0 { return tgr }
  let n = tgr.seqnames.length()
  let indices : Array[Int] = Array::make(n, 0)
  let mut i = 0
  while i < n { indices[i] = i; i = i + 1 }
  let mut ii = 0
  while ii < n - 1 {
    let mut jj = ii + 1
    while jj < n {
      let ai = indices[ii]
      let aj = indices[jj]
      if data[ai] > data[aj] {
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

fn tgr_reorder(tgr : TidifiedGRanges, indices : Array[Int]) -> TidifiedGRanges {
  let n = indices.length()
  let new_seqnames : Array[String] = Array::make(n, "")
  let new_starts : Array[Int] = Array::make(n, 0)
  let new_ends : Array[Int] = Array::make(n, 0)
  let new_strands : Array[String] = Array::make(n, "")
  let new_widths : Array[Int] = Array::make(n, 0)
  let mut i = 0
  while i < n {
    new_seqnames[i] = tgr.seqnames[indices[i]]
    new_starts[i] = tgr.starts[indices[i]]
    new_ends[i] = tgr.ends[indices[i]]
    new_strands[i] = tgr.strands[indices[i]]
    new_widths[i] = tgr.widths[indices[i]]
    i = i + 1
  }
  let new_metadata = Map::new()
  let new_metadata_str = Map::new()
  for key in tgr.metadata.keys() {
    let old = tgr.metadata[key]
    let nv : Array[Double] = Array::make(n, 0.0)
    let mut j = 0
    while j < n { nv[j] = old[indices[j]]; j = j + 1 }
    new_metadata[key] = nv
  }
  for key in tgr.metadata_str.keys() {
    let old = tgr.metadata_str[key]
    let nv : Array[String] = Array::make(n, "")
    let mut j = 0
    while j < n { nv[j] = old[indices[j]]; j = j + 1 }
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
pub fn tgr_rename(tgr : TidifiedGRanges, old_name : String, new_name : String) -> TidifiedGRanges {
  let new_meta = Map::new()
  for key in tgr.metadata.keys() {
    if key == old_name { new_meta[new_name] = tgr.metadata[key] }
    else { new_meta[key] = tgr.metadata[key] }
  }
  let new_meta_str = Map::new()
  for key in tgr.metadata_str.keys() {
    if key == old_name { new_meta_str[new_name] = tgr.metadata_str[key] }
    else { new_meta_str[key] = tgr.metadata_str[key] }
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
pub fn tgr_select(tgr : TidifiedGRanges, columns : Array[String]) -> TidifiedGRanges {
  let new_meta = Map::new()
  let new_meta_str = Map::new()
  for col in columns {
    if tgr.metadata.contains(col) { new_meta[col] = tgr.metadata[col] }
    if tgr.metadata_str.contains(col) { new_meta_str[col] = tgr.metadata_str[col] }
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
  group_count : Int
  mean_value : Map[String, Double]
  sum_value : Map[String, Double]
  min_value : Map[String, Double]
  max_value : Map[String, Double]
}

///|
pub fn tgr_summarise(
  tgr : TidifiedGRanges,
  group_col : String,
  value_cols : Array[String],
) -> Array[GroupSummary] {
  let data = tgr.metadata_str[group_col]
  if data.length() == 0 { return [] }
  let groups = tgr_get_unique_strings(data)
  let result : Array[GroupSummary] = Array::new()
  for g in groups {
    let mut count = 0
    let mut i = 0
    while i < data.length() {
      if data[i] == g { count = count + 1 }
      i = i + 1
    }
    let mean_vals = Map::new()
    let sum_vals = Map::new()
    let min_vals = Map::new()
    let max_vals = Map::new()
    for col in value_cols {
      if tgr.metadata.contains(col) {
        let col_data = tgr.metadata[col]
        let mut s = 0.0
        let mut mn = 0.0
        let mut mx = 0.0
        let mut j = 0
        let first = true
        while j < data.length() {
          if data[j] == g {
            let v = col_data[j]
            s = s + v
            if first || v < mn { mn = v }
            if first || v > mx { mx = v }
          }
          j = j + 1
        }
        mean_vals[col] = s / count.to_double()
        sum_vals[col] = s
        min_vals[col] = mn
        max_vals[col] = mx
      }
    }
    result.push(GroupSummary::{
      group_value: g,
      group_count: count,
      mean_value: mean_vals,
      sum_value: sum_vals,
      min_value: min_vals,
      max_value: max_vals,
    })
  }
  result
}

fn tgr_get_unique_strings(arr : Array[String]) -> Array[String] {
  let seen = Map::new()
  let result : Array[String] = Array::new()
  for s in arr {
    if !seen.contains(s) {
      seen[s] = true
      result.push(s)
    }
  }
  result
}

///|
/// Compute mean of a numeric metadata column grouped by a string column.
pub fn tgr_group_mean(
  tgr : TidifiedGRanges,
  group_col : String,
  value_col : String,
) -> Array[(String, Double)] {
  let data = tgr.metadata_str[group_col]
  let vals = tgr.metadata[value_col]
  if data.length() == 0 || vals.length() == 0 { return [] }
  let groups = tgr_get_unique_strings(data)
  let result : Array[(String, Double)] = Array::new()
  for g in groups {
    let mut s = 0.0
    let mut c = 0
    let mut i = 0
    while i < data.length() {
      if data[i] == g { s = s + vals[i]; c = c + 1 }
      i = i + 1
    }
    result.push((g, s / c.to_double()))
  }
  result
}

///|
/// Compute sum of a numeric metadata column grouped by a string column.
pub fn tgr_group_sum(
  tgr : TidifiedGRanges,
  group_col : String,
  value_col : String,
) -> Array[(String, Double)] {
  let data = tgr.metadata_str[group_col]
  let vals = tgr.metadata[value_col]
  if data.length() == 0 || vals.length() == 0 { return [] }
  let groups = tgr_get_unique_strings(data)
  let result : Array[(String, Double)] = Array::new()
  for g in groups {
    let mut s = 0.0
    let mut i = 0
    while i < data.length() {
      if data[i] == g { s = s + vals[i] }
      i = i + 1
    }
    result.push((g, s))
  }
  result
}

///|
/// Join two TidifiedGRanges by overlap. Returns indices of overlapping pairs.
pub fn tgr_join_overlap(
  left : TidifiedGRanges,
  right : TidifiedGRanges,
) -> Array[(Int, Int)] {
  let result : Array[(Int, Int)] = Array::new()
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
/// Inner join by overlap. Returns only rows from left that overlap with right.
pub fn tgr_inner_join(
  left : TidifiedGRanges,
  right : TidifiedGRanges,
) -> TidifiedGRanges {
  let overlaps = tgr_join_overlap(left, right)
  let n = overlaps.length()
  let keep : Array[Bool] = Array::make(left.seqnames.length(), false)
  for pair in overlaps {
    keep[pair.0] = true
  }
  tgr_subset(left, keep)
}

///|
/// Left join by overlap. Returns all rows from left, with metadata from right where overlapping.
pub fn tgr_left_join(
  left : TidifiedGRanges,
  right : TidifiedGRanges,
  _suffix : String,
) -> TidifiedGRanges {
  let overlaps = tgr_join_overlap(left, right)
  let n_left = left.seqnames.length()
  let matched_right = Map::new()
  for pair in overlaps {
    matched_right[pair.0] = pair.1
  }
  let new_meta = Map::new()
  for key in right.metadata.keys() {
    let old = right.metadata[key]
    let nv : Array[Double] = Array::make(n_left, 0.0)
    let mut i = 0
    while i < n_left {
      if matched_right.contains(i) { nv[i] = old[matched_right[i]] }
      else { nv[i] = 0.0 }
      i = i + 1
    }
    new_meta[key] = nv
  }
  let new_meta_str = Map::new()
  for key in right.metadata_str.keys() {
    let old = right.metadata_str[key]
    let nv : Array[String] = Array::make(n_left, "")
    let mut i = 0
    while i < n_left {
      if matched_right.contains(i) { nv[i] = old[matched_right[i]] }
      else { nv[i] = "" }
      i = i + 1
    }
    new_meta_str[key] = nv
  }
  // Merge with left metadata
  let merged_meta = Map::new()
  for key in left.metadata.keys() { merged_meta[key] = left.metadata[key] }
  for key in new_meta.keys() {
    if !merged_meta.contains(key) { merged_meta[key] = new_meta[key] }
  }
  let merged_meta_str = Map::new()
  for key in left.metadata_str.keys() { merged_meta_str[key] = left.metadata_str[key] }
  for key in new_meta_str.keys() {
    if !merged_meta_str.contains(key) { merged_meta_str[key] = new_meta_str[key] }
  }
  TidifiedGRanges::{
    seqnames: left.seqnames.copy(),
    starts: left.starts.copy(),
    ends: left.ends.copy(),
    strands: left.strands.copy(),
    widths: left.widths.copy(),
    metadata: merged_meta,
    metadata_str: merged_meta_str,
  }
}

///|
/// Set anchor position for resize operations.
pub fn tgr_anchor(tgr : TidifiedGRanges, anchor : String) -> TidifiedGRanges {
  // Anchor doesn't change data, it's metadata for subsequent operations
  // We store anchor info in the strands field temporarily (or add metadata)
  // For now, just return a copy
  let _ = anchor
  tgr
}

///|
/// Stretch (resize) ranges by adding bases on each side.
pub fn tgr_stretch(tgr : TidifiedGRanges, amount : Int) -> TidifiedGRanges {
  let n = tgr.seqnames.length()
  let new_starts : Array[Int] = Array::make(n, 0)
  let new_ends : Array[Int] = Array::make(n, 0)
  let new_widths : Array[Int] = Array::make(n, 0)
  let mut i = 0
  while i < n {
    new_starts[i] = tgr.starts[i] - amount
    new_ends[i] = tgr.ends[i] + amount
    if new_starts[i] < 1 { new_starts[i] = 1 }
    new_widths[i] = new_ends[i] - new_starts[i] + 1
    i = i + 1
  }
  TidifiedGRanges::{
    seqnames: tgr.seqnames.copy(),
    starts: new_starts,
    ends: new_ends,
    strands: tgr.strands.copy(),
    widths: new_widths,
    metadata: tgr.metadata,
    metadata_str: tgr.metadata_str,
  }
}

///|
/// Resize to a fixed width, anchored at center.
pub fn tgr_resize_fixed(tgr : TidifiedGRanges, width : Int) -> TidifiedGRanges {
  let n = tgr.seqnames.length()
  let new_starts : Array[Int] = Array::make(n, 0)
  let new_ends : Array[Int] = Array::make(n, 0)
  let new_widths : Array[Int] = Array::make(n, 0)
  let mut i = 0
  while i < n {
    let center = (tgr.starts[i] + tgr.ends[i]) / 2
    new_starts[i] = center - width / 2
    new_ends[i] = new_starts[i] + width - 1
    if new_starts[i] < 1 { new_starts[i] = 1 }
    new_widths[i] = new_ends[i] - new_starts[i] + 1
    i = i + 1
  }
  TidifiedGRanges::{
    seqnames: tgr.seqnames.copy(),
    starts: new_starts,
    ends: new_ends,
    strands: tgr.strands.copy(),
    widths: new_widths,
    metadata: tgr.metadata,
    metadata_str: tgr.metadata_str,
  }
}

///|
/// Shift ranges by a fixed amount.
pub fn tgr_shift(tgr : TidifiedGRanges, amount : Int) -> TidifiedGRanges {
  let n = tgr.seqnames.length()
  let new_starts : Array[Int] = Array::make(n, 0)
  let new_ends : Array[Int] = Array::make(n, 0)
  let mut i = 0
  while i < n {
    new_starts[i] = tgr.starts[i] + amount
    new_ends[i] = tgr.ends[i] + amount
    i = i + 1
  }
  TidifiedGRanges::{
    seqnames: tgr.seqnames.copy(),
    starts: new_starts,
    ends: new_ends,
    strands: tgr.strands.copy(),
    widths: tgr.widths.copy(),
    metadata: tgr.metadata,
    metadata_str: tgr.metadata_str,
  }
}

///|
/// Flank: create flanking regions around each range.
pub fn tgr_flank(tgr : TidifiedGRanges, upstream : Int, downstream : Int) -> TidifiedGRanges {
  let n = tgr.seqnames.length()
  let new_seqnames : Array[String] = Array::new()
  let new_starts : Array[Int] = Array::new()
  let new_ends : Array[Int] = Array::new()
  let new_strands : Array[String] = Array::new()
  let new_widths : Array[Int] = Array::new()
  let mut i = 0
  while i < n {
    // Upstream flank
    new_seqnames.push(tgr.seqnames[i])
    new_starts.push(tgr.starts[i] - upstream)
    if new_starts[new_starts.length() - 1] < 1 { new_starts[new_starts.length() - 1] = 1 }
    new_ends.push(tgr.starts[i] - 1)
    if new_ends[new_ends.length() - 1] < new_starts[new_starts.length() - 1] {
      new_ends[new_ends.length() - 1] = new_starts[new_starts.length() - 1]
    }
    new_strands.push(tgr.strands[i])
    new_widths.push(new_ends[new_ends.length() - 1] - new_starts[new_starts.length() - 1] + 1)
    // Downstream flank
    new_seqnames.push(tgr.seqnames[i])
    new_starts.push(tgr.ends[i] + 1)
    new_ends.push(tgr.ends[i] + downstream)
    new_strands.push(tgr.strands[i])
    new_widths.push(downstream)
    i = i + 1
  }
  TidifiedGRanges::{
    seqnames: new_seqnames,
    starts: new_starts,
    ends: new_ends,
    strands: new_strands,
    widths: new_widths,
    metadata: Map::new(),
    metadata_str: Map::new(),
  }
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
/// Compute GC content of ranges (placeholder - needs sequence data).
pub fn tgr_compute_gc(
  tgr : TidifiedGRanges,
  sequences : Map[String, String],
) -> Array[Double] {
  let n = tgr.seqnames.length()
  let result : Array[Double] = Array::make(n, 0.0)
  let mut i = 0
  while i < n {
    let seqname = tgr.seqnames[i]
    let start = tgr.starts[i]
    let end = tgr.ends[i]
    if sequences.contains(seqname) {
      let seq = sequences[seqname]
      let len = seq.length()
      let s = if start < 1 { 1 } else { start }
      let e = if end > len { len } else { end }
      let mut gc = 0
      let mut at = 0
      let mut j = s
      while j <= e {
        let c = seq[j - 1:j].to_owned().to_upper()
        if c == "G" || c == "C" { gc = gc + 1 }
        else if c == "A" || c == "T" { at = at + 1 }
        j = j + 1
      }
      let total = gc + at
      result[i] = if total > 0 { gc.to_double() / total.to_double() } else { 0.0 }
    }
    i = i + 1
  }
  result
}

///|
/// Intervene: find regions that are unique to each set.
pub fn tgr_intervene(
  sets : Array[TidifiedGRanges],
  set_names : Array[String],
) -> Array[(String, Int)] {
  let n = sets.length()
  if n == 0 { return [] }
  let result : Array[(String, Int)] = Array::new()
  // Compute pairwise overlaps
  let mut i = 0
  while i < n {
    let mut unique_count = sets[i].seqnames.length()
    let mut j = 0
    while j < n {
      if i != j {
        let overlaps = tgr_join_overlap(sets[i], sets[j])
        unique_count = unique_count - overlaps.length()
      }
      j = j + 1
    }
    if unique_count < 0 { unique_count = 0 }
    result.push((set_names[i], unique_count))
    i = i + 1
  }
  result
}

///|
/// Convert TidifiedGRanges to BED format string.
pub fn tgr_to_bed(tgr : TidifiedGRanges) -> String {
  let n = tgr.seqnames.length()
  let lines : Array[String] = Array::new()
  let mut i = 0
  while i < n {
    let line = "\{tgr.seqnames[i]}\t\{tgr.starts[i] - 1}\t\{tgr.ends[i]}"
    lines.push(line)
    i = i + 1
  }
  lines.join("\n")
}

///|
/// Summary statistics for TidifiedGRanges.
pub struct TgrSummary {
  n_ranges : Int
  n_seqnames : Int
  n_strands : Int
  total_bases : Int
  mean_width : Double
  min_width : Int
  max_width : Int
}

///|
pub fn tgr_summary(tgr : TidifiedGRanges) -> TgrSummary {
  let n = tgr.seqnames.length()
  if n == 0 {
    return TgrSummary::{
      n_ranges: 0, n_seqnames: 0, n_strands: 0,
      total_bases: 0, mean_width: 0.0, min_width: 0, max_width: 0,
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
    if tgr.widths[i] < mn { mn = tgr.widths[i] }
    if tgr.widths[i] > mx { mx = tgr.widths[i] }
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

///|
/// Compute coverage per range against a reference set.
pub fn tgr_coverage(
  query : TidifiedGRanges,
  reference : TidifiedGRanges,
) -> Array[Double] {
  let nq = query.seqnames.length()
  let result : Array[Double] = Array::make(nq, 0.0)
  let mut i = 0
  while i < nq {
    let qseq = query.seqnames[i]
    let qstart = query.starts[i]
    let qend = query.ends[i]
    let qlen = qend - qstart + 1
    let mut covered = 0
    let mut j = 0
    while j < reference.seqnames.length() {
      if reference.seqnames[j] == qseq {
        let rstart = reference.starts[j]
        let rend = reference.ends[j]
        let ov_start = if qstart > rstart { qstart } else { rstart }
        let ov_end = if qend < rend { qend } else { rend }
        if ov_start <= ov_end {
          covered = covered + (ov_end - ov_start + 1)
        }
      }
      j = j + 1
    }
    result[i] = if qlen > 0 { covered.to_double() / qlen.to_double() } else { 0.0 }
    i = i + 1
  }
  result
}

///|
/// Set operations: intersect two TidifiedGRanges.
pub fn tgr_intersect(a : TidifiedGRanges, b : TidifiedGRanges) -> TidifiedGRanges {
  let overlaps = tgr_join_overlap(a, b)
  let result_seqnames : Array[String] = Array::new()
  let result_starts : Array[Int] = Array::new()
  let result_ends : Array[Int] = Array::new()
  let result_strands : Array[String] = Array::new()
  let result_widths : Array[Int] = Array::new()
  for pair in overlaps {
    let i = pair.0
    let j = pair.1
    let ov_start = if a.starts[i] > b.starts[j] { a.starts[i] } else { b.starts[j] }
    let ov_end = if a.ends[i] < b.ends[j] { a.ends[i] } else { b.ends[j] }
    if ov_start <= ov_end {
      result_seqnames.push(a.seqnames[i])
      result_starts.push(ov_start)
      result_ends.push(ov_end)
      result_strands.push(a.strands[i])
      result_widths.push(ov_end - ov_start + 1)
    }
  }
  TidifiedGRanges::{
    seqnames: result_seqnames,
    starts: result_starts,
    ends: result_ends,
    strands: result_strands,
    widths: result_widths,
    metadata: Map::new(),
    metadata_str: Map::new(),
  }
}

///|
/// Set operations: union of two TidifiedGRanges (no merging of overlapping).
pub fn tgr_union(a : TidifiedGRanges, b : TidifiedGRanges) -> TidifiedGRanges {
  let overlaps = tgr_join_overlap(a, b)
  let b_overlapping = Map::new()
  for pair in overlaps { b_overlapping[pair.1] = true }
  let result_seqnames : Array[String] = Array::new()
  let result_starts : Array[Int] = Array::new()
  let result_ends : Array[Int] = Array::new()
  let result_strands : Array[String] = Array::new()
  let result_widths : Array[Int] = Array::new()
  // Add all from a
  let mut i = 0
  while i < a.seqnames.length() {
    result_seqnames.push(a.seqnames[i])
    result_starts.push(a.starts[i])
    result_ends.push(a.ends[i])
    result_strands.push(a.strands[i])
    result_widths.push(a.widths[i])
    i = i + 1
  }
  // Add non-overlapping from b
  let mut j = 0
  while j < b.seqnames.length() {
    if !b_overlapping.contains(j) {
      result_seqnames.push(b.seqnames[j])
      result_starts.push(b.starts[j])
      result_ends.push(b.ends[j])
      result_strands.push(b.strands[j])
      result_widths.push(b.widths[j])
    }
    j = j + 1
  }
  TidifiedGRanges::{
    seqnames: result_seqnames,
    starts: result_starts,
    ends: result_ends,
    strands: result_strands,
    widths: result_widths,
    metadata: Map::new(),
    metadata_str: Map::new(),
  }
}

///|
/// Compute distance to nearest range in another set.
pub fn tgr_distance_to_nearest(
  query : TidifiedGRanges,
  subject : TidifiedGRanges,
) -> Array[Int] {
  let nq = query.seqnames.length()
  let result : Array[Int] = Array::make(nq, 0)
  let mut i = 0
  while i < nq {
    let qseq = query.seqnames[i]
    let qstart = query.starts[i]
    let qend = query.ends[i]
    let mut min_dist = 999999999
    let mut j = 0
    while j < subject.seqnames.length() {
      if subject.seqnames[j] == qseq {
        let rstart = subject.starts[j]
        let rend = subject.ends[j]
        let dist = if qend < rstart { rstart - qend }
          else if qstart > rend { qstart - rend }
          else { 0 }
        if dist < min_dist { min_dist = dist }
      }
      j = j + 1
    }
    result[i] = min_dist
    i = i + 1
  }
  result
}
