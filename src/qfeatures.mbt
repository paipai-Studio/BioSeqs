///|
/// QFeatures - Quantitative proteomics data management container.
/// Inspired by Bioconductor QFeatures package.
///
/// Provides a unified container for quantitative proteomics data, managing
/// multiple assay levels (e.g., MS1, PSM, peptide, protein) with linking
/// between levels.
///
/// References:
///   - https://www.bioconductor.org/packages/release/bioc/html/QFeatures.html
///   - Gatto L et al. (2020) QFeatures: Quantitative features for mass spectrometry data.
///     bioRxiv 2020.10.27.358273.

// ============================================================================
// Struct definitions
// ============================================================================

///|
/// A single quantitative assay (matrix of features x samples).
pub struct QfAssay {
  name : String
  row_names : Array[String]
  col_names : Array[String]
  data : Array[Array[Double]]
  n_rows : Int
  n_cols : Int
}

///|
/// A link between two features across assay levels.
pub struct QfFeatureLink {
  from_assay : String
  from_row : Int
  to_assay : String
  to_row : Int
}

///|
/// Collection of feature links connecting assay levels.
pub struct QfFeatureLinks {
  links : Array[QfFeatureLink]
}

///|
/// QFeatures container holding multiple assays and their links.
pub struct QFeatures {
  assays : Array[QfAssay]
  links : QfFeatureLinks
  mut n_assays : Int
}

// ============================================================================
// QfAssay methods
// ============================================================================

///|
/// Construct a QfAssay from name, row names, column names, and data matrix.
pub fn QfAssay::new(
  name : String,
  row_names : Array[String],
  col_names : Array[String],
  data : Array[Array[Double]],
) -> QfAssay {
  let nr = data.length()
  let nc = if nr > 0 { data[0].length() } else { 0 }
  QfAssay::{ name, row_names, col_names, data, n_rows: nr, n_cols: nc }
}

///|
/// Construct a QfAssay filled with zeros.
pub fn QfAssay::zeros(name : String, n_rows : Int, n_cols : Int) -> QfAssay {
  let data : Array[Array[Double]] = Array::new()
  for _ in 0..<n_rows {
    let row : Array[Double] = Array::make(n_cols, 0.0)
    data.push(row)
  }
  let row_names : Array[String] = Array::new()
  for i in 0..<n_rows {
    row_names.push("row\{i}")
  }
  let col_names : Array[String] = Array::new()
  for j in 0..<n_cols {
    col_names.push("col\{j}")
  }
  QfAssay::new(name, row_names, col_names, data)
}

///|
/// Get assay name.
pub fn QfAssay::name(self : QfAssay) -> String {
  self.name
}

///|
/// Get row names.
pub fn QfAssay::row_names(self : QfAssay) -> Array[String] {
  self.row_names
}

///|
/// Get column names.
pub fn QfAssay::col_names(self : QfAssay) -> Array[String] {
  self.col_names
}

///|
/// Get data matrix.
pub fn QfAssay::data(self : QfAssay) -> Array[Array[Double]] {
  self.data
}

///|
/// Get number of rows.
pub fn QfAssay::n_rows(self : QfAssay) -> Int {
  self.n_rows
}

///|
/// Get number of columns.
pub fn QfAssay::n_cols(self : QfAssay) -> Int {
  self.n_cols
}

///|
/// Get value at row i, column j.
pub fn QfAssay::get(self : QfAssay, i : Int, j : Int) -> Double {
  self.data[i][j]
}

///|
/// Set value at row i, column j.
pub fn QfAssay::set(self : QfAssay, i : Int, j : Int, val : Double) -> Unit {
  self.data[i][j] = val
}

///|
/// Get a copy of row i.
pub fn QfAssay::get_row(self : QfAssay, i : Int) -> Array[Double] {
  let row : Array[Double] = Array::new()
  for j in 0..<self.n_cols {
    row.push(self.data[i][j])
  }
  row
}

///|
/// Get a copy of column j.
pub fn QfAssay::get_col(self : QfAssay, j : Int) -> Array[Double] {
  let col : Array[Double] = Array::new()
  for i in 0..<self.n_rows {
    col.push(self.data[i][j])
  }
  col
}

///|
/// Compute mean of each row (skipping NaN values).
pub fn QfAssay::row_means(self : QfAssay) -> Array[Double] {
  let means : Array[Double] = Array::new()
  for i in 0..<self.n_rows {
    let mut s = 0.0
    let mut c = 0
    for j in 0..<self.n_cols {
      let v = self.data[i][j]
      if !qf_is_na(v) {
        s = s + v
        c = c + 1
      }
    }
    if c == 0 {
      means.push(0.0 / 0.0)
    } else {
      means.push(s / c.to_double())
    }
  }
  means
}

///|
/// Compute mean of each column (skipping NaN values).
pub fn QfAssay::col_means(self : QfAssay) -> Array[Double] {
  let means : Array[Double] = Array::new()
  for j in 0..<self.n_cols {
    let mut s = 0.0
    let mut c = 0
    for i in 0..<self.n_rows {
      let v = self.data[i][j]
      if !qf_is_na(v) {
        s = s + v
        c = c + 1
      }
    }
    if c == 0 {
      means.push(0.0 / 0.0)
    } else {
      means.push(s / c.to_double())
    }
  }
  means
}

///|
/// Compute sum of each row (skipping NaN values).
pub fn QfAssay::row_sums(self : QfAssay) -> Array[Double] {
  let sums : Array[Double] = Array::new()
  for i in 0..<self.n_rows {
    let mut s = 0.0
    for j in 0..<self.n_cols {
      let v = self.data[i][j]
      if !qf_is_na(v) {
        s = s + v
      }
    }
    sums.push(s)
  }
  sums
}

///|
/// Compute sum of each column (skipping NaN values).
pub fn QfAssay::col_sums(self : QfAssay) -> Array[Double] {
  let sums : Array[Double] = Array::new()
  for j in 0..<self.n_cols {
    let mut s = 0.0
    for i in 0..<self.n_rows {
      let v = self.data[i][j]
      if !qf_is_na(v) {
        s = s + v
      }
    }
    sums.push(s)
  }
  sums
}

///|
/// String representation of the assay.
pub fn QfAssay::to_string(self : QfAssay) -> String {
  let buf = StringBuilder::new()
  buf.write_string("QfAssay '\{self.name}': \{self.n_rows} rows x \{self.n_cols} cols")
  if self.n_rows > 0 && self.n_cols > 0 {
    buf.write_string("\n  data[0][0] = \{self.data[0][0]}")
  }
  buf.to_string()
}

// ============================================================================
// QfFeatureLink methods
// ============================================================================

///|
/// Construct a feature link.
pub fn QfFeatureLink::new(
  from_assay : String,
  from_row : Int,
  to_assay : String,
  to_row : Int,
) -> QfFeatureLink {
  QfFeatureLink::{ from_assay, from_row, to_assay, to_row }
}

///|
/// Get source assay name.
pub fn QfFeatureLink::from_assay(self : QfFeatureLink) -> String {
  self.from_assay
}

///|
/// Get source row index.
pub fn QfFeatureLink::from_row(self : QfFeatureLink) -> Int {
  self.from_row
}

///|
/// Get target assay name.
pub fn QfFeatureLink::to_assay(self : QfFeatureLink) -> String {
  self.to_assay
}

///|
/// Get target row index.
pub fn QfFeatureLink::to_row(self : QfFeatureLink) -> Int {
  self.to_row
}

// ============================================================================
// QfFeatureLinks methods
// ============================================================================

///|
/// Construct an empty QfFeatureLinks.
pub fn QfFeatureLinks::new() -> QfFeatureLinks {
  QfFeatureLinks::{ links: Array::new() }
}

///|
/// Add a link to the collection.
pub fn QfFeatureLinks::add_link(self : QfFeatureLinks, link : QfFeatureLink) -> Unit {
  self.links.push(link)
}

///|
/// Get the number of links.
pub fn QfFeatureLinks::n_links(self : QfFeatureLinks) -> Int {
  self.links.length()
}

///|
/// Get all links originating from (assay, row).
pub fn QfFeatureLinks::get_links_from(
  self : QfFeatureLinks,
  assay : String,
  row : Int,
) -> Array[QfFeatureLink] {
  let result : Array[QfFeatureLink] = Array::new()
  for l in self.links {
    if l.from_assay == assay && l.from_row == row {
      result.push(l)
    }
  }
  result
}

///|
/// Get all links pointing to (assay, row).
pub fn QfFeatureLinks::get_links_to(
  self : QfFeatureLinks,
  assay : String,
  row : Int,
) -> Array[QfFeatureLink] {
  let result : Array[QfFeatureLink] = Array::new()
  for l in self.links {
    if l.to_assay == assay && l.to_row == row {
      result.push(l)
    }
  }
  result
}

// ============================================================================
// QFeatures methods
// ============================================================================

///|
/// Construct an empty QFeatures.
pub fn QFeatures::new() -> QFeatures {
  QFeatures::{
    assays: Array::new(),
    links: QfFeatureLinks::new(),
    n_assays: 0,
  }
}

///|
/// Get all assays.
pub fn QFeatures::assays(self : QFeatures) -> Array[QfAssay] {
  self.assays
}

///|
/// Get feature links.
pub fn QFeatures::links(self : QFeatures) -> QfFeatureLinks {
  self.links
}

///|
/// Get number of assays.
pub fn QFeatures::n_assays(self : QFeatures) -> Int {
  self.n_assays
}

///|
/// Get an assay by name, returning None if not found.
pub fn QFeatures::get_assay(self : QFeatures, name : String) -> QfAssay? {
  for a in self.assays {
    if a.name == name {
      return Some(a)
    }
  }
  None
}

///|
/// Get names of all assays.
pub fn QFeatures::assay_names(self : QFeatures) -> Array[String] {
  let names : Array[String] = Array::new()
  for a in self.assays {
    names.push(a.name)
  }
  names
}

///|
/// Add an assay to the container.
pub fn QFeatures::add_assay(self : QFeatures, assay : QfAssay) -> Unit {
  self.assays.push(assay)
  self.n_assays = self.assays.length()
}

///|
/// Add a link between two assay features.
pub fn QFeatures::add_link(
  self : QFeatures,
  from_assay : String,
  from_row : Int,
  to_assay : String,
  to_row : Int,
) -> Unit {
  self.links.add_link(QfFeatureLink::new(from_assay, from_row, to_assay, to_row))
}

///|
/// Get number of features (rows) in an assay.
pub fn QFeatures::n_features(self : QFeatures, assay_name : String) -> Int {
  match self.get_assay(assay_name) {
    Some(a) => a.n_rows
    None => 0
  }
}

///|
/// Get number of samples (columns) in an assay.
pub fn QFeatures::n_samples(self : QFeatures, assay_name : String) -> Int {
  match self.get_assay(assay_name) {
    Some(a) => a.n_cols
    None => 0
  }
}

// ============================================================================
// Internal helper functions
// ============================================================================

///|
/// Check if a value is NaN.
pub fn qf_is_na(val : Double) -> Bool {
  val != val
}

fn qf_infinity() -> Double {
  1.7976931348623157e+308
}

/// Compute the median of an array (sorts a copy, does not modify input).
fn qf_median(values : Array[Double]) -> Double {
  let n = values.length()
  if n == 0 {
    return 0.0 / 0.0
  }
  let sorted : Array[Double] = Array::new()
  for v in values {
    sorted.push(v)
  }
  sorted.sort_by(fn(a : Double, b : Double) -> Int {
    if a < b { -1 } else if a > b { 1 } else { 0 }
  })
  let mid = n / 2
  if n % 2 == 0 {
    (sorted[mid - 1] + sorted[mid]) / 2.0
  } else {
    sorted[mid]
  }
}

/// Aggregate a list of values using the specified method (skipping NaN).
fn qf_aggregate_values(values : Array[Double], method : String) -> Double {
  let clean : Array[Double] = Array::new()
  for v in values {
    if !qf_is_na(v) {
      clean.push(v)
    }
  }
  if clean.length() == 0 {
    return 0.0 / 0.0
  }
  match method {
    "sum" => {
      let mut s = 0.0
      for v in clean {
        s = s + v
      }
      s
    }
    "mean" => {
      let mut s = 0.0
      for v in clean {
        s = s + v
      }
      s / clean.length().to_double()
    }
    "median" => qf_median(clean)
    "max" => {
      let mut m = clean[0]
      for v in clean {
        if v > m {
          m = v
        }
      }
      m
    }
    _ => {
      let mut s = 0.0
      for v in clean {
        s = s + v
      }
      s / clean.length().to_double()
    }
  }
}

/// Euclidean distance between two rows on shared non-NaN columns.
fn qf_row_distance(a : Array[Double], b : Array[Double]) -> Double {
  let n = a.length()
  let mut sum_sq = 0.0
  let mut count = 0
  for j in 0..<n {
    if !qf_is_na(a[j]) && !qf_is_na(b[j]) {
      let d = a[j] - b[j]
      sum_sq = sum_sq + d * d
      count = count + 1
    }
  }
  if count == 0 {
    qf_infinity()
  } else {
    (sum_sq / count.to_double()).sqrt()
  }
}

/// Deep-copy a data matrix.
fn qf_copy_matrix(data : Array[Array[Double]]) -> Array[Array[Double]] {
  let nr = data.length()
  let out : Array[Array[Double]] = Array::new()
  for i in 0..<nr {
    let row : Array[Double] = Array::new()
    for j in 0..<data[i].length() {
      row.push(data[i][j])
    }
    out.push(row)
  }
  out
}

// ============================================================================
// Aggregation
// ============================================================================

///|
/// Aggregate features from one assay to create a new assay.
///
/// mapping: array of (row_index, group_name) pairs mapping each source row
/// to a group. The new assay has one row per unique group name.
/// method: "sum", "mean", "median", or "max".
/// Links are automatically created between source rows and new assay rows.
pub fn qf_aggregate(
  qf : QFeatures,
  from_assay : String,
  to_assay : String,
  mapping : Array[(Int, String)],
  method : String,
) -> QFeatures {
  let src = match qf.get_assay(from_assay) {
    Some(a) => a
    None => return qf
  }
  // Determine unique group names in order of first appearance.
  let group_names : Array[String] = Array::new()
  for m in mapping {
    let gname = m.1
    let mut found = false
    for g in group_names {
      if g == gname {
        found = true
        break
      }
    }
    if !found {
      group_names.push(gname)
    }
  }
  let n_groups = group_names.length()
  let nc = src.n_cols
  // Collect source row indices for each group.
  let group_rows : Array[Array[Int]] = Array::new()
  for _ in 0..<n_groups {
    group_rows.push(Array::new())
  }
  for m in mapping {
    let row_idx = m.0
    let gname = m.1
    for g in 0..<n_groups {
      if group_names[g] == gname {
        group_rows[g].push(row_idx)
        break
      }
    }
  }
  // Compute aggregated data.
  let new_data : Array[Array[Double]] = Array::new()
  for g in 0..<n_groups {
    let row : Array[Double] = Array::make(nc, 0.0)
    for j in 0..<nc {
      let values : Array[Double] = Array::new()
      for r in group_rows[g] {
        if r >= 0 && r < src.n_rows {
          values.push(src.data[r][j])
        }
      }
      row[j] = qf_aggregate_values(values, method)
    }
    new_data.push(row)
  }
  let new_assay = QfAssay::new(to_assay, group_names, src.col_names, new_data)
  // Build new QFeatures with the assay and links added.
  let new_assays : Array[QfAssay] = Array::new()
  for a in qf.assays {
    new_assays.push(a)
  }
  new_assays.push(new_assay)
  let new_links = QfFeatureLinks::new()
  for l in qf.links.links {
    new_links.add_link(l)
  }
  for g in 0..<n_groups {
    for r in group_rows[g] {
      new_links.add_link(QfFeatureLink::new(from_assay, r, to_assay, g))
    }
  }
  QFeatures::{
    assays: new_assays,
    links: new_links,
    n_assays: new_assays.length(),
  }
}

///|
/// Aggregate columns (e.g., combine technical replicates).
///
/// group_labels: one label per column in from_assay. Columns sharing a label
/// are combined into a single column in the new assay.
/// method: "sum" or "mean".
pub fn qf_aggregate_by_col(
  qf : QFeatures,
  from_assay : String,
  to_assay : String,
  group_labels : Array[String],
  method : String,
) -> QFeatures {
  let src = match qf.get_assay(from_assay) {
    Some(a) => a
    None => return qf
  }
  let nc = src.n_cols
  // Determine unique labels in order of first appearance.
  let unique_labels : Array[String] = Array::new()
  for j in 0..<nc {
    if j < group_labels.length() {
      let label = group_labels[j]
      let mut found = false
      for l in unique_labels {
        if l == label {
          found = true
          break
        }
      }
      if !found {
        unique_labels.push(label)
      }
    }
  }
  // Collect column indices for each unique label.
  let group_cols : Array[Array[Int]] = Array::new()
  for _ in 0..<unique_labels.length() {
    group_cols.push(Array::new())
  }
  for j in 0..<nc {
    if j < group_labels.length() {
      let label = group_labels[j]
      for g in 0..<unique_labels.length() {
        if unique_labels[g] == label {
          group_cols[g].push(j)
          break
        }
      }
    }
  }
  // Compute aggregated data.
  let nr = src.n_rows
  let new_data : Array[Array[Double]] = Array::new()
  for i in 0..<nr {
    let row : Array[Double] = Array::new()
    for g in 0..<unique_labels.length() {
      let values : Array[Double] = Array::new()
      for j in group_cols[g] {
        values.push(src.data[i][j])
      }
      match method {
        "sum" => {
          let mut s = 0.0
          for v in values {
            if !qf_is_na(v) {
              s = s + v
            }
          }
          row.push(s)
        }
        _ => {
          let mut s = 0.0
          let mut c = 0
          for v in values {
            if !qf_is_na(v) {
              s = s + v
              c = c + 1
            }
          }
          if c == 0 {
            row.push(0.0 / 0.0)
          } else {
            row.push(s / c.to_double())
          }
        }
      }
    }
    new_data.push(row)
  }
  let new_assay = QfAssay::new(to_assay, src.row_names, unique_labels, new_data)
  let new_assays : Array[QfAssay] = Array::new()
  for a in qf.assays {
    new_assays.push(a)
  }
  new_assays.push(new_assay)
  let new_links = QfFeatureLinks::new()
  for l in qf.links.links {
    new_links.add_link(l)
  }
  QFeatures::{
    assays: new_assays,
    links: new_links,
    n_assays: new_assays.length(),
  }
}

// ============================================================================
// Filtering
// ============================================================================

///|
/// Filter rows in an assay, keeping only specified rows. Updates links
/// accordingly (removed rows drop their links; kept rows are reindexed).
pub fn qf_filter_features(
  qf : QFeatures,
  assay_name : String,
  keep_rows : Array[Int],
) -> QFeatures {
  // Find source assay dimension.
  let mut src_n_rows = 0
  let mut found = false
  for a in qf.assays {
    if a.name == assay_name {
      src_n_rows = a.n_rows
      found = true
      break
    }
  }
  if !found {
    return qf
  }
  // Build old-to-new index map.
  let old_to_new : Array[Int] = Array::make(src_n_rows, -1)
  for new_idx in 0..<keep_rows.length() {
    let old_idx = keep_rows[new_idx]
    if old_idx >= 0 && old_idx < src_n_rows {
      old_to_new[old_idx] = new_idx
    }
  }
  // Build new assays.
  let new_assays : Array[QfAssay] = Array::new()
  for a in qf.assays {
    if a.name == assay_name {
      let new_row_names : Array[String] = Array::new()
      let new_data : Array[Array[Double]] = Array::new()
      for r in keep_rows {
        if r >= 0 && r < a.n_rows {
          new_row_names.push(a.row_names[r])
          let row_copy : Array[Double] = Array::new()
          for j in 0..<a.n_cols {
            row_copy.push(a.data[r][j])
          }
          new_data.push(row_copy)
        }
      }
      let new_assay = QfAssay::new(a.name, new_row_names, a.col_names, new_data)
      new_assays.push(new_assay)
    } else {
      new_assays.push(a)
    }
  }
  // Update links.
  let new_links = QfFeatureLinks::new()
  for l in qf.links.links {
    let mut keep_link = true
    let mut new_from_row = l.from_row
    let mut new_to_row = l.to_row
    if l.from_assay == assay_name {
      if l.from_row >= 0 && l.from_row < src_n_rows {
        new_from_row = old_to_new[l.from_row]
        if new_from_row < 0 {
          keep_link = false
        }
      } else {
        keep_link = false
      }
    }
    if l.to_assay == assay_name {
      if l.to_row >= 0 && l.to_row < src_n_rows {
        new_to_row = old_to_new[l.to_row]
        if new_to_row < 0 {
          keep_link = false
        }
      } else {
        keep_link = false
      }
    }
    if keep_link {
      new_links.add_link(QfFeatureLink::new(
        l.from_assay,
        new_from_row,
        l.to_assay,
        new_to_row,
      ))
    }
  }
  QFeatures::{
    assays: new_assays,
    links: new_links,
    n_assays: new_assays.length(),
  }
}

///|
/// Filter columns in an assay, keeping only specified columns.
pub fn qf_filter_samples(
  qf : QFeatures,
  assay_name : String,
  keep_cols : Array[Int],
) -> QFeatures {
  let new_assays : Array[QfAssay] = Array::new()
  for a in qf.assays {
    if a.name == assay_name {
      let new_col_names : Array[String] = Array::new()
      for c in keep_cols {
        if c >= 0 && c < a.n_cols {
          new_col_names.push(a.col_names[c])
        }
      }
      let new_data : Array[Array[Double]] = Array::new()
      for i in 0..<a.n_rows {
        let row : Array[Double] = Array::new()
        for c in keep_cols {
          if c >= 0 && c < a.n_cols {
            row.push(a.data[i][c])
          }
        }
        new_data.push(row)
      }
      let new_assay = QfAssay::new(a.name, a.row_names, new_col_names, new_data)
      new_assays.push(new_assay)
    } else {
      new_assays.push(a)
    }
  }
  // Links are row-based and unaffected by column filtering.
  let new_links = QfFeatureLinks::new()
  for l in qf.links.links {
    new_links.add_link(l)
  }
  QFeatures::{
    assays: new_assays,
    links: new_links,
    n_assays: new_assays.length(),
  }
}

///|
/// Remove rows with more than max_na_frac fraction of NaN values.
pub fn qf_filter_na(
  qf : QFeatures,
  assay_name : String,
  max_na_frac : Double,
) -> QFeatures {
  let src = match qf.get_assay(assay_name) {
    Some(a) => a
    None => return qf
  }
  let nc = src.n_cols
  let keep_rows : Array[Int] = Array::new()
  for i in 0..<src.n_rows {
    let mut na_count = 0
    for j in 0..<nc {
      if qf_is_na(src.data[i][j]) {
        na_count = na_count + 1
      }
    }
    let frac = if nc > 0 {
      na_count.to_double() / nc.to_double()
    } else {
      0.0
    }
    if frac <= max_na_frac {
      keep_rows.push(i)
    }
  }
  qf_filter_features(qf, assay_name, keep_rows)
}

///|
/// Remove rows where all values are below threshold.
pub fn qf_filter_low_abundance(
  qf : QFeatures,
  assay_name : String,
  threshold : Double,
) -> QFeatures {
  let src = match qf.get_assay(assay_name) {
    Some(a) => a
    None => return qf
  }
  let keep_rows : Array[Int] = Array::new()
  for i in 0..<src.n_rows {
    let mut all_below = true
    for j in 0..<src.n_cols {
      let v = src.data[i][j]
      if !qf_is_na(v) && v >= threshold {
        all_below = false
        break
      }
    }
    if !all_below {
      keep_rows.push(i)
    }
  }
  qf_filter_features(qf, assay_name, keep_rows)
}

// ============================================================================
// Normalization
// ============================================================================

///|
/// Quantile normalization across columns.
///
/// Sorts each column independently, computes the mean of each rank across
/// columns, then replaces each value with the reference mean at its rank.
pub fn qf_normalize_quantiles(assay : QfAssay) -> QfAssay {
  let nr = assay.n_rows
  let nc = assay.n_cols
  if nr == 0 || nc == 0 {
    return assay
  }
  // Sort each column independently.
  let sorted_cols : Array[Array[Double]] = Array::new()
  for j in 0..<nc {
    let col = assay.get_col(j)
    col.sort_by(fn(a : Double, b : Double) -> Int {
      if a < b { -1 } else if a > b { 1 } else { 0 }
    })
    sorted_cols.push(col)
  }
  // Compute rank means (mean across columns for each rank position).
  let rank_means : Array[Double] = Array::make(nr, 0.0)
  for i in 0..<nr {
    let mut s = 0.0
    for j in 0..<nc {
      s = s + sorted_cols[j][i]
    }
    rank_means[i] = s / nc.to_double()
  }
  // For each column, assign rank means based on the rank of each value.
  let new_data : Array[Array[Double]] = Array::new()
  for _ in 0..<nr {
    new_data.push(Array::make(nc, 0.0))
  }
  for j in 0..<nc {
    // Build (original_index, value) pairs and sort by value to determine ranks.
    let pairs : Array[(Int, Double)] = Array::new()
    for i in 0..<nr {
      pairs.push((i, assay.data[i][j]))
    }
    pairs.sort_by(fn(a : (Int, Double), b : (Int, Double)) -> Int {
      if a.1 < b.1 { -1 } else if a.1 > b.1 { 1 } else { 0 }
    })
    for rank in 0..<nr {
      let orig_idx = pairs[rank].0
      new_data[orig_idx][j] = rank_means[rank]
    }
  }
  QfAssay::new(assay.name, assay.row_names, assay.col_names, new_data)
}

///|
/// Center (subtract mean) and/or scale (divide by std) each column.
///
/// NaN values are skipped when computing mean and std, and remain NaN
/// in the output.
pub fn qf_normalize_center_scale(
  assay : QfAssay,
  center : Bool,
  scale : Bool,
) -> QfAssay {
  let nr = assay.n_rows
  let nc = assay.n_cols
  if nr == 0 || nc == 0 {
    return assay
  }
  let new_data : Array[Array[Double]] = Array::new()
  for _ in 0..<nr {
    new_data.push(Array::make(nc, 0.0))
  }
  for j in 0..<nc {
    // Compute column mean (skipping NaN).
    let mut s = 0.0
    let mut c = 0
    for i in 0..<nr {
      let v = assay.data[i][j]
      if !qf_is_na(v) {
        s = s + v
        c = c + 1
      }
    }
    let mean = if c > 0 { s / c.to_double() } else { 0.0 }
    // Compute column std (skipping NaN).
    let mut ss = 0.0
    if c > 0 {
      for i in 0..<nr {
        let v = assay.data[i][j]
        if !qf_is_na(v) {
          let d = v - mean
          ss = ss + d * d
        }
      }
    }
    let std = if c > 1 {
      (ss / (c - 1).to_double()).sqrt()
    } else {
      1.0
    }
    let denom = if scale && std > 0.0 { std } else { 1.0 }
    for i in 0..<nr {
      let v = assay.data[i][j]
      if qf_is_na(v) {
        new_data[i][j] = v
      } else {
        let centered = if center { v - mean } else { v }
        new_data[i][j] = centered / denom
      }
    }
  }
  QfAssay::new(assay.name, assay.row_names, assay.col_names, new_data)
}

///|
/// Log transform: log_base(x + offset).
///
/// Uses @math.ln(x + offset) / @math.ln(base) for arbitrary base.
/// NaN values remain NaN; non-positive values after offset become NaN.
pub fn qf_log_transform(
  assay : QfAssay,
  base : Double,
  offset : Double,
) -> QfAssay {
  let nr = assay.n_rows
  let nc = assay.n_cols
  let log_base = @math.ln(base)
  let new_data : Array[Array[Double]] = Array::new()
  for i in 0..<nr {
    let row : Array[Double] = Array::new()
    for j in 0..<nc {
      let v = assay.data[i][j]
      if qf_is_na(v) {
        row.push(v)
      } else {
        let val = v + offset
        if val > 0.0 {
          row.push(@math.ln(val) / log_base)
        } else {
          row.push(0.0 / 0.0)
        }
      }
    }
    new_data.push(row)
  }
  QfAssay::new(assay.name, assay.row_names, assay.col_names, new_data)
}

// ============================================================================
// Imputation
// ============================================================================

///|
/// K-nearest neighbor imputation for NaN values.
///
/// For each row with NaN, finds k nearest rows by Euclidean distance on
/// shared non-NaN columns, then imputes NaN positions with the
/// inverse-distance-weighted average of the neighbors' values.
pub fn qf_impute_knn(assay : QfAssay, k : Int) -> QfAssay {
  let nr = assay.n_rows
  let nc = assay.n_cols
  if nr == 0 || nc == 0 {
    return assay
  }
  let data = qf_copy_matrix(assay.data)
  let eps = 1.0e-8
  for i in 0..<nr {
    // Check if row i has any NaN.
    let mut has_na = false
    for j in 0..<nc {
      if qf_is_na(data[i][j]) {
        has_na = true
        break
      }
    }
    if !has_na {
      continue
    }
    // Compute distances to all other rows.
    let pairs : Array[(Int, Double)] = Array::new()
    for other in 0..<nr {
      if other == i {
        continue
      }
      let dist = qf_row_distance(data[i], data[other])
      pairs.push((other, dist))
    }
    // Sort by distance (ascending).
    pairs.sort_by(fn(a : (Int, Double), b : (Int, Double)) -> Int {
      if a.1 < b.1 { -1 } else if a.1 > b.1 { 1 } else { 0 }
    })
    let actual_k = if k < pairs.length() { k } else { pairs.length() }
    // Impute each NaN position.
    for j in 0..<nc {
      if !qf_is_na(data[i][j]) {
        continue
      }
      let mut num = 0.0
      let mut den = 0.0
      for idx in 0..<actual_k {
        let ni = pairs[idx].0
        let d = pairs[idx].1
        let v = data[ni][j]
        if qf_is_na(v) {
          continue
        }
        let w = 1.0 / (d + eps)
        num = num + w * v
        den = den + w
      }
      if den > 0.0 {
        data[i][j] = num / den
      } else {
        data[i][j] = 0.0
      }
    }
  }
  QfAssay::new(assay.name, assay.row_names, assay.col_names, data)
}

///|
/// Replace NaN with column mean (computed from non-NaN values).
pub fn qf_impute_mean(assay : QfAssay) -> QfAssay {
  let nr = assay.n_rows
  let nc = assay.n_cols
  if nr == 0 || nc == 0 {
    return assay
  }
  // Compute column means (skipping NaN).
  let col_means : Array[Double] = Array::make(nc, 0.0)
  for j in 0..<nc {
    let mut s = 0.0
    let mut c = 0
    for i in 0..<nr {
      let v = assay.data[i][j]
      if !qf_is_na(v) {
        s = s + v
        c = c + 1
      }
    }
    col_means[j] = if c > 0 { s / c.to_double() } else { 0.0 }
  }
  let new_data : Array[Array[Double]] = Array::new()
  for i in 0..<nr {
    let row : Array[Double] = Array::new()
    for j in 0..<nc {
      let v = assay.data[i][j]
      if qf_is_na(v) {
        row.push(col_means[j])
      } else {
        row.push(v)
      }
    }
    new_data.push(row)
  }
  QfAssay::new(assay.name, assay.row_names, assay.col_names, new_data)
}

///|
/// Replace NaN with 0.0.
pub fn qf_impute_zero(assay : QfAssay) -> QfAssay {
  let nr = assay.n_rows
  let nc = assay.n_cols
  let new_data : Array[Array[Double]] = Array::new()
  for i in 0..<nr {
    let row : Array[Double] = Array::new()
    for j in 0..<nc {
      let v = assay.data[i][j]
      if qf_is_na(v) {
        row.push(0.0)
      } else {
        row.push(v)
      }
    }
    new_data.push(row)
  }
  QfAssay::new(assay.name, assay.row_names, assay.col_names, new_data)
}

// ============================================================================
// Missing value handling
// ============================================================================

///|
/// Count total NaN values in an assay.
pub fn qf_count_na(assay : QfAssay) -> Int {
  let mut count = 0
  for i in 0..<assay.n_rows {
    for j in 0..<assay.n_cols {
      if qf_is_na(assay.data[i][j]) {
        count = count + 1
      }
    }
  }
  count
}

///|
/// Count NaN values per row.
pub fn qf_count_na_by_row(assay : QfAssay) -> Array[Int] {
  let counts : Array[Int] = Array::new()
  for i in 0..<assay.n_rows {
    let mut c = 0
    for j in 0..<assay.n_cols {
      if qf_is_na(assay.data[i][j]) {
        c = c + 1
      }
    }
    counts.push(c)
  }
  counts
}

// ============================================================================
// Summary and utility
// ============================================================================

///|
/// Print a summary of the QFeatures container: number of assays, names,
/// dimensions, and link count.
pub fn qf_summary(qf : QFeatures) -> String {
  let buf = StringBuilder::new()
  buf.write_string("QFeatures Summary\n")
  buf.write_string("================\n")
  buf.write_string("Number of assays: \{qf.n_assays}\n")
  buf.write_string("Number of links: \{qf.links.n_links()}\n")
  buf.write_string("\nAssays:\n")
  for a in qf.assays {
    buf.write_string(
      "  \{a.name}: \{a.n_rows} features x \{a.n_cols} samples\n",
    )
  }
  buf.to_string()
}

///|
/// Print an assay summary: dimensions, mean, range, and NA count.
pub fn qf_assay_summary(assay : QfAssay) -> String {
  let buf = StringBuilder::new()
  buf.write_string("QfAssay Summary: '\{assay.name}'\n")
  buf.write_string("  Dimensions: \{assay.n_rows} rows x \{assay.n_cols} cols\n")
  // Collect non-NaN values for statistics.
  let all_vals : Array[Double] = Array::new()
  for i in 0..<assay.n_rows {
    for j in 0..<assay.n_cols {
      let v = assay.data[i][j]
      if !qf_is_na(v) {
        all_vals.push(v)
      }
    }
  }
  let na_count = qf_count_na(assay)
  let total = assay.n_rows * assay.n_cols
  if all_vals.length() > 0 {
    let mut s = 0.0
    let mut mn = all_vals[0]
    let mut mx = all_vals[0]
    for v in all_vals {
      s = s + v
      if v < mn {
        mn = v
      }
      if v > mx {
        mx = v
      }
    }
    let mean = s / all_vals.length().to_double()
    buf.write_string("  Mean: \{mean}\n")
    buf.write_string("  Range: [\{mn}, \{mx}]\n")
  } else {
    buf.write_string("  Mean: NA (all values missing)\n")
    buf.write_string("  Range: NA\n")
  }
  buf.write_string("  NA count: \{na_count} / \{total}\n")
  buf.to_string()
}

///|
/// Convert an assay to long format: array of (row_name, col_name, value).
pub fn qf_to_long_format(
  qf : QFeatures,
  assay_name : String,
) -> Array[(String, String, Double)] {
  let result : Array[(String, String, Double)] = Array::new()
  let src = match qf.get_assay(assay_name) {
    Some(a) => a
    None => return result
  }
  for i in 0..<src.n_rows {
    for j in 0..<src.n_cols {
      result.push((src.row_names[i], src.col_names[j], src.data[i][j]))
    }
  }
  result
}

// ============================================================================
// Sample data
// ============================================================================

///|
/// Create a sample QFeatures with three assay levels (psm, peptide, protein)
/// and links between all levels.
///
/// - "psm": 6 PSMs x 4 samples
/// - "peptide": 4 peptides x 4 samples (aggregated from PSM by mean)
/// - "protein": 3 proteins x 4 samples (aggregated from peptide by mean)
pub fn qf_sample_data() -> QFeatures {
  // Create PSM assay (6 PSMs x 4 samples).
  let psm_row_names = ["PSM1", "PSM2", "PSM3", "PSM4", "PSM5", "PSM6"]
  let psm_col_names = ["S1", "S2", "S3", "S4"]
  let psm_data : Array[Array[Double]] = [
    [10.0, 12.0, 11.0, 13.0],
    [20.0, 22.0, 21.0, 23.0],
    [30.0, 32.0, 31.0, 33.0],
    [40.0, 42.0, 41.0, 43.0],
    [50.0, 52.0, 51.0, 53.0],
    [60.0, 62.0, 61.0, 63.0],
  ]
  let psm_assay = QfAssay::new("psm", psm_row_names, psm_col_names, psm_data)
  let mut qf = QFeatures::new()
  qf.add_assay(psm_assay)
  // Aggregate PSM -> peptide.
  // PSM1, PSM2 -> PepA
  // PSM3, PSM4 -> PepB
  // PSM5 -> PepC
  // PSM6 -> PepD
  let psm_to_peptide : Array[(Int, String)] = [
    (0, "PepA"),
    (1, "PepA"),
    (2, "PepB"),
    (3, "PepB"),
    (4, "PepC"),
    (5, "PepD"),
  ]
  qf = qf_aggregate(qf, "psm", "peptide", psm_to_peptide, "mean")
  // Aggregate peptide -> protein.
  // PepA, PepB -> Prot1
  // PepC -> Prot2
  // PepD -> Prot3
  let peptide_to_protein : Array[(Int, String)] = [
    (0, "Prot1"),
    (1, "Prot1"),
    (2, "Prot2"),
    (3, "Prot3"),
  ]
  qf = qf_aggregate(qf, "peptide", "protein", peptide_to_protein, "mean")
  qf
}
