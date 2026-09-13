///|
/// Affymetrix CEL file parser, ported from Biopython's `Bio.Affy.CelFile`
/// module.
///
/// Only the text-based ASCII CEL format (versions 3.x and 4 as produced by
/// GCOS / AGCC "make-cel-file" tools) is supported, mirroring the scope of
/// the reference implementation. The parser walks the `[HEADER]`,
/// `[INTENSITY]`, `[MASK]`, `[OUTLIERS]`, and `[MODIFIED]` sections and
/// reconstructs the row-major intensity matrices.

///|
/// Error raised for malformed CEL input.
pub suberror CelError {
  CelError(String)
}

///|
/// Parsed CEL record.
///
/// `intensities[y][x]`, `stdevs[y][x]`, and `npix[y][x]` follow the
/// Biopython layout: the first index is the row (`y`), the second the
/// column (`x`).
pub struct CelRecord {
  version : String
  headers : Array[(String, String)]
  n_rows : Int
  n_cols : Int
  n_cells : Int
  intensities : Array[Array[Double]]
  stdevs : Array[Array[Double]]
  npix : Array[Array[Int]]
  mask : Array[Int]
  outliers : Array[Int]
  modified : Array[Int]
  cells_read : Int
} derive(Eq)

///|
fn cel_fail(message : String) -> Unit raise CelError {
  raise CelError(message)
}

///|
fn cel_parse_int(text : String) -> Int? {
  let s = text.trim().to_owned()
  if s.length() == 0 {
    return None
  }
  let mut neg = false
  let mut start = 0
  if s[0] == '-' {
    neg = true
    start = 1
  } else if s[0] == '+' {
    start = 1
  }
  if start >= s.length() {
    return None
  }
  let mut value = 0
  let mut i = start
  while i < s.length() {
    let c = s[i]
    if c < '0' || c > '9' {
      return None
    }
    value = value * 10 + (c.to_int() - '0'.to_int())
    i = i + 1
  }
  Some(if neg { -value } else { value })
}

///|
fn cel_parse_double(text : String) -> Double? {
  let s = text.trim().to_owned()
  if s.length() == 0 {
    return None
  }
  let mut neg = false
  let mut start = 0
  if s[0] == '-' {
    neg = true
    start = 1
  } else if s[0] == '+' {
    start = 1
  }
  if start >= s.length() {
    return None
  }
  let mut int_part = 0.0
  let mut frac_part = 0.0
  let mut frac_div = 1.0
  let mut in_frac = false
  let mut seen_digit = false
  let mut i = start
  while i < s.length() {
    let c = s[i]
    if c >= '0' && c <= '9' {
      seen_digit = true
      let d = (c.to_int() - '0'.to_int()).to_double()
      if in_frac {
        frac_div = frac_div * 10.0
        frac_part = frac_part + d / frac_div
      } else {
        int_part = int_part * 10.0 + d
      }
    } else if c == '.' && !in_frac {
      in_frac = true
    } else {
      return None
    }
    i = i + 1
  }
  if !seen_digit {
    return None
  }
  let v = int_part + frac_part
  Some(if neg { 0.0 - v } else { v })
}

///|
/// Look up a header value by key (case-sensitive, as in the file).
pub fn CelRecord::get_header(self : CelRecord, key : String) -> String? {
  for kv in self.headers {
    if kv.0 == key {
      return Some(kv.1)
    }
  }
  None
}

///|
/// Mean intensity over all cells (probe-level global average).
pub fn CelRecord::mean_intensity(self : CelRecord) -> Double {
  let mut total = 0.0
  let mut count = 0
  for row in self.intensities {
    for value in row {
      total += value
      count += 1
    }
  }
  if count == 0 {
    0.0
  } else {
    total / count.to_double()
  }
}

///|
/// Median intensity over all cells.
pub fn CelRecord::median_intensity(self : CelRecord) -> Double {
  let flat : Array[Double] = Array::new()
  for row in self.intensities {
    for value in row {
      flat.push(value)
    }
  }
  if flat.is_empty() {
    return 0.0
  }
  flat.sort()
  let mid = flat.length() / 2
  if flat.length() % 2 == 1 {
    flat[mid]
  } else {
    (flat[mid - 1] + flat[mid]) / 2.0
  }
}

///|
/// Minimum and maximum intensity, as `(min, max)`.
pub fn CelRecord::intensity_range(self : CelRecord) -> (Double, Double) {
  let mut lo = 0.0
  let mut hi = 0.0
  let mut seen = false
  for row in self.intensities {
    for value in row {
      if !seen {
        lo = value
        hi = value
        seen = true
      } else {
        if value < lo {
          lo = value
        }
        if value > hi {
          hi = value
        }
      }
    }
  }
  (lo, hi)
}

///|
/// Fraction of cells flagged as outliers.
pub fn CelRecord::percent_outliers(self : CelRecord) -> Double {
  if self.n_cells == 0 {
    0.0
  } else {
    self.outliers.length().to_double() / self.n_cells.to_double()
  }
}

///|
/// Fraction of cells flagged as masked.
pub fn CelRecord::percent_masked(self : CelRecord) -> Double {
  if self.n_cells == 0 {
    0.0
  } else {
    self.mask.length().to_double() / self.n_cells.to_double()
  }
}

///|
/// All cell intensities flattened in row-major order.
pub fn CelRecord::intensity_values(self : CelRecord) -> Array[Double] {
  let flat : Array[Double] = Array::new()
  for row in self.intensities {
    for value in row {
      flat.push(value)
    }
  }
  flat
}

///|
/// Intensity percentile (0-100) using linear interpolation between ranks.
pub fn CelRecord::intensity_percentile(
  self : CelRecord,
  percentile : Double,
) -> Double {
  let flat = self.intensity_values()
  if flat.is_empty() {
    return 0.0
  }
  flat.sort()
  let n = flat.length()
  let clamped = if percentile < 0.0 {
    0.0
  } else if percentile > 100.0 {
    100.0
  } else {
    percentile
  }
  let position = (clamped / 100.0) * (n - 1).to_double()
  let lower = position.to_int()
  let upper = if lower + 1 < n { lower + 1 } else { lower }
  let frac = position - lower.to_double()
  flat[lower] + (flat[upper] - flat[lower]) * frac
}

///|
/// Background estimate: the 5th-percentile intensity, the convention used
/// by Affymetrix MAS 5.0 background correction.
pub fn CelRecord::background(self : CelRecord) -> Double {
  self.intensity_percentile(5.0)
}

///|
/// Robust scale estimate: the interquartile range of the intensities.
pub fn CelRecord::intensity_iqr(self : CelRecord) -> Double {
  self.intensity_percentile(75.0) - self.intensity_percentile(25.0)
}

///|
fn cel_lookup(rows : Array[Array[Double]], y : Int, x : Int) -> Double {
  if y < rows.length() && x < rows[y].length() {
    rows[y][x]
  } else {
    0.0
  }
}

///|
/// Probe-level intensity at column `x`, row `y`.
pub fn CelRecord::intensity_at(self : CelRecord, x : Int, y : Int) -> Double {
  cel_lookup(self.intensities, y, x)
}

///|
/// Standard deviation at column `x`, row `y`.
pub fn CelRecord::stdev_at(self : CelRecord, x : Int, y : Int) -> Double {
  cel_lookup(self.stdevs, y, x)
}

///|
fn cel_matrix(rows : Int, cols : Int) -> Array[Array[Double]] {
  let grid : Array[Array[Double]] = Array::new()
  for _y in 0..<rows {
    let row : Array[Double] = Array::new()
    for _x in 0..<cols {
      row.push(0.0)
    }
    grid.push(row)
  }
  grid
}

///|
fn cel_int_matrix(rows : Int, cols : Int) -> Array[Array[Int]] {
  let grid : Array[Array[Int]] = Array::new()
  for _y in 0..<rows {
    let row : Array[Int] = Array::new()
    for _x in 0..<cols {
      row.push(0)
    }
    grid.push(row)
  }
  grid
}

///|
/// Split a header/dat line into whitespace-separated fields
/// (both spaces and tabs act as separators).
fn cel_fields(line : String) -> Array[String] {
  let fields : Array[String] = Array::new()
  for tab_chunk_view in line.split("\t") {
    for field_view in tab_chunk_view.split(" ") {
      let field = field_view.trim().to_owned()
      if field.length() > 0 {
        fields.push(field)
      }
    }
  }
  fields
}

///|
/// Parse grid dimensions from a DatHeader value.
///
/// The leading field encodes the pixel grid as `COLS:ROWS` (e.g.
/// `712:712`); when absent the function returns `(0, 0)`.
pub fn cel_parse_dat_header(value : String) -> (Int, Int) {
  for field_view in value.split(" ") {
    let field = field_view.trim().to_owned()
    if field.length() == 0 || !(field.contains(":")) {
      continue
    }
    let mut cols_text = ""
    let mut rows_text = ""
    let mut in_rows = false
    for ch_view in field.split(":") {
      if in_rows {
        rows_text = ch_view.trim().to_owned()
      } else {
        cols_text = ch_view.trim().to_owned()
        in_rows = true
      }
    }
    match (cel_parse_int(cols_text), cel_parse_int(rows_text)) {
      (Some(cols), Some(rows)) => return (cols, rows)
      _ => continue
    }
  }
  (0, 0)
}

///|
/// Parse a text CEL file, equivalent to `Bio.Affy.CelFile.read`.
pub fn cel_read(content : String) -> CelRecord raise CelError {
  let mut version = ""
  let headers : Array[(String, String)] = Array::new()
  let mut n_rows = 0
  let mut n_cols = 0
  let mut n_cells = 0
  let mut intensities : Array[Array[Double]] = []
  let mut stdevs : Array[Array[Double]] = []
  let mut npix : Array[Array[Int]] = []
  let mask : Array[Int] = Array::new()
  let outliers : Array[Int] = Array::new()
  let modified : Array[Int] = Array::new()
  let mut cells_read = 0
  let mut section = ""
  let mut intensity_meta_seen = false
  let mut grid_set = false
  for line_view in content.split("\n") {
    let raw = line_view.trim().to_owned()
    if raw.has_prefix("[") && raw.has_suffix("]") {
      section = raw
      intensity_meta_seen = false
      continue
    }
    if raw.length() == 0 || raw.has_prefix("#") {
      continue
    }
    if section == "[HEADER]" {
      match raw.find("=") {
        None => continue
        Some(idx) => {
          let key = raw[0:idx].trim().to_owned()
          let value = raw[idx + 1:].trim().to_owned()
          if key == "Version" && version == "" {
            version = value
          }
          if key == "DatHeader" {
            let (cols, rows) = cel_parse_dat_header(value)
            if cols > 0 && rows > 0 {
              n_cols = cols
              n_rows = rows
              grid_set = true
            }
          }
          if key == "Rows" && !grid_set {
            n_rows = cel_parse_int(value).unwrap_or(0)
          }
          if key == "Cols" && !grid_set {
            n_cols = cel_parse_int(value).unwrap_or(0)
          }
          headers.push((key, value))
        }
      }
    } else if section == "[INTENSITY]" {
      if !intensity_meta_seen {
        if raw.has_prefix("NumberCells=") {
          n_cells = cel_parse_int(raw[12:].to_owned()).unwrap_or(0)
          intensities = cel_matrix(n_rows, n_cols)
          stdevs = cel_matrix(n_rows, n_cols)
          npix = cel_int_matrix(n_rows, n_cols)
          intensity_meta_seen = true
        }
        continue
      }
      let parts = cel_fields(raw)
      if parts.length() < 5 {
        continue
      }
      let x = cel_parse_int(parts[0]).unwrap_or(-1)
      let y = cel_parse_int(parts[1]).unwrap_or(-1)
      let mean = cel_parse_double(parts[2]).unwrap_or(0.0)
      let stdv = cel_parse_double(parts[3]).unwrap_or(0.0)
      let npix_value = cel_parse_int(parts[4]).unwrap_or(0)
      if x < 0 || y < 0 || y >= n_rows || x >= n_cols {
        cel_fail(
          "CEL intensity cell (\{parts[0]}, \{parts[1]}) outside \{n_cols}x\{n_rows} grid",
        )
      }
      intensities[y][x] = mean
      stdevs[y][x] = stdv
      npix[y][x] = npix_value
      cells_read += 1
    } else if section == "[MASK]" {
      if raw.has_prefix("NumberCells=") {
        continue
      }
      match cel_parse_int(raw) {
        Some(index) => mask.push(index)
        None => ()
      }
    } else if section == "[OUTLIERS]" {
      if raw.has_prefix("NumberCells=") {
        continue
      }
      match cel_parse_int(raw) {
        Some(index) => outliers.push(index)
        None => ()
      }
    } else if section == "[MODIFIED]" {
      if raw.has_prefix("NumberCells=") {
        continue
      }
      match cel_parse_int(raw) {
        Some(index) => modified.push(index)
        None => ()
      }
    }
  }
  if version == "" && headers.is_empty() && cells_read == 0 {
    cel_fail("no CEL sections found (missing [HEADER]/[INTENSITY])")
  }
  CelRecord::{
    version,
    headers,
    n_rows,
    n_cols,
    n_cells,
    intensities,
    stdevs,
    npix,
    mask,
    outliers,
    modified,
    cells_read,
  }
}

///|
/// Sample CEL text in the version 3.x ASCII format, 3 columns x 2 rows.
pub fn cel_sample_text() -> String {
  "[HEADER]\n" +
  "Version=3.4\n" +
  "Cols=3\n" +
  "Rows=2\n" +
  "MaxX=3\n" +
  "MaxY=2\n" +
  "DatHeader=[0..3]  3:2:  3  2  Expressions\n" +
  "Algorithm=Percentile\n" +
  "AlgorithmParameters=Percentile:75;CellMargin:2;\n" +
  "\n" +
  "[INTENSITY]\n" +
  "NumberCells=6\n" +
  "CellFile=C:\\data\\sample.cel\n" +
  "1 1 130.0 7.0 25\n" +
  "0 0 91.5 4.5 25\n" +
  "1 0 88.0 5.0 25\n" +
  "2 0 102.3 6.1 25\n" +
  "0 1 141.2 8.3 25\n" +
  "2 1 120.9 6.7 25\n" +
  "\n" +
  "[MASK]\n" +
  "NumberCells=1\n" +
  "3\n" +
  "\n" +
  "[OUTLIERS]\n" +
  "NumberCells=1\n" +
  "5\n" +
  "\n" +
  "[MODIFIED]\n" +
  "NumberCells=0\n"
}
