///|
/// BIOM 1.0 (Biological Observation Matrix, JSON) table — parse and write.
/// Mirrors scikit-bio's `skbio.io.biom` `Table`. BIOM 2.1 (HDF5) is out of
/// scope (no HDF5 dependency).

///|
pub struct BiomTable {
  table_id : String
  matrix_type : String
  matrix_element_type : String
  observation_ids : Array[String]
  sample_ids : Array[String]
  observation_metadata : Array[Map[String, JsonValue]]
  sample_metadata : Array[Map[String, JsonValue]]
  data : Array[Array[Double]]
}

///|
fn json_as_str(v : JsonValue) -> String? {
  match v {
    JsonStr(s) => Some(s)
    _ => None
  }
}

///|
fn json_as_num(v : JsonValue) -> Double? {
  match v {
    JsonNum(d) => Some(d)
    _ => None
  }
}

///|
fn json_as_obj(v : JsonValue) -> Map[String, JsonValue]? {
  match v {
    JsonObj(m) => Some(m)
    _ => None
  }
}

///|
fn json_as_arr(v : JsonValue) -> Array[JsonValue]? {
  match v {
    JsonArr(a) => Some(a)
    _ => None
  }
}

///|
/// Extract (id, metadata) from a BIOM row/column entry.
fn biom_id_metadata(
  entry : JsonValue,
) -> (String, Map[String, JsonValue])? {
  let m = match json_as_obj(entry) {
    Some(x) => x
    None => return None
  }
  let id = match m.get("id") {
    Some(x) => {
      match json_as_str(x) {
        Some(s) => s
        None => return None
      }
    }
    None => return None
  }
  let meta = match m.get("metadata") {
    Some(JsonObj(mm)) => mm
    _ => Map([], capacity=0)
  }
  Some((id, meta))
}

///|
/// Parse a BIOM 1.0 JSON document into a `BiomTable`.
pub fn parse_biom(text : String) -> BiomTable? {
  let root = match parse_json(text) {
    Some(x) => x
    None => return None
  }
  let obj = match json_as_obj(root) {
    Some(x) => x
    None => return None
  }

  let table_id = match obj.get("id") {
    Some(x) => json_as_str(x).unwrap_or("")
    None => ""
  }
  let matrix_type = match obj.get("matrix_type") {
    Some(x) => json_as_str(x).unwrap_or("dense")
    None => "dense"
  }
  let matrix_element_type = match obj.get("matrix_element_type") {
    Some(x) => json_as_str(x).unwrap_or("float")
    None => "float"
  }

  let shape = match json_as_arr(obj.get("shape").unwrap_or(JsonNull)) {
    Some(x) => x
    None => return None
  }
  let n_rows = json_as_num(shape[0]).unwrap_or(0.0).to_int()
  let n_cols = json_as_num(shape[1]).unwrap_or(0.0).to_int()

  let rows_arr = match json_as_arr(obj.get("rows").unwrap_or(JsonNull)) {
    Some(x) => x
    None => return None
  }
  let cols_arr = match json_as_arr(obj.get("columns").unwrap_or(JsonNull)) {
    Some(x) => x
    None => return None
  }

  let observation_ids : Array[String] = []
  let observation_metadata : Array[Map[String, JsonValue]] = []
  for r in rows_arr {
    match biom_id_metadata(r) {
      Some((id, meta)) => {
        observation_ids.push(id)
        observation_metadata.push(meta)
      }
      None => return None
    }
  }
  let sample_ids : Array[String] = []
  let sample_metadata : Array[Map[String, JsonValue]] = []
  for c in cols_arr {
    match biom_id_metadata(c) {
      Some((id, meta)) => {
        sample_ids.push(id)
        sample_metadata.push(meta)
      }
      None => return None
    }
  }

  // data: 稠密或稀疏
  let dens : Array[Array[Double]] = Array(capacity=n_rows)
  for i = 0; i < n_rows; i = i + 1 {
    dens.push(Array::make(n_cols, 0.0))
  }
  let data_node = obj.get("data").unwrap_or(JsonNull)
  if matrix_type == "sparse" {
    let entries = match json_as_arr(data_node) {
      Some(x) => x
      None => return None
    }
    for e in entries {
      let t = match json_as_arr(e) {
        Some(x) => x
        None => return None
      }
      let r = json_as_num(t[0]).unwrap_or(0.0).to_int()
      let c = json_as_num(t[1]).unwrap_or(0.0).to_int()
      let v = json_as_num(t[2]).unwrap_or(0.0)
      if r >= 0 && r < n_rows && c >= 0 && c < n_cols {
        dens[r][c] = v
      }
    }
  } else {
    let rows_data = match json_as_arr(data_node) {
      Some(x) => x
      None => return None
    }
    for i = 0; i < n_rows; i = i + 1 {
      let vals = match json_as_arr(rows_data[i]) {
        Some(x) => x
        None => return None
      }
      for j = 0; j < n_cols; j = j + 1 {
        dens[i][j] = json_as_num(vals[j]).unwrap_or(0.0)
      }
    }
  }

  Some(BiomTable::{
    table_id,
    matrix_type,
    matrix_element_type,
    observation_ids,
    sample_ids,
    observation_metadata,
    sample_metadata,
    data: dens,
  })
}

///|
/// Serialize a `BiomTable` back to BIOM 1.0 JSON (compact).
pub fn BiomTable::to_json(self : BiomTable) -> String {
  let n_rows = self.observation_ids.length()
  let n_cols = self.sample_ids.length()

  let rows : Array[JsonValue] = []
  for i = 0; i < n_rows; i = i + 1 {
    let m : Map[String, JsonValue] = Map([], capacity=2)
    m.set("id", JsonStr(self.observation_ids[i]))
    let meta = self.observation_metadata[i]
    if meta.length() > 0 {
      m.set("metadata", JsonObj(meta))
    }
    rows.push(JsonObj(m))
  }
  let cols : Array[JsonValue] = []
  for j = 0; j < n_cols; j = j + 1 {
    let m : Map[String, JsonValue] = Map([], capacity=2)
    m.set("id", JsonStr(self.sample_ids[j]))
    let meta = self.sample_metadata[j]
    if meta.length() > 0 {
      m.set("metadata", JsonObj(meta))
    }
    cols.push(JsonObj(m))
  }

  let data_node : JsonValue =
    if self.matrix_type == "sparse" {
      let entries : Array[JsonValue] = []
      for i = 0; i < n_rows; i = i + 1 {
        for j = 0; j < n_cols; j = j + 1 {
          if self.data[i][j] != 0.0 {
            entries.push(
              JsonArr([
                JsonNum(i.to_double()),
                JsonNum(j.to_double()),
                JsonNum(self.data[i][j]),
              ]),
            )
          }
        }
      }
      JsonArr(entries)
    } else {
      let rows_d : Array[JsonValue] = []
      for i = 0; i < n_rows; i = i + 1 {
        let vals : Array[JsonValue] = []
        for j = 0; j < n_cols; j = j + 1 {
          vals.push(JsonNum(self.data[i][j]))
        }
        rows_d.push(JsonArr(vals))
      }
      JsonArr(rows_d)
    }

  let root : Map[String, JsonValue] = Map([], capacity=8)
  root.set("id", JsonStr(self.table_id))
  root.set("format", JsonStr("Biological Observation Matrix 1.0.0"))
  root.set("matrix_type", JsonStr(self.matrix_type))
  root.set("matrix_element_type", JsonStr(self.matrix_element_type))
  root.set(
    "shape",
    JsonArr([JsonNum(n_rows.to_double()), JsonNum(n_cols.to_double())]),
  )
  root.set("rows", JsonArr(rows))
  root.set("columns", JsonArr(cols))
  root.set("data", data_node)

  json_to_string(JsonObj(root))
}

// ============================================================
// Convenience constructors and accessors
// ============================================================

///|
/// Construct a `BiomTable` from raw arrays (no JSON needed).
/// Observation metadata and sample metadata default to empty maps.
pub fn BiomTable::new(
  observation_ids : Array[String],
  sample_ids : Array[String],
  data : Array[Array[Double]],
  table_id~ : String,
) -> BiomTable {
  let n_rows = observation_ids.length()
  let n_cols = sample_ids.length()
  let obs_meta : Array[Map[String, JsonValue]] = []
  let samp_meta : Array[Map[String, JsonValue]] = []
  for _ in 0..<n_rows {
    obs_meta.push(Map([], capacity=0))
  }
  for _ in 0..<n_cols {
    samp_meta.push(Map([], capacity=0))
  }
  BiomTable::{
    table_id,
    matrix_type: "dense",
    matrix_element_type: "float",
    observation_ids,
    sample_ids,
    observation_metadata: obs_meta,
    sample_metadata: samp_meta,
    data,
  }
}

///|
/// Number of observations (rows).
pub fn BiomTable::n_observations(self : BiomTable) -> Int {
  self.observation_ids.length()
}

///|
/// Number of samples (columns).
pub fn BiomTable::n_samples(self : BiomTable) -> Int {
  self.sample_ids.length()
}

///|
/// Get a cell value by observation index and sample index.
pub fn BiomTable::get_value(
  self : BiomTable,
  obs_idx : Int,
  samp_idx : Int
) -> Double? {
  if obs_idx >= 0 && obs_idx < self.data.length() {
    let row = self.data[obs_idx]
    if samp_idx >= 0 && samp_idx < row.length() {
      return Some(row[samp_idx])
    }
  }
  None
}

///|
/// Find the index of an observation ID, or -1 if not found.
pub fn BiomTable::observation_index(self : BiomTable, id : String) -> Int {
  for i in 0..<self.observation_ids.length() {
    if self.observation_ids[i] == id {
      return i
    }
  }
  -1
}

///|
/// Find the index of a sample ID, or -1 if not found.
pub fn BiomTable::sample_index(self : BiomTable, id : String) -> Int {
  for i in 0..<self.sample_ids.length() {
    if self.sample_ids[i] == id {
      return i
    }
  }
  -1
}

// ============================================================
// Row/column sums
// ============================================================

///|
/// Sum of each sample (column) — total counts per sample.
pub fn BiomTable::sample_sums(self : BiomTable) -> Array[Double] {
  let n_cols = self.sample_ids.length()
  let sums = Array::make(n_cols, 0.0)
  for i in 0..<self.data.length() {
    let row = self.data[i]
    for j in 0..<row.length() {
      sums[j] = sums[j] + row[j]
    }
  }
  sums
}

///|
/// Sum of each observation (row) — total counts per OTU/taxon.
pub fn BiomTable::observation_sums(self : BiomTable) -> Array[Double] {
  let sums : Array[Double] = []
  for i in 0..<self.data.length() {
    let mut s = 0.0
    let row = self.data[i]
    for j in 0..<row.length() {
      s = s + row[j]
    }
    sums.push(s)
  }
  sums
}

// ============================================================
// Transform & normalize
// ============================================================

///|
/// Apply a function to every cell in the data matrix (in-place mutation
/// of a copy). Returns a new BiomTable; the original is unchanged.
///
/// Example: `table.transform(fn(v) { if v > 0.0 { @math.log2(v + 1.0) } else { 0.0 } })`
pub fn BiomTable::transform(
  self : BiomTable,
  f : (Double, Int) -> Double
) -> BiomTable {
  let new_data : Array[Array[Double]] = []
  for i in 0..<self.data.length() {
    let row = self.data[i]
    let new_row : Array[Double] = []
    for j in 0..<row.length() {
      new_row.push(f(row[j], j))
    }
    new_data.push(new_row)
  }
  BiomTable::{
    table_id: self.table_id,
    matrix_type: self.matrix_type,
    matrix_element_type: self.matrix_element_type,
    observation_ids: self.observation_ids,
    sample_ids: self.sample_ids,
    observation_metadata: self.observation_metadata,
    sample_metadata: self.sample_metadata,
    data: new_data,
  }
}

///|
/// Normalize to relative abundance: divide each cell by its column (sample)
/// sum. Samples with zero total are left as zeros.
/// Equivalent to QIIME's `biom normalize-table`.
pub fn BiomTable::normalize_samples(self : BiomTable) -> BiomTable {
  let sums = self.sample_sums()
  self.transform(fn(v, j) {
    if j < sums.length() && sums[j] > 0.0 {
      v / sums[j]
    } else {
      0.0
    }
  })
}

// ============================================================
// Filtering
// ============================================================

///|
/// Filter observations (rows) by a list of IDs to keep.
/// Returns a new BiomTable containing only the specified observations.
pub fn BiomTable::filter_observations(
  self : BiomTable,
  keep_ids : Array[String]
) -> BiomTable {
  let new_ids : Array[String] = []
  let new_meta : Array[Map[String, JsonValue]] = []
  let new_data : Array[Array[Double]] = []
  for id in keep_ids {
    let idx = self.observation_index(id)
    if idx >= 0 {
      new_ids.push(id)
      new_meta.push(self.observation_metadata[idx])
      new_data.push(self.data[idx])
    }
  }
  BiomTable::{
    table_id: self.table_id,
    matrix_type: self.matrix_type,
    matrix_element_type: self.matrix_element_type,
    observation_ids: new_ids,
    sample_ids: self.sample_ids,
    observation_metadata: new_meta,
    sample_metadata: self.sample_metadata,
    data: new_data,
  }
}

///|
/// Filter samples (columns) by a list of IDs to keep.
/// Returns a new BiomTable containing only the specified samples.
pub fn BiomTable::filter_samples(
  self : BiomTable,
  keep_ids : Array[String]
) -> BiomTable {
  let keep_indices : Array[Int] = []
  for id in keep_ids {
    let idx = self.sample_index(id)
    if idx >= 0 {
      keep_indices.push(idx)
    }
  }
  let new_ids : Array[String] = []
  let new_meta : Array[Map[String, JsonValue]] = []
  for idx in keep_indices {
    new_ids.push(self.sample_ids[idx])
    new_meta.push(self.sample_metadata[idx])
  }
  let new_data : Array[Array[Double]] = []
  for i in 0..<self.data.length() {
    let row = self.data[i]
    let new_row : Array[Double] = []
    for idx in keep_indices {
      new_row.push(row[idx])
    }
    new_data.push(new_row)
  }
  BiomTable::{
    table_id: self.table_id,
    matrix_type: self.matrix_type,
    matrix_element_type: self.matrix_element_type,
    observation_ids: self.observation_ids,
    sample_ids: new_ids,
    observation_metadata: self.observation_metadata,
    sample_metadata: new_meta,
    data: new_data,
  }
}

// ============================================================
// TSV conversion (biom convert --to-tsv / --from-tsv)
// ============================================================

///|
/// Convert to tab-separated values format (TSV), matching `biom convert
/// --to-tsv`. First column header is `#OTU ID`, followed by sample IDs.
/// Each row is an observation with its counts. A final line `# Observation
/// count: N; Sample count: M` is appended (matching biom format).
pub fn BiomTable::to_tsv(self : BiomTable) -> String {
  let sb = StringBuilder()
  // Header
  sb.write_string("#OTU ID")
  for j in 0..<self.sample_ids.length() {
    sb.write_char('\t')
    sb.write_string(self.sample_ids[j])
  }
  sb.write_char('\n')
  // Data rows
  for i in 0..<self.observation_ids.length() {
    sb.write_string(self.observation_ids[i])
    let row = self.data[i]
    for j in 0..<row.length() {
      sb.write_char('\t')
      sb.write_string(row[j].to_string())
    }
    sb.write_char('\n')
  }
  // Footer comment
  sb.write_string("# Observation count: ")
  sb.write_string(self.observation_ids.length().to_string())
  sb.write_string("; Sample count: ")
  sb.write_string(self.sample_ids.length().to_string())
  sb.write_char('\n')
  sb.to_string()
}

///|
/// Parse a string to Double (simple parser for TSV cell values).
/// Returns 0.0 for non-numeric strings.
fn biom_parse_double(s : String) -> Double {
  let s = s.trim()
  if s.length() == 0 {
    return 0.0
  }
  let mut neg = false
  let mut start = 0
  if s[0] == 45 { // '-'
    neg = true
    start = 1
  } else if s[0] == 43 { // '+'
    start = 1
  }
  if start >= s.length() {
    return 0.0
  }
  let mut int_part = 0.0
  let mut frac_part = 0.0
  let mut frac_div = 1.0
  let mut in_frac = false
  let mut seen_digit = false
  let mut i = start
  while i < s.length() {
    let c = s.unsafe_get(i)
    if c >= 48 && c <= 57 { // '0'..'9'
      seen_digit = true
      let digit = (c - 48).to_int().to_double()
      if in_frac {
        frac_div = frac_div * 10.0
        frac_part = frac_part + digit / frac_div
      } else {
        int_part = int_part * 10.0 + digit
      }
    } else if c == 46 && !(in_frac) { // '.'
      in_frac = true
    } else if c == 101 || c == 69 { // 'e' or 'E'
      // Simple exponent handling
      break
    } else {
      break
    }
    i = i + 1
  }
  if !(seen_digit) {
    return 0.0
  }
  let result = int_part + frac_part
  if neg { -result } else { result }
}

///|
/// Parse a TSV table into a `BiomTable`. Inverse of `to_tsv`.
/// Expects first column header to be `#OTU ID` (case-insensitive prefix `#`).
/// Lines starting with `#` (except the header) are treated as comments and
/// ignored. Numeric cells are parsed as Double; non-numeric cells become 0.0.
pub fn parse_biom_tsv(text : String) -> BiomTable? {
  let lines = text.split("\n").to_array()
  if lines.length() == 0 {
    return None
  }
  // Find header line (starts with #OTU)
  let mut header_idx = -1
  let mut k = 0
  while k < lines.length() {
    let line = lines[k].trim()
    if line.has_prefix("#OTU") || line.has_prefix("#otu") {
      header_idx = k
      break
    }
    k = k + 1
  }
  if header_idx < 0 {
    // No explicit header — assume first non-comment, non-empty line is header
    k = 0
    while k < lines.length() {
      let line = lines[k].trim()
      if line.length() > 0 && !(line.has_prefix("#")) {
        header_idx = k
        break
      }
      k = k + 1
    }
  }
  if header_idx < 0 {
    return None
  }
  // Parse header
  let header_parts = lines[header_idx].trim().split("\t").to_array()
  let sample_ids : Array[String] = []
  let mut j = 1
  while j < header_parts.length() {
    sample_ids.push(header_parts[j].to_owned())
    j = j + 1
  }
  let n_cols = sample_ids.length()
  // Parse data rows
  let observation_ids : Array[String] = []
  let data : Array[Array[Double]] = []
  let mut i = header_idx + 1
  while i < lines.length() {
    let line = lines[i].trim()
    if line.length() == 0 {
      i = i + 1
      continue
    }
    if line.has_prefix("#") {
      i = i + 1
      continue // comment/footer
    }
    let parts = line.split("\t").to_array()
    if parts.length() < 1 {
      i = i + 1
      continue
    }
    observation_ids.push(parts[0].to_owned())
    let row : Array[Double] = Array::make(n_cols, 0.0)
    let mut jj = 1
    while jj < parts.length() {
      if jj - 1 < n_cols {
        let v = biom_parse_double(parts[jj].trim().to_owned())
        row[jj - 1] = v
      }
      jj = jj + 1
    }
    data.push(row)
    i = i + 1
  }
  Some(BiomTable::new(observation_ids, sample_ids, data, table_id="from_tsv"))
}

// ============================================================
// Merge two tables
// ============================================================

///|
/// Merge two BiomTables: union of observations and samples.
/// Missing cells (present in one table but not the other) are filled
/// with 0.0. Observation/sample metadata is taken from whichever table
/// has it (table `a` takes priority).
pub fn BiomTable::merge(
  self : BiomTable,
  other : BiomTable
) -> BiomTable {
  // Union of sample IDs
  let all_sample_ids : Array[String] = []
  let samp_seen : Map[String, Int] = {}
  for id in self.sample_ids {
    if !(samp_seen.contains(id)) {
      samp_seen[id] = all_sample_ids.length()
      all_sample_ids.push(id)
    }
  }
  for id in other.sample_ids {
    if !(samp_seen.contains(id)) {
      samp_seen[id] = all_sample_ids.length()
      all_sample_ids.push(id)
    }
  }
  // Union of observation IDs
  let all_obs_ids : Array[String] = []
  let obs_seen : Map[String, Int] = {}
  for id in self.observation_ids {
    if !(obs_seen.contains(id)) {
      obs_seen[id] = all_obs_ids.length()
      all_obs_ids.push(id)
    }
  }
  for id in other.observation_ids {
    if !(obs_seen.contains(id)) {
      obs_seen[id] = all_obs_ids.length()
      all_obs_ids.push(id)
    }
  }
  // Build merged data matrix
  let n_rows = all_obs_ids.length()
  let n_cols = all_sample_ids.length()
  let merged_data : Array[Array[Double]] = Array(capacity=n_rows)
  for _ in 0..<n_rows {
    merged_data.push(Array::make(n_cols, 0.0))
  }
  // Fill from self
  for i in 0..<self.observation_ids.length() {
    let row_idx = obs_seen[self.observation_ids[i]]
    let row = self.data[i]
    for j in 0..<self.sample_ids.length() {
      let col_idx = samp_seen[self.sample_ids[j]]
      merged_data[row_idx][col_idx] = merged_data[row_idx][col_idx] + row[j]
    }
  }
  // Fill from other
  for i in 0..<other.observation_ids.length() {
    let row_idx = obs_seen[other.observation_ids[i]]
    let row = other.data[i]
    for j in 0..<other.sample_ids.length() {
      let col_idx = samp_seen[other.sample_ids[j]]
      merged_data[row_idx][col_idx] = merged_data[row_idx][col_idx] + row[j]
    }
  }
  // Metadata: empty for merged
  let obs_meta : Array[Map[String, JsonValue]] = []
  let samp_meta : Array[Map[String, JsonValue]] = []
  for _ in 0..<n_rows {
    obs_meta.push(Map([], capacity=0))
  }
  for _ in 0..<n_cols {
    samp_meta.push(Map([], capacity=0))
  }
  BiomTable::{
    table_id: self.table_id + "_merged",
    matrix_type: "dense",
    matrix_element_type: "float",
    observation_ids: all_obs_ids,
    sample_ids: all_sample_ids,
    observation_metadata: obs_meta,
    sample_metadata: samp_meta,
    data: merged_data,
  }
}
