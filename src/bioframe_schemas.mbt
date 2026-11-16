///|
/// bioframe_schemas.mbt — Schema registry and tabular I/O for Bioframe.
///
/// Bioframe equivalent: `bioframe.SCHEMAS` dict + `bioframe.read_table(schema=)`.
///
/// Bioframe's killer feature: give a file path + a schema name, get back
/// a Pandas DataFrame with columns already named.  Here we can't use
/// Pandas, but we can:
///
///   1. enumerate the most common schema names as String constants
///   2. map each schema name → array of canonical column names (via fn)
///   3. supply a `bf_read_table` function that parses tab-delimited
///      text into rows + a column-name array, applying the schema to
///      fill column names when the file itself doesn't have a header
///   4. supply helpers to extract BedInterval arrays from parsed rows
///
/// All coordinates remain in the canonical 0-based half-open convention
/// that BioSeqs uses internally.

// ================================================================
// 1. Schema NAME constants (primitive String — MoonBit allows these)
// ================================================================

const BF_SCHEMA_BED3 : String = "bed3"

///|
const BF_SCHEMA_BED4 : String = "bed4"

///|
const BF_SCHEMA_BED5 : String = "bed5"

///|
const BF_SCHEMA_BED6 : String = "bed6"

///|
const BF_SCHEMA_BED7 : String = "bed7"

///|
const BF_SCHEMA_BED8 : String = "bed8"

///|
const BF_SCHEMA_BED9 : String = "bed9"

///|
const BF_SCHEMA_BED12 : String = "bed12"

///|
const BF_SCHEMA_NARROWPEAK : String = "narrowPeak"

///|
const BF_SCHEMA_BROADPEAK : String = "broadPeak"

///|
const BF_SCHEMA_GAPPEDPEAK : String = "gappedPeak"

///|
const BF_SCHEMA_PAIRS : String = "pairs"

///|
const BF_SCHEMA_PAIRS_PX : String = "pairs_px"

///|
const BF_SCHEMA_TAD : String = "tad"

///|
const BF_SCHEMA_REFSEQ : String = "refSeq"

///|
const BF_SCHEMA_SIMPLEREPEAT : String = "simpleRepeat"

///|
const BF_SCHEMA_JASPAR : String = "jaspar"

///|
const BF_SCHEMA_BEDMETHYL : String = "bedMethyl"

// ================================================================
// 2. Column-list FNs  (MoonBit const disallows Array types, so we
//                     expose each schema's columns through a fn)
// ================================================================

///|
fn bf_cols_bed3() -> Array[String] {
  ["chrom", "start", "end"]
}

///|
fn bf_cols_bed4() -> Array[String] {
  ["chrom", "start", "end", "name"]
}

///|
fn bf_cols_bed5() -> Array[String] {
  ["chrom", "start", "end", "name", "score"]
}

///|
fn bf_cols_bed6() -> Array[String] {
  ["chrom", "start", "end", "name", "score", "strand"]
}

///|
fn bf_cols_bed7() -> Array[String] {
  ["chrom", "start", "end", "name", "score", "strand", "thickStart", "thickEnd"]
}

///|
fn bf_cols_bed8() -> Array[String] {
  [
    "chrom", "start", "end", "name", "score", "strand", "thickStart", "thickEnd",
    "itemRgb",
  ]
}

///|
fn bf_cols_bed9() -> Array[String] {
  [
    "chrom", "start", "end", "name", "score", "strand", "thickStart", "thickEnd",
    "itemRgb", "blockCount", "blockSizes", "blockStarts",
  ]
}

///|
fn bf_cols_bed12() -> Array[String] {
  [
    "chrom", "start", "end", "name", "score", "strand", "thickStart", "thickEnd",
    "itemRgb", "blockCount", "blockSizes", "blockStarts",
  ]
}

///|
fn bf_cols_narrowpeak() -> Array[String] {
  [
    "chrom", "start", "end", "name", "score", "strand", "fc", "negLog10p", "negLog10q",
    "relSummit",
  ]
}

///|
fn bf_cols_broadpeak() -> Array[String] {
  [
    "chrom", "start", "end", "name", "score", "strand", "fc", "negLog10p", "negLog10q",
  ]
}

///|
fn bf_cols_gappedpeak() -> Array[String] {
  [
    "chrom", "start", "end", "name", "score", "strand", "fc", "negLog10p", "negLog10q",
    "relSummit", "blockCount", "blockSizes", "blockStarts", "blockAq", "blockP",
    "blockQ",
  ]
}

///|
fn bf_cols_pairs() -> Array[String] {
  ["chrom1", "pos1", "chrom2", "pos2"]
}

///|
fn bf_cols_pairs_px() -> Array[String] {
  ["read_id", "chrom1", "pos1", "chrom2", "pos2", "strand1", "strand2"]
}

///|
fn bf_cols_tad() -> Array[String] {
  ["chrom", "start", "end", "name"]
}

///|
fn bf_cols_refseq() -> Array[String] {
  [
    "bin", "name", "chrom", "strand", "txStart", "txEnd", "cdsStart", "cdsEnd", "exonCount",
    "exonStarts", "exonEnds", "score", "name2", "cdsStartStat", "cdsEndStat", "exonFrames",
  ]
}

///|
fn bf_cols_simplerepeat() -> Array[String] {
  [
    "bin", "chrom", "chromStart", "chromEnd", "name", "period", "copyNum", "consensusSize",
    "perMatch", "perIndel", "score", "A", "C", "G", "T", "entropy", "sequence",
  ]
}

///|
fn bf_cols_jaspar() -> Array[String] {
  ["chrom", "start", "end", "name", "score", "strand"]
}

///|
fn bf_cols_bedmethyl() -> Array[String] {
  [
    "chrom", "start", "end", "name", "score", "strand", "thickStart", "thickEnd",
    "itemRgb", "coverage", "methylationPct",
  ]
}

// ================================================================
// 3. Public resolver
// ================================================================

///|
/// Look up the canonical column list for a known schema name.
/// Returns None if the schema is unknown.
pub fn bf_schema_cols(name : String) -> Array[String]? {
  match name {
    BF_SCHEMA_BED3 => Some(bf_cols_bed3())
    BF_SCHEMA_BED4 => Some(bf_cols_bed4())
    BF_SCHEMA_BED5 => Some(bf_cols_bed5())
    BF_SCHEMA_BED6 => Some(bf_cols_bed6())
    BF_SCHEMA_BED7 => Some(bf_cols_bed7())
    BF_SCHEMA_BED8 => Some(bf_cols_bed8())
    BF_SCHEMA_BED9 => Some(bf_cols_bed9())
    BF_SCHEMA_BED12 => Some(bf_cols_bed12())
    BF_SCHEMA_NARROWPEAK => Some(bf_cols_narrowpeak())
    BF_SCHEMA_BROADPEAK => Some(bf_cols_broadpeak())
    BF_SCHEMA_GAPPEDPEAK => Some(bf_cols_gappedpeak())
    BF_SCHEMA_PAIRS => Some(bf_cols_pairs())
    BF_SCHEMA_PAIRS_PX => Some(bf_cols_pairs_px())
    BF_SCHEMA_TAD => Some(bf_cols_tad())
    BF_SCHEMA_REFSEQ => Some(bf_cols_refseq())
    BF_SCHEMA_SIMPLEREPEAT => Some(bf_cols_simplerepeat())
    BF_SCHEMA_JASPAR => Some(bf_cols_jaspar())
    BF_SCHEMA_BEDMETHYL => Some(bf_cols_bedmethyl())
    _ => None
  }
}

///|
/// List all known schema names.
pub fn bf_schemas_available() -> Array[String] {
  [
    BF_SCHEMA_BED3,
    BF_SCHEMA_BED4,
    BF_SCHEMA_BED5,
    BF_SCHEMA_BED6,
    BF_SCHEMA_BED7,
    BF_SCHEMA_BED8,
    BF_SCHEMA_BED9,
    BF_SCHEMA_BED12,
    BF_SCHEMA_NARROWPEAK,
    BF_SCHEMA_BROADPEAK,
    BF_SCHEMA_GAPPEDPEAK,
    BF_SCHEMA_PAIRS,
    BF_SCHEMA_PAIRS_PX,
    BF_SCHEMA_TAD,
    BF_SCHEMA_REFSEQ,
    BF_SCHEMA_SIMPLEREPEAT,
    BF_SCHEMA_JASPAR,
    BF_SCHEMA_BEDMETHYL,
  ]
}

// ================================================================
// 4. Row + I/O types
// ================================================================

///|
/// A parsed row from a tabular schema-driven file.
pub struct BfTableRow {
  fields : Array[String]
}

///|
/// Get the raw field by column index.
pub fn BfTableRow::field(self : BfTableRow, idx : Int32) -> String? {
  if idx >= 0 && idx < self.fields.length() {
    Some(self.fields[idx])
  } else {
    None
  }
}

// ================================================================
// 5. bf_read_table
// ================================================================

///|
/// Parse tab-delimited text into (header_cols, rows).
///
/// - `schema`: if provided, column names come from the schema and the
///   first data line is NOT treated as a header.  Otherwise: first
///   non-comment line is used as header when it looks like text, or
///   default column names are invented if it looks numeric.
pub fn bf_read_table(
  text : String,
  schema~ : String? = None,
) -> (Array[String], Array[BfTableRow]) {
  let mut header : Array[String] = []
  let rows : Array[BfTableRow] = Array(capacity=0)
  let mut current = StringBuilder()
  let mut is_first_line = true
  for c in text {
    if c == '\n' {
      let line = current.to_string()
      current = StringBuilder()
      let trimmed = bf_trim(line)
      if trimmed.length() == 0 {
        is_first_line = false
        continue
      }
      if trimmed.has_prefix("#") ||
        trimmed.has_prefix("track ") ||
        trimmed.has_prefix("browser ") ||
        trimmed.has_prefix("##") {
        is_first_line = false
        continue
      }
      let cols = bf_sch_split_tab(line)
      if is_first_line {
        match schema {
          Some(sch) => {
            match bf_schema_cols(sch) {
              Some(h) => header = h
              None => header = cols
            }
            rows.push(BfTableRow::{ fields: cols, })
          }
          None => {
            let looks_like_header = !bf_all_numeric(cols)
            if looks_like_header {
              header = cols
            } else {
              header = bf_default_header(cols.length())
              rows.push(BfTableRow::{ fields: cols, })
            }
          }
        }
        is_first_line = false
      } else {
        rows.push(BfTableRow::{ fields: cols, })
      }
    } else {
      current.write_char(c)
    }
  }
  (header, rows)
}

///|
/// Get a column from a row by column name.
pub fn bf_row_get(
  header : Array[String],
  row : BfTableRow,
  col : String,
) -> String? {
  let mut i : Int32 = 0
  while i < header.length() {
    if header[i] == col {
      return row.field(i)
    }
    i = i + 1
  }
  None
}

///|
/// Extract BedInterval array from parsed rows using schema header names.
pub fn bf_rows_to_intervals(
  header : Array[String],
  rows : Array[BfTableRow],
) -> Array[BedInterval] {
  let out : Array[BedInterval] = Array(capacity=0)
  let mut ci : Int32 = 0
  let mut si : Int32 = 1
  let mut ei : Int32 = 2
  let mut i : Int32 = 0
  while i < header.length() {
    let h = header[i].to_lower()
    if h == "chrom" || h == "chr" {
      ci = i
    } else if h == "start" || h == "chromstart" || h == "pos1" {
      si = i
    } else if h == "end" || h == "chromend" || h == "pos2" {
      ei = i
    }
    i = i + 1
  }
  for r in rows {
    let chrom = match r.field(ci) {
      Some(c) => c
      None => continue
    }
    let s_txt = match r.field(si) {
      Some(v) => v
      None => continue
    }
    let e_txt = match r.field(ei) {
      Some(v) => v
      None => continue
    }
    let s = parse_int(s_txt)
    let e = parse_int(e_txt)
    if chrom.length() > 0 {
      out.push(BedInterval::new(chrom, s, e))
    }
  }
  out
}

// ================================================================
// 6. Internal helpers
// ================================================================

///|
fn bf_trim(s : String) -> String {
  let mut start = 0
  let mut end = s.length()
  while start < end {
    let c = s.unsafe_get(start).to_int()
    if c == ' '.to_int() ||
      c == '\t'.to_int() ||
      c == '\n'.to_int() ||
      c == '\r'.to_int() {
      start = start + 1
    } else {
      break
    }
  }
  while end > start {
    let c = s.unsafe_get(end - 1).to_int()
    if c == ' '.to_int() ||
      c == '\t'.to_int() ||
      c == '\n'.to_int() ||
      c == '\r'.to_int() {
      end = end - 1
    } else {
      break
    }
  }
  s[start:end].to_owned()
}

///|
/// Tab split helper, unique name (assemblies.mbt has bf_split_tab too).
fn bf_sch_split_tab(line : String) -> Array[String] {
  let out : Array[String] = []
  let mut sb = StringBuilder()
  for c in line {
    if c == '\t' {
      out.push(sb.to_string())
      sb = StringBuilder()
    } else {
      sb.write_char(c)
    }
  }
  out.push(sb.to_string())
  out
}

///|
/// Heuristic: true when every element looks like a number
/// (digits, optional dot/scientific, optional +/-).
fn bf_all_numeric(cols : Array[String]) -> Bool {
  if cols.length() == 0 {
    return false
  }
  for c in cols {
    if c.length() == 0 {
      return false
    }
    let mut i : Int32 = 0
    let first = c.unsafe_get(0).to_int()
    if first == '-'.to_int() || first == '+'.to_int() {
      i = 1
    }
    let mut saw_digit = false
    let mut saw_dot = false
    while i < c.length() {
      let ch = c.unsafe_get(i).to_int()
      if ch >= '0'.to_int() && ch <= '9'.to_int() {
        saw_digit = true
      } else if ch == '.'.to_int() && !saw_dot {
        saw_dot = true
      } else {
        return false
      }
      i = i + 1
    }
    if !saw_digit {
      return false
    }
  }
  true
}

///|
fn bf_default_header(n : Int32) -> Array[String] {
  let out : Array[String] = []
  let mut i : Int32 = 0
  while i < n {
    out.push("col" + i.to_string())
    i = i + 1
  }
  out
}
