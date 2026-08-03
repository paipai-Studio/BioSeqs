// INSDC (GenBank/EMBL) feature table parser and writer, mirroring
// BioPython's Bio.SeqIO.InsdcIO module.
//
// The INSDC Feature Table format is shared by GenBank and EMBL. This module
// implements feature-table parsing/writing plus higher-level record
// parsing/writing that complements the existing genbank_io.mbt and
// embl_io.mbt parsers.
//
// Feature table layout (GenBank):
//   FEATURES    Location/Qualifiers
//       CDS       join(1031..1186,1188..1770)
//                 /gene="rpoB"
//                 /product="DNA-directed RNA polymerase subunit beta"
//                 /translation="MSITK..."
//       gene      1031..1770
//                 /gene="rpoB"
//                 /locus_tag="ECK_1234"
//
// Feature table layout (EMBL, with "FT   " prefix on every line):
//   FH   Key             Location/Qualifiers
//   FH
//   FT   CDS             join(1031..1186,1188..1770)
//   FT                   /gene="rpoB"
//   FT   gene            1031..1770
//   FT                   /gene="rpoB"

///|
/// INSDC feature table feature record.
pub struct InsdcFeature {
  feature_type : String
  location_str : String
  qualifiers : Map[String, Array[String]]
} derive(Eq)

///|
/// INSDC record combining header metadata, features, and sequence.
pub struct InsdcRecord {
  locus_name : String
  sequence_length : Int
  molecule_type : String
  topology : String
  division : String
  date : String
  definition : String
  accession : String
  version : String
  organism : String
  features : Array[InsdcFeature]
  sequence : String
} derive(Eq)

// ---------------------------------------------------------------------------
// Internal helpers
// ---------------------------------------------------------------------------

///|
/// Push a qualifier value into a mutable qualifiers map, appending to the
/// existing list of values for that key when present.
fn insdc_push_qualifier(
  quals : Map[String, Array[String]],
  key : String,
  value : String,
) -> Unit {
  match quals.get(key) {
    Some(arr) => {
      let new_arr = arr.copy()
      new_arr.push(value)
      quals[key] = new_arr
    }
    None => quals[key] = [value]
  }
}

///|
/// Copy a qualifiers map (deep-copying each value array).
fn insdc_copy_quals(
  src : Map[String, Array[String]],
) -> Map[String, Array[String]] {
  let dst : Map[String, Array[String]] = Map([], capacity=src.length())
  for k, v in src {
    dst[k] = v.copy()
  }
  dst
}

///|
/// Split a string on ASCII whitespace, returning non-empty tokens.
fn insdc_split_whitespace(s : String) -> Array[String] {
  let result : Array[String] = Array::new()
  let len = s.length()
  let sp : UInt16 = ' '.to_int().to_uint16()
  let tab : UInt16 = '\t'.to_int().to_uint16()
  let mut i = 0
  while i < len {
    while i < len {
      let c = s.unsafe_get(i)
      if c == sp || c == tab {
        i = i + 1
      } else {
        break
      }
    }
    if i >= len {
      break
    }
    let start = i
    while i < len {
      let c = s.unsafe_get(i)
      if c == sp || c == tab {
        break
      }
      i = i + 1
    }
    result.push(s[start:i].to_owned())
  }
  result
}

///|
/// Split a string on a single delimiter character code.
fn insdc_split_char(s : String, delim : Int) -> Array[String] {
  let result : Array[String] = Array::new()
  let n = s.length()
  let mut start = 0
  let mut i = 0
  while i < n {
    if s.unsafe_get(i).to_int() == delim {
      result.push(s[start:i].to_owned())
      start = i + 1
    }
    i = i + 1
  }
  if start <= n {
    result.push(s[start:n].to_owned())
  }
  result
}

///|
/// Strip leading and trailing ASCII whitespace from a string.
fn insdc_trim(s : String) -> String {
  let n = s.length()
  let sp : UInt16 = ' '.to_int().to_uint16()
  let tab : UInt16 = '\t'.to_int().to_uint16()
  let cr : UInt16 = '\r'.to_int().to_uint16()
  let nl : UInt16 = '\n'.to_int().to_uint16()
  let mut start = 0
  while start < n {
    let c = s.unsafe_get(start)
    if c == sp || c == tab || c == cr || c == nl {
      start = start + 1
    } else {
      break
    }
  }
  let mut end = n
  while end > start {
    let c = s.unsafe_get(end - 1)
    if c == sp || c == tab || c == cr || c == nl {
      end = end - 1
    } else {
      break
    }
  }
  if start >= end {
    ""
  } else {
    s[start:end].to_owned()
  }
}

// ---------------------------------------------------------------------------
// Feature table parsing
// ---------------------------------------------------------------------------

///|
/// Parse INSDC feature table lines and return a list of features.
///
/// Accepts both GenBank-style lines (5-space indented feature type) and
/// EMBL-style lines (each line prefixed with `FT   `). The `FEATURES` header
/// line and `FH` lines are tolerated and skipped.
pub fn parse_insdc_feature_table(lines : Array[String]) -> Array[InsdcFeature] {
  let features : Array[InsdcFeature] = Array::new()
  let n = lines.length()
  let sp : UInt16 = ' '.to_int().to_uint16()
  let cr : UInt16 = '\r'.to_int().to_uint16()
  let slash : UInt16 = '/'.to_int().to_uint16()
  let eq : UInt16 = '='.to_int().to_uint16()
  let quote : UInt16 = '"'.to_int().to_uint16()
  let dot : UInt16 = '.'.to_int().to_uint16()
  let cap_f : UInt16 = 'F'.to_int().to_uint16()
  let cap_t : UInt16 = 'T'.to_int().to_uint16()

  let mut cur_type = ""
  let cur_loc = StringBuilder::new()
  let mut cur_quals : Map[String, Array[String]] = Map([], capacity=0)
  let mut cur_qual_key = ""
  let cur_qual_val = StringBuilder::new()
  let mut in_quote = false
  let mut have_feature = false

  for i = 0; i < n; i = i + 1 {
    let raw = lines[i]
    let rl = raw.length()
    if rl == 0 {
      continue
    }
    let mut le = rl
    if raw.unsafe_get(rl - 1) == cr {
      le = rl - 1
    }
    if le == 0 {
      continue
    }

    // Strip "FT   " prefix from EMBL lines when present.
    let line = if le >= 5 &&
      raw.unsafe_get(0) == cap_f &&
      raw.unsafe_get(1) == cap_t &&
      raw.unsafe_get(2) == sp &&
      raw.unsafe_get(3) == sp &&
      raw.unsafe_get(4) == sp {
      raw[5:le].to_owned()
    } else {
      raw[0:le].to_owned()
    }
    let ll = line.length()
    if ll == 0 {
      continue
    }

    // Count leading spaces.
    let mut leading = 0
    while leading < ll && line.unsafe_get(leading) == sp {
      leading = leading + 1
    }
    if leading == ll {
      // whitespace-only line
      continue
    }

    if leading >= 21 {
      // Qualifier line (or continuation).
      let content = line[21:ll].to_owned()
      let cl = content.length()
      if cl == 0 {
        continue
      }
      let first = content.unsafe_get(0)

      if in_quote {
        // Continuation of a quoted qualifier value.
        cur_qual_val.write_string(content)
        if content.unsafe_get(cl - 1) == quote {
          // Closing quote reached: strip it and finalize.
          let full = cur_qual_val.to_string()
          let fl = full.length()
          let val_clean = if fl > 0 && full.unsafe_get(fl - 1) == quote {
            full[0:fl - 1].to_owned()
          } else {
            full
          }
          insdc_push_qualifier(cur_quals, cur_qual_key, val_clean)
          in_quote = false
          cur_qual_key = ""
          cur_qual_val.reset()
        }
      } else if first == slash {
        // Start of a new qualifier: "/key=value" or "/key".
        let rest = content[1:cl].to_owned()
        let rl2 = rest.length()
        let mut eq_pos = -1
        for k = 0; k < rl2; k = k + 1 {
          if rest.unsafe_get(k) == eq {
            eq_pos = k
            break
          }
        }
        if eq_pos == -1 {
          // "/key" with no value.
          insdc_push_qualifier(cur_quals, rest, "")
        } else {
          let key = rest[0:eq_pos].to_owned()
          let val_part = rest[eq_pos + 1:rl2].to_owned()
          let vp_len = val_part.length()
          if vp_len > 0 && val_part.unsafe_get(0) == quote {
            // Quoted value: may span multiple lines.
            let inner = val_part[1:vp_len].to_owned()
            let inner_len = inner.length()
            if inner_len > 0 && inner.unsafe_get(inner_len - 1) == quote {
              // Single-line quoted value.
              let val_clean = inner[0:inner_len - 1].to_owned()
              insdc_push_qualifier(cur_quals, key, val_clean)
            } else {
              // Multi-line quoted value: keep accumulating.
              cur_qual_key = key
              cur_qual_val.reset()
              cur_qual_val.write_string(inner)
              in_quote = true
            }
          } else {
            // Unquoted value (e.g., "/kb=123").
            insdc_push_qualifier(cur_quals, key, val_part)
          }
        }
      } else if have_feature && cur_qual_key.length() == 0 && !in_quote {
        // Continuation of a multi-line location (rare in practice).
        cur_loc.write_string(content)
      }
      // else: stray content; ignore.
    } else if leading >= 5 && leading < 21 {
      // New feature line: "     TYPE    LOCATION".
      if have_feature {
        features.push(InsdcFeature::{
          feature_type: cur_type,
          location_str: cur_loc.to_string(),
          qualifiers: insdc_copy_quals(cur_quals),
        })
        cur_quals = Map([], capacity=0)
        cur_loc.reset()
        cur_qual_key = ""
        cur_qual_val.reset()
        in_quote = false
      }

      // Parse "TYPE LOCATION" starting at column 5.
      let rest = line[5:ll].to_owned()
      let rl2 = rest.length()
      let mut type_end = 0
      while type_end < rl2 && rest.unsafe_get(type_end) != sp {
        type_end = type_end + 1
      }
      cur_type = rest[0:type_end].to_owned()
      let mut loc_start = type_end
      while loc_start < rl2 && rest.unsafe_get(loc_start) == sp {
        loc_start = loc_start + 1
      }
      if loc_start < rl2 {
        cur_loc.write_string(rest[loc_start:rl2].to_owned())
      }
      have_feature = true
    } else {
      // Could be the "FEATURES" header line or "FH" line; ignore.
      ()
    }
  }

  // Flush the last feature in progress.
  if have_feature {
    features.push(InsdcFeature::{
      feature_type: cur_type,
      location_str: cur_loc.to_string(),
      qualifiers: insdc_copy_quals(cur_quals),
    })
  }
  // Reference dot/quote locals so unused warnings stay quiet.
  let _ = dot
  features
}

// ---------------------------------------------------------------------------
// Location parsing
// ---------------------------------------------------------------------------

///|
/// Parse and normalize an INSDC location string.
///
/// Returns the location string with all internal whitespace removed so that
/// callers see a canonical form regardless of how the original file wrapped
/// long expressions. Recognized operators (`complement`, `join`, `order`)
/// and fuzzy markers (`<`, `>`, `?`) are preserved verbatim.
pub fn parse_insdc_location(loc_str : String) -> String {
  let sp : UInt16 = ' '.to_int().to_uint16()
  let tab : UInt16 = '\t'.to_int().to_uint16()
  let nl : UInt16 = '\n'.to_int().to_uint16()
  let cr : UInt16 = '\r'.to_int().to_uint16()
  let buf = StringBuilder::new()
  let n = loc_str.length()
  for i = 0; i < n; i = i + 1 {
    let c = loc_str.unsafe_get(i)
    if c != sp && c != tab && c != nl && c != cr {
      buf.write_char(c.unsafe_to_char())
    }
  }
  buf.to_string()
}

// ---------------------------------------------------------------------------
// Feature formatting
// ---------------------------------------------------------------------------

///|
/// Format an INSDC feature as a list of feature-table lines (GenBank style).
///
/// `indent` is the number of leading spaces before the feature type (5 for
/// GenBank, 5 for EMBL where the leading "FT" is added separately). The
/// feature type is left-justified in a 16-column field so that the location
/// starts at column 21, and qualifier lines begin at column 21. Long quoted
/// values (e.g., `/translation`) are wrapped at column 79 with the opening
/// quote on the first line and the closing quote on the last line.
pub fn format_insdc_feature(
  feature : InsdcFeature,
  indent : Int,
) -> Array[String] {
  let out : Array[String] = Array::new()

  // First line: indent + feature_type (padded to 16) + location.
  let first = StringBuilder::new()
  for _j = 0; _j < indent; _j = _j + 1 {
    first.write_char(' ')
  }
  first.write_string(feature.feature_type)
  let target_col = indent + 16
  let cur_len = indent + feature.feature_type.length()
  let pad = target_col - cur_len
  for _j = 0; pad > 0 && _j < pad; _j = _j + 1 {
    first.write_char(' ')
  }
  first.write_string(feature.location_str)
  out.push(first.to_string())

  // Qualifier lines, each at column 21 indent.
  for key, values in feature.qualifiers {
    for val in values {
      let q_lines = format_insdc_qualifier(key, val)
      for ql in q_lines {
        out.push(ql)
      }
    }
  }
  out
}

///|
/// Format a single qualifier with line wrapping for long values.
/// - Empty value: `                     /key`
/// - Short value: `                     /key="value"`
/// - Long value: opening quote on first line, closing quote on last line,
///   continuation lines indented to column 21.
fn format_insdc_qualifier(key : String, value : String) -> Array[String] {
  let out : Array[String] = Array::new()
  let indent_str = "                     " // 21 spaces

  if value.length() == 0 {
    out.push(indent_str + "/" + key)
    return out
  }

  let prefix = indent_str + "/" + key + "=\""
  let prefix_len = prefix.length()
  // Wrap so each line fits within 79 columns.
  let max_first = if 79 - prefix_len > 0 { 79 - prefix_len } else { 1 }
  let max_cont = if 79 - 21 > 0 { 79 - 21 } else { 1 }
  let total = value.length()

  if total <= max_first {
    out.push(prefix + value + "\"")
    return out
  }

  // Multi-line: opening quote on first line, closing quote on last.
  let mut pos = 0
  let first_end = if pos + max_first < total { pos + max_first } else { total }
  out.push(prefix + value[pos:first_end].to_owned())
  pos = first_end
  while pos < total {
    let end = if pos + max_cont < total { pos + max_cont } else { total }
    let chunk = value[pos:end].to_owned()
    if end >= total {
      out.push(indent_str + chunk + "\"")
    } else {
      out.push(indent_str + chunk)
    }
    pos = end
  }
  out
}

// ---------------------------------------------------------------------------
// Qualifier extraction helpers
// ---------------------------------------------------------------------------

///|
/// Extract the first value of a given qualifier key from a feature, or None.
pub fn extract_feature_qualifier(
  feature : InsdcFeature,
  key : String,
) -> String? {
  match feature.qualifiers.get(key) {
    Some(arr) => if arr.length() > 0 { Some(arr[0]) } else { None }
    None => None
  }
}

///|
/// Convenience helper: extract the `/translation` value of a CDS feature.
pub fn extract_translation(feature : InsdcFeature) -> String? {
  extract_feature_qualifier(feature, "translation")
}

///|
/// Convenience helper: extract the `/gene` value of a feature.
pub fn extract_gene_name(feature : InsdcFeature) -> String? {
  extract_feature_qualifier(feature, "gene")
}

// ---------------------------------------------------------------------------
// GenBank record parsing/writing
// ---------------------------------------------------------------------------

///|
/// Parse a single GenBank record into an InsdcRecord.
///
/// Extracts the LOCUS line metadata (name, length, molecule type, topology,
/// division, date), the DEFINITION, ACCESSION, VERSION, ORGANISM fields,
/// the FEATURES table, and the ORIGIN sequence. The record terminator is
/// `//`.
pub fn parse_genbank_record(content : String) -> InsdcRecord {
  let lines = content.split("\n").to_array()
  let n = lines.length()

  let mut locus_name = ""
  let mut sequence_length = 0
  let mut molecule_type = ""
  let mut topology = "linear"
  let mut division = "UNK"
  let mut date = ""
  let mut definition = ""
  let mut accession = ""
  let mut version = ""
  let mut organism = ""
  let features : Array[InsdcFeature] = Array::new()
  let seq_buf = StringBuilder::new()

  let feature_lines : Array[String] = Array::new()
  let mut in_features = false
  let mut in_origin = false

  let sp : UInt16 = ' '.to_int().to_uint16()
  let cr : UInt16 = '\r'.to_int().to_uint16()
  let dot : UInt16 = '.'.to_int().to_uint16()
  let slash : UInt16 = '/'.to_int().to_uint16()

  let mut i = 0
  while i < n {
    let raw = lines[i]
    let rl = raw.length()
    let mut le = rl
    if rl > 0 && raw.unsafe_get(rl - 1) == cr {
      le = rl - 1
    }
    let line = raw[0:le].to_owned()
    let ll = line.length()

    // Record terminator.
    if ll >= 2 &&
      line.unsafe_get(0) == slash &&
      line.unsafe_get(1) == slash {
      i = i + 1
      break
    }

    // LOCUS line.
    if ll >= 5 && line.has_prefix("LOCUS") {
      let rest = insdc_trim(line[5:ll].to_owned())
      let tokens = insdc_split_whitespace(rest)
      if tokens.length() >= 1 {
        locus_name = tokens[0]
      }
      if tokens.length() >= 2 {
        sequence_length = parse_int(tokens[1])
      }
      let mut bp_idx = -1
      for t = 0; t < tokens.length(); t = t + 1 {
        if tokens[t] == "bp" || tokens[t] == "bp," {
          bp_idx = t
          break
        }
      }
      if bp_idx >= 0 {
        let mut idx = bp_idx + 1
        if idx < tokens.length() {
          molecule_type = tokens[idx]
          idx = idx + 1
        }
        // Topology is optional ("linear" or "circular").
        if idx < tokens.length() &&
          (tokens[idx] == "linear" || tokens[idx] == "circular") {
          topology = tokens[idx]
          idx = idx + 1
        }
        if idx < tokens.length() {
          division = tokens[idx]
          idx = idx + 1
        }
        if idx < tokens.length() {
          let date_buf = StringBuilder::new()
          for d = idx; d < tokens.length(); d = d + 1 {
            if d > idx {
              date_buf.write_char(' ')
            }
            date_buf.write_string(tokens[d])
          }
          date = date_buf.to_string()
        }
      }
      i = i + 1
      continue
    }

    // DEFINITION line (and continuations indented by spaces).
    if ll >= 10 && line.has_prefix("DEFINITION") {
      let def_buf = StringBuilder::new()
      if ll > 10 {
        def_buf.write_string(insdc_trim(line[10:ll].to_owned()))
      }
      i = i + 1
      while i < n {
        let cont = lines[i]
        let cl = cont.length()
        let mut ct = cl
        if cl > 0 && cont.unsafe_get(cl - 1) == cr {
          ct = cl - 1
        }
        if ct == 0 || cont.unsafe_get(0) != sp {
          break
        }
        def_buf.write_char(' ')
        def_buf.write_string(insdc_trim(cont[0:ct].to_owned()))
        i = i + 1
      }
      definition = insdc_trim(def_buf.to_string())
      if definition.length() > 0 &&
        definition.unsafe_get(definition.length() - 1) == dot {
        definition = definition[0:definition.length() - 1].to_owned()
      }
      continue
    }

    // ACCESSION line.
    if ll >= 9 && line.has_prefix("ACCESSION") {
      if ll > 9 {
        let acc = insdc_trim(line[9:ll].to_owned())
        let acc_tokens = insdc_split_whitespace(acc)
        if acc_tokens.length() > 0 {
          accession = acc_tokens[0]
        }
      }
      i = i + 1
      continue
    }

    // VERSION line.
    if ll >= 7 && line.has_prefix("VERSION") {
      if ll > 7 {
        let ver = insdc_trim(line[7:ll].to_owned())
        let ver_tokens = insdc_split_whitespace(ver)
        if ver_tokens.length() > 0 {
          version = ver_tokens[0]
        }
      }
      i = i + 1
      continue
    }

    // SOURCE/ORGANISM lines.
    if ll >= 8 && line.has_prefix("ORGANISM") {
      if ll > 8 {
        organism = insdc_trim(line[8:ll].to_owned())
      }
      i = i + 1
      continue
    }

    // FEATURES header line.
    if ll >= 8 && line.has_prefix("FEATURES") {
      in_features = true
      in_origin = false
      i = i + 1
      continue
    }
    if in_features {
      // End of features section.
      if ll >= 6 && line.has_prefix("ORIGIN") {
        in_features = false
        in_origin = true
        i = i + 1
        continue
      }
      // Collect feature-table lines (those that start with at least 5 spaces).
      if ll >= 5 && line.unsafe_get(0) == sp {
        feature_lines.push(line)
      }
      i = i + 1
      continue
    }

    // ORIGIN section start.
    if ll >= 6 && line.has_prefix("ORIGIN") {
      in_origin = true
      i = i + 1
      continue
    }
    if in_origin {
      // Extract sequence letters from the line.
      for j = 0; j < ll; j = j + 1 {
        let c = line.unsafe_get(j)
        let ci = c.to_int()
        if (ci >= 'A'.to_int() && ci <= 'Z'.to_int()) ||
          (ci >= 'a'.to_int() && ci <= 'z'.to_int()) {
          seq_buf.write_char(c.unsafe_to_char())
        }
      }
    }

    i = i + 1
  }

  // Parse the collected feature-table lines.
  if feature_lines.length() > 0 {
    let parsed = parse_insdc_feature_table(feature_lines)
    for f in parsed {
      features.push(f)
    }
  }

  InsdcRecord::{
    locus_name,
    sequence_length,
    molecule_type,
    topology,
    division,
    date,
    definition,
    accession,
    version,
    organism,
    features,
    sequence: seq_buf.to_string(),
  }
}

///|
/// Write an InsdcRecord to GenBank format text.
///
/// Emits a LOCUS line, DEFINITION, ACCESSION, VERSION, a FEATURES section
/// (one feature per record), and an ORIGIN section with the sequence wrapped
/// at 60 bases per line in groups of 10. The record ends with `//`.
pub fn write_genbank_record(record : InsdcRecord) -> String {
  let buf = StringBuilder::new()

  // --- LOCUS line ---
  // LOCUS       <name>      <len> bp    <mol>    <topo>   <div> <date>
  buf.write_string("LOCUS       ")
  buf.write_string(record.locus_name)
  let pad1 = 28 - (12 + record.locus_name.length())
  for _j = 0; pad1 > 0 && _j < pad1; _j = _j + 1 {
    buf.write_char(' ')
  }
  let len_str = record.sequence_length.to_string() + " bp"
  let pad_len = 11 - len_str.length()
  for _j = 0; pad_len > 0 && _j < pad_len; _j = _j + 1 {
    buf.write_char(' ')
  }
  buf.write_string(len_str)
  buf.write_string(" ")
  buf.write_string(if record.molecule_type.length() > 0 {
    record.molecule_type
  } else {
    "DNA"
  })
  buf.write_string(" ")
  buf.write_string(if record.topology.length() > 0 {
    record.topology
  } else {
    "linear"
  })
  buf.write_string("  ")
  buf.write_string(if record.division.length() > 0 {
    record.division
  } else {
    "UNK"
  })
  buf.write_string(" ")
  buf.write_string(if record.date.length() > 0 {
    record.date
  } else {
    "01-JAN-2024"
  })
  buf.write_char('\n')

  // --- DEFINITION ---
  buf.write_string("DEFINITION  ")
  buf.write_string(record.definition)
  buf.write_char('.')
  buf.write_char('\n')

  // --- ACCESSION ---
  buf.write_string("ACCESSION   ")
  buf.write_string(record.accession)
  buf.write_char('\n')

  // --- VERSION ---
  buf.write_string("VERSION     ")
  buf.write_string(record.version)
  buf.write_char('\n')

  // --- FEATURES ---
  buf.write_string("FEATURES         Location/Qualifiers\n")
  for f in record.features {
    let lines = format_insdc_feature(f, 5)
    for l in lines {
      buf.write_string(l)
      buf.write_char('\n')
    }
  }

  // --- ORIGIN ---
  buf.write_string("ORIGIN      \n")
  if record.sequence.length() > 0 {
    buf.write_string(insdc_build_genbank_sequence(record.sequence))
  }

  // End record.
  buf.write_string("//\n")
  buf.to_string()
}

///|
/// Build the ORIGIN sequence section (GenBank style):
/// "      1 atgcatgcat gcatgcatgc ..." with a 1-indexed position label
/// right-aligned in 9 columns and 6 groups of 10 bases per line.
fn insdc_build_genbank_sequence(seq : String) -> String {
  let buf = StringBuilder::new()
  let total = seq.length()
  let mut i = 0
  let mut line_start = 0
  while i < total {
    let pos_label = (line_start + 1).to_string()
    let pad_pos = 9 - pos_label.length()
    for _j = 0; pad_pos > 0 && _j < pad_pos; _j = _j + 1 {
      buf.write_char(' ')
    }
    buf.write_string(pos_label)
    buf.write_char(' ')
    let mut group = 0
    while group < 6 && i < total {
      let group_end = if i + 10 < total { i + 10 } else { total }
      buf.write_string(seq[i:group_end].to_owned())
      i = group_end
      buf.write_char(' ')
      group = group + 1
    }
    buf.write_char('\n')
    line_start = i
  }
  buf.to_string()
}

// ---------------------------------------------------------------------------
// EMBL record parsing/writing
// ---------------------------------------------------------------------------

///|
/// Parse a single EMBL record into an InsdcRecord.
///
/// Extracts the ID line (locus name, sequence version, topology, molecule
/// type, division, length), AC (accession), DE (definition), OS (organism),
/// DT (date), the FH/FT feature table, and the SQ sequence section.
pub fn parse_embl_record(content : String) -> InsdcRecord {
  let lines = content.split("\n").to_array()
  let n = lines.length()

  let mut locus_name = ""
  let mut sequence_length = 0
  let mut molecule_type = ""
  let mut topology = "linear"
  let mut division = "UNK"
  let mut date = ""
  let mut definition = ""
  let mut accession = ""
  let mut version = ""
  let mut organism = ""
  let features : Array[InsdcFeature] = Array::new()
  let seq_buf = StringBuilder::new()

  let feature_lines : Array[String] = Array::new()
  let mut in_ft = false
  let mut in_sq = false

  let sp : UInt16 = ' '.to_int().to_uint16()
  let cr : UInt16 = '\r'.to_int().to_uint16()
  let dot : UInt16 = '.'.to_int().to_uint16()
  let slash : UInt16 = '/'.to_int().to_uint16()

  let mut sv_num = ""

  let mut i = 0
  while i < n {
    let raw = lines[i]
    let rl = raw.length()
    let mut le = rl
    if rl > 0 && raw.unsafe_get(rl - 1) == cr {
      le = rl - 1
    }
    let line = raw[0:le].to_owned()
    let ll = line.length()

    // Record terminator.
    if ll >= 2 &&
      line.unsafe_get(0) == slash &&
      line.unsafe_get(1) == slash {
      i = i + 1
      break
    }

    // Skip XX separator lines.
    if ll >= 2 &&
      line.unsafe_get(0) == 'X'.to_int().to_uint16() &&
      line.unsafe_get(1) == 'X'.to_int().to_uint16() &&
      (ll == 2 || line.unsafe_get(2) == sp) {
      i = i + 1
      continue
    }

    // ID line: "ID   name; SV ver; topology; mol_type; class; division; length BP."
    if ll >= 2 &&
      line.unsafe_get(0) == 'I'.to_int().to_uint16() &&
      line.unsafe_get(1) == 'D'.to_int().to_uint16() &&
      (ll == 2 || line.unsafe_get(2) == sp) {
      let rest = if ll > 2 { insdc_trim(line[2:ll].to_owned()) } else { "" }
      let parts = insdc_split_char(rest, ';'.to_int())
      // parts[0]=name parts[1]=SV ver parts[2]=topology parts[3]=mol_type
      // parts[4]=class parts[5]=division parts[6]=length BP.
      if parts.length() >= 1 {
        locus_name = insdc_trim(parts[0])
      }
      if parts.length() >= 2 {
        let sv_str = insdc_trim(parts[1])
        let sv_tokens = insdc_split_whitespace(sv_str)
        if sv_tokens.length() >= 2 && sv_tokens[0] == "SV" {
          sv_num = sv_tokens[1]
        }
      }
      if parts.length() >= 3 {
        let topo = insdc_trim(parts[2])
        if topo == "linear" || topo == "circular" {
          topology = topo
        }
      }
      if parts.length() >= 4 {
        molecule_type = insdc_trim(parts[3])
      }
      if parts.length() >= 6 {
        division = insdc_trim(parts[5])
      }
      if parts.length() >= 7 {
        let len_str = insdc_trim(parts[6])
        let len_tokens = insdc_split_whitespace(len_str)
        if len_tokens.length() >= 1 {
          sequence_length = parse_int(len_tokens[0])
        }
      }
      i = i + 1
      continue
    }

    // DT line: take the first date token.
    if ll >= 2 &&
      line.unsafe_get(0) == 'D'.to_int().to_uint16() &&
      line.unsafe_get(1) == 'T'.to_int().to_uint16() &&
      (ll == 2 || line.unsafe_get(2) == sp) {
      if date.length() == 0 && ll > 2 {
        let rest = insdc_trim(line[2:ll].to_owned())
        let dt_tokens = insdc_split_whitespace(rest)
        if dt_tokens.length() > 0 {
          date = dt_tokens[0]
        }
      }
      i = i + 1
      continue
    }

    // AC line: "AC   ACC1; ACC2;"
    if ll >= 2 &&
      line.unsafe_get(0) == 'A'.to_int().to_uint16() &&
      line.unsafe_get(1) == 'C'.to_int().to_uint16() &&
      (ll == 2 || line.unsafe_get(2) == sp) {
      if ll > 2 {
        let rest = insdc_trim(line[2:ll].to_owned())
        let acc_parts = insdc_split_char(rest, ';'.to_int())
        if acc_parts.length() > 0 {
          let first = insdc_trim(acc_parts[0])
          if first.length() > 0 {
            accession = first
          }
        }
      }
      i = i + 1
      continue
    }

    // DE line (definition, may span multiple lines).
    if ll >= 2 &&
      line.unsafe_get(0) == 'D'.to_int().to_uint16() &&
      line.unsafe_get(1) == 'E'.to_int().to_uint16() &&
      (ll == 2 || line.unsafe_get(2) == sp) {
      let def_buf = StringBuilder::new()
      if ll > 2 {
        def_buf.write_string(insdc_trim(line[2:ll].to_owned()))
      }
      i = i + 1
      while i < n {
        let cont = lines[i]
        let cl = cont.length()
        let mut ct = cl
        if cl > 0 && cont.unsafe_get(cl - 1) == cr {
          ct = cl - 1
        }
        if ct < 2 ||
          !(cont.unsafe_get(0) == 'D'.to_int().to_uint16() &&
            cont.unsafe_get(1) == 'E'.to_int().to_uint16() &&
            (ct == 2 || cont.unsafe_get(2) == sp)) {
          break
        }
        def_buf.write_char(' ')
        if ct > 2 {
          def_buf.write_string(insdc_trim(cont[2:ct].to_owned()))
        }
        i = i + 1
      }
      definition = insdc_trim(def_buf.to_string())
      if definition.length() > 0 &&
        definition.unsafe_get(definition.length() - 1) == dot {
        definition = definition[0:definition.length() - 1].to_owned()
      }
      continue
    }

    // OS line (organism, may span multiple OS lines).
    if ll >= 2 &&
      line.unsafe_get(0) == 'O'.to_int().to_uint16() &&
      line.unsafe_get(1) == 'S'.to_int().to_uint16() &&
      (ll == 2 || line.unsafe_get(2) == sp) {
      let org_buf = StringBuilder::new()
      if ll > 2 {
        org_buf.write_string(insdc_trim(line[2:ll].to_owned()))
      }
      i = i + 1
      while i < n {
        let cont = lines[i]
        let cl = cont.length()
        let mut ct = cl
        if cl > 0 && cont.unsafe_get(cl - 1) == cr {
          ct = cl - 1
        }
        if ct < 2 ||
          !(cont.unsafe_get(0) == 'O'.to_int().to_uint16() &&
            cont.unsafe_get(1) == 'S'.to_int().to_uint16() &&
            (ct == 2 || cont.unsafe_get(2) == sp)) {
          break
        }
        org_buf.write_char(' ')
        if ct > 2 {
          org_buf.write_string(insdc_trim(cont[2:ct].to_owned()))
        }
        i = i + 1
      }
      organism = insdc_trim(org_buf.to_string())
      continue
    }

    // FH header lines are skipped.
    if ll >= 2 &&
      line.unsafe_get(0) == 'F'.to_int().to_uint16() &&
      line.unsafe_get(1) == 'H'.to_int().to_uint16() &&
      (ll == 2 || line.unsafe_get(2) == sp) {
      i = i + 1
      continue
    }

    // FT feature-table lines: collect them verbatim.
    if ll >= 2 &&
      line.unsafe_get(0) == 'F'.to_int().to_uint16() &&
      line.unsafe_get(1) == 'T'.to_int().to_uint16() &&
      (ll == 2 || line.unsafe_get(2) == sp) {
      in_ft = true
      in_sq = false
      feature_lines.push(line)
      i = i + 1
      continue
    }
    if in_ft {
      // A non-FT line ends the feature table.
      if ll < 2 ||
        !(line.unsafe_get(0) == 'F'.to_int().to_uint16() &&
          line.unsafe_get(1) == 'T'.to_int().to_uint16()) {
        in_ft = false
        // Fall through to other section checks below.
      }
    }

    // SQ line marks the start of the sequence section.
    if ll >= 2 &&
      line.unsafe_get(0) == 'S'.to_int().to_uint16() &&
      line.unsafe_get(1) == 'Q'.to_int().to_uint16() &&
      (ll == 2 || line.unsafe_get(2) == sp) {
      in_sq = true
      in_ft = false
      i = i + 1
      continue
    }

    if in_sq {
      // Extract sequence letters from the line.
      for j = 0; j < ll; j = j + 1 {
        let c = line.unsafe_get(j)
        let ci = c.to_int()
        if (ci >= 'A'.to_int() && ci <= 'Z'.to_int()) ||
          (ci >= 'a'.to_int() && ci <= 'z'.to_int()) {
          seq_buf.write_char(c.unsafe_to_char())
        }
      }
    }

    i = i + 1
  }

  // Compose the version string from accession and SV number.
  if version.length() == 0 {
    if accession.length() > 0 && sv_num.length() > 0 {
      version = accession + "." + sv_num
    } else if accession.length() > 0 {
      version = accession
    }
  }

  // Parse the collected feature-table lines.
  if feature_lines.length() > 0 {
    let parsed = parse_insdc_feature_table(feature_lines)
    for f in parsed {
      features.push(f)
    }
  }

  InsdcRecord::{
    locus_name,
    sequence_length,
    molecule_type,
    topology,
    division,
    date,
    definition,
    accession,
    version,
    organism,
    features,
    sequence: seq_buf.to_string(),
  }
}

///|
/// Write an InsdcRecord to EMBL format text.
///
/// Emits an ID line, AC, DE, OS, an FH/FT feature table, and an SQ section
/// with the sequence wrapped at 60 bases per line in groups of 10. The
/// record ends with `//`.
pub fn write_embl_record(record : InsdcRecord) -> String {
  let buf = StringBuilder::new()

  // --- ID line ---
  // ID   <name>; SV <ver>; <topology>; <mol_type>; STD; <division>; <len> BP.
  buf.write_string("ID   ")
  buf.write_string(record.locus_name)
  buf.write_string("; SV ")
  // Extract SV number from version (after the '.') if available.
  let sv_num = if record.version.length() > 0 {
    let dot_pos = insdc_find_char(record.version, '.'.to_int())
    if dot_pos >= 0 && dot_pos + 1 < record.version.length() {
      record.version[dot_pos + 1:].to_owned()
    } else {
      "1"
    }
  } else {
    "1"
  }
  buf.write_string(sv_num)
  buf.write_string("; ")
  buf.write_string(if record.topology.length() > 0 {
    record.topology
  } else {
    "linear"
  })
  buf.write_string("; ")
  buf.write_string(if record.molecule_type.length() > 0 {
    record.molecule_type
  } else {
    "DNA"
  })
  buf.write_string("; STD; ")
  buf.write_string(if record.division.length() > 0 {
    record.division
  } else {
    "UNK"
  })
  buf.write_string("; ")
  buf.write_string(record.sequence_length.to_string())
  buf.write_string(" BP.")
  buf.write_char('\n')
  buf.write_string("XX\n")

  // --- AC ---
  buf.write_string("AC   ")
  buf.write_string(record.accession)
  buf.write_string(";\n")
  buf.write_string("XX\n")

  // --- DE ---
  buf.write_string("DE   ")
  buf.write_string(record.definition)
  buf.write_char('\n')
  buf.write_string("XX\n")

  // --- OS ---
  if record.organism.length() > 0 {
    buf.write_string("OS   ")
    buf.write_string(record.organism)
    buf.write_char('\n')
    buf.write_string("XX\n")
  }

  // --- FH/FT feature table ---
  buf.write_string("FH   Key             Location/Qualifiers\n")
  buf.write_string("FH\n")
  for f in record.features {
    let lines = format_insdc_feature(f, 5)
    for l in lines {
      // Convert the GenBank-style 5-space indent into "FT   " by replacing
      // the first two spaces with the literal "FT".
      let l_to_write = if l.length() >= 2 {
        "FT" + l[2:].to_owned()
      } else {
        l
      }
      buf.write_string(l_to_write)
      buf.write_char('\n')
    }
  }
  buf.write_string("XX\n")

  // --- SQ ---
  buf.write_string("SQ   Sequence ")
  buf.write_string(record.sequence_length.to_string())
  buf.write_string(" BP; 0 A; 0 C; 0 G; 0 T; 0 other;\n")
  if record.sequence.length() > 0 {
    buf.write_string(insdc_build_embl_sequence(record.sequence))
  }

  // End record.
  buf.write_string("//\n")
  buf.to_string()
}

///|
/// Build the EMBL sequence section: 5-space indent, 6 groups of 10 bases.
fn insdc_build_embl_sequence(seq : String) -> String {
  let buf = StringBuilder::new()
  let total = seq.length()
  let mut i = 0
  while i < total {
    buf.write_string("     ")
    let mut group = 0
    while group < 6 && i < total {
      let group_end = if i + 10 < total { i + 10 } else { total }
      buf.write_string(seq[i:group_end].to_owned())
      i = group_end
      buf.write_char(' ')
      group = group + 1
    }
    buf.write_char('\n')
  }
  buf.to_string()
}

///|
/// Find the first index of a character code in a string, or -1 if absent.
fn insdc_find_char(s : String, c : Int) -> Int {
  let n = s.length()
  let mut i = 0
  while i < n {
    if s.unsafe_get(i).to_int() == c {
      return i
    }
    i = i + 1
  }
  -1
}
