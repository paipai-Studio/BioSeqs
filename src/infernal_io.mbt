///|
/// Bio.SearchIO.InfernalIO support for Infernal cmscan/cmsearch output.
///
/// Infernal searches RNA sequence databases with covariance models. Its
/// terminology calls the model the query and the sequence the target. This
/// module follows Biopython SearchIO and exposes model coordinates as query
/// coordinates and sequence coordinates as hit coordinates.

///|
pub suberror InfernalError {
  InfernalError(String)
}

///|
pub enum InfernalTabularFormat {
  Format1
  Format2
  Format3
} derive(Eq, Debug)
pub extend InfernalTabularFormat with Eq::{equal, not_equal}
pub extend InfernalTabularFormat with @moonbitlang/core/debug.Debug::{to_repr}

///|
pub struct InfernalFragment {
  query_start : Int
  query_end : Int
  hit_start : Int
  hit_end : Int
  hit_strand : Strand
  query_sequence : String
  hit_sequence : String
  consensus_structure : String
  noncanonical_pairs : String
  similarity : String
  posterior_probability : String
  model_omission_before : Int
  sequence_omission_before : Int
} derive(Debug)
pub extend InfernalFragment with @moonbitlang/core/debug.Debug::{to_repr}

///|
pub struct InfernalHSP {
  fragments : Array[InfernalFragment]
  model : String
  truncated : String
  pipeline_pass : Int
  gc : Double
  bias : Double
  bitscore : Double
  evalue : Double
  is_included : Bool
  query_end_type : String
  hit_end_type : String
  average_accuracy : Double
  overlap : String
  any_index : Int?
  any_overlap_fraction : Double?
  reciprocal_overlap_fraction : Double?
  winner_index : Int?
  winner_overlap_fraction : Double?
  winner_reciprocal_fraction : Double?
} derive(Debug)
pub extend InfernalHSP with @moonbitlang/core/debug.Debug::{to_repr}

///|
pub struct InfernalHit {
  id : String
  accession : String
  description : String
  model_length : Int
  hsps : Array[InfernalHSP]
} derive(Debug)
pub extend InfernalHit with @moonbitlang/core/debug.Debug::{to_repr}

///|
pub struct InfernalQueryResult {
  id : String
  accession : String
  description : String
  clan : String
  sequence_length : Int
  program : String
  version : String
  target : String
  hits : Array[InfernalHit]
} derive(Debug)
pub extend InfernalQueryResult with @moonbitlang/core/debug.Debug::{to_repr}

///|
priv struct InfernalTabRow {
  query_id : String
  query_accession : String
  query_clan : String
  query_length : Int
  hit_id : String
  hit_accession : String
  hit_description : String
  hit_model_length : Int
  hsp : InfernalHSP
}

///|
priv struct InfernalMarker {
  start : Int
  end : Int
  model_size : Int
  sequence_size : Int
}

///|
fn infernal_split_whitespace(s : String) -> Array[String] {
  let values : Array[String] = []
  let mut start = 0
  let mut in_token = false
  for index in 0..<s.length() {
    let value = s.unsafe_get(index).to_int()
    let whitespace = value == ' '.to_int() ||
      value == '\t'.to_int() ||
      value == '\r'.to_int()
    if whitespace {
      if in_token {
        values.push(s[start:index].to_owned())
        in_token = false
      }
    } else if !in_token {
      start = index
      in_token = true
    }
  }
  if in_token {
    values.push(s[start:s.length()].to_owned())
  }
  values
}

///|
fn infernal_join(values : Array[String], start : Int) -> String {
  let mut result = ""
  for index in start..<values.length() {
    if index > start {
      result = result + " "
    }
    result = result + values[index]
  }
  result
}

///|
fn infernal_parse_int(
  value : String,
  label : String,
) -> Int raise InfernalError {
  if value.length() == 0 {
    raise InfernalError("missing integer field " + label)
  }
  let mut start = 0
  if value.unsafe_get(0).to_int() == '-'.to_int() {
    if value.length() == 1 {
      raise InfernalError("invalid integer in " + label + ": " + value)
    }
    start = 1
  }
  for index in start..<value.length() {
    let digit = value.unsafe_get(index).to_int()
    if digit < '0'.to_int() || digit > '9'.to_int() {
      raise InfernalError("invalid integer in " + label + ": " + value)
    }
  }
  parse_int(value)
}

///|
fn infernal_parse_double(
  value : String,
  label : String,
) -> Double raise InfernalError {
  match parse_double(value) {
    Some(number) =>
      if number.is_nan() {
        raise InfernalError("NaN is not valid in " + label)
      } else {
        number
      }
    None => raise InfernalError("invalid number in " + label + ": " + value)
  }
}

///|
fn infernal_optional_int(value : String) -> Int? raise InfernalError {
  if value == "-" || value == "\"" {
    None
  } else {
    Some(infernal_parse_int(value, "overlap index"))
  }
}

///|
fn infernal_optional_double(value : String) -> Double? raise InfernalError {
  if value == "-" || value == "\"" {
    None
  } else {
    Some(infernal_parse_double(value, "overlap fraction"))
  }
}

///|
fn infernal_digit_count(value : String) -> Int {
  let mut result = 0
  for index in 0..<value.length() {
    let digit = value.unsafe_get(index).to_int()
    if digit >= '0'.to_int() && digit <= '9'.to_int() {
      result = result * 10 + digit - '0'.to_int()
    }
  }
  result
}

///|
fn infernal_is_integer(value : String) -> Bool {
  if value.length() == 0 {
    return false
  }
  let mut start = 0
  if value.unsafe_get(0).to_int() == '-'.to_int() {
    start = 1
  }
  if start == value.length() {
    return false
  }
  for index in start..<value.length() {
    let digit = value.unsafe_get(index).to_int()
    if digit < '0'.to_int() || digit > '9'.to_int() {
      return false
    }
  }
  true
}

///|
fn infernal_coordinate_fragment(
  model_from : Int,
  model_to : Int,
  sequence_from : Int,
  sequence_to : Int,
  strand : String,
) -> InfernalFragment raise InfernalError {
  if model_from < 1 || model_to < 1 || sequence_from < 1 || sequence_to < 1 {
    raise InfernalError("Infernal coordinates must be positive")
  }
  let (query_start, query_end) = if model_from <= model_to {
    (model_from - 1, model_to)
  } else {
    (model_to - 1, model_from)
  }
  let (hit_start, hit_end, hit_strand) = match strand {
    "+" => (sequence_from - 1, sequence_to, Strand::Plus)
    "-" => (sequence_to - 1, sequence_from, Strand::Minus)
    _ => raise InfernalError("invalid Infernal strand: " + strand)
  }
  InfernalFragment::{
    query_start,
    query_end,
    hit_start,
    hit_end,
    hit_strand,
    query_sequence: "",
    hit_sequence: "",
    consensus_structure: "",
    noncanonical_pairs: "",
    similarity: "",
    posterior_probability: "",
    model_omission_before: 0,
    sequence_omission_before: 0,
  }
}

///|
fn infernal_detect_tabular_format(
  content : String,
) -> InfernalTabularFormat raise InfernalError {
  let lines = split_by_char(content, '\n'.to_int())
  for line in lines {
    let trimmed = line.trim().to_owned()
    if trimmed.has_prefix("#") && trimmed.contains("target name") {
      if trimmed.contains("clan name") {
        return InfernalTabularFormat::Format2
      }
      if trimmed.contains("mdl len") {
        return InfernalTabularFormat::Format3
      }
      return InfernalTabularFormat::Format1
    }
  }
  raise InfernalError(
    "cannot determine Infernal tabular format; supply format=1, 2, or 3",
  )
}

///|
pub fn infernal_tabular_format(
  content : String,
) -> InfernalTabularFormat raise InfernalError {
  infernal_detect_tabular_format(content)
}

///|
fn infernal_format_number(
  content : String,
  requested : Int,
) -> Int raise InfernalError {
  if requested != 0 {
    if requested < 1 || requested > 3 {
      raise InfernalError("Infernal tabular format must be 1, 2, or 3")
    }
    return requested
  }
  match infernal_detect_tabular_format(content) {
    Format1 => 1
    Format2 => 2
    Format3 => 3
  }
}

///|
fn infernal_parse_tab_row(
  values : Array[String],
  format : Int,
) -> InfernalTabRow raise InfernalError {
  let minimum = if format == 1 { 18 } else if format == 2 { 29 } else { 20 }
  if values.length() < minimum {
    raise InfernalError(
      "Infernal tabular format " +
      format.to_string() +
      " expects at least " +
      minimum.to_string() +
      " fields, found " +
      values.length().to_string(),
    )
  }
  let mut offset = 0
  let mut clan = ""
  let mut model_length = 0
  let mut sequence_length = 0
  let mut overlap = ""
  let mut any_index : Int? = None
  let mut any_fraction : Double? = None
  let mut reciprocal_fraction : Double? = None
  let mut winner_index : Int? = None
  let mut winner_fraction : Double? = None
  let mut winner_reciprocal : Double? = None
  let description_start = if format == 1 {
    17
  } else if format == 2 {
    offset = 1
    clan = values[5]
    overlap = values[19]
    any_index = infernal_optional_int(values[20])
    any_fraction = infernal_optional_double(values[21])
    reciprocal_fraction = infernal_optional_double(values[22])
    winner_index = infernal_optional_int(values[23])
    winner_fraction = infernal_optional_double(values[24])
    winner_reciprocal = infernal_optional_double(values[25])
    model_length = infernal_parse_int(values[26], "mdl len")
    sequence_length = infernal_parse_int(values[27], "seq len")
    28
  } else {
    model_length = infernal_parse_int(values[17], "mdl len")
    sequence_length = infernal_parse_int(values[18], "seq len")
    19
  }
  let target_name = values[offset]
  let target_accession = values[offset + 1]
  let query_name = values[offset + 2]
  let query_accession = values[offset + 3]
  let model_offset = if format == 2 { 6 } else { 4 }
  let model = values[model_offset]
  let model_from = infernal_parse_int(values[model_offset + 1], "mdl from")
  let model_to = infernal_parse_int(values[model_offset + 2], "mdl to")
  let sequence_from = infernal_parse_int(values[model_offset + 3], "seq from")
  let sequence_to = infernal_parse_int(values[model_offset + 4], "seq to")
  let strand = values[model_offset + 5]
  let fragment = infernal_coordinate_fragment(
    model_from, model_to, sequence_from, sequence_to, strand,
  )
  let hsp = InfernalHSP::{
    fragments: [fragment],
    model,
    truncated: values[model_offset + 6],
    pipeline_pass: infernal_parse_int(values[model_offset + 7], "pipeline pass"),
    gc: infernal_parse_double(values[model_offset + 8], "gc"),
    bias: infernal_parse_double(values[model_offset + 9], "bias"),
    bitscore: infernal_parse_double(values[model_offset + 10], "score"),
    evalue: infernal_parse_double(values[model_offset + 11], "E-value"),
    is_included: values[model_offset + 12] == "!",
    query_end_type: "",
    hit_end_type: "",
    average_accuracy: 0.0,
    overlap,
    any_index,
    any_overlap_fraction: any_fraction,
    reciprocal_overlap_fraction: reciprocal_fraction,
    winner_index,
    winner_overlap_fraction: winner_fraction,
    winner_reciprocal_fraction: winner_reciprocal,
  }
  InfernalTabRow::{
    query_id: query_name,
    query_accession,
    query_clan: clan,
    query_length: sequence_length,
    hit_id: target_name,
    hit_accession: target_accession,
    hit_description: infernal_join(values, description_start),
    hit_model_length: model_length,
    hsp,
  }
}

///|
fn infernal_build_results(
  rows : Array[InfernalTabRow],
  program : String,
  version : String,
  target : String,
) -> Array[InfernalQueryResult] {
  let query_order : Array[String] = []
  let by_query : Map[String, Array[InfernalTabRow]] = Map([], capacity=0)
  for row in rows {
    match by_query.get(row.query_id) {
      Some(values) => values.push(row)
      None => {
        query_order.push(row.query_id)
        by_query.set(row.query_id, [row])
      }
    }
  }
  let results : Array[InfernalQueryResult] = []
  for query_id in query_order {
    let query_rows = by_query.get(query_id).unwrap()
    let hit_order : Array[String] = []
    let by_hit : Map[String, Array[InfernalTabRow]] = Map([], capacity=0)
    for row in query_rows {
      match by_hit.get(row.hit_id) {
        Some(values) => values.push(row)
        None => {
          hit_order.push(row.hit_id)
          by_hit.set(row.hit_id, [row])
        }
      }
    }
    let hits : Array[InfernalHit] = []
    for hit_id in hit_order {
      let hit_rows = by_hit.get(hit_id).unwrap()
      let hsps : Array[InfernalHSP] = []
      for row in hit_rows {
        hsps.push(row.hsp)
      }
      let first = hit_rows[0]
      hits.push(InfernalHit::{
        id: hit_id,
        accession: first.hit_accession,
        description: first.hit_description,
        model_length: first.hit_model_length,
        hsps,
      })
    }
    let first = query_rows[0]
    results.push(InfernalQueryResult::{
      id: query_id,
      accession: first.query_accession,
      description: "",
      clan: first.query_clan,
      sequence_length: first.query_length,
      program,
      version,
      target,
      hits,
    })
  }
  results
}

///|
/// Parse Infernal `--tblout` formats 1, 2, and 3.
///
/// The format is inferred from the header unless `format` is supplied.
/// Coordinates are normalized to zero-based, half-open intervals. Reverse
/// strand sequence coordinates are reordered while retaining `Strand::Minus`.
pub fn infernal_parse_tabular(
  content : String,
  format? : Int = 0,
) -> Array[InfernalQueryResult] raise InfernalError {
  let selected = infernal_format_number(content, format)
  let rows : Array[InfernalTabRow] = []
  let lines = split_by_char(content, '\n'.to_int())
  for line in lines {
    let trimmed = line.trim().to_owned()
    if trimmed.length() == 0 || trimmed.has_prefix("#") {
      continue
    }
    rows.push(
      infernal_parse_tab_row(infernal_split_whitespace(trimmed), selected),
    )
  }
  infernal_build_results(rows, "", "", "")
}

///|
fn infernal_marker_at(sequence : String, start : Int) -> (Int, Int)? {
  if start + 4 >= sequence.length() ||
    sequence.unsafe_get(start).to_int() != '*'.to_int() ||
    sequence.unsafe_get(start + 1).to_int() != '['.to_int() {
    return None
  }
  let mut end = start + 2
  let mut size = 0
  let mut has_digit = false
  while end < sequence.length() {
    let value = sequence.unsafe_get(end).to_int()
    if value == ']'.to_int() {
      if end + 1 < sequence.length() &&
        sequence.unsafe_get(end + 1).to_int() == '*'.to_int() &&
        has_digit {
        return Some((end + 2, size))
      }
      return None
    }
    if value == ' '.to_int() {
      end = end + 1
      continue
    }
    if value < '0'.to_int() || value > '9'.to_int() {
      return None
    }
    has_digit = true
    size = size * 10 + value - '0'.to_int()
    end = end + 1
  }
  None
}

///|
fn infernal_find_markers(
  model_sequence : String,
  hit_sequence : String,
) -> Array[InfernalMarker] {
  let markers : Array[InfernalMarker] = []
  let length = model_sequence.length().max(hit_sequence.length())
  let mut index = 0
  while index < length {
    let model_marker = infernal_marker_at(model_sequence, index)
    let hit_marker = infernal_marker_at(hit_sequence, index)
    match (model_marker, hit_marker) {
      (None, None) => index = index + 1
      _ => {
        let model_end = match model_marker {
          Some(value) => value.0
          None => index
        }
        let hit_end = match hit_marker {
          Some(value) => value.0
          None => index
        }
        let end = model_end.max(hit_end)
        markers.push(InfernalMarker::{
          start: index,
          end,
          model_size: match model_marker {
            Some(value) => value.1
            None => 0
          },
          sequence_size: match hit_marker {
            Some(value) => value.1
            None => 0
          },
        })
        index = if end > index { end } else { index + 1 }
      }
    }
  }
  markers
}

///|
fn infernal_residue_count(sequence : String, start : Int, end : Int) -> Int {
  let bounded_end = end.min(sequence.length())
  if start >= bounded_end {
    return 0
  }
  let mut count = 0
  for index in start..<bounded_end {
    let value = sequence.unsafe_get(index).to_int()
    if (value >= 'A'.to_int() && value <= 'Z'.to_int()) ||
      (value >= 'a'.to_int() && value <= 'z'.to_int()) {
      count = count + 1
    }
  }
  count
}

///|
fn infernal_slice(sequence : String, start : Int, end : Int) -> String {
  let bounded_start = start.min(sequence.length())
  let bounded_end = end.min(sequence.length())
  if bounded_start >= bounded_end {
    ""
  } else {
    sequence[bounded_start:bounded_end].to_owned()
  }
}

///|
fn infernal_add_fragment(
  fragments : Array[InfernalFragment],
  model_sequence : String,
  hit_sequence : String,
  consensus_structure : String,
  noncanonical_pairs : String,
  similarity : String,
  posterior_probability : String,
  start : Int,
  end : Int,
  model_position : Int,
  hit_position : Int,
  strand : Strand,
  model_omission : Int,
  sequence_omission : Int,
) -> (Int, Int) {
  let model_count = infernal_residue_count(model_sequence, start, end)
  let hit_count = infernal_residue_count(hit_sequence, start, end)
  if model_count == 0 && hit_count == 0 {
    return (model_position, hit_position)
  }
  let (hit_start, hit_end, next_hit) = if strand == Strand::Minus {
    (hit_position - hit_count, hit_position, hit_position - hit_count)
  } else {
    (hit_position, hit_position + hit_count, hit_position + hit_count)
  }
  fragments.push(InfernalFragment::{
    query_start: model_position,
    query_end: model_position + model_count,
    hit_start,
    hit_end,
    hit_strand: strand,
    query_sequence: infernal_slice(model_sequence, start, end),
    hit_sequence: infernal_slice(hit_sequence, start, end),
    consensus_structure: infernal_slice(consensus_structure, start, end),
    noncanonical_pairs: infernal_slice(noncanonical_pairs, start, end),
    similarity: infernal_slice(similarity, start, end),
    posterior_probability: infernal_slice(posterior_probability, start, end),
    model_omission_before: model_omission,
    sequence_omission_before: sequence_omission,
  })
  (model_position + model_count, next_hit)
}

///|
fn infernal_alignment_fragments(
  base : InfernalFragment,
  model_sequence : String,
  hit_sequence : String,
  consensus_structure : String,
  noncanonical_pairs : String,
  similarity : String,
  posterior_probability : String,
) -> Array[InfernalFragment] {
  let markers = infernal_find_markers(model_sequence, hit_sequence)
  let fragments : Array[InfernalFragment] = []
  let mut model_position = base.query_start
  let mut hit_position = if base.hit_strand == Strand::Minus {
    base.hit_end
  } else {
    base.hit_start
  }
  let mut display_position = 0
  let mut model_omission = 0
  let mut sequence_omission = 0
  for marker in markers {
    let next = infernal_add_fragment(
      fragments,
      model_sequence,
      hit_sequence,
      consensus_structure,
      noncanonical_pairs,
      similarity,
      posterior_probability,
      display_position,
      marker.start,
      model_position,
      hit_position,
      base.hit_strand,
      model_omission,
      sequence_omission,
    )
    model_position = next.0
    hit_position = next.1
    let skipped_model = if marker.model_size > 0 {
      marker.model_size
    } else {
      infernal_residue_count(model_sequence, marker.start, marker.end)
    }
    let skipped_hit = if marker.sequence_size > 0 {
      marker.sequence_size
    } else {
      infernal_residue_count(hit_sequence, marker.start, marker.end)
    }
    model_position = model_position + skipped_model
    hit_position = if base.hit_strand == Strand::Minus {
      hit_position - skipped_hit
    } else {
      hit_position + skipped_hit
    }
    model_omission = skipped_model
    sequence_omission = skipped_hit
    display_position = marker.end
  }
  let display_end = model_sequence.length().max(hit_sequence.length())
  ignore(
    infernal_add_fragment(
      fragments,
      model_sequence,
      hit_sequence,
      consensus_structure,
      noncanonical_pairs,
      similarity,
      posterior_probability,
      display_position,
      display_end,
      model_position,
      hit_position,
      base.hit_strand,
      model_omission,
      sequence_omission,
    ),
  )
  if fragments.length() == 0 {
    [base]
  } else {
    fragments
  }
}

///|
fn infernal_text_metadata(
  lines : Array[String],
) -> (String, String, String, Bool) {
  let mut program = ""
  let mut version = ""
  let mut target = ""
  let mut show_alignments = true
  for line in lines {
    let trimmed = line.trim().to_owned()
    if trimmed.has_prefix("# cmsearch") || trimmed.has_prefix("# cmscan") {
      let values = infernal_split_whitespace(trimmed)
      if values.length() > 1 {
        program = values[1]
      }
    } else if trimmed.has_prefix("# INFERNAL") {
      let values = infernal_split_whitespace(trimmed)
      if values.length() > 2 {
        version = values[2]
      }
    } else if trimmed.has_prefix("#") && trimmed.contains(":") {
      if trimmed.contains("target") {
        let parts = split_by_char(trimmed, ':'.to_int())
        if parts.length() > 1 {
          target = parts[1].trim().to_owned()
        }
      }
      if trimmed.contains("show alignments in output") &&
        trimmed.has_suffix("no") {
        show_alignments = false
      }
    }
  }
  (program, version, target, show_alignments)
}

///|
fn infernal_text_hit_header(line : String) -> (String, String) {
  let values = infernal_split_whitespace(line.trim().to_owned())
  if values.length() < 2 {
    return ("", "")
  }
  (values[1], infernal_join(values, 2))
}

///|
fn infernal_parse_text_hsp(
  row : Array[String],
  model_sequence : String,
  hit_sequence : String,
  consensus_structure : String,
  noncanonical_pairs : String,
  similarity : String,
  posterior_probability : String,
) -> InfernalHSP raise InfernalError {
  if row.length() != 16 {
    raise InfernalError(
      "Infernal text alignment score row must contain 16 fields",
    )
  }
  let base = infernal_coordinate_fragment(
    infernal_parse_int(row[6], "mdl from"),
    infernal_parse_int(row[7], "mdl to"),
    infernal_parse_int(row[9], "seq from"),
    infernal_parse_int(row[10], "seq to"),
    row[11],
  )
  let fragments = if model_sequence.length() == 0 && hit_sequence.length() == 0 {
    [base]
  } else {
    infernal_alignment_fragments(
      base, model_sequence, hit_sequence, consensus_structure, noncanonical_pairs,
      similarity, posterior_probability,
    )
  }
  InfernalHSP::{
    fragments,
    model: row[5],
    truncated: row[14],
    pipeline_pass: 0,
    gc: infernal_parse_double(row[15], "gc"),
    bias: infernal_parse_double(row[4], "bias"),
    bitscore: infernal_parse_double(row[3], "score"),
    evalue: infernal_parse_double(row[2], "E-value"),
    is_included: row[1] == "!",
    query_end_type: row[8],
    hit_end_type: row[12],
    average_accuracy: infernal_parse_double(row[13], "accuracy"),
    overlap: "",
    any_index: None,
    any_overlap_fraction: None,
    reciprocal_overlap_fraction: None,
    winner_index: None,
    winner_overlap_fraction: None,
    winner_reciprocal_fraction: None,
  }
}

///|
fn infernal_alignment_strings(
  lines : Array[String],
  start : Int,
  end : Int,
  query_id : String,
  hit_id : String,
) -> (String, String, String, String, String, String) {
  let mut model_sequence = ""
  let mut hit_sequence = ""
  let mut consensus = ""
  let mut noncanonical = ""
  let mut similarity = ""
  let mut posterior = ""
  let mut pending_similarity = ""
  for index in start..<end {
    let trimmed = lines[index].trim().to_owned()
    let values = infernal_split_whitespace(trimmed)
    if values.length() == 0 {
      continue
    }
    if values.length() >= 4 &&
      infernal_is_integer(values[1]) &&
      infernal_is_integer(values[values.length() - 1]) {
      let sequence = values[2]
      if values[0] == hit_id {
        hit_sequence = hit_sequence + sequence
        if pending_similarity.length() > 0 {
          similarity = similarity + pending_similarity
          pending_similarity = ""
        }
      } else if values[0] == query_id || model_sequence.length() == 0 {
        model_sequence = model_sequence + sequence
      }
    } else {
      match values[0] {
        "CS" | "SS_cons" =>
          if values.length() > 1 {
            consensus = consensus + values[1]
          }
        "NC" =>
          if values.length() > 1 {
            noncanonical = noncanonical + values[1]
          }
        "PP" => if values.length() > 1 { posterior = posterior + values[1] }
        _ =>
          if model_sequence.length() > hit_sequence.length() {
            pending_similarity = values[0]
          }
      }
    }
  }
  (model_sequence, hit_sequence, consensus, noncanonical, similarity, posterior)
}

///|
fn infernal_parse_text_hit(
  lines : Array[String],
  start : Int,
  end : Int,
  query_id : String,
) -> InfernalHit raise InfernalError {
  let (hit_id, description) = infernal_text_hit_header(lines[start])
  if hit_id.length() == 0 {
    raise InfernalError("invalid Infernal text hit header")
  }
  let mut row_index = -1
  for index in (start + 1)..<end {
    let trimmed = lines[index].trim().to_owned()
    if trimmed.has_prefix("----") {
      let mut candidate = index + 1
      while candidate < end && lines[candidate].trim().length() == 0 {
        candidate = candidate + 1
      }
      if candidate < end {
        row_index = candidate
      }
      break
    }
  }
  if row_index < 0 {
    raise InfernalError("missing Infernal text alignment score row")
  }
  let row = infernal_split_whitespace(lines[row_index].trim().to_owned())
  let strings = infernal_alignment_strings(
    lines,
    row_index + 1,
    end,
    query_id,
    hit_id,
  )
  let hsp = infernal_parse_text_hsp(
    row,
    strings.0,
    strings.1,
    strings.2,
    strings.3,
    strings.4,
    strings.5,
  )
  InfernalHit::{
    id: hit_id,
    accession: "",
    description,
    model_length: 0,
    hsps: [hsp],
  }
}

///|
fn infernal_parse_noali_rows(
  lines : Array[String],
  start : Int,
  end : Int,
) -> Array[InfernalHit] raise InfernalError {
  let hits : Array[InfernalHit] = []
  let mut table_started = false
  for index in start..<end {
    let trimmed = lines[index].trim().to_owned()
    if trimmed.has_prefix("----") {
      table_started = true
      continue
    }
    if !table_started ||
      trimmed.length() == 0 ||
      trimmed.has_prefix("------") ||
      trimmed.has_prefix("Hit alignments:") {
      continue
    }
    let values = infernal_split_whitespace(trimmed)
    if values.length() < 12 || !infernal_is_integer(values[0]) {
      continue
    }
    let hit_id = values[5]
    let sequence_from = infernal_parse_int(values[6], "seq from")
    let sequence_to = infernal_parse_int(values[7], "seq to")
    let strand = values[8]
    let (hit_start, hit_end, hit_strand) = match strand {
      "+" => (sequence_from - 1, sequence_to, Strand::Plus)
      "-" => (sequence_to - 1, sequence_from, Strand::Minus)
      _ => raise InfernalError("invalid Infernal strand: " + strand)
    }
    let fragment = InfernalFragment::{
      query_start: -1,
      query_end: -1,
      hit_start,
      hit_end,
      hit_strand,
      query_sequence: "",
      hit_sequence: "",
      consensus_structure: "",
      noncanonical_pairs: "",
      similarity: "",
      posterior_probability: "",
      model_omission_before: 0,
      sequence_omission_before: 0,
    }
    let hsp = InfernalHSP::{
      fragments: [fragment],
      model: values[9],
      truncated: values[10],
      pipeline_pass: 0,
      gc: infernal_parse_double(values[11], "gc"),
      bias: infernal_parse_double(values[4], "bias"),
      bitscore: infernal_parse_double(values[3], "score"),
      evalue: infernal_parse_double(values[2], "E-value"),
      is_included: values[1] == "!",
      query_end_type: "",
      hit_end_type: "",
      average_accuracy: 0.0,
      overlap: "",
      any_index: None,
      any_overlap_fraction: None,
      reciprocal_overlap_fraction: None,
      winner_index: None,
      winner_overlap_fraction: None,
      winner_reciprocal_fraction: None,
    }
    hits.push(InfernalHit::{
      id: hit_id,
      accession: "",
      description: infernal_join(values, 12),
      model_length: 0,
      hsps: [hsp],
    })
  }
  hits
}

///|
fn infernal_merge_hits(source : Array[InfernalHit]) -> Array[InfernalHit] {
  let order : Array[String] = []
  let grouped : Map[String, Array[InfernalHit]] = Map([], capacity=0)
  for hit in source {
    match grouped.get(hit.id) {
      Some(values) => values.push(hit)
      None => {
        order.push(hit.id)
        grouped.set(hit.id, [hit])
      }
    }
  }
  let results : Array[InfernalHit] = []
  for hit_id in order {
    let matches = grouped.get(hit_id).unwrap()
    let hsps : Array[InfernalHSP] = []
    for hit in matches {
      for hsp in hit.hsps {
        hsps.push(hsp)
      }
    }
    results.push(InfernalHit::{ ..matches[0], hsps, })
  }
  results
}

///|
/// Parse non-verbose Infernal plain-text output.
///
/// Both alignment output and `--noali` score tables are supported. Alignment
/// local ends (`*[NN]*`) split one HSP into coordinate-aware fragments.
pub fn infernal_parse_text(
  content : String,
) -> Array[InfernalQueryResult] raise InfernalError {
  let lines = split_by_char(content, '\n'.to_int())
  let metadata = infernal_text_metadata(lines)
  let results : Array[InfernalQueryResult] = []
  let mut index = 0
  while index < lines.length() {
    let trimmed = lines[index].trim().to_owned()
    if !trimmed.has_prefix("Query:") {
      index = index + 1
      continue
    }
    let query_values = infernal_split_whitespace(trimmed)
    if query_values.length() < 2 {
      raise InfernalError("invalid Infernal Query line")
    }
    let query_id = query_values[1]
    let mut query_length = 0
    for value in query_values {
      if value.contains("M=") || value.contains("CLEN=") || value.contains("L=") {
        query_length = infernal_digit_count(value)
      }
    }
    let query_start = index
    let mut query_end = index + 1
    while query_end < lines.length() &&
          !lines[query_end].trim().has_prefix("//") {
      query_end = query_end + 1
    }
    let mut accession = ""
    let mut description = ""
    for cursor in (query_start + 1)..<query_end {
      let current = lines[cursor].trim().to_owned()
      if current.has_prefix("Accession:") {
        let values = infernal_split_whitespace(current)
        if values.length() > 1 {
          accession = values[1]
        }
      } else if current.has_prefix("Description:") {
        let values = infernal_split_whitespace(current)
        description = infernal_join(values, 1)
      }
    }
    let hits : Array[InfernalHit] = []
    if metadata.3 {
      let hit_starts : Array[Int] = []
      for cursor in (query_start + 1)..<query_end {
        if lines[cursor].trim().has_prefix(">> ") {
          hit_starts.push(cursor)
        }
      }
      for hit_index in 0..<hit_starts.length() {
        let hit_end = if hit_index + 1 < hit_starts.length() {
          hit_starts[hit_index + 1]
        } else {
          query_end
        }
        hits.push(
          infernal_parse_text_hit(
            lines,
            hit_starts[hit_index],
            hit_end,
            query_id,
          ),
        )
      }
    } else {
      for hit in infernal_parse_noali_rows(lines, query_start + 1, query_end) {
        hits.push(hit)
      }
    }
    results.push(InfernalQueryResult::{
      id: query_id,
      accession,
      description,
      clan: "",
      sequence_length: query_length,
      program: metadata.0,
      version: metadata.1,
      target: metadata.2,
      hits: infernal_merge_hits(hits),
    })
    index = if query_end < lines.length() { query_end + 1 } else { query_end }
  }
  results
}

///|
pub fn InfernalHSP::query_start(self : InfernalHSP) -> Int {
  if self.fragments.length() == 0 {
    -1
  } else {
    let mut result = self.fragments[0].query_start
    for fragment in self.fragments {
      if fragment.query_start >= 0 && fragment.query_start < result {
        result = fragment.query_start
      }
    }
    result
  }
}

///|
pub fn InfernalHSP::query_end(self : InfernalHSP) -> Int {
  let mut result = -1
  for fragment in self.fragments {
    result = result.max(fragment.query_end)
  }
  result
}

///|
pub fn InfernalHSP::hit_start(self : InfernalHSP) -> Int {
  if self.fragments.length() == 0 {
    -1
  } else {
    let mut result = self.fragments[0].hit_start
    for fragment in self.fragments {
      result = result.min(fragment.hit_start)
    }
    result
  }
}

///|
pub fn InfernalHSP::hit_end(self : InfernalHSP) -> Int {
  let mut result = -1
  for fragment in self.fragments {
    result = result.max(fragment.hit_end)
  }
  result
}

///|
pub fn InfernalHSP::alignment_span(self : InfernalHSP) -> Int {
  let mut result = 0
  for fragment in self.fragments {
    result = result + fragment.hit_end - fragment.hit_start
  }
  result
}

///|
pub fn InfernalHit::best_hsp(self : InfernalHit) -> InfernalHSP? {
  if self.hsps.length() == 0 {
    return None
  }
  let mut best = self.hsps[0]
  for hsp in self.hsps {
    if hsp.evalue < best.evalue ||
      (hsp.evalue == best.evalue && hsp.bitscore > best.bitscore) {
      best = hsp
    }
  }
  Some(best)
}

///|
pub fn InfernalQueryResult::count_hsps(self : InfernalQueryResult) -> Int {
  let mut count = 0
  for hit in self.hits {
    count = count + hit.hsps.length()
  }
  count
}

///|
pub fn InfernalQueryResult::filter(
  self : InfernalQueryResult,
  maximum_evalue? : Double = 1.0,
  included_only? : Bool = false,
) -> InfernalQueryResult raise InfernalError {
  if maximum_evalue < 0.0 || maximum_evalue.is_nan() {
    raise InfernalError("maximum_evalue must be non-negative")
  }
  let hits : Array[InfernalHit] = []
  for hit in self.hits {
    let hsps : Array[InfernalHSP] = []
    for hsp in hit.hsps {
      if hsp.evalue <= maximum_evalue && (!included_only || hsp.is_included) {
        hsps.push(hsp)
      }
    }
    if hsps.length() > 0 {
      hits.push(InfernalHit::{ ..hit, hsps, })
    }
  }
  InfernalQueryResult::{ ..self, hits, }
}

///|
/// Convert the Infernal-specific model to the repository SearchIO hierarchy.
pub fn InfernalQueryResult::to_searchio(
  self : InfernalQueryResult,
) -> QueryResult {
  let hits : Array[Hit] = []
  for source_hit in self.hits {
    let hsps : Array[HSP] = []
    let mut hit_score = 0.0
    let mut hit_evalue = 1.0
    for source_hsp in source_hit.hsps {
      let fragments : Array[HSPFragment] = []
      for source_fragment in source_hsp.fragments {
        fragments.push(
          HSPFragment::new(
            source_fragment.query_start,
            source_fragment.query_end,
            source_fragment.hit_start,
            source_fragment.hit_end,
            Strand::Star,
            source_fragment.hit_strand,
            (source_fragment.hit_end - source_fragment.hit_start).max(0),
            0,
            source_hsp.bitscore,
            source_hsp.evalue,
            source_fragment.query_sequence,
            source_fragment.hit_sequence,
            source_fragment.similarity,
          ),
        )
      }
      hsps.push(
        HSP::new(
          source_hit.id,
          source_hit.description,
          fragments,
          source_hsp.bitscore,
          source_hsp.evalue,
          source_hsp.bitscore,
        ),
      )
      hit_score = hit_score + source_hsp.bitscore
      if source_hsp.evalue < hit_evalue {
        hit_evalue = source_hsp.evalue
      }
    }
    hits.push(
      Hit::new(
        source_hit.id,
        source_hit.description,
        source_hit.accession,
        hsps,
        hit_score,
        hit_evalue,
      ),
    )
  }
  QueryResult::new(
    self.id,
    self.description,
    self.sequence_length,
    hits,
    self.program,
    self.version,
    self.target,
  )
}

///|
pub fn InfernalQueryResult::summary(self : InfernalQueryResult) -> String {
  "InfernalQueryResult(query=" +
  self.id +
  ", hits=" +
  self.hits.length().to_string() +
  ", hsps=" +
  self.count_hsps().to_string() +
  ", program=" +
  self.program +
  ")"
}

///|
pub fn infernal_tabular_sample() -> String {
  "#target name accession query name accession mdl mdl from mdl to seq from seq to strand trunc pass gc bias score E-value inc mdl len seq len description of target\n" +
  "#----------- --------- ---------- --------- --- -------- -------- -------- -------- ------ ----- ---- -- ---- ----- ------- --- ------- ------- ---------------------\n" +
  "chrA RFSEQ1 RF00001 RF00001 cm 1 71 101 171 + no 1 0.52 0.0 48.5 2.0e-14 ! 71 1000 bacterial 5S RNA locus\n" +
  "chrA RFSEQ1 RF00001 RF00001 cm 2 70 480 412 - 3' 2 0.47 0.1 39.2 4.5e-09 ? 71 1000 reverse strand copy\n" +
  "chrB RFSEQ2 RF00002 RF00002 hmm 3 66 20 83 + no 1 0.61 0.0 31.0 7.0e-06 ! 68 800 riboswitch candidate\n"
}

///|
pub fn infernal_text_sample() -> String {
  "# cmsearch :: search a covariance model against a sequence database\n" +
  "# INFERNAL 1.1.5 (Sep 2023)\n" +
  "# - - -\n" +
  "# target sequence database: transcripts.fa\n" +
  "# show alignments in output: yes\n" +
  "# - - -\n" +
  "Query: RFTEST [M=19]\n" +
  "Accession: RFTEST.1\n" +
  "Description: synthetic structured RNA\n" +
  "Hit scores:\n" +
  "Hit alignments:\n" +
  ">> targetA  synthetic RNA target\n" +
  " ----   --------- ------ ----- ---- --- -------- -------- --- -------- -------- --- --- ---- ----- --\n" +
  " 1 ! 2.0e-10 45.0 0.1 cm 1 19 [] 100 118 + .. 0.96 no 0.50\n" +
  "CS (((....))).........\n" +
  "RFTEST 1 ACGUACGU*[5]*AACCGG 19\n" +
  "targetA 100 ACGUACGUGGGGGAACCGG 118\n" +
  "PP 9999999999999999999\n" +
  "Internal CM pipeline statistics summary:\n" +
  "//\n" +
  "[ok]\n"
}
