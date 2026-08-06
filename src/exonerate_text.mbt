// Biopython-compatible Bio.SearchIO.ExonerateIO.exonerate_text support.
//
// This module parses Exonerate's human-readable C4 alignments. It is separate
// from exonerate.mbt (vulgar/cigar summaries) and align_exonerate.mbt (modern
// operation paths) because SearchIO text reports preserve aligned sequences,
// annotations, spliced fragments, split codons, and frame shifts.

///|
/// Error raised for malformed or internally inconsistent Exonerate text.
pub suberror ExonerateTextError {
  ExonerateTextError(String)
}

///|
/// Strand semantics used by Exonerate SearchIO fragments.
pub(all) enum ExonerateTextStrand {
  ExonerateTextForward
  ExonerateTextReverse
  ExonerateTextProtein
} derive(Eq, Debug)

///|
/// A normalized zero-based half-open coordinate range.
pub struct ExonerateTextRange {
  start : Int
  end : Int
} derive(Eq, Debug)

///|
/// File-level metadata from an Exonerate report.
pub struct ExonerateTextMetadata {
  program : String
  command_line : String
  hostname : String
} derive(Eq, Debug)

///|
/// One uninterrupted aligned fragment.
///
/// Translated alignments expose one-letter amino-acid sequences on both rows.
/// `query_step` and `hit_step` retain the coordinate width represented by one
/// displayed residue (one for protein, usually three for translated DNA).
pub struct ExonerateTextFragment {
  query_sequence : String
  hit_sequence : String
  similarity : String
  query_annotation : String?
  hit_annotation : String?
  query_start : Int
  query_end : Int
  hit_start : Int
  hit_end : Int
  query_strand : ExonerateTextStrand
  hit_strand : ExonerateTextStrand
  query_frame : Int
  hit_frame : Int
  phase : Int
  query_step : Int
  hit_step : Int
} derive(Eq, Debug)

///|
/// One Exonerate high-scoring pair, potentially split into many fragments.
pub struct ExonerateTextHsp {
  score : Int
  model : String
  query_start : Int
  query_end : Int
  hit_start : Int
  hit_end : Int
  query_strand : ExonerateTextStrand
  hit_strand : ExonerateTextStrand
  fragments : Array[ExonerateTextFragment]
  query_split_codons : Array[ExonerateTextRange]
  hit_split_codons : Array[ExonerateTextRange]
} derive(Eq, Debug)

///|
/// One target and all HSPs associated with it.
pub struct ExonerateTextHit {
  id : String
  description : String
  hsps : Array[ExonerateTextHsp]
} derive(Eq, Debug)

///|
/// One query result and all target hits associated with it.
pub struct ExonerateTextQueryResult {
  id : String
  description : String
  model : String
  hits : Array[ExonerateTextHit]
} derive(Eq, Debug)

///|
/// A complete Exonerate human-readable report.
pub struct ExonerateTextDocument {
  metadata : ExonerateTextMetadata
  queries : Array[ExonerateTextQueryResult]
} derive(Eq, Debug)

///|
/// Aggregate statistics for one HSP.
pub struct ExonerateTextCounts {
  fragments : Int
  alignment_columns : Int
  identities : Int
  positives : Int
  mismatches : Int
  gap_columns : Int
  query_gap_opens : Int
  hit_gap_opens : Int
  introns_or_ner_gaps : Int
  frameshift_gaps : Int
} derive(Eq, Debug)

///|
/// Construct a normalized coordinate range.
pub fn ExonerateTextRange::create(
  start : Int,
  end : Int,
) -> ExonerateTextRange raise ExonerateTextError {
  if start < 0 || end < start {
    raise ExonerateTextError(
      "Exonerate coordinate ranges must satisfy 0 <= start <= end",
    )
  }
  ExonerateTextRange::{ start, end }
}

///|
/// Parse a complete Exonerate C4 text report.
///
/// Coordinates are returned as zero-based half-open ranges. Empty reports are
/// valid, but metadata and `-- completed exonerate analysis` are required.
pub fn exonerate_text_parse(
  text : String,
) -> ExonerateTextDocument raise ExonerateTextError {
  let lines = exonerate_text_lines(text)
  if lines.length() == 0 {
    raise ExonerateTextError("Empty Exonerate text input")
  }
  let command_line = exonerate_text_metadata_value(lines, "Command line:")
  let hostname = exonerate_text_metadata_value(lines, "Hostname:")
  let mut completed = false
  for line in lines {
    if line.trim().has_prefix("-- completed exonerate analysis") {
      completed = true
    }
  }
  if !completed {
    raise ExonerateTextError(
      "Exonerate text report is missing its completion marker",
    )
  }
  let parsed : Array[ExonerateTextParsedAlignment] = []
  let mut index = 0
  while index < lines.length() {
    if lines[index].trim() == "C4 Alignment:" {
      let mut stop = index + 1
      while stop < lines.length() &&
            lines[stop].trim() != "C4 Alignment:" &&
            !lines[stop].trim().has_prefix("-- completed exonerate analysis") {
        stop = stop + 1
      }
      parsed.push(exonerate_text_parse_alignment(lines, index, stop))
      index = stop
    } else {
      index = index + 1
    }
  }
  let queries = exonerate_text_aggregate(parsed)
  ExonerateTextDocument::{
    metadata: ExonerateTextMetadata::{
      program: "exonerate",
      command_line,
      hostname,
    },
    queries,
  }
}

///|
pub fn ExonerateTextDocument::num_queries(self : ExonerateTextDocument) -> Int {
  self.queries.length()
}

///|
pub fn ExonerateTextDocument::num_hits(self : ExonerateTextDocument) -> Int {
  let mut count = 0
  for query in self.queries {
    count = count + query.hits.length()
  }
  count
}

///|
pub fn ExonerateTextDocument::num_hsps(self : ExonerateTextDocument) -> Int {
  let mut count = 0
  for query in self.queries {
    for hit in query.hits {
      count = count + hit.hsps.length()
    }
  }
  count
}

///|
pub fn ExonerateTextDocument::find_query(
  self : ExonerateTextDocument,
  id : String,
) -> ExonerateTextQueryResult? {
  for query in self.queries {
    if query.id == id {
      return Some(exonerate_text_copy_query(query))
    }
  }
  None
}

///|
pub fn ExonerateTextDocument::best_hsp(
  self : ExonerateTextDocument,
) -> ExonerateTextHsp? {
  let mut found = false
  let mut best = ExonerateTextHsp::{
    score: 0,
    model: "",
    query_start: 0,
    query_end: 0,
    hit_start: 0,
    hit_end: 0,
    query_strand: ExonerateTextForward,
    hit_strand: ExonerateTextForward,
    fragments: [],
    query_split_codons: [],
    hit_split_codons: [],
  }
  for query in self.queries {
    for hit in query.hits {
      for hsp in hit.hsps {
        if !found || hsp.score > best.score {
          best = hsp
          found = true
        }
      }
    }
  }
  if found {
    Some(exonerate_text_copy_hsp(best))
  } else {
    None
  }
}

///|
pub fn ExonerateTextDocument::summary(self : ExonerateTextDocument) -> String {
  "Exonerate text report(host=" +
  self.metadata.hostname +
  ", queries=" +
  self.num_queries().to_string() +
  ", hits=" +
  self.num_hits().to_string() +
  ", hsps=" +
  self.num_hsps().to_string() +
  ")"
}

///|
pub fn ExonerateTextQueryResult::find_hit(
  self : ExonerateTextQueryResult,
  id : String,
) -> ExonerateTextHit? {
  for hit in self.hits {
    if hit.id == id {
      return Some(exonerate_text_copy_hit(hit))
    }
  }
  None
}

///|
pub fn ExonerateTextHsp::num_fragments(self : ExonerateTextHsp) -> Int {
  self.fragments.length()
}

///|
pub fn ExonerateTextHsp::query_ranges(
  self : ExonerateTextHsp,
) -> Array[ExonerateTextRange] {
  let result : Array[ExonerateTextRange] = []
  for fragment in self.fragments {
    result.push(ExonerateTextRange::{
      start: fragment.query_start,
      end: fragment.query_end,
    })
  }
  result
}

///|
pub fn ExonerateTextHsp::hit_ranges(
  self : ExonerateTextHsp,
) -> Array[ExonerateTextRange] {
  let result : Array[ExonerateTextRange] = []
  for fragment in self.fragments {
    result.push(ExonerateTextRange::{
      start: fragment.hit_start,
      end: fragment.hit_end,
    })
  }
  result
}

///|
pub fn ExonerateTextHsp::query_inter_ranges(
  self : ExonerateTextHsp,
) -> Array[ExonerateTextRange] {
  exonerate_text_inter_ranges(self.fragments, true)
}

///|
pub fn ExonerateTextHsp::hit_inter_ranges(
  self : ExonerateTextHsp,
) -> Array[ExonerateTextRange] {
  exonerate_text_inter_ranges(self.fragments, false)
}

///|
pub fn ExonerateTextHsp::counts(self : ExonerateTextHsp) -> ExonerateTextCounts {
  let mut alignment_columns = 0
  let mut identities = 0
  let mut positives = 0
  let mut mismatches = 0
  let mut gap_columns = 0
  let mut query_gap_opens = 0
  let mut hit_gap_opens = 0
  let mut query_in_gap = false
  let mut hit_in_gap = false
  for fragment in self.fragments {
    let width = if fragment.query_sequence.length() <
      fragment.hit_sequence.length() {
      fragment.query_sequence.length()
    } else {
      fragment.hit_sequence.length()
    }
    alignment_columns = alignment_columns + width
    for index = 0; index < width; index = index + 1 {
      let qgap = fragment.query_sequence.unsafe_get(index).to_int() ==
        '-'.to_int()
      let hgap = fragment.hit_sequence.unsafe_get(index).to_int() ==
        '-'.to_int()
      if qgap || hgap {
        gap_columns = gap_columns + 1
      } else {
        let qcode = exonerate_text_upper(
          fragment.query_sequence.unsafe_get(index).to_int(),
        )
        let hcode = exonerate_text_upper(
          fragment.hit_sequence.unsafe_get(index).to_int(),
        )
        if qcode == hcode {
          identities = identities + 1
        } else {
          mismatches = mismatches + 1
        }
        if qcode == hcode ||
          exonerate_text_similarity_positive(fragment.similarity, index, width) {
          positives = positives + 1
        }
      }
      if qgap && !query_in_gap {
        query_gap_opens = query_gap_opens + 1
      }
      if hgap && !hit_in_gap {
        hit_gap_opens = hit_gap_opens + 1
      }
      query_in_gap = qgap
      hit_in_gap = hgap
    }
    query_in_gap = false
    hit_in_gap = false
  }
  let query_inter = self.query_inter_ranges()
  let hit_inter = self.hit_inter_ranges()
  let mut intervening = 0
  let mut frameshift_gaps = 0
  for index = 0; index < query_inter.length(); index = index + 1 {
    let qlen = query_inter[index].end - query_inter[index].start
    let hlen = hit_inter[index].end - hit_inter[index].start
    if qlen > 0 || hlen > 0 {
      intervening = intervening + 1
    }
    if (qlen > 0 && qlen <= 2) || (hlen > 0 && hlen <= 2) {
      frameshift_gaps = frameshift_gaps + 1
    }
  }
  ExonerateTextCounts::{
    fragments: self.fragments.length(),
    alignment_columns,
    identities,
    positives,
    mismatches,
    gap_columns,
    query_gap_opens,
    hit_gap_opens,
    introns_or_ner_gaps: intervening,
    frameshift_gaps,
  }
}

///|
/// Map an aligned residue index to its query coordinate.
pub fn ExonerateTextFragment::query_position(
  self : ExonerateTextFragment,
  aligned_index : Int,
) -> Int? {
  exonerate_text_fragment_position(
    self.query_sequence,
    aligned_index,
    self.query_start,
    self.query_end,
    self.query_strand,
    self.query_step,
  )
}

///|
/// Map an aligned residue index to its hit coordinate.
pub fn ExonerateTextFragment::hit_position(
  self : ExonerateTextFragment,
  aligned_index : Int,
) -> Int? {
  exonerate_text_fragment_position(
    self.hit_sequence,
    aligned_index,
    self.hit_start,
    self.hit_end,
    self.hit_strand,
    self.hit_step,
  )
}

///|
/// Small official-style report used by the offline example.
pub fn exonerate_text_example() -> String {
  let marker = "  >>>> Target Intron 1 >>>>  "
  let splice = "gt" + ".".repeat(marker.length() - 4) + "ag"
  let label = "7 bp"
  let left = (marker.length() - label.length()) / 2
  let intron_similarity = " ".repeat(left) +
    label +
    " ".repeat(marker.length() - label.length() - left)
  "Command line: [exonerate -m est2genome query.fa genome.fa]\n" +
  "Hostname: [offline]\n\n" +
  "C4 Alignment:\n" +
  "------------\n" +
  "         Query: transcript demo transcript\n" +
  "        Target: chromosome demo chromosome\n" +
  "         Model: est2genome\n" +
  "     Raw score: 42\n" +
  "   Query range: 0 -> 12\n" +
  "  Target range: 100 -> 119\n\n" +
  "   1 : ACGTAC" +
  marker +
  "GTACGT :  12\n" +
  "       ||||||" +
  intron_similarity +
  "||||||\n" +
  " 101 : ACGTAC" +
  splice +
  "GTACGT : 119\n\n" +
  "-- completed exonerate analysis\n"
}

///|
priv struct ExonerateTextParsedAlignment {
  query_id : String
  query_description : String
  hit_id : String
  hit_description : String
  model : String
  hsp : ExonerateTextHsp
}

///|
priv struct ExonerateTextHeader {
  query_id : String
  query_description : String
  hit_id : String
  hit_description : String
  model : String
  score : Int
  query_raw_start : Int
  query_raw_end : Int
  hit_raw_start : Int
  hit_raw_end : Int
  query_strand : ExonerateTextStrand
  hit_strand : ExonerateTextStrand
  body_start : Int
}

///|
priv struct ExonerateTextRowMap {
  query : Int
  similarity : Int
  hit : Int
  query_annotation : Int
  hit_annotation : Int
}

///|
priv struct ExonerateTextSpan {
  start : Int
  end : Int
}

///|
priv struct ExonerateTextRawFragment {
  query : String
  hit : String
  similarity : String
  query_annotation : String?
  hit_annotation : String?
  query_start : Int
  query_end : Int
  hit_start : Int
  hit_end : Int
  phase : Int
}

///|
fn exonerate_text_parse_alignment(
  lines : Array[String],
  start : Int,
  stop : Int,
) -> ExonerateTextParsedAlignment raise ExonerateTextError {
  let header = exonerate_text_parse_header(lines, start, stop)
  let raw_blocks = exonerate_text_read_physical_blocks(
    lines,
    header.body_start,
    stop,
  )
  let rows = exonerate_text_stitch_rows(raw_blocks)
  let row_map = exonerate_text_row_map(rows.length(), header.model)
  let has_ner = header.model.to_upper().contains("NER")
  let separators = if has_ner {
    exonerate_text_find_ner_spans(rows[row_map.query])
  } else {
    exonerate_text_find_intron_spans(rows[row_map.query])
  }
  let blocks = exonerate_text_complement_spans(
    rows[row_map.query].length(),
    separators,
  )
  if blocks.length() == 0 {
    raise ExonerateTextError(
      "Exonerate alignment contains no sequence fragment",
    )
  }
  let query_gaps : Array[Int] = []
  let hit_gaps : Array[Int] = []
  for separator in separators {
    if has_ner {
      query_gaps.push(
        exonerate_text_parse_ner_gap(
          rows[row_map.query][separator.start:separator.end].to_owned(),
          "query",
        ),
      )
      hit_gaps.push(
        exonerate_text_parse_ner_gap(
          rows[row_map.hit][separator.start:separator.end].to_owned(),
          "hit",
        ),
      )
    } else {
      let query_inter = rows[row_map.query][separator.start:separator.end].to_owned()
      let hit_inter = rows[row_map.hit][separator.start:separator.end].to_owned()
      let similarity = rows[row_map.similarity][separator.start:separator.end].to_owned()
      let lengths = exonerate_text_bp_lengths(similarity)
      let query_has_intron = hit_inter.contains("Intron")
      let hit_has_intron = query_inter.contains("Intron")
      if !query_has_intron && !hit_has_intron {
        raise ExonerateTextError(
          "Exonerate intron separator has no Intron annotation",
        )
      }
      if lengths.length() == 2 {
        query_gaps.push(if query_has_intron { lengths[0] } else { 0 })
        hit_gaps.push(if hit_has_intron { lengths[1] } else { 0 })
      } else if lengths.length() == 1 {
        query_gaps.push(if query_has_intron { lengths[0] } else { 0 })
        hit_gaps.push(if hit_has_intron { lengths[0] } else { 0 })
      } else {
        raise ExonerateTextError(
          "Exonerate intron separator must contain one or two bp lengths",
        )
      }
    }
  }
  if query_gaps.length() + 1 != blocks.length() ||
    hit_gaps.length() + 1 != blocks.length() {
    raise ExonerateTextError(
      "Exonerate intervening-gap count does not match fragment count",
    )
  }
  let built = exonerate_text_build_fragments(
    rows, row_map, blocks, query_gaps, hit_gaps, header,
  )
  ExonerateTextParsedAlignment::{
    query_id: header.query_id,
    query_description: header.query_description,
    hit_id: header.hit_id,
    hit_description: header.hit_description,
    model: header.model,
    hsp: built,
  }
}

///|
fn exonerate_text_parse_header(
  lines : Array[String],
  start : Int,
  stop : Int,
) -> ExonerateTextHeader raise ExonerateTextError {
  let mut query_value = ""
  let mut hit_value = ""
  let mut model = ""
  let mut score_value = ""
  let mut query_range = ""
  let mut hit_range = ""
  let mut body_start = -1
  let mut index = start + 1
  while index < stop {
    let trimmed = lines[index].trim().to_owned()
    if trimmed.has_prefix("Query:") {
      query_value = exonerate_text_after_prefix(trimmed, "Query:")
    } else if trimmed.has_prefix("Target:") {
      hit_value = exonerate_text_after_prefix(trimmed, "Target:")
    } else if trimmed.has_prefix("Model:") {
      model = exonerate_text_after_prefix(trimmed, "Model:")
    } else if trimmed.has_prefix("Raw score:") {
      score_value = exonerate_text_after_prefix(trimmed, "Raw score:")
    } else if trimmed.has_prefix("Query range:") {
      query_range = exonerate_text_after_prefix(trimmed, "Query range:")
    } else if trimmed.has_prefix("Target range:") {
      hit_range = exonerate_text_after_prefix(trimmed, "Target range:")
    } else if query_value.length() > 0 &&
      hit_value.length() > 0 &&
      model.length() > 0 &&
      score_value.length() > 0 &&
      query_range.length() > 0 &&
      hit_range.length() > 0 &&
      trimmed.length() == 0 {
      body_start = index + 1
      break
    }
    index = index + 1
  }
  if query_value.length() == 0 ||
    hit_value.length() == 0 ||
    model.length() == 0 ||
    score_value.length() == 0 ||
    query_range.length() == 0 ||
    hit_range.length() == 0 ||
    body_start < 0 {
    raise ExonerateTextError(
      "Exonerate C4 header is incomplete or lacks a separating blank line",
    )
  }
  let query_pair = exonerate_text_id_description(query_value, "query")
  let hit_pair = exonerate_text_id_description(hit_value, "target")
  let query_coords = exonerate_text_parse_range(query_range, "query")
  let hit_coords = exonerate_text_parse_range(hit_range, "target")
  let query_is_protein = model.to_lower().contains("protein2")
  let hit_is_protein = model.to_lower().contains("2protein")
  let query_strand_pair = exonerate_text_strand(query_pair.1, query_is_protein)
  let hit_strand_pair = exonerate_text_strand(hit_pair.1, hit_is_protein)
  exonerate_text_validate_orientation(
    query_coords.0,
    query_coords.1,
    query_strand_pair.0,
    "query",
  )
  exonerate_text_validate_orientation(
    hit_coords.0,
    hit_coords.1,
    hit_strand_pair.0,
    "target",
  )
  ExonerateTextHeader::{
    query_id: query_pair.0,
    query_description: query_strand_pair.1,
    hit_id: hit_pair.0,
    hit_description: hit_strand_pair.1,
    model,
    score: exonerate_text_parse_nonnegative_int(score_value, "Raw score"),
    query_raw_start: query_coords.0,
    query_raw_end: query_coords.1,
    hit_raw_start: hit_coords.0,
    hit_raw_end: hit_coords.1,
    query_strand: query_strand_pair.0,
    hit_strand: hit_strand_pair.0,
    body_start,
  }
}

///|
fn exonerate_text_read_physical_blocks(
  lines : Array[String],
  start : Int,
  stop : Int,
) -> Array[Array[String]] raise ExonerateTextError {
  let result : Array[Array[String]] = []
  let mut index = start
  while index < stop {
    while index < stop && exonerate_text_coordinate_layout(lines[index]) is None {
      let trimmed = lines[index].trim().to_owned()
      if trimmed.has_prefix("vulgar:") || trimmed.has_prefix("cigar:") {
        index = stop
        break
      }
      if trimmed.length() > 0 {
        raise ExonerateTextError(
          "Unexpected line before an Exonerate alignment row: " + trimmed,
        )
      }
      index = index + 1
    }
    if index >= stop {
      break
    }
    let layout = match exonerate_text_coordinate_layout(lines[index]) {
      Some(value) => value
      None => raise ExonerateTextError("Internal Exonerate row-layout failure")
    }
    let block : Array[String] = []
    block.push(exonerate_text_fixed_slice(lines[index], layout.0, layout.1))
    index = index + 1
    let mut closed = false
    while index < stop {
      let row = exonerate_text_fixed_slice(lines[index], layout.0, layout.1)
      block.push(row)
      if exonerate_text_coordinate_layout(lines[index]) is Some(_) {
        closed = true
        index = index + 1
        break
      }
      index = index + 1
    }
    if !closed {
      raise ExonerateTextError(
        "Exonerate physical alignment block lacks its closing coordinate row",
      )
    }
    if block.length() < 3 || block.length() > 5 {
      raise ExonerateTextError(
        "Exonerate physical alignment blocks must contain 3, 4, or 5 rows",
      )
    }
    result.push(block)
  }
  if result.length() == 0 {
    raise ExonerateTextError("Exonerate C4 block has no alignment rows")
  }
  result
}

///|
fn exonerate_text_stitch_rows(
  blocks : Array[Array[String]],
) -> Array[String] raise ExonerateTextError {
  let mut row_count = 0
  for block in blocks {
    if block.length() > row_count {
      row_count = block.length()
    }
  }
  if row_count < 3 || row_count > 5 {
    raise ExonerateTextError("Unexpected Exonerate stitched row count")
  }
  let builders : Array[StringBuilder] = []
  for index = 0; index < row_count; index = index + 1 {
    builders.push(StringBuilder::new())
  }
  for block in blocks {
    if block.length() == row_count {
      for row = 0; row < row_count; row = row + 1 {
        builders[row].write_string(block[row])
      }
    } else if block.length() + 2 == row_count {
      let blank = " ".repeat(block[0].length())
      builders[0].write_string(blank)
      for row = 0; row < block.length(); row = row + 1 {
        builders[row + 1].write_string(block[row])
      }
      builders[row_count - 1].write_string(blank)
    } else {
      raise ExonerateTextError(
        "Inconsistent Exonerate physical row counts while stitching",
      )
    }
  }
  let rows : Array[String] = []
  for builder in builders {
    rows.push(builder.to_string())
  }
  if row_count == 5 {
    let query_flip = exonerate_text_flip_codons(rows[0], rows[1])
    rows[0] = query_flip.0
    rows[1] = query_flip.1
    let hit_flip = exonerate_text_flip_codons(rows[4], rows[3])
    rows[4] = hit_flip.0
    rows[3] = hit_flip.1
  }
  rows
}

///|
fn exonerate_text_row_map(
  count : Int,
  model : String,
) -> ExonerateTextRowMap raise ExonerateTextError {
  if count == 3 {
    ExonerateTextRowMap::{
      query: 0,
      similarity: 1,
      hit: 2,
      query_annotation: -1,
      hit_annotation: -1,
    }
  } else if count == 4 {
    let lower = model.to_lower()
    if lower.contains("protein2") {
      ExonerateTextRowMap::{
        query: 0,
        similarity: 1,
        hit: 2,
        query_annotation: -1,
        hit_annotation: 3,
      }
    } else if lower.contains("2protein") {
      ExonerateTextRowMap::{
        query: 1,
        similarity: 2,
        hit: 3,
        query_annotation: 0,
        hit_annotation: -1,
      }
    } else {
      raise ExonerateTextError(
        "Four-row Exonerate alignment requires a protein2* or *2protein model",
      )
    }
  } else if count == 5 {
    ExonerateTextRowMap::{
      query: 1,
      similarity: 2,
      hit: 3,
      query_annotation: 0,
      hit_annotation: 4,
    }
  } else {
    raise ExonerateTextError(
      "Exonerate stitched alignment must contain 3, 4, or 5 rows",
    )
  }
}

///|
fn exonerate_text_build_fragments(
  rows : Array[String],
  row_map : ExonerateTextRowMap,
  blocks : Array[ExonerateTextSpan],
  query_gaps : Array[Int],
  hit_gaps : Array[Int],
  header : ExonerateTextHeader,
) -> ExonerateTextHsp raise ExonerateTextError {
  let raw_fragments : Array[ExonerateTextRawFragment] = []
  let query_split : Array[ExonerateTextRange] = []
  let hit_split : Array[ExonerateTextRange] = []
  let query_direction = exonerate_text_direction(header.query_strand)
  let hit_direction = exonerate_text_direction(header.hit_strand)
  let mut query_cursor = header.query_raw_start
  let mut hit_cursor = header.hit_raw_start
  let mut raw_phase_offset = 0
  for block_index = 0
      block_index < blocks.length()
      block_index = block_index + 1 {
    let span = blocks[block_index]
    let raw_query = rows[row_map.query][span.start:span.end].to_owned()
    let raw_hit = rows[row_map.hit][span.start:span.end].to_owned()
    let raw_similarity = rows[row_map.similarity][span.start:span.end].to_owned()
    let raw_query_annotation = if row_map.query_annotation >= 0 {
      Some(rows[row_map.query_annotation][span.start:span.end].to_owned())
    } else {
      None
    }
    let raw_hit_annotation = if row_map.hit_annotation >= 0 {
      Some(rows[row_map.hit_annotation][span.start:span.end].to_owned())
    } else {
      None
    }
    let query_scodon = exonerate_text_split_codon_moves(raw_query)
    let hit_scodon = exonerate_text_split_codon_moves(raw_hit)
    let query = exonerate_text_remove_braces(raw_query)
    let hit = exonerate_text_remove_braces(raw_hit)
    let similarity = exonerate_text_remove_braces(raw_similarity)
    let query_annotation = match raw_query_annotation {
      Some(value) => Some(exonerate_text_remove_braces(value))
      None => None
    }
    let hit_annotation = match raw_hit_annotation {
      Some(value) => Some(exonerate_text_remove_braces(value))
      None => None
    }
    if query.length() != hit.length() || query.length() != similarity.length() {
      raise ExonerateTextError(
        "Exonerate fragment rows have unequal widths after split-codon cleanup",
      )
    }
    match query_annotation {
      Some(value) =>
        if value.length() != query.length() {
          raise ExonerateTextError(
            "Exonerate query annotation width does not match its fragment",
          )
        }
      None => ()
    }
    match hit_annotation {
      Some(value) =>
        if value.length() != hit.length() {
          raise ExonerateTextError(
            "Exonerate hit annotation width does not match its fragment",
          )
        }
      None => ()
    }
    let block_query_start = query_cursor
    let block_hit_start = hit_cursor
    let shifts = exonerate_text_hash_runs(similarity)
    let mut column = 0
    for shift in shifts {
      if shift.start > column {
        let built = exonerate_text_make_raw_fragment(
          query[column:shift.start].to_owned(),
          hit[column:shift.start].to_owned(),
          similarity[column:shift.start].to_owned(),
          exonerate_text_option_slice(query_annotation, column, shift.start),
          exonerate_text_option_slice(hit_annotation, column, shift.start),
          query_cursor,
          hit_cursor,
          query_direction,
          hit_direction,
          raw_phase_offset,
        )
        query_cursor = built.1
        hit_cursor = built.2
        raw_phase_offset = raw_phase_offset + shift.start - column
        if built.0.query.length() > 0 || built.0.hit.length() > 0 {
          raw_fragments.push(built.0)
        }
      }
      let use_annotations = query_annotation is Some(_) &&
        hit_annotation is Some(_)
      let qshift_source = if use_annotations {
        match query_annotation {
          Some(value) => value[shift.start:shift.end].to_owned()
          None => ""
        }
      } else {
        query[shift.start:shift.end].to_owned()
      }
      let hshift_source = if use_annotations {
        match hit_annotation {
          Some(value) => value[shift.start:shift.end].to_owned()
          None => ""
        }
      } else {
        hit[shift.start:shift.end].to_owned()
      }
      if qshift_source.contains("#") {
        query_cursor = query_cursor +
          (shift.end - shift.start) * query_direction
      } else if hshift_source.contains("#") {
        hit_cursor = hit_cursor + (shift.end - shift.start) * hit_direction
      } else {
        raise ExonerateTextError(
          "Exonerate frameshift marker has no query or hit owner",
        )
      }
      column = shift.end
    }
    if column < query.length() {
      let built = exonerate_text_make_raw_fragment(
        query[column:query.length()].to_owned(),
        hit[column:hit.length()].to_owned(),
        similarity[column:similarity.length()].to_owned(),
        exonerate_text_option_slice(query_annotation, column, query.length()),
        exonerate_text_option_slice(hit_annotation, column, hit.length()),
        query_cursor,
        hit_cursor,
        query_direction,
        hit_direction,
        raw_phase_offset,
      )
      query_cursor = built.1
      hit_cursor = built.2
      raw_phase_offset = raw_phase_offset + query.length() - column
      if built.0.query.length() > 0 || built.0.hit.length() > 0 {
        raw_fragments.push(built.0)
      }
    }
    let block_query_end = query_cursor
    let block_hit_end = hit_cursor
    exonerate_text_append_split_codons(
      query_split,
      query_scodon,
      block_query_start,
      block_query_end,
      header.query_strand,
    )
    exonerate_text_append_split_codons(
      hit_split,
      hit_scodon,
      block_hit_start,
      block_hit_end,
      header.hit_strand,
    )
    if block_index < query_gaps.length() {
      query_cursor = query_cursor + query_gaps[block_index] * query_direction
      hit_cursor = hit_cursor + hit_gaps[block_index] * hit_direction
    }
  }
  let expected_raw_query_end = if header.query_strand == ExonerateTextProtein {
    header.query_raw_start + (header.query_raw_end - header.query_raw_start) * 3
  } else {
    header.query_raw_end
  }
  let expected_raw_hit_end = if header.hit_strand == ExonerateTextProtein {
    header.hit_raw_start + (header.hit_raw_end - header.hit_raw_start) * 3
  } else {
    header.hit_raw_end
  }
  if query_cursor != expected_raw_query_end {
    raise ExonerateTextError(
      "Exonerate query fragment coordinates do not reach the header range",
    )
  }
  if hit_cursor != expected_raw_hit_end {
    raise ExonerateTextError(
      "Exonerate hit fragment coordinates do not reach the header range",
    )
  }
  let fragments = exonerate_text_finalize_fragments(raw_fragments, header)
  if fragments.length() == 0 {
    raise ExonerateTextError("Exonerate HSP contains no non-empty fragment")
  }
  ExonerateTextHsp::{
    score: header.score,
    model: header.model,
    query_start: exonerate_text_min(
      header.query_raw_start,
      header.query_raw_end,
    ),
    query_end: exonerate_text_max(header.query_raw_start, header.query_raw_end),
    hit_start: exonerate_text_min(header.hit_raw_start, header.hit_raw_end),
    hit_end: exonerate_text_max(header.hit_raw_start, header.hit_raw_end),
    query_strand: header.query_strand,
    hit_strand: header.hit_strand,
    fragments,
    query_split_codons: query_split,
    hit_split_codons: hit_split,
  }
}

///|
fn exonerate_text_make_raw_fragment(
  query : String,
  hit : String,
  similarity : String,
  query_annotation : String?,
  hit_annotation : String?,
  query_cursor : Int,
  hit_cursor : Int,
  query_direction : Int,
  hit_direction : Int,
  phase_offset : Int,
) -> (ExonerateTextRawFragment, Int, Int) {
  let query_units = exonerate_text_raw_units(query)
  let hit_units = exonerate_text_raw_units(hit)
  let query_next = query_cursor + query_units * query_direction
  let hit_next = hit_cursor + hit_units * hit_direction
  (
    ExonerateTextRawFragment::{
      query,
      hit,
      similarity,
      query_annotation,
      hit_annotation,
      query_start: exonerate_text_min(query_cursor, query_next),
      query_end: exonerate_text_max(query_cursor, query_next),
      hit_start: exonerate_text_min(hit_cursor, hit_next),
      hit_end: exonerate_text_max(hit_cursor, hit_next),
      phase: (3 - phase_offset % 3) % 3,
    },
    query_next,
    hit_next,
  )
}

///|
fn exonerate_text_finalize_fragments(
  raw : Array[ExonerateTextRawFragment],
  header : ExonerateTextHeader,
) -> Array[ExonerateTextFragment] raise ExonerateTextError {
  let translated = header.query_strand == ExonerateTextProtein ||
    header.hit_strand == ExonerateTextProtein
  let result : Array[ExonerateTextFragment] = []
  let mut protein_cursor = if header.query_strand == ExonerateTextProtein {
    header.query_raw_start
  } else {
    header.hit_raw_start
  }
  let protein_direction = 1
  for item in raw {
    let query_sequence = if translated {
      exonerate_text_aa_sequence(item.query, item.phase)
    } else {
      item.query
    }
    let hit_sequence = if translated {
      exonerate_text_aa_sequence(item.hit, item.phase)
    } else {
      item.hit
    }
    let similarity = if translated {
      exonerate_text_triplet_annotation(item.similarity, item.phase)
    } else {
      item.similarity
    }
    let protein_raw = if header.query_strand == ExonerateTextProtein {
      item.query
    } else {
      item.hit
    }
    let protein_leading_partial = if item.phase > 0 && protein_raw.length() > 0 {
      1
    } else {
      0
    }
    let protein_start_cursor = protein_cursor +
      protein_leading_partial * protein_direction
    let protein_units = exonerate_text_complete_residue_units(
      protein_raw,
      item.phase,
    )
    let protein_next = protein_start_cursor + protein_units * protein_direction
    let query_start = if header.query_strand == ExonerateTextProtein {
      exonerate_text_min(protein_start_cursor, protein_next)
    } else {
      item.query_start
    }
    let query_end = if header.query_strand == ExonerateTextProtein {
      exonerate_text_max(protein_start_cursor, protein_next)
    } else {
      item.query_end
    }
    let hit_start = if header.hit_strand == ExonerateTextProtein {
      exonerate_text_min(protein_start_cursor, protein_next)
    } else {
      item.hit_start
    }
    let hit_end = if header.hit_strand == ExonerateTextProtein {
      exonerate_text_max(protein_start_cursor, protein_next)
    } else {
      item.hit_end
    }
    let query_step = if header.query_strand == ExonerateTextProtein {
      1
    } else if translated {
      3
    } else {
      1
    }
    let hit_step = if header.hit_strand == ExonerateTextProtein {
      1
    } else if translated {
      3
    } else {
      1
    }
    result.push(ExonerateTextFragment::{
      query_sequence,
      hit_sequence,
      similarity,
      query_annotation: item.query_annotation,
      hit_annotation: item.hit_annotation,
      query_start,
      query_end,
      hit_start,
      hit_end,
      query_strand: header.query_strand,
      hit_strand: header.hit_strand,
      query_frame: exonerate_text_frame(query_start, header.query_strand),
      hit_frame: exonerate_text_frame(hit_start, header.hit_strand),
      phase: item.phase,
      query_step,
      hit_step,
    })
    protein_cursor = protein_next
  }
  let expected_protein_end = if header.query_strand == ExonerateTextProtein {
    header.query_raw_end
  } else {
    header.hit_raw_end
  }
  if translated && protein_cursor != expected_protein_end {
    raise ExonerateTextError(
      "Translated Exonerate residues do not reach the protein header range",
    )
  }
  result
}

///|
fn exonerate_text_aggregate(
  parsed : Array[ExonerateTextParsedAlignment],
) -> Array[ExonerateTextQueryResult] raise ExonerateTextError {
  let queries : Array[ExonerateTextQueryResult] = []
  for alignment in parsed {
    let mut query_index = -1
    for index = 0; index < queries.length(); index = index + 1 {
      if queries[index].id == alignment.query_id {
        query_index = index
        break
      }
    }
    if query_index < 0 {
      queries.push(ExonerateTextQueryResult::{
        id: alignment.query_id,
        description: alignment.query_description,
        model: alignment.model,
        hits: [
          ExonerateTextHit::{
            id: alignment.hit_id,
            description: alignment.hit_description,
            hsps: [alignment.hsp],
          },
        ],
      })
    } else {
      if queries[query_index].description != alignment.query_description ||
        queries[query_index].model != alignment.model {
        raise ExonerateTextError(
          "Repeated Exonerate query has inconsistent description or model",
        )
      }
      let mut hit_index = -1
      for index = 0
          index < queries[query_index].hits.length()
          index = index + 1 {
        if queries[query_index].hits[index].id == alignment.hit_id {
          hit_index = index
          break
        }
      }
      if hit_index < 0 {
        queries[query_index].hits.push(ExonerateTextHit::{
          id: alignment.hit_id,
          description: alignment.hit_description,
          hsps: [alignment.hsp],
        })
      } else {
        if queries[query_index].hits[hit_index].description !=
          alignment.hit_description {
          raise ExonerateTextError(
            "Repeated Exonerate target has an inconsistent description",
          )
        }
        queries[query_index].hits[hit_index].hsps.push(alignment.hsp)
      }
    }
  }
  let copied : Array[ExonerateTextQueryResult] = []
  for query in queries {
    copied.push(exonerate_text_copy_query(query))
  }
  copied
}

///|
fn exonerate_text_find_ner_spans(
  sequence : String,
) -> Array[ExonerateTextSpan] raise ExonerateTextError {
  let result : Array[ExonerateTextSpan] = []
  let mut cursor = 0
  while cursor + 2 < sequence.length() {
    let start = exonerate_text_find_substring(sequence, "--<", cursor)
    if start < 0 {
      break
    }
    let close = exonerate_text_find_substring(sequence, ">--", start + 3)
    if close < 0 {
      raise ExonerateTextError("Unterminated Exonerate NER marker")
    }
    result.push(ExonerateTextSpan::{ start, end: close + 3 })
    cursor = close + 3
  }
  result
}

///|
fn exonerate_text_find_intron_spans(
  sequence : String,
) -> Array[ExonerateTextSpan] raise ExonerateTextError {
  let result : Array[ExonerateTextSpan] = []
  let mut cursor = 0
  while cursor < sequence.length() {
    let intron = exonerate_text_find_substring(sequence, "Intron", cursor)
    let dots = exonerate_text_find_dot_run(sequence, cursor)
    if intron < 0 && dots.0 < 0 {
      break
    }
    if intron >= 0 && (dots.0 < 0 || intron < dots.0) {
      let mut start = intron
      while start > 0 &&
            sequence.unsafe_get(start).to_int() != '<'.to_int() &&
            sequence.unsafe_get(start).to_int() != '>'.to_int() {
        start = start - 1
      }
      while start > 0 &&
            sequence.unsafe_get(start - 1).to_int() ==
            sequence.unsafe_get(start).to_int() {
        start = start - 1
      }
      let mut end = intron + 6
      while end < sequence.length() &&
            sequence.unsafe_get(end).to_int() != '<'.to_int() &&
            sequence.unsafe_get(end).to_int() != '>'.to_int() {
        end = end + 1
      }
      if end >= sequence.length() {
        raise ExonerateTextError("Unterminated Exonerate intron marker")
      }
      let arrow = sequence.unsafe_get(end).to_int()
      while end < sequence.length() &&
            sequence.unsafe_get(end).to_int() == arrow {
        end = end + 1
      }
      let expanded = exonerate_text_expand_splice_flanks(sequence, start, end)
      result.push(expanded)
      cursor = expanded.end
    } else {
      let expanded = exonerate_text_expand_splice_flanks(
        sequence,
        dots.0,
        dots.1,
      )
      result.push(expanded)
      cursor = expanded.end
    }
  }
  result
}

///|
fn exonerate_text_expand_splice_flanks(
  sequence : String,
  raw_start : Int,
  raw_end : Int,
) -> ExonerateTextSpan {
  let mut start = raw_start
  let mut moved = 0
  while start > 0 &&
        moved < 2 &&
        exonerate_text_is_splice_flank(sequence.unsafe_get(start - 1).to_int()) {
    start = start - 1
    moved = moved + 1
  }
  let mut end = raw_end
  moved = 0
  while end < sequence.length() &&
        moved < 2 &&
        exonerate_text_is_splice_flank(sequence.unsafe_get(end).to_int()) {
    end = end + 1
    moved = moved + 1
  }
  ExonerateTextSpan::{ start, end }
}

///|
fn exonerate_text_complement_spans(
  length : Int,
  separators : Array[ExonerateTextSpan],
) -> Array[ExonerateTextSpan] raise ExonerateTextError {
  let result : Array[ExonerateTextSpan] = []
  let mut cursor = 0
  for separator in separators {
    if separator.start < cursor ||
      separator.end <= separator.start ||
      separator.end > length {
      raise ExonerateTextError("Overlapping or invalid Exonerate separator")
    }
    if separator.start > cursor {
      result.push(ExonerateTextSpan::{ start: cursor, end: separator.start })
    }
    cursor = separator.end
  }
  if cursor < length {
    result.push(ExonerateTextSpan::{ start: cursor, end: length })
  }
  result
}

///|
fn exonerate_text_bp_lengths(
  text : String,
) -> Array[Int] raise ExonerateTextError {
  let fields = exonerate_text_split_whitespace(text)
  let result : Array[Int] = []
  for index = 0; index + 1 < fields.length(); index = index + 1 {
    if fields[index + 1] == "bp" {
      result.push(
        exonerate_text_parse_nonnegative_int(fields[index], "intron length"),
      )
    }
  }
  result
}

///|
fn exonerate_text_parse_ner_gap(
  text : String,
  axis : String,
) -> Int raise ExonerateTextError {
  let mut start = -1
  let mut end = -1
  for index = 0; index < text.length(); index = index + 1 {
    let code = text.unsafe_get(index).to_int()
    if code >= '0'.to_int() && code <= '9'.to_int() {
      if start < 0 {
        start = index
      }
      end = index + 1
    } else if start >= 0 {
      break
    }
  }
  if start < 0 || end <= start {
    raise ExonerateTextError(
      "Exonerate " + axis + " NER marker lacks a gap length",
    )
  }
  exonerate_text_parse_nonnegative_int(
    text[start:end].to_owned(),
    axis + " NER length",
  )
}

///|
fn exonerate_text_hash_runs(text : String) -> Array[ExonerateTextSpan] {
  let result : Array[ExonerateTextSpan] = []
  let mut index = 0
  while index < text.length() {
    if text.unsafe_get(index).to_int() == '#'.to_int() {
      let start = index
      while index < text.length() &&
            text.unsafe_get(index).to_int() == '#'.to_int() {
        index = index + 1
      }
      result.push(ExonerateTextSpan::{ start, end: index })
    } else {
      index = index + 1
    }
  }
  result
}

///|
fn exonerate_text_split_codon_moves(sequence : String) -> (Int, Int) {
  let mut suffix = 0
  if sequence.length() >= 3 &&
    sequence.unsafe_get(sequence.length() - 1).to_int() == '}'.to_int() {
    let mut index = sequence.length() - 2
    while index >= 0 && sequence.unsafe_get(index).to_int() != '{'.to_int() {
      suffix = suffix + 1
      index = index - 1
    }
    if index < 0 || suffix < 1 || suffix > 2 {
      suffix = 0
    }
  }
  let mut prefix = 0
  if sequence.length() >= 3 && sequence.unsafe_get(0).to_int() == '{'.to_int() {
    let mut index = 1
    while index < sequence.length() &&
          sequence.unsafe_get(index).to_int() != '}'.to_int() {
      prefix = prefix + 1
      index = index + 1
    }
    if index >= sequence.length() || prefix < 1 || prefix > 2 {
      prefix = 0
    }
  }
  (suffix, prefix)
}

///|
fn exonerate_text_append_split_codons(
  output : Array[ExonerateTextRange],
  moves : (Int, Int),
  raw_start : Int,
  raw_end : Int,
  strand : ExonerateTextStrand,
) -> Unit {
  let direction = exonerate_text_direction(strand)
  if moves.0 > 0 {
    let anchor = raw_end
    let other = anchor - moves.0 * direction
    output.push(ExonerateTextRange::{
      start: exonerate_text_min(anchor, other),
      end: exonerate_text_max(anchor, other),
    })
  }
  if moves.1 > 0 {
    let anchor = raw_start
    let other = anchor + moves.1 * direction
    output.push(ExonerateTextRange::{
      start: exonerate_text_min(anchor, other),
      end: exonerate_text_max(anchor, other),
    })
  }
}

///|
fn exonerate_text_aa_sequence(sequence : String, phase : Int) -> String {
  let output = StringBuilder::new()
  let mut index = 0
  if phase > 0 && sequence.length() > 0 {
    output.write_char('X')
    index = if phase < sequence.length() { phase } else { sequence.length() }
  }
  while index + 3 <= sequence.length() {
    let triplet = sequence[index:index + 3].to_owned()
    output.write_char(exonerate_text_aa_code(triplet))
    index = index + 3
  }
  if index < sequence.length() {
    output.write_char('X')
  }
  output.to_string()
}

///|
fn exonerate_text_triplet_annotation(
  annotation : String,
  phase : Int,
) -> String {
  let output = StringBuilder::new()
  let mut index = 0
  if phase > 0 && annotation.length() > 0 {
    output.write_char(
      exonerate_text_annotation_code(
        annotation[0:exonerate_text_min(phase, annotation.length())].to_owned(),
      ),
    )
    index = exonerate_text_min(phase, annotation.length())
  }
  while index + 3 <= annotation.length() {
    output.write_char(
      exonerate_text_annotation_code(annotation[index:index + 3].to_owned()),
    )
    index = index + 3
  }
  if index < annotation.length() {
    output.write_char(
      exonerate_text_annotation_code(
        annotation[index:annotation.length()].to_owned(),
      ),
    )
  }
  output.to_string()
}

///|
fn exonerate_text_aa_code(value : String) -> Char {
  match value {
    "Ala" => 'A'
    "Arg" => 'R'
    "Asn" => 'N'
    "Asp" => 'D'
    "Cys" => 'C'
    "Gln" => 'Q'
    "Glu" => 'E'
    "Gly" => 'G'
    "His" => 'H'
    "Ile" => 'I'
    "Leu" => 'L'
    "Lys" => 'K'
    "Met" => 'M'
    "Phe" => 'F'
    "Pro" => 'P'
    "Ser" => 'S'
    "Thr" => 'T'
    "Trp" => 'W'
    "Tyr" => 'Y'
    "Val" => 'V'
    "Sec" => 'U'
    "Pyl" => 'O'
    "Asx" => 'B'
    "Glx" => 'Z'
    "Xle" => 'J'
    "Unk" => 'X'
    "***" => '*'
    "<->" => '-'
    "---" => '-'
    _ => 'X'
  }
}

///|
fn exonerate_text_annotation_code(value : String) -> Char {
  if value.contains("|") {
    '|'
  } else if value.contains(":") {
    ':'
  } else if value.contains(".") {
    '.'
  } else if value.contains("!") {
    '!'
  } else if value.contains("#") {
    '#'
  } else {
    ' '
  }
}

///|
fn exonerate_text_frame(start : Int, strand : ExonerateTextStrand) -> Int {
  match strand {
    ExonerateTextProtein => 0
    ExonerateTextForward => start % 3 + 1
    ExonerateTextReverse => -(start % 3 + 1)
  }
}

///|
fn exonerate_text_raw_units(sequence : String) -> Int {
  let mut count = 0
  for index = 0; index < sequence.length(); index = index + 1 {
    let code = sequence.unsafe_get(index).to_int()
    if code != '-'.to_int() &&
      code != '<'.to_int() &&
      code != '>'.to_int() &&
      code != ' '.to_int() &&
      code != '\t'.to_int() {
      count = count + 1
    }
  }
  count
}

///|
fn exonerate_text_complete_residue_units(sequence : String, phase : Int) -> Int {
  let mut index = exonerate_text_min(phase, sequence.length())
  let mut count = 0
  while index + 3 <= sequence.length() {
    let triplet = sequence[index:index + 3].to_owned()
    if exonerate_text_aa_code(triplet) != '-' {
      count = count + 1
    }
    index = index + 3
  }
  count
}

///|
fn exonerate_text_inter_ranges(
  fragments : Array[ExonerateTextFragment],
  query : Bool,
) -> Array[ExonerateTextRange] {
  let result : Array[ExonerateTextRange] = []
  if fragments.length() < 2 {
    return result
  }
  for index = 0; index + 1 < fragments.length(); index = index + 1 {
    let strand = if query {
      fragments[index].query_strand
    } else {
      fragments[index].hit_strand
    }
    let first_start = if query {
      fragments[index].query_start
    } else {
      fragments[index].hit_start
    }
    let first_end = if query {
      fragments[index].query_end
    } else {
      fragments[index].hit_end
    }
    let second_start = if query {
      fragments[index + 1].query_start
    } else {
      fragments[index + 1].hit_start
    }
    let second_end = if query {
      fragments[index + 1].query_end
    } else {
      fragments[index + 1].hit_end
    }
    if strand == ExonerateTextReverse {
      result.push(ExonerateTextRange::{ start: second_end, end: first_start })
    } else {
      result.push(ExonerateTextRange::{ start: first_end, end: second_start })
    }
  }
  result
}

///|
fn exonerate_text_fragment_position(
  sequence : String,
  aligned_index : Int,
  start : Int,
  end : Int,
  strand : ExonerateTextStrand,
  step : Int,
) -> Int? {
  if aligned_index < 0 || aligned_index >= sequence.length() {
    return None
  }
  if sequence.unsafe_get(aligned_index).to_int() == '-'.to_int() {
    return None
  }
  let mut residues = 0
  for index = 0; index < aligned_index; index = index + 1 {
    if sequence.unsafe_get(index).to_int() != '-'.to_int() {
      residues = residues + 1
    }
  }
  if strand == ExonerateTextReverse {
    Some(end - (residues + 1) * step)
  } else {
    Some(start + residues * step)
  }
}

///|
fn exonerate_text_similarity_positive(
  similarity : String,
  index : Int,
  sequence_width : Int,
) -> Bool {
  if similarity.length() == sequence_width && index < similarity.length() {
    let code = similarity.unsafe_get(index).to_int()
    code == '|'.to_int() || code == ':'.to_int() || code == '.'.to_int()
  } else {
    false
  }
}

///|
fn exonerate_text_coordinate_layout(line : String) -> (Int, Int)? {
  let mut first = -1
  let mut last = -1
  for index = 0; index < line.length(); index = index + 1 {
    if line.unsafe_get(index).to_int() == ':'.to_int() {
      if first < 0 {
        first = index
      }
      last = index
    }
  }
  if first <= 0 || last <= first {
    return None
  }
  let left = line[0:first].to_owned().trim().to_owned()
  let right = line[last + 1:line.length()].to_owned().trim().to_owned()
  if !exonerate_text_is_digits(left) || !exonerate_text_is_digits(right) {
    return None
  }
  let mut start = first + 1
  if start < last && line.unsafe_get(start).to_int() == ' '.to_int() {
    start = start + 1
  }
  let mut end = last
  if end > start && line.unsafe_get(end - 1).to_int() == ' '.to_int() {
    end = end - 1
  }
  if end <= start {
    None
  } else {
    Some((start, end - start))
  }
}

///|
fn exonerate_text_fixed_slice(
  line : String,
  start : Int,
  width : Int,
) -> String raise ExonerateTextError {
  if start < 0 || width <= 0 || start + width > line.length() {
    raise ExonerateTextError(
      "Exonerate alignment row is shorter than its coordinate-row width",
    )
  }
  line[start:start + width].to_owned()
}

///|
fn exonerate_text_flip_codons(
  codon : String,
  target : String,
) -> (String, String) raise ExonerateTextError {
  if codon.length() != target.length() {
    raise ExonerateTextError("Exonerate codon rows have unequal widths")
  }
  let first = StringBuilder::new(size_hint=codon.length())
  let second = StringBuilder::new(size_hint=codon.length())
  for index = 0; index < codon.length(); index = index + 1 {
    let c = codon.unsafe_get(index).to_int()
    let t = target.unsafe_get(index).to_int()
    if c == ' '.to_int() {
      first.write_char(c.unsafe_to_char())
      second.write_char(t.unsafe_to_char())
    } else {
      first.write_char(t.unsafe_to_char())
      second.write_char(c.unsafe_to_char())
    }
  }
  (first.to_string(), second.to_string())
}

///|
fn exonerate_text_metadata_value(
  lines : Array[String],
  prefix : String,
) -> String raise ExonerateTextError {
  for line in lines {
    let trimmed = line.trim().to_owned()
    if trimmed.has_prefix(prefix) {
      let value = exonerate_text_after_prefix(trimmed, prefix)
      if value.length() >= 2 &&
        value.unsafe_get(0).to_int() == '['.to_int() &&
        value.unsafe_get(value.length() - 1).to_int() == ']'.to_int() {
        return value[1:value.length() - 1].to_owned()
      }
      if value.length() == 0 {
        raise ExonerateTextError(
          "Exonerate " + prefix + " metadata value is empty",
        )
      }
      return value
    }
    if trimmed == "C4 Alignment:" {
      break
    }
  }
  raise ExonerateTextError("Exonerate report is missing " + prefix)
}

///|
fn exonerate_text_id_description(
  value : String,
  label : String,
) -> (String, String) raise ExonerateTextError {
  let trimmed = value.trim().to_owned()
  if trimmed.length() == 0 {
    raise ExonerateTextError("Exonerate " + label + " identifier is empty")
  }
  let mut split = trimmed.length()
  for index = 0; index < trimmed.length(); index = index + 1 {
    if exonerate_text_is_whitespace(trimmed.unsafe_get(index).to_int()) {
      split = index
      break
    }
  }
  let id = trimmed[0:split].to_owned()
  let description = if split < trimmed.length() {
    trimmed[split:trimmed.length()].to_owned().trim().to_owned()
  } else {
    ""
  }
  (id, description)
}

///|
fn exonerate_text_strand(
  description : String,
  protein : Bool,
) -> (ExonerateTextStrand, String) {
  if protein {
    return (ExonerateTextProtein, description)
  }
  if description.has_suffix(":[revcomp]") {
    (
      ExonerateTextReverse,
      description[0:description.length() - 10].to_owned().trim().to_owned(),
    )
  } else if description.has_suffix("[revcomp]") {
    (
      ExonerateTextReverse,
      description[0:description.length() - 9].to_owned().trim().to_owned(),
    )
  } else {
    (ExonerateTextForward, description)
  }
}

///|
fn exonerate_text_parse_range(
  value : String,
  label : String,
) -> (Int, Int) raise ExonerateTextError {
  let arrow = exonerate_text_find_substring(value, "->", 0)
  if arrow < 0 {
    raise ExonerateTextError("Exonerate " + label + " range must contain ->")
  }
  let left = value[0:arrow].to_owned().trim().to_owned()
  let right = value[arrow + 2:value.length()].to_owned().trim().to_owned()
  (
    exonerate_text_parse_nonnegative_int(left, label + " range start"),
    exonerate_text_parse_nonnegative_int(right, label + " range end"),
  )
}

///|
fn exonerate_text_validate_orientation(
  start : Int,
  end : Int,
  strand : ExonerateTextStrand,
  label : String,
) -> Unit raise ExonerateTextError {
  if start == end {
    raise ExonerateTextError(
      "Exonerate " + label + " range must have positive length",
    )
  }
  match strand {
    ExonerateTextReverse =>
      if end >= start {
        raise ExonerateTextError(
          "Reverse Exonerate " + label + " range must decrease",
        )
      }
    ExonerateTextForward | ExonerateTextProtein =>
      if end <= start {
        raise ExonerateTextError(
          "Forward/protein Exonerate " + label + " range must increase",
        )
      }
  }
}

///|
fn exonerate_text_parse_nonnegative_int(
  text : String,
  label : String,
) -> Int raise ExonerateTextError {
  if text.length() == 0 {
    raise ExonerateTextError(label + " is empty")
  }
  let mut value = 0
  for index = 0; index < text.length(); index = index + 1 {
    let code = text.unsafe_get(index).to_int()
    if code < '0'.to_int() || code > '9'.to_int() {
      raise ExonerateTextError(label + " must be a non-negative integer")
    }
    let digit = code - '0'.to_int()
    if value > 214748364 || (value == 214748364 && digit > 7) {
      raise ExonerateTextError(
        label + " is outside the supported integer range",
      )
    }
    value = value * 10 + digit
  }
  value
}

///|
fn exonerate_text_find_substring(
  text : String,
  pattern : String,
  from : Int,
) -> Int {
  if pattern.length() == 0 {
    return from
  }
  let mut index = if from < 0 { 0 } else { from }
  while index + pattern.length() <= text.length() {
    let mut matches = true
    for offset = 0; offset < pattern.length(); offset = offset + 1 {
      if text.unsafe_get(index + offset).to_int() !=
        pattern.unsafe_get(offset).to_int() {
        matches = false
        break
      }
    }
    if matches {
      return index
    }
    index = index + 1
  }
  -1
}

///|
fn exonerate_text_find_dot_run(text : String, from : Int) -> (Int, Int) {
  let mut index = if from < 0 { 0 } else { from }
  while index < text.length() {
    if text.unsafe_get(index).to_int() == '.'.to_int() {
      let start = index
      while index < text.length() &&
            text.unsafe_get(index).to_int() == '.'.to_int() {
        index = index + 1
      }
      if index - start >= 2 {
        return (start, index)
      }
    } else {
      index = index + 1
    }
  }
  (-1, -1)
}

///|
fn exonerate_text_remove_braces(text : String) -> String {
  let output = StringBuilder::new(size_hint=text.length())
  for index = 0; index < text.length(); index = index + 1 {
    let code = text.unsafe_get(index).to_int()
    if code != '{'.to_int() && code != '}'.to_int() {
      output.write_char(code.unsafe_to_char())
    }
  }
  output.to_string()
}

///|
fn exonerate_text_option_slice(
  value : String?,
  start : Int,
  end : Int,
) -> String? {
  match value {
    Some(text) => Some(text[start:end].to_owned())
    None => None
  }
}

///|
fn exonerate_text_after_prefix(value : String, prefix : String) -> String {
  value[prefix.length():value.length()].to_owned().trim().to_owned()
}

///|
fn exonerate_text_direction(strand : ExonerateTextStrand) -> Int {
  if strand == ExonerateTextReverse {
    -1
  } else {
    1
  }
}

///|
fn exonerate_text_is_splice_flank(code : Int) -> Bool {
  code == ' '.to_int() ||
  code == 'a'.to_int() ||
  code == 't'.to_int() ||
  code == 'g'.to_int() ||
  code == 'c'.to_int()
}

///|
fn exonerate_text_is_digits(value : String) -> Bool {
  if value.length() == 0 {
    return false
  }
  for index = 0; index < value.length(); index = index + 1 {
    let code = value.unsafe_get(index).to_int()
    if code < '0'.to_int() || code > '9'.to_int() {
      return false
    }
  }
  true
}

///|
fn exonerate_text_split_whitespace(value : String) -> Array[String] {
  let result : Array[String] = []
  let mut start = 0
  let mut active = false
  for index = 0; index < value.length(); index = index + 1 {
    if exonerate_text_is_whitespace(value.unsafe_get(index).to_int()) {
      if active {
        result.push(value[start:index].to_owned())
        active = false
      }
    } else if !active {
      start = index
      active = true
    }
  }
  if active {
    result.push(value[start:value.length()].to_owned())
  }
  result
}

///|
fn exonerate_text_is_whitespace(code : Int) -> Bool {
  code == ' '.to_int() ||
  code == '\t'.to_int() ||
  code == '\r'.to_int() ||
  code == '\n'.to_int()
}

///|
fn exonerate_text_upper(code : Int) -> Int {
  if code >= 'a'.to_int() && code <= 'z'.to_int() {
    code - 'a'.to_int() + 'A'.to_int()
  } else {
    code
  }
}

///|
fn exonerate_text_min(first : Int, second : Int) -> Int {
  if first < second {
    first
  } else {
    second
  }
}

///|
fn exonerate_text_max(first : Int, second : Int) -> Int {
  if first > second {
    first
  } else {
    second
  }
}

///|
fn exonerate_text_lines(text : String) -> Array[String] {
  let raw = text.split("\n").to_array()
  let lines : Array[String] = []
  for line in raw {
    if line.length() > 0 &&
      line.unsafe_get(line.length() - 1).to_int() == '\r'.to_int() {
      lines.push(line[0:line.length() - 1].to_owned())
    } else {
      lines.push(line.to_owned())
    }
  }
  lines
}

///|
fn exonerate_text_copy_hsp(hsp : ExonerateTextHsp) -> ExonerateTextHsp {
  let fragments : Array[ExonerateTextFragment] = []
  let query_split_codons : Array[ExonerateTextRange] = []
  let hit_split_codons : Array[ExonerateTextRange] = []
  for fragment in hsp.fragments {
    fragments.push(fragment)
  }
  for range in hsp.query_split_codons {
    query_split_codons.push(range)
  }
  for range in hsp.hit_split_codons {
    hit_split_codons.push(range)
  }
  ExonerateTextHsp::{
    score: hsp.score,
    model: hsp.model,
    query_start: hsp.query_start,
    query_end: hsp.query_end,
    hit_start: hsp.hit_start,
    hit_end: hsp.hit_end,
    query_strand: hsp.query_strand,
    hit_strand: hsp.hit_strand,
    fragments,
    query_split_codons,
    hit_split_codons,
  }
}

///|
fn exonerate_text_copy_hit(hit : ExonerateTextHit) -> ExonerateTextHit {
  let hsps : Array[ExonerateTextHsp] = []
  for hsp in hit.hsps {
    hsps.push(exonerate_text_copy_hsp(hsp))
  }
  ExonerateTextHit::{ id: hit.id, description: hit.description, hsps }
}

///|
fn exonerate_text_copy_query(
  query : ExonerateTextQueryResult,
) -> ExonerateTextQueryResult {
  let hits : Array[ExonerateTextHit] = []
  for hit in query.hits {
    hits.push(exonerate_text_copy_hit(hit))
  }
  ExonerateTextQueryResult::{
    id: query.id,
    description: query.description,
    model: query.model,
    hits,
  }
}
