///|
/// Parser and data model for HH-suite HHR output.
///
/// Mirrors the semantics of Biopython's `Bio.Align.hhr` module. An HHR file
/// describes one query profile and a ranked collection of pairwise
/// profile-profile alignments produced by HHsearch or HHblits.

///|
/// Error raised for malformed or inconsistent HHR input.
pub suberror HhrError {
  HhrError(String)
}

///|
/// File-level metadata from the HHR header.
pub struct HhrMetadata {
  query_name : String
  match_columns : Int
  sequence_count : Int
  sequence_count_total : Int
  neff : Double
  template_neff : Double?
  searched_hmms : Int
  run_date : String
  command_line : String
} derive(Eq)
pub extend HhrMetadata with Eq::{equal, not_equal}

///|
/// Create HHR file-level metadata.
pub fn HhrMetadata::new(
  query_name~ : String,
  match_columns~ : Int,
  sequence_count? : Int = 0,
  sequence_count_total? : Int = 0,
  neff? : Double = 0.0,
  template_neff? : Double? = None,
  searched_hmms? : Int = 0,
  run_date? : String = "",
  command_line? : String = "",
) -> HhrMetadata {
  HhrMetadata::{
    query_name,
    match_columns,
    sequence_count,
    sequence_count_total,
    neff,
    template_neff,
    searched_hmms,
    run_date,
    command_line,
  }
}

///|
/// A row in the ranked HHR hit summary table.
///
/// Coordinates are 0-based and end-exclusive, matching Biopython's alignment
/// coordinate convention.
pub struct HhrHitSummary {
  rank : Int
  target_name : String
  description : String
  probability : Double
  evalue : Double
  pvalue : Double
  score : Double
  secondary_structure_score : Double
  aligned_columns : Int
  query_start : Int
  query_end : Int
  target_start : Int
  target_end : Int
  target_length : Int
} derive(Eq)
pub extend HhrHitSummary with Eq::{equal, not_equal}

///|
/// One detailed pairwise profile alignment from an HHR file.
///
/// Sequence and annotation strings contain alignment gaps and therefore all
/// non-empty annotation strings have the same length as `query_sequence`.
pub struct HhrAlignment {
  rank : Int
  query_name : String
  query_length : Int
  target_name : String
  target_description : String
  target_length : Int
  probability : Double
  evalue : Double
  pvalue : Double
  score : Double
  secondary_structure_score : Double
  identities : Double
  similarity : Double
  sum_probabilities : Double
  aligned_columns : Int
  query_start : Int
  query_end : Int
  target_start : Int
  target_end : Int
  query_sequence : String
  target_sequence : String
  query_consensus : String
  target_consensus : String
  query_secondary_structure : String
  target_secondary_structure : String
  target_dssp : String
  column_score : String
  confidence : String
} derive(Eq)
pub extend HhrAlignment with Eq::{equal, not_equal}

///|
/// A complete HHR document.
pub struct HhrRecord {
  metadata : HhrMetadata
  hit_summaries : Array[HhrHitSummary]
  alignments : Array[HhrAlignment]
} derive(Eq)
pub extend HhrRecord with Eq::{equal, not_equal}

///|
/// Return the number of detailed hit alignments.
pub fn HhrRecord::num_hits(self : HhrRecord) -> Int {
  self.alignments.length()
}

///|
/// Return an alignment by zero-based result index.
pub fn HhrRecord::get(self : HhrRecord, index : Int) -> HhrAlignment? {
  if index < 0 || index >= self.alignments.length() {
    None
  } else {
    Some(self.alignments[index])
  }
}

///|
/// Return all alignments for a target identifier.
///
/// HH-suite can report multiple local alignments for the same target.
pub fn HhrRecord::alignments_for_target(
  self : HhrRecord,
  target_name : String,
) -> Array[HhrAlignment] {
  let result : Array[HhrAlignment] = []
  for alignment in self.alignments {
    if alignment.target_name == target_name {
      result.push(alignment)
    }
  }
  result
}

///|
/// Return the first ranked alignment for a target identifier.
pub fn HhrRecord::find_target(
  self : HhrRecord,
  target_name : String,
) -> HhrAlignment? {
  for alignment in self.alignments {
    if alignment.target_name == target_name {
      return Some(alignment)
    }
  }
  None
}

///|
/// Return the best alignment by probability, then E-value, then rank.
pub fn HhrRecord::best_alignment(self : HhrRecord) -> HhrAlignment? {
  if self.alignments.length() == 0 {
    return None
  }
  let mut best = self.alignments[0]
  for i = 1; i < self.alignments.length(); i = i + 1 {
    let candidate = self.alignments[i]
    if candidate.probability > best.probability ||
      (
        candidate.probability == best.probability &&
        candidate.evalue < best.evalue
      ) ||
      (
        candidate.probability == best.probability &&
        candidate.evalue == best.evalue &&
        candidate.rank < best.rank
      ) {
      best = candidate
    }
  }
  Some(best)
}

///|
/// Return alignments whose HH-suite probability meets a minimum percentage.
pub fn HhrRecord::filter_by_probability(
  self : HhrRecord,
  minimum : Double,
) -> Array[HhrAlignment] {
  let result : Array[HhrAlignment] = []
  for alignment in self.alignments {
    if alignment.probability >= minimum {
      result.push(alignment)
    }
  }
  result
}

///|
/// Return alignments whose E-value does not exceed a threshold.
pub fn HhrRecord::filter_by_evalue(
  self : HhrRecord,
  maximum : Double,
) -> Array[HhrAlignment] {
  let result : Array[HhrAlignment] = []
  for alignment in self.alignments {
    if alignment.evalue <= maximum {
      result.push(alignment)
    }
  }
  result
}

///|
/// Return the number of alignment columns, including gap columns.
pub fn HhrAlignment::alignment_length(self : HhrAlignment) -> Int {
  self.query_sequence.length()
}

///|
/// Return the aligned query sequence without gaps.
pub fn HhrAlignment::ungapped_query(self : HhrAlignment) -> String {
  hhr_remove_gaps(self.query_sequence)
}

///|
/// Return the aligned target sequence without gaps.
pub fn HhrAlignment::ungapped_target(self : HhrAlignment) -> String {
  hhr_remove_gaps(self.target_sequence)
}

///|
/// Return exact identical residue columns, excluding gaps.
pub fn HhrAlignment::identity_count(self : HhrAlignment) -> Int {
  let mut count = 0
  for i = 0; i < self.query_sequence.length(); i = i + 1 {
    let query = self.query_sequence.unsafe_get(i).to_int()
    let target = self.target_sequence.unsafe_get(i).to_int()
    if query != '-'.to_int() && target != '-'.to_int() && query == target {
      count = count + 1
    }
  }
  count
}

///|
/// Return query coverage as aligned query residues divided by query length.
pub fn HhrAlignment::query_coverage(self : HhrAlignment) -> Double {
  if self.query_length <= 0 {
    0.0
  } else {
    (self.query_end - self.query_start).to_double() /
    self.query_length.to_double()
  }
}

///|
/// Return target coverage as aligned target residues divided by target length.
pub fn HhrAlignment::target_coverage(self : HhrAlignment) -> Double {
  if self.target_length <= 0 {
    0.0
  } else {
    (self.target_end - self.target_start).to_double() /
    self.target_length.to_double()
  }
}

///|
/// Map an absolute zero-based query position to its aligned target position.
///
/// Returns `None` when the query position is outside this local alignment or
/// is aligned to a target gap.
pub fn HhrAlignment::query_to_target(
  self : HhrAlignment,
  query_position : Int,
) -> Int? {
  if query_position < self.query_start || query_position >= self.query_end {
    return None
  }
  let mut query_coordinate = self.query_start
  let mut target_coordinate = self.target_start
  for i = 0; i < self.query_sequence.length(); i = i + 1 {
    let query_gap = self.query_sequence.unsafe_get(i).to_int() == '-'.to_int()
    let target_gap = self.target_sequence.unsafe_get(i).to_int() == '-'.to_int()
    if !query_gap && query_coordinate == query_position {
      if target_gap {
        return None
      }
      return Some(target_coordinate)
    }
    if !query_gap {
      query_coordinate = query_coordinate + 1
    }
    if !target_gap {
      target_coordinate = target_coordinate + 1
    }
  }
  None
}

///|
/// Return per-column absolute coordinates.
///
/// A gap is represented by `None`; residue coordinates are zero-based.
pub fn HhrAlignment::aligned_pairs(self : HhrAlignment) -> Array[(Int?, Int?)] {
  let pairs : Array[(Int?, Int?)] = []
  let mut query_coordinate = self.query_start
  let mut target_coordinate = self.target_start
  for i = 0; i < self.query_sequence.length(); i = i + 1 {
    let query_gap = self.query_sequence.unsafe_get(i).to_int() == '-'.to_int()
    let target_gap = self.target_sequence.unsafe_get(i).to_int() == '-'.to_int()
    let query_value : Int? = if query_gap {
      None
    } else {
      Some(query_coordinate)
    }
    let target_value : Int? = if target_gap {
      None
    } else {
      Some(target_coordinate)
    }
    pairs.push((query_value, target_value))
    if !query_gap {
      query_coordinate = query_coordinate + 1
    }
    if !target_gap {
      target_coordinate = target_coordinate + 1
    }
  }
  pairs
}

///|
/// Return a compact alignment summary.
pub fn HhrAlignment::summary(self : HhrAlignment) -> String {
  "HhrAlignment(rank=" +
  self.rank.to_string() +
  ", target=" +
  self.target_name +
  ", probability=" +
  self.probability.to_string() +
  ", evalue=" +
  self.evalue.to_string() +
  ", columns=" +
  self.aligned_columns.to_string() +
  ")"
}

///|
/// Return a compact document summary.
pub fn HhrRecord::summary(self : HhrRecord) -> String {
  "HhrRecord(query=" +
  self.metadata.query_name +
  ", columns=" +
  self.metadata.match_columns.to_string() +
  ", hits=" +
  self.alignments.length().to_string() +
  ")"
}

///|
/// Parse one HHsearch or HHblits HHR document.
///
/// The parser validates header keys, contiguous ranks, summary/detail
/// agreement, block widths, coordinates, full sequence lengths, and all
/// per-column annotation lengths.
pub fn hhr_parse(text : String) -> HhrRecord raise HhrError {
  let lines = hhr_normalize_lines(text)
  if lines.length() == 0 {
    raise HhrError("Empty HHR input")
  }
  let expected_header = [
    "No", "Hit", "Prob", "E-value", "P-value", "Score", "SS", "Cols", "Query", "HMM",
    "Template", "HMM",
  ]
  let mut index = 0
  let mut query_name = ""
  let mut match_columns = -1
  let mut sequence_count = 0
  let mut sequence_count_total = 0
  let mut neff = 0.0
  let mut template_neff : Double? = None
  let mut searched_hmms = 0
  let mut run_date = ""
  let mut command_line = ""
  while index < lines.length() && lines[index].trim().length() > 0 {
    if hhr_split_whitespace(lines[index]) == expected_header {
      break
    }
    let line = lines[index].trim().to_owned()
    let split = hhr_first_whitespace(line)
    if split is None {
      raise HhrError("Malformed HHR header at line " + (index + 1).to_string())
    }
    let (key, value) = split.unwrap()
    match key {
      "Query" => query_name = value
      "Match_columns" =>
        match_columns = hhr_parse_nonnegative_int(
          value,
          "Match_columns",
          index + 1,
        )
      "No_of_seqs" => {
        let marker = " out of "
        let marker_index = hhr_find_substring(value, marker)
        if marker_index < 0 {
          raise HhrError(
            "Malformed No_of_seqs at line " + (index + 1).to_string(),
          )
        }
        sequence_count = hhr_parse_nonnegative_int(
          value[0:marker_index].to_owned(),
          "No_of_seqs",
          index + 1,
        )
        sequence_count_total = hhr_parse_nonnegative_int(
          value[marker_index + marker.length():].to_owned(),
          "No_of_seqs total",
          index + 1,
        )
      }
      "Neff" => neff = hhr_parse_number(value, "Neff", index + 1)
      "Template_Neff" =>
        template_neff = Some(
          hhr_parse_number(value, "Template_Neff", index + 1),
        )
      "Searched_HMMs" =>
        searched_hmms = hhr_parse_nonnegative_int(
          value,
          "Searched_HMMs",
          index + 1,
        )
      "Date" => run_date = value
      "Command" => command_line = value
      _ => raise HhrError("Unknown HHR header key '" + key + "'")
    }
    index = index + 1
  }
  if query_name.length() == 0 {
    raise HhrError("HHR header is missing Query")
  }
  if match_columns < 0 {
    raise HhrError("HHR header is missing Match_columns")
  }
  while index < lines.length() && lines[index].trim().length() == 0 {
    index = index + 1
  }
  if index >= lines.length() {
    raise HhrError("Truncated HHR file before hit table")
  }
  let table_header = hhr_split_whitespace(lines[index])
  if table_header != expected_header {
    raise HhrError("Invalid HHR hit table header")
  }
  index = index + 1
  let hit_summaries : Array[HhrHitSummary] = []
  while index < lines.length() {
    let trimmed = lines[index].trim().to_owned()
    if trimmed.length() == 0 {
      index = index + 1
      break
    }
    if hhr_is_rank_line(trimmed) || trimmed == "Done!" {
      break
    }
    let summary = hhr_parse_hit_summary(
      lines[index],
      hit_summaries.length() + 1,
      index + 1,
    )
    hit_summaries.push(summary)
    index = index + 1
  }
  while index < lines.length() && lines[index].trim().length() == 0 {
    index = index + 1
  }
  let alignments : Array[HhrAlignment] = []
  let mut current : HhrAlignmentBuilder? = None
  let mut found_done = false
  while index < lines.length() {
    let line = lines[index]
    let trimmed = line.trim().to_owned()
    if trimmed.length() == 0 {
      index = index + 1
      continue
    }
    if trimmed == "Done!" {
      match current {
        Some(builder) => {
          alignments.push(
            hhr_finish_alignment(
              builder, query_name, match_columns, hit_summaries,
            ),
          )
          current = None
        }
        None => ()
      }
      found_done = true
      index = index + 1
      while index < lines.length() {
        if lines[index].trim().length() > 0 {
          raise HhrError("Additional data found after Done!")
        }
        index = index + 1
      }
      break
    }
    if hhr_has_prefix(trimmed, "No ") {
      match current {
        Some(builder) =>
          alignments.push(
            hhr_finish_alignment(
              builder, query_name, match_columns, hit_summaries,
            ),
          )
        None => ()
      }
      let words = hhr_split_whitespace(trimmed)
      if words.length() != 2 || words[0] != "No" {
        raise HhrError(
          "Malformed alignment rank at line " + (index + 1).to_string(),
        )
      }
      let rank = hhr_parse_nonnegative_int(
        words[1],
        "alignment rank",
        index + 1,
      )
      if rank != alignments.length() + 1 {
        raise HhrError(
          "Non-contiguous alignment rank at line " + (index + 1).to_string(),
        )
      }
      current = Some(HhrAlignmentBuilder::new(rank))
      index = index + 1
      continue
    }
    let builder = match current {
      Some(value) => value
      None =>
        raise HhrError(
          "Alignment content before a No line at line " +
          (index + 1).to_string(),
        )
    }
    if hhr_has_prefix(line, ">") {
      if builder.target_name.length() > 0 {
        raise HhrError(
          "Duplicate target header at line " + (index + 1).to_string(),
        )
      }
      let target_header = line[1:].trim().to_owned()
      let split = hhr_first_whitespace(target_header)
      match split {
        Some((name, description)) => {
          builder.target_name = name
          builder.target_description = description
        }
        None => builder.target_name = target_header
      }
    } else if hhr_has_prefix(trimmed, "Probab=") {
      hhr_parse_alignment_statistics(builder, trimmed, index + 1)
    } else if hhr_has_prefix(line, "Q ss_pred ") {
      builder.query_secondary_structure = builder.query_secondary_structure +
        hhr_last_token(line, "query ss_pred", index + 1)
    } else if hhr_has_prefix(line, "Q Consensus ") {
      let part = hhr_parse_sequence_line(line, "Q Consensus", index + 1)
      builder.query_consensus = builder.query_consensus + part.sequence
      builder.block_width = part.sequence.length()
    } else if hhr_has_prefix(line, "Q ") {
      let part = hhr_parse_sequence_line(line, "Q", index + 1)
      if builder.query_start < 0 {
        builder.query_start = part.start
      } else if part.start != builder.query_end {
        raise HhrError(
          "Discontinuous query coordinates at line " + (index + 1).to_string(),
        )
      }
      builder.query_end = part.end
      if builder.query_length == 0 {
        builder.query_length = part.total
      } else if builder.query_length != part.total {
        raise HhrError(
          "Inconsistent query length at line " + (index + 1).to_string(),
        )
      }
      builder.query_sequence = builder.query_sequence + part.sequence
      builder.block_width = part.sequence.length()
    } else if hhr_has_prefix(line, "T ss_pred ") {
      builder.target_secondary_structure = builder.target_secondary_structure +
        hhr_last_token(line, "target ss_pred", index + 1)
    } else if hhr_has_prefix(line, "T ss_dssp ") {
      builder.target_dssp = builder.target_dssp +
        hhr_last_token(line, "target ss_dssp", index + 1)
    } else if hhr_has_prefix(line, "T Consensus ") {
      let part = hhr_parse_sequence_line(line, "T Consensus", index + 1)
      builder.target_consensus = builder.target_consensus + part.sequence
      builder.block_width = part.sequence.length()
    } else if hhr_has_prefix(line, "T ") {
      let part = hhr_parse_sequence_line(line, "T", index + 1)
      if builder.detail_target_name.length() == 0 {
        builder.detail_target_name = part.name
      } else if builder.detail_target_name != part.name {
        raise HhrError(
          "Inconsistent target name at line " + (index + 1).to_string(),
        )
      }
      if builder.target_start < 0 {
        builder.target_start = part.start
      } else if part.start != builder.target_end {
        raise HhrError(
          "Discontinuous target coordinates at line " + (index + 1).to_string(),
        )
      }
      builder.target_end = part.end
      if builder.target_length == 0 {
        builder.target_length = part.total
      } else if builder.target_length != part.total {
        raise HhrError(
          "Inconsistent target length at line " + (index + 1).to_string(),
        )
      }
      builder.target_sequence = builder.target_sequence + part.sequence
      builder.block_width = part.sequence.length()
    } else if hhr_has_prefix(line, "Confidence") {
      if builder.block_width <= 0 {
        raise HhrError(
          "Confidence before alignment block at line " + (index + 1).to_string(),
        )
      }
      builder.confidence = builder.confidence +
        hhr_suffix(line, builder.block_width, "Confidence", index + 1)
    } else if hhr_is_indented(line) {
      if builder.block_width <= 0 {
        raise HhrError(
          "Column score before alignment block at line " +
          (index + 1).to_string(),
        )
      }
      builder.column_score = builder.column_score +
        hhr_suffix(line, builder.block_width, "column score", index + 1)
    } else {
      raise HhrError(
        "Failed to parse HHR line " + (index + 1).to_string() + ": " + trimmed,
      )
    }
    index = index + 1
  }
  if !found_done {
    match current {
      Some(builder) =>
        alignments.push(
          hhr_finish_alignment(
            builder, query_name, match_columns, hit_summaries,
          ),
        )
      None => ()
    }
  }
  if alignments.length() != hit_summaries.length() {
    raise HhrError(
      "Expected " +
      hit_summaries.length().to_string() +
      " alignments, found " +
      alignments.length().to_string(),
    )
  }
  let metadata = HhrMetadata::new(
    query_name~,
    match_columns~,
    sequence_count~,
    sequence_count_total~,
    neff~,
    template_neff~,
    searched_hmms~,
    run_date~,
    command_line~,
  )
  HhrRecord::{ metadata, hit_summaries, alignments, }
}

///|
/// Serialize a parsed HHR document to a canonical, parseable HHR form.
pub fn HhrRecord::to_string(self : HhrRecord) -> String {
  let output = StringBuilder::new()
  output.write_string("Query         " + self.metadata.query_name + "\n")
  output.write_string(
    "Match_columns " + self.metadata.match_columns.to_string() + "\n",
  )
  output.write_string(
    "No_of_seqs    " +
    self.metadata.sequence_count.to_string() +
    " out of " +
    self.metadata.sequence_count_total.to_string() +
    "\n",
  )
  output.write_string("Neff          " + self.metadata.neff.to_string() + "\n")
  match self.metadata.template_neff {
    Some(value) =>
      output.write_string("Template_Neff " + value.to_string() + "\n")
    None => ()
  }
  output.write_string(
    "Searched_HMMs " + self.metadata.searched_hmms.to_string() + "\n",
  )
  if self.metadata.run_date.length() > 0 {
    output.write_string("Date          " + self.metadata.run_date + "\n")
  }
  if self.metadata.command_line.length() > 0 {
    output.write_string("Command       " + self.metadata.command_line + "\n")
  }
  output.write_string("\n")
  output.write_string(
    " No Hit                             Prob E-value P-value  Score    SS Cols Query HMM  Template HMM\n",
  )
  for summary in self.hit_summaries {
    let description = if summary.description.length() > 0 {
      " " + summary.description
    } else {
      ""
    }
    output.write_string(
      " " +
      summary.rank.to_string() +
      " " +
      summary.target_name +
      description +
      " " +
      summary.probability.to_string() +
      " " +
      summary.evalue.to_string() +
      " " +
      summary.pvalue.to_string() +
      " " +
      summary.score.to_string() +
      " " +
      summary.secondary_structure_score.to_string() +
      " " +
      summary.aligned_columns.to_string() +
      " " +
      (summary.query_start + 1).to_string() +
      "-" +
      summary.query_end.to_string() +
      " " +
      (summary.target_start + 1).to_string() +
      "-" +
      summary.target_end.to_string() +
      " (" +
      summary.target_length.to_string() +
      ")\n",
    )
  }
  output.write_string("\n")
  for alignment in self.alignments {
    output.write_string("No " + alignment.rank.to_string() + "\n")
    output.write_string(">" + alignment.target_name)
    if alignment.target_description.length() > 0 {
      output.write_string(" " + alignment.target_description)
    }
    output.write_string("\n")
    output.write_string(
      "Probab=" +
      alignment.probability.to_string() +
      " E-value=" +
      alignment.evalue.to_string() +
      " Score=" +
      alignment.score.to_string() +
      " Aligned_cols=" +
      alignment.aligned_columns.to_string() +
      " Identities=" +
      alignment.identities.to_string() +
      "% Similarity=" +
      alignment.similarity.to_string() +
      " Sum_probs=" +
      alignment.sum_probabilities.to_string() +
      "\n",
    )
    if alignment.query_secondary_structure.length() > 0 {
      output.write_string(
        "Q ss_pred " + alignment.query_secondary_structure + "\n",
      )
    }
    output.write_string(
      "Q " +
      alignment.query_name +
      " " +
      (alignment.query_start + 1).to_string() +
      " " +
      alignment.query_sequence +
      " " +
      alignment.query_end.to_string() +
      " (" +
      alignment.query_length.to_string() +
      ")\n",
    )
    if alignment.query_consensus.length() > 0 {
      output.write_string(
        "Q Consensus " +
        (alignment.query_start + 1).to_string() +
        " " +
        alignment.query_consensus +
        " " +
        alignment.query_end.to_string() +
        " (" +
        alignment.query_length.to_string() +
        ")\n",
      )
    }
    if alignment.column_score.length() > 0 {
      output.write_string(
        "                      " + alignment.column_score + "\n",
      )
    }
    if alignment.target_consensus.length() > 0 {
      output.write_string(
        "T Consensus " +
        (alignment.target_start + 1).to_string() +
        " " +
        alignment.target_consensus +
        " " +
        alignment.target_end.to_string() +
        " (" +
        alignment.target_length.to_string() +
        ")\n",
      )
    }
    output.write_string(
      "T " +
      alignment.target_name +
      " " +
      (alignment.target_start + 1).to_string() +
      " " +
      alignment.target_sequence +
      " " +
      alignment.target_end.to_string() +
      " (" +
      alignment.target_length.to_string() +
      ")\n",
    )
    if alignment.target_dssp.length() > 0 {
      output.write_string("T ss_dssp " + alignment.target_dssp + "\n")
    }
    if alignment.target_secondary_structure.length() > 0 {
      output.write_string(
        "T ss_pred " + alignment.target_secondary_structure + "\n",
      )
    }
    if alignment.confidence.length() > 0 {
      output.write_string("Confidence " + alignment.confidence + "\n")
    }
  }
  output.write_string("Done!\n")
  output.to_string()
}

///|
/// Return a compact two-hit HHR document for examples and tests.
pub fn hhr_sample_text() -> String {
  "Query         demo_query\n" +
  "Match_columns 12\n" +
  "No_of_seqs    24 out of 40\n" +
  "Neff          3.5\n" +
  "Template_Neff 2.1\n" +
  "Searched_HMMs 2\n" +
  "Date          Tue Aug 4 12:00:00 2026\n" +
  "Command       hhsearch -i demo.a3m -d pdb70\n" +
  " No Hit                             Prob E-value P-value  Score    SS Cols Query HMM  Template HMM\n" +
  " 1 target_A Alpha beta enzyme        99.8 1.2E-20 4.0E-25 85.4 2.1 10 1-12 5-14 (40)\n" +
  " 2 target_B Membrane protein          72.5 0.004 2.0E-6 31.2 0.0 8 3-10 2-9 (20)\n" +
  "No 1\n" +
  ">target_A Alpha beta enzyme\n" +
  "Probab=99.8 E-value=1.2e-20 Score=85.4 Aligned_cols=10 Identities=72.7% Similarity=1.45 Sum_probs=10.8\n" +
  "Q ss_pred          CCHHHH\n" +
  "Q demo_query     1 ACDEFG   6 (12)\n" +
  "Q Consensus      1 AcDEfG   6 (12)\n" +
  "                     ||+|.|\n" +
  "T Consensus      5 AC-EfG   9 (40)\n" +
  "T target_A       5 AC-EFG   9 (40)\n" +
  "T ss_dssp          CCEEEE\n" +
  "T ss_pred          CCHHHH\n" +
  "Confidence         998879\n" +
  "Q ss_pred          HHEECC\n" +
  "Q demo_query     7 HIKLMN  12 (12)\n" +
  "Q Consensus      7 HiKLMn  12 (12)\n" +
  "                     |.||||\n" +
  "T Consensus     10 H-KLMN  14 (40)\n" +
  "T target_A      10 H-KLMN  14 (40)\n" +
  "T ss_dssp          E-EECC\n" +
  "T ss_pred          H-EECC\n" +
  "Confidence         976899\n" +
  "No 2\n" +
  ">target_B Membrane protein\n" +
  "Probab=72.5 E-value=0.004 Score=31.2 Aligned_cols=8 Identities=50% Similarity=0.72 Sum_probs=6.2\n" +
  "Q demo_query     3 DEFGHIKL  10 (12)\n" +
  "Q Consensus      3 dEfGHiKl  10 (12)\n" +
  "                     +||..||+\n" +
  "T Consensus      2 TEFGYIKV   9 (20)\n" +
  "T target_B       2 TEFGYIKV   9 (20)\n" +
  "Confidence         67887765\n" +
  "Done!\n"
}

///|
priv struct HhrAlignmentBuilder {
  rank : Int
  mut target_name : String
  mut target_description : String
  mut detail_target_name : String
  mut probability : Double
  mut evalue : Double
  mut score : Double
  mut identities : Double
  mut similarity : Double
  mut sum_probabilities : Double
  mut aligned_columns : Int
  mut has_statistics : Bool
  mut query_length : Int
  mut target_length : Int
  mut query_start : Int
  mut query_end : Int
  mut target_start : Int
  mut target_end : Int
  mut query_sequence : String
  mut target_sequence : String
  mut query_consensus : String
  mut target_consensus : String
  mut query_secondary_structure : String
  mut target_secondary_structure : String
  mut target_dssp : String
  mut column_score : String
  mut confidence : String
  mut block_width : Int
}

///|
fn HhrAlignmentBuilder::new(rank : Int) -> HhrAlignmentBuilder {
  HhrAlignmentBuilder::{
    rank,
    target_name: "",
    target_description: "",
    detail_target_name: "",
    probability: 0.0,
    evalue: 0.0,
    score: 0.0,
    identities: 0.0,
    similarity: 0.0,
    sum_probabilities: 0.0,
    aligned_columns: 0,
    has_statistics: false,
    query_length: 0,
    target_length: 0,
    query_start: -1,
    query_end: -1,
    target_start: -1,
    target_end: -1,
    query_sequence: "",
    target_sequence: "",
    query_consensus: "",
    target_consensus: "",
    query_secondary_structure: "",
    target_secondary_structure: "",
    target_dssp: "",
    column_score: "",
    confidence: "",
    block_width: 0,
  }
}

///|
priv struct HhrSequencePart {
  name : String
  start : Int
  sequence : String
  end : Int
  total : Int
}

///|
fn hhr_finish_alignment(
  builder : HhrAlignmentBuilder,
  query_name : String,
  match_columns : Int,
  summaries : Array[HhrHitSummary],
) -> HhrAlignment raise HhrError {
  if builder.rank < 1 || builder.rank > summaries.length() {
    raise HhrError("Detailed alignment rank has no summary row")
  }
  if builder.target_name.length() == 0 {
    raise HhrError(
      "Alignment " + builder.rank.to_string() + " is missing target header",
    )
  }
  if !builder.has_statistics {
    raise HhrError(
      "Alignment " + builder.rank.to_string() + " is missing statistics",
    )
  }
  if builder.query_sequence.length() == 0 ||
    builder.target_sequence.length() == 0 {
    raise HhrError("Alignment " + builder.rank.to_string() + " has no sequence")
  }
  let width = builder.query_sequence.length()
  if builder.target_sequence.length() != width {
    raise HhrError(
      "Alignment " + builder.rank.to_string() + " has unequal row widths",
    )
  }
  if builder.query_length != match_columns {
    raise HhrError(
      "Alignment " +
      builder.rank.to_string() +
      " query length disagrees with Match_columns",
    )
  }
  if builder.query_start < 0 || builder.target_start < 0 {
    raise HhrError(
      "Alignment " + builder.rank.to_string() + " has missing coordinates",
    )
  }
  let query_residues = hhr_count_residues(builder.query_sequence)
  let target_residues = hhr_count_residues(builder.target_sequence)
  if query_residues != builder.query_end - builder.query_start {
    raise HhrError(
      "Alignment " +
      builder.rank.to_string() +
      " query coordinates disagree with sequence",
    )
  }
  if target_residues != builder.target_end - builder.target_start {
    raise HhrError(
      "Alignment " +
      builder.rank.to_string() +
      " target coordinates disagree with sequence",
    )
  }
  hhr_validate_annotation(
    builder.query_consensus,
    width,
    "query consensus",
    builder.rank,
  )
  hhr_validate_annotation(
    builder.target_consensus,
    width,
    "target consensus",
    builder.rank,
  )
  hhr_validate_annotation(
    builder.query_secondary_structure,
    width,
    "query secondary structure",
    builder.rank,
  )
  hhr_validate_annotation(
    builder.target_secondary_structure,
    width,
    "target secondary structure",
    builder.rank,
  )
  hhr_validate_annotation(
    builder.target_dssp,
    width,
    "target DSSP",
    builder.rank,
  )
  hhr_validate_annotation(
    builder.column_score,
    width,
    "column score",
    builder.rank,
  )
  hhr_validate_annotation(builder.confidence, width, "confidence", builder.rank)
  let summary = summaries[builder.rank - 1]
  if summary.rank != builder.rank || summary.target_name != builder.target_name {
    raise HhrError(
      "Alignment " + builder.rank.to_string() + " disagrees with summary target",
    )
  }
  if builder.detail_target_name.length() > 0 &&
    builder.detail_target_name != builder.target_name {
    raise HhrError(
      "Alignment " +
      builder.rank.to_string() +
      " has inconsistent target identifiers",
    )
  }
  if builder.aligned_columns != summary.aligned_columns {
    raise HhrError(
      "Alignment " +
      builder.rank.to_string() +
      " aligned column count disagrees with summary",
    )
  }
  if builder.aligned_columns !=
    hhr_count_paired_residues(builder.query_sequence, builder.target_sequence) {
    raise HhrError(
      "Alignment " +
      builder.rank.to_string() +
      " Aligned_cols disagrees with sequence rows",
    )
  }
  if builder.query_start != summary.query_start ||
    builder.query_end != summary.query_end ||
    builder.target_start != summary.target_start ||
    builder.target_end != summary.target_end ||
    builder.target_length != summary.target_length {
    raise HhrError(
      "Alignment " +
      builder.rank.to_string() +
      " coordinates disagree with summary",
    )
  }
  HhrAlignment::{
    rank: builder.rank,
    query_name,
    query_length: builder.query_length,
    target_name: builder.target_name,
    target_description: builder.target_description,
    target_length: builder.target_length,
    probability: builder.probability,
    evalue: builder.evalue,
    pvalue: summary.pvalue,
    score: builder.score,
    secondary_structure_score: summary.secondary_structure_score,
    identities: builder.identities,
    similarity: builder.similarity,
    sum_probabilities: builder.sum_probabilities,
    aligned_columns: builder.aligned_columns,
    query_start: builder.query_start,
    query_end: builder.query_end,
    target_start: builder.target_start,
    target_end: builder.target_end,
    query_sequence: builder.query_sequence,
    target_sequence: builder.target_sequence,
    query_consensus: builder.query_consensus,
    target_consensus: builder.target_consensus,
    query_secondary_structure: builder.query_secondary_structure,
    target_secondary_structure: builder.target_secondary_structure,
    target_dssp: builder.target_dssp,
    column_score: builder.column_score,
    confidence: builder.confidence,
  }
}

///|
fn hhr_parse_hit_summary(
  line : String,
  expected_rank : Int,
  line_number : Int,
) -> HhrHitSummary raise HhrError {
  let words = hhr_split_whitespace(line)
  if words.length() < 11 {
    raise HhrError(
      "Malformed HHR hit summary at line " + line_number.to_string(),
    )
  }
  let rank = hhr_parse_nonnegative_int(words[0], "hit rank", line_number)
  if rank != expected_rank {
    raise HhrError("Non-contiguous hit rank at line " + line_number.to_string())
  }
  let numeric_start = words.length() - 9
  if numeric_start < 2 {
    raise HhrError(
      "Malformed HHR hit summary at line " + line_number.to_string(),
    )
  }
  let description = hhr_join_words(words, 2, numeric_start)
  let probability = hhr_parse_number(
    words[numeric_start],
    "hit probability",
    line_number,
  )
  let evalue = hhr_parse_number(
    words[numeric_start + 1],
    "hit E-value",
    line_number,
  )
  let pvalue = hhr_parse_number(
    words[numeric_start + 2],
    "hit P-value",
    line_number,
  )
  let score = hhr_parse_number(
    words[numeric_start + 3],
    "hit score",
    line_number,
  )
  let secondary_structure_score = hhr_parse_number(
    words[numeric_start + 4],
    "hit SS score",
    line_number,
  )
  let aligned_columns = hhr_parse_nonnegative_int(
    words[numeric_start + 5],
    "hit columns",
    line_number,
  )
  let (query_start, query_end) = hhr_parse_range(
    words[numeric_start + 6],
    "query range",
    line_number,
  )
  let (target_start, target_end) = hhr_parse_range(
    words[numeric_start + 7],
    "target range",
    line_number,
  )
  let target_length = hhr_parse_parenthesized_int(
    words[numeric_start + 8],
    "target length",
    line_number,
  )
  HhrHitSummary::{
    rank,
    target_name: words[1],
    description,
    probability,
    evalue,
    pvalue,
    score,
    secondary_structure_score,
    aligned_columns,
    query_start,
    query_end,
    target_start,
    target_end,
    target_length,
  }
}

///|
fn hhr_parse_alignment_statistics(
  builder : HhrAlignmentBuilder,
  line : String,
  line_number : Int,
) -> Unit raise HhrError {
  if builder.has_statistics {
    raise HhrError(
      "Duplicate alignment statistics at line " + line_number.to_string(),
    )
  }
  let values : Map[String, String] = Map([])
  for word in hhr_split_whitespace(line) {
    let equal = hhr_find_char(word, '='.to_int())
    if equal <= 0 || equal == word.length() - 1 {
      raise HhrError(
        "Malformed alignment statistic at line " + line_number.to_string(),
      )
    }
    let key = word[0:equal].to_owned()
    let mut value = word[equal + 1:].to_owned()
    if key == "Identities" &&
      value.length() > 0 &&
      value.unsafe_get(value.length() - 1).to_int() == '%'.to_int() {
      value = value[0:value.length() - 1].to_owned()
    }
    values[key] = value
  }
  let required = [
    "Probab", "E-value", "Score", "Aligned_cols", "Identities", "Similarity", "Sum_probs",
  ]
  for key in required {
    if values.get(key) is None {
      raise HhrError("Missing " + key + " at line " + line_number.to_string())
    }
  }
  builder.probability = hhr_parse_number(
    values["Probab"],
    "Probab",
    line_number,
  )
  builder.evalue = hhr_parse_number(values["E-value"], "E-value", line_number)
  builder.score = hhr_parse_number(values["Score"], "Score", line_number)
  builder.aligned_columns = hhr_parse_nonnegative_int(
    values["Aligned_cols"],
    "Aligned_cols",
    line_number,
  )
  builder.identities = hhr_parse_number(
    values["Identities"],
    "Identities",
    line_number,
  )
  builder.similarity = hhr_parse_number(
    values["Similarity"],
    "Similarity",
    line_number,
  )
  builder.sum_probabilities = hhr_parse_number(
    values["Sum_probs"],
    "Sum_probs",
    line_number,
  )
  builder.has_statistics = true
}

///|
fn hhr_parse_sequence_line(
  line : String,
  kind : String,
  line_number : Int,
) -> HhrSequencePart raise HhrError {
  let words = hhr_split_whitespace(line)
  if words.length() != 6 {
    raise HhrError("Malformed " + kind + " line at " + line_number.to_string())
  }
  let name = words[1]
  let start = hhr_parse_nonnegative_int(words[2], kind + " start", line_number)
  if start <= 0 {
    raise HhrError(
      kind + " start must be positive at line " + line_number.to_string(),
    )
  }
  let sequence = words[3]
  let end = hhr_parse_nonnegative_int(words[4], kind + " end", line_number)
  let total = hhr_parse_parenthesized_int(
    words[5],
    kind + " total",
    line_number,
  )
  let converted_start = start - 1
  if end < converted_start || end > total {
    raise HhrError(
      "Invalid " + kind + " coordinates at line " + line_number.to_string(),
    )
  }
  if hhr_count_residues(sequence) != end - converted_start {
    raise HhrError(
      kind +
      " coordinates disagree with sequence at line " +
      line_number.to_string(),
    )
  }
  HhrSequencePart::{ name, start: converted_start, sequence, end, total, }
}

///|
fn hhr_validate_annotation(
  annotation : String,
  width : Int,
  name : String,
  rank : Int,
) -> Unit raise HhrError {
  if annotation.length() > 0 && annotation.length() != width {
    raise HhrError(
      "Alignment " +
      rank.to_string() +
      " " +
      name +
      " width " +
      annotation.length().to_string() +
      " does not match alignment width " +
      width.to_string(),
    )
  }
}

///|
fn hhr_normalize_lines(text : String) -> Array[String] {
  let result : Array[String] = []
  for view in text.split("\n") {
    let mut line = view.to_owned()
    if line.length() > 0 &&
      line.unsafe_get(line.length() - 1).to_int() == '\r'.to_int() {
      line = line[0:line.length() - 1].to_owned()
    }
    result.push(line)
  }
  result
}

///|
fn hhr_split_whitespace(text : String) -> Array[String] {
  let result : Array[String] = []
  let mut start = 0
  let mut in_word = false
  for i = 0; i < text.length(); i = i + 1 {
    let code = text.unsafe_get(i).to_int()
    let whitespace = code == ' '.to_int() ||
      code == '\t'.to_int() ||
      code == '\r'.to_int()
    if whitespace {
      if in_word {
        result.push(text[start:i].to_owned())
        in_word = false
      }
    } else if !in_word {
      start = i
      in_word = true
    }
  }
  if in_word {
    result.push(text[start:].to_owned())
  }
  result
}

///|
fn hhr_first_whitespace(text : String) -> (String, String)? {
  for i = 0; i < text.length(); i = i + 1 {
    let code = text.unsafe_get(i).to_int()
    if code == ' '.to_int() || code == '\t'.to_int() {
      let key = text[0:i].to_owned()
      let value = text[i:].trim().to_owned()
      if key.length() > 0 && value.length() > 0 {
        return Some((key, value))
      }
      return None
    }
  }
  None
}

///|
fn hhr_has_prefix(text : String, prefix : String) -> Bool {
  if prefix.length() > text.length() {
    return false
  }
  for i = 0; i < prefix.length(); i = i + 1 {
    if text.unsafe_get(i).to_int() != prefix.unsafe_get(i).to_int() {
      return false
    }
  }
  true
}

///|
fn hhr_is_indented(text : String) -> Bool {
  text.length() > 0 &&
  (
    text.unsafe_get(0).to_int() == ' '.to_int() ||
    text.unsafe_get(0).to_int() == '\t'.to_int()
  )
}

///|
fn hhr_is_rank_line(text : String) -> Bool {
  let words = hhr_split_whitespace(text)
  words.length() == 2 && words[0] == "No"
}

///|
fn hhr_find_char(text : String, needle : Int) -> Int {
  for i = 0; i < text.length(); i = i + 1 {
    if text.unsafe_get(i).to_int() == needle {
      return i
    }
  }
  -1
}

///|
fn hhr_find_substring(text : String, needle : String) -> Int {
  if needle.length() == 0 {
    return 0
  }
  if needle.length() > text.length() {
    return -1
  }
  for i = 0; i <= text.length() - needle.length(); i = i + 1 {
    let mut matches = true
    for j = 0; j < needle.length(); j = j + 1 {
      if text.unsafe_get(i + j).to_int() != needle.unsafe_get(j).to_int() {
        matches = false
        break
      }
    }
    if matches {
      return i
    }
  }
  -1
}

///|
fn hhr_parse_nonnegative_int(
  text : String,
  field : String,
  line_number : Int,
) -> Int raise HhrError {
  if text.length() == 0 {
    raise HhrError("Empty " + field + " at line " + line_number.to_string())
  }
  let mut value = 0
  for i = 0; i < text.length(); i = i + 1 {
    let code = text.unsafe_get(i).to_int()
    if code < '0'.to_int() || code > '9'.to_int() {
      raise HhrError("Invalid " + field + " at line " + line_number.to_string())
    }
    value = value * 10 + code - '0'.to_int()
  }
  value
}

///|
fn hhr_parse_number(
  text : String,
  field : String,
  line_number : Int,
) -> Double raise HhrError {
  if text.length() == 0 {
    raise HhrError("Empty " + field + " at line " + line_number.to_string())
  }
  let mut has_digit = false
  for i = 0; i < text.length(); i = i + 1 {
    let code = text.unsafe_get(i).to_int()
    if code >= '0'.to_int() && code <= '9'.to_int() {
      has_digit = true
    } else if code != '-'.to_int() &&
      code != '+'.to_int() &&
      code != '.'.to_int() &&
      code != 'e'.to_int() &&
      code != 'E'.to_int() {
      raise HhrError("Invalid " + field + " at line " + line_number.to_string())
    }
  }
  if !has_digit {
    raise HhrError("Invalid " + field + " at line " + line_number.to_string())
  }
  match parse_double(text) {
    Some(value) => value
    None =>
      raise HhrError("Invalid " + field + " at line " + line_number.to_string())
  }
}

///|
fn hhr_parse_range(
  text : String,
  field : String,
  line_number : Int,
) -> (Int, Int) raise HhrError {
  let separator = hhr_find_char(text, '-'.to_int())
  if separator <= 0 || separator >= text.length() - 1 {
    raise HhrError("Invalid " + field + " at line " + line_number.to_string())
  }
  let start = hhr_parse_nonnegative_int(
    text[0:separator].to_owned(),
    field,
    line_number,
  )
  let end = hhr_parse_nonnegative_int(
    text[separator + 1:].to_owned(),
    field,
    line_number,
  )
  if start <= 0 || end < start {
    raise HhrError("Invalid " + field + " at line " + line_number.to_string())
  }
  (start - 1, end)
}

///|
fn hhr_parse_parenthesized_int(
  text : String,
  field : String,
  line_number : Int,
) -> Int raise HhrError {
  if text.length() < 3 ||
    text.unsafe_get(0).to_int() != '('.to_int() ||
    text.unsafe_get(text.length() - 1).to_int() != ')'.to_int() {
    raise HhrError("Invalid " + field + " at line " + line_number.to_string())
  }
  hhr_parse_nonnegative_int(
    text[1:text.length() - 1].to_owned(),
    field,
    line_number,
  )
}

///|
fn hhr_last_token(
  line : String,
  field : String,
  line_number : Int,
) -> String raise HhrError {
  let words = hhr_split_whitespace(line)
  if words.length() < 3 {
    raise HhrError("Malformed " + field + " at line " + line_number.to_string())
  }
  words[words.length() - 1]
}

///|
fn hhr_suffix(
  line : String,
  width : Int,
  field : String,
  line_number : Int,
) -> String raise HhrError {
  if width < 0 || line.length() < width {
    raise HhrError("Short " + field + " at line " + line_number.to_string())
  }
  line[line.length() - width:].to_owned()
}

///|
fn hhr_join_words(words : Array[String], start : Int, end : Int) -> String {
  let output = StringBuilder::new()
  for i = start; i < end; i = i + 1 {
    if i > start {
      output.write_char(' ')
    }
    output.write_string(words[i])
  }
  output.to_string()
}

///|
fn hhr_remove_gaps(sequence : String) -> String {
  let output = StringBuilder::new()
  for i = 0; i < sequence.length(); i = i + 1 {
    let code = sequence.unsafe_get(i).to_int()
    if code != '-'.to_int() {
      output.write_char(code.unsafe_to_char())
    }
  }
  output.to_string()
}

///|
fn hhr_count_residues(sequence : String) -> Int {
  let mut count = 0
  for i = 0; i < sequence.length(); i = i + 1 {
    if sequence.unsafe_get(i).to_int() != '-'.to_int() {
      count = count + 1
    }
  }
  count
}

///|
fn hhr_count_paired_residues(query : String, target : String) -> Int {
  let mut count = 0
  let width = if query.length() < target.length() {
    query.length()
  } else {
    target.length()
  }
  for i = 0; i < width; i = i + 1 {
    if query.unsafe_get(i).to_int() != '-'.to_int() &&
      target.unsafe_get(i).to_int() != '-'.to_int() {
      count = count + 1
    }
  }
  count
}
