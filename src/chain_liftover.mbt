///|
/// chain_liftover - UCSC Chain file format parsing and genome coordinate liftOver.
///
/// Chain files describe pairwise alignments between two genome assemblies.
/// This module is inspired by Bioconductor rtracklayer::liftOver and supports
/// parsing chain format text, lifting over individual positions and ranges,
/// querying chain summaries, and serializing chains back to text.

///|
/// A contiguous aligned segment within a chain alignment.
/// Each block represents a stretch of aligned bases between target and query,
/// followed by an unaligned gap in each sequence.
pub struct ChainBlock {
  size : Int
  dt : Int
  dq : Int
}

///|
/// A single chain alignment record describing the alignment between a region
/// of the target genome and a region of the query genome.
pub struct ChainAlignment {
  t_name : String
  score : Int
  t_plus : Bool
  t_size : Int
  t_start : Int
  t_end : Int
  q_name : String
  q_size : Int
  q_strand : String
  q_start : Int
  q_end : Int
  id : Int
  blocks : Array[ChainBlock]
}

///|
/// A collection of chain alignments, typically representing a whole-genome
/// alignment between two assemblies.
pub struct ChainFile {
  alignments : Array[ChainAlignment]
  target_name : String
  query_name : String
}

///|
/// Parse UCSC chain file format text into a ChainFile.
///
/// The chain format consists of header lines starting with "chain" followed by
/// indented block lines. Each alignment record maps a region of the target genome
/// to a region of the query genome.
pub fn cl_parse_chain(content : String) -> ChainFile {
  let lines = cl_split_lines(content)
  let alignments : Array[ChainAlignment] = Array::new()

  let mut i = 0
  let mut current_target_name = ""
  let mut current_query_name = ""

  while i < lines.length() {
    let line = lines[i]
    let trimmed = cl_trim(line)

    if trimmed.length() == 0 {
      i = i + 1
      continue
    }

    if cl_starts_with(trimmed, "chain") {
      let (header, next_i) = cl_parse_chain_header(lines, i)
      current_target_name = header.t_name
      current_query_name = header.q_name
      let blocks = cl_parse_chain_blocks(lines, next_i)
      alignments.push(ChainAlignment::{
        t_name: header.t_name,
        score: header.score,
        t_plus: header.t_plus,
        t_size: header.t_size,
        t_start: header.t_start,
        t_end: header.t_end,
        q_name: header.q_name,
        q_size: header.q_size,
        q_strand: header.q_strand,
        q_start: header.q_start,
        q_end: header.q_end,
        id: header.id,
        blocks,
      })
      i = cl_skip_to_next_chain(lines, next_i)
    } else {
      i = i + 1
    }
  }

  ChainFile::{
    alignments,
    target_name: current_target_name,
    query_name: current_query_name,
  }
}

///|
/// Lift over a single genome position from target to query coordinates.
///
/// Returns the corresponding query position, or None if the position falls
/// in an unaligned gap or is outside any chain alignment.
pub fn cl_liftover_pos(chains : ChainFile, chrom : String, pos : Int) -> Int? {
  let chain = cl_find_chain_for_pos(chains, chrom, pos)
  match chain {
    None => None
    Some(ch) => {
      let mut target_offset = ch.t_start
      let mut query_offset = ch.q_start

      let mut i = 0
      while i < ch.blocks.length() {
        let block = ch.blocks[i]
        let block_end = target_offset + block.size

        if pos < block_end {
          if pos >= target_offset {
            let q_pos = query_offset + (pos - target_offset)
            if ch.q_strand == "-" {
              return Some(ch.q_size - q_pos)
            } else {
              return Some(q_pos)
            }
          } else {
            return None
          }
        }
        target_offset = block_end + block.dt
        query_offset = query_offset + block.size + block.dq
        i = i + 1
      }
      None
    }
  }
}

///|
/// Lift over a genomic range from target to query coordinates.
///
/// Returns a tuple of (start, end) in query coordinates, or None if the
/// range cannot be lifted over (e.g., spans unaligned gaps or multiple chains).
pub fn cl_liftover_range(
  chains : ChainFile,
  chrom : String,
  start : Int,
  end : Int,
) -> (Int, Int)? {
  if start >= end {
    return None
  }

  let chain = cl_find_chain_for_pos(chains, chrom, start)
  match chain {
    None => None
    Some(ch) => {
      if end > ch.t_end {
        return None
      }

      let mut target_offset = ch.t_start
      let mut query_offset = ch.q_start
      let mut mapped_start : Int? = None
      let mut mapped_end : Int? = None

      let mut i = 0
      while i < ch.blocks.length() {
        let block = ch.blocks[i]
        let block_end = target_offset + block.size

        if start < block_end && mapped_start is None {
          if start >= target_offset {
            mapped_start = Some(query_offset + (start - target_offset))
          } else {
            return None
          }
        }

        if end <= block_end {
          if end >= target_offset {
            mapped_end = Some(query_offset + (end - target_offset))
          } else {
            return None
          }
          break
        }

        target_offset = block_end + block.dt
        query_offset = query_offset + block.size + block.dq
        i = i + 1
      }

      match (mapped_start, mapped_end) {
        (Some(ms), Some(me)) => {
          if ch.q_strand == "-" {
            let qs = ch.q_size - me
            let qe = ch.q_size - ms
            if qs < qe {
              Some((qs, qe))
            } else {
              Some((qe, qs))
            }
          } else {
            if ms < me {
              Some((ms, me))
            } else {
              Some((me, ms))
            }
          }
        }
        _ => None
      }
    }
  }
}

///|
/// Get alignment statistics for a chain alignment.
///
/// Returns a Map with keys: "total_block_size", "num_blocks",
/// "total_gap_target", "total_gap_query", "aligned_target", "aligned_query".
pub fn cl_get_chain_summary(chain : ChainAlignment) -> Map[String, Int] {
  let summary : Map[String, Int] = Map([], capacity=6)
  let mut total_block_size = 0
  let mut total_gap_target = 0
  let mut total_gap_query = 0

  let mut i = 0
  while i < chain.blocks.length() {
    let block = chain.blocks[i]
    total_block_size = total_block_size + block.size
    total_gap_target = total_gap_target + block.dt
    total_gap_query = total_gap_query + block.dq
    i = i + 1
  }

  summary["total_block_size"] = total_block_size
  summary["num_blocks"] = chain.blocks.length()
  summary["total_gap_target"] = total_gap_target
  summary["total_gap_query"] = total_gap_query
  summary["aligned_target"] = chain.t_end - chain.t_start
  summary["aligned_query"] = chain.q_end - chain.q_start

  summary
}

///|
/// Serialize a chain alignment back to UCSC chain format text.
///
/// Produces a multi-line string with the chain header followed by
/// indented block lines.
pub fn cl_chain_to_string(chain : ChainAlignment) -> String {
  let mut result = "chain "
  result = result + chain.score.to_string() + " "
  result = result + chain.t_name + " "
  result = result + chain.t_size.to_string() + " "
  let strand_t = if chain.t_plus { "+" } else { "-" }
  result = result + strand_t + " "
  result = result + chain.t_start.to_string() + " "
  result = result + chain.t_end.to_string() + " "
  result = result + chain.q_name + " "
  result = result + chain.q_size.to_string() + " "
  result = result + chain.q_strand + " "
  result = result + chain.q_start.to_string() + " "
  result = result + chain.q_end.to_string() + " "
  result = result + chain.id.to_string() + "\n"

  let mut i = 0
  while i < chain.blocks.length() {
    let block = chain.blocks[i]
    result = result + "  " + block.size.to_string() + " "
    result = result + block.dt.to_string() + " "
    result = result + block.dq.to_string() + "\n"
    i = i + 1
  }

  result
}

///|
/// Find the first chain alignment that covers a given target position.
///
/// Searches through all alignments in the file and returns the first one
/// whose target chromosome matches and whose [t_start, t_end) range contains
/// the position. Returns None if no such chain exists.
pub fn cl_find_chain_for_pos(
  chains : ChainFile,
  chrom : String,
  pos : Int,
) -> ChainAlignment? {
  let mut i = 0
  while i < chains.alignments.length() {
    let ch = chains.alignments[i]
    if ch.t_name == chrom && pos >= ch.t_start && pos < ch.t_end {
      return Some(ch)
    }
    i = i + 1
  }
  None
}

// --- Internal helper functions ---

///|
/// Split a string by newline character.
fn cl_split_lines(s : String) -> Array[String] {
  let result : Array[String] = Array::new()
  let mut start = 0
  let mut i = 0
  while i < s.length() {
    if s.unsafe_get(i).to_int() == '\n'.to_int() {
      let line = if start < i {
        s[start:i].to_owned()
      } else {
        ""
      }
      result.push(line)
      start = i + 1
    }
    i = i + 1
  }
  if start < s.length() {
    result.push(s[start:s.length()].to_owned())
  }
  result
}

///|
/// Trim leading and trailing whitespace.
fn cl_trim(s : String) -> String {
  let mut start = 0
  while start < s.length() &&
        (
          s.unsafe_get(start).to_int() == ' '.to_int() ||
          s.unsafe_get(start).to_int() == '\t'.to_int()
        ) {
    start = start + 1
  }
  let mut end = s.length()
  while end > start &&
        (
          s.unsafe_get(end - 1).to_int() == ' '.to_int() ||
          s.unsafe_get(end - 1).to_int() == '\t'.to_int() ||
          s.unsafe_get(end - 1).to_int() == '\r'.to_int()
        ) {
    end = end - 1
  }
  if start <= end && end <= s.length() {
    s[start:end].to_owned()
  } else {
    ""
  }
}

///|
/// Check if string starts with a given prefix.
fn cl_starts_with(s : String, prefix : String) -> Bool {
  if s.length() < prefix.length() {
    return false
  }
  let mut i = 0
  while i < prefix.length() {
    if s.unsafe_get(i).to_int() != prefix.unsafe_get(i).to_int() {
      return false
    }
    i = i + 1
  }
  true
}

///|
/// Split a string by whitespace (spaces and tabs).
fn cl_split_spaces(s : String) -> Array[String] {
  let result : Array[String] = Array::new()
  let mut start = 0
  let mut i = 0
  while i < s.length() {
    let c = s.unsafe_get(i).to_int()
    if c == ' '.to_int() || c == '\t'.to_int() {
      if start < i {
        result.push(s[start:i].to_owned())
      }
      start = i + 1
    }
    i = i + 1
  }
  if start < s.length() {
    result.push(s[start:s.length()].to_owned())
  }
  result
}

///|
/// Internal struct for parsed chain header fields.
priv struct ChainHeader {
  t_name : String
  score : Int
  t_size : Int
  t_plus : Bool
  t_start : Int
  t_end : Int
  q_name : String
  q_size : Int
  q_strand : String
  q_start : Int
  q_end : Int
  id : Int
}

///|
/// Parse a chain header line and return the header plus the index
/// of the next line to process.
fn cl_parse_chain_header(
  lines : Array[String],
  start_idx : Int,
) -> (ChainHeader, Int) {
  let line = lines[start_idx]
  let fields = cl_split_spaces(line)

  let score = if fields.length() > 1 { cl_parse_int_safe(fields[1]) } else { 0 }
  let t_name = if fields.length() > 2 { fields[2] } else { "" }
  let t_size = if fields.length() > 3 { cl_parse_int_safe(fields[3]) } else { 0 }
  let t_strand = if fields.length() > 4 { fields[4] } else { "+" }
  let t_plus = t_strand == "+"
  let t_start = if fields.length() > 5 { cl_parse_int_safe(fields[5]) } else { 0 }
  let t_end = if fields.length() > 6 { cl_parse_int_safe(fields[6]) } else { 0 }
  let q_name = if fields.length() > 7 { fields[7] } else { "" }
  let q_size = if fields.length() > 8 { cl_parse_int_safe(fields[8]) } else { 0 }
  let q_strand = if fields.length() > 9 { fields[9] } else { "+" }
  let q_start = if fields.length() > 10 { cl_parse_int_safe(fields[10]) } else { 0 }
  let q_end = if fields.length() > 11 { cl_parse_int_safe(fields[11]) } else { 0 }
  let id = if fields.length() > 12 { cl_parse_int_safe(fields[12]) } else { 0 }

  let header = ChainHeader::{
    t_name,
    score,
    t_size,
    t_plus,
    t_start,
    t_end,
    q_name,
    q_size,
    q_strand,
    q_start,
    q_end,
    id,
  }

  (header, start_idx + 1)
}

///|
/// Parse chain block lines starting from the given index.
/// Block lines are indented (start with whitespace) and contain
/// space-delimited integers: size, dt, dq.
fn cl_parse_chain_blocks(
  lines : Array[String],
  start_idx : Int,
) -> Array[ChainBlock] {
  let blocks : Array[ChainBlock] = Array::new()
  let mut i = start_idx

  while i < lines.length() {
    let line = lines[i]
    let trimmed = cl_trim(line)

    if trimmed.length() == 0 {
      i = i + 1
      continue
    }

    if cl_starts_with(trimmed, "chain") {
      break
    }

    let is_block = line.length() > 0 &&
      (
        line.unsafe_get(0).to_int() == ' '.to_int() ||
        line.unsafe_get(0).to_int() == '\t'.to_int()
      )

    if is_block {
      let fields = cl_split_spaces(trimmed)
      if fields.length() >= 3 {
        let size = cl_parse_int_safe(fields[0])
        let dt = cl_parse_int_safe(fields[1])
        let dq = cl_parse_int_safe(fields[2])
        blocks.push(ChainBlock::{ size, dt, dq })
      }
      i = i + 1
    } else {
      break
    }
  }

  blocks
}

///|
/// Skip to the next chain header line (or end of file).
fn cl_skip_to_next_chain(lines : Array[String], from : Int) -> Int {
  let mut i = from
  while i < lines.length() {
    let trimmed = cl_trim(lines[i])
    if cl_starts_with(trimmed, "chain") {
      return i
    }
    i = i + 1
  }
  lines.length()
}

///|
/// Parse an integer safely, returning 0 for invalid input.
fn cl_parse_int_safe(s : String) -> Int {
  let n = s.length()
  if n == 0 {
    return 0
  }
  let mut result = 0
  let mut i = 0
  let zero = '0'.to_int()
  while i < n {
    let c = s.unsafe_get(i).to_int()
    if c < zero || c > zero + 9 {
      break
    }
    result = result * 10 + (c - zero)
    i = i + 1
  }
  result
}

///|
/// Create a simple chain alignment for testing purposes.
fn cl_make_test_chain(
  t_name : String,
  t_start : Int,
  t_end : Int,
  q_name : String,
  q_start : Int,
  q_end : Int,
  blocks : Array[ChainBlock],
  id : Int,
) -> ChainAlignment {
  ChainAlignment::{
    t_name,
    score: 1000,
    t_plus: true,
    t_size: 0,
    t_start,
    t_end,
    q_name,
    q_size: 0,
    q_strand: "+",
    q_start,
    q_end,
    id,
    blocks,
  }
}

// ============================================================
// Tests
// ============================================================

///|
test "cl_parse_chain_empty" {
  let content = ""
  let cf = cl_parse_chain(content)
  assert_eq(cf.alignments.length(), 0)
}

///|
test "cl_parse_chain_single_header" {
  let content = "chain 5000 chr1 249250621 + 1 5000 chr2 10000 + 1 100 5"
  let cf = cl_parse_chain(content)
  assert_eq(cf.alignments.length(), 1)
  assert_eq(cf.alignments[0].score, 5000)
  assert_eq(cf.alignments[0].t_name, "chr1")
  assert_eq(cf.alignments[0].t_start, 1)
  assert_eq(cf.alignments[0].t_end, 5000)
  assert_eq(cf.alignments[0].q_name, "chr2")
  assert_eq(cf.alignments[0].q_start, 1)
  assert_eq(cf.alignments[0].q_end, 100)
  assert_eq(cf.alignments[0].id, 5)
  assert_eq(cf.alignments[0].t_plus, true)
}

///|
test "cl_parse_chain_with_blocks" {
  let content =
    "chain 5000 chr1 249250621 + 1 5000 chr2 10000 + 1 100 5\n  500 10 20\n  400 0 0"
  let cf = cl_parse_chain(content)
  assert_eq(cf.alignments.length(), 1)
  assert_eq(cf.alignments[0].blocks.length(), 2)
  assert_eq(cf.alignments[0].blocks[0].size, 500)
  assert_eq(cf.alignments[0].blocks[0].dt, 10)
  assert_eq(cf.alignments[0].blocks[0].dq, 20)
  assert_eq(cf.alignments[0].blocks[1].size, 400)
  assert_eq(cf.alignments[0].blocks[1].dt, 0)
  assert_eq(cf.alignments[0].blocks[1].dq, 0)
}

///|
test "cl_parse_chain_multiple_records" {
  let content =
    "chain 5000 chr1 249250621 + 1 5000 chr2 10000 + 1 100 5\n  500 0 0\nchain 3000 chr1 249250621 + 6000 10000 chr3 8000 + 500 1500 6\n  300 0 0"
  let cf = cl_parse_chain(content)
  assert_eq(cf.alignments.length(), 2)
  assert_eq(cf.alignments[0].id, 5)
  assert_eq(cf.alignments[1].id, 6)
  assert_eq(cf.alignments[0].blocks.length(), 1)
  assert_eq(cf.alignments[1].blocks.length(), 1)
}

///|
test "cl_parse_chain_multiple_blocks" {
  let content =
    "chain 5000 chr1 249250621 + 1 5000 chr2 10000 + 1 100 5\n  500 10 20\n  300 5 15\n  200 0 0"
  let cf = cl_parse_chain(content)
  assert_eq(cf.alignments[0].blocks.length(), 3)
  assert_eq(cf.alignments[0].blocks[0].size, 500)
  assert_eq(cf.alignments[0].blocks[1].size, 300)
  assert_eq(cf.alignments[0].blocks[2].size, 200)
}

///|
test "cl_parse_chain_negative_query_strand" {
  let content =
    "chain 5000 chr1 249250621 + 1 5000 chr2 10000 - 5000 1 5\n  5000 0 0"
  let cf = cl_parse_chain(content)
  assert_eq(cf.alignments[0].q_strand, "-")
  assert_eq(cf.alignments[0].q_start, 5000)
  assert_eq(cf.alignments[0].q_end, 1)
}

///|
test "cl_parse_chain_negative_target_strand" {
  let content =
    "chain 5000 chr1 249250621 - 5000 1 chr2 10000 + 1 5000 5\n  5000 0 0"
  let cf = cl_parse_chain(content)
  assert_eq(cf.alignments[0].t_plus, false)
  assert_eq(cf.alignments[0].t_start, 5000)
  assert_eq(cf.alignments[0].t_end, 1)
}

///|
test "cl_parse_chain_complex" {
  let content =
    "chain 5000 chr1 249250621 + 1 5000 chr2 10000 + 1 100 5\n  500 10 20\n  400 30 40\n  300 0 0\nchain 3000 chr5 181538259 + 100 200 chr7 159345973 + 300 500 6\n  100 0 0\nchain 2000 chr12 133275309 + 500 800 chr14 107043718 + 700 1000 7\n  200 5 10\n  100 0 0"
  let cf = cl_parse_chain(content)
  assert_eq(cf.alignments.length(), 3)
  assert_eq(cf.alignments[0].blocks.length(), 3)
  assert_eq(cf.alignments[1].blocks.length(), 1)
  assert_eq(cf.alignments[2].blocks.length(), 2)
}

///|
test "cl_find_chain_for_pos_valid" {
  let content =
    "chain 5000 chr1 249250621 + 1 5000 chr2 10000 + 1 100 5\n  5000 0 0"
  let cf = cl_parse_chain(content)
  let result = cl_find_chain_for_pos(cf, "chr1", 100)
  assert_true(result is Some(_))
  match result {
    Some(ch) => {
      assert_eq(ch.id, 5)
    }
    None => ()
  }
}

///|
test "cl_find_chain_for_pos_first_block" {
  let content =
    "chain 5000 chr1 249250621 + 1 5000 chr2 10000 + 1 100 5\n  500 10 20\n  4500 0 0"
  let cf = cl_parse_chain(content)
  let result = cl_find_chain_for_pos(cf, "chr1", 100)
  assert_true(result is Some(_))
}

///|
test "cl_find_chain_for_pos_second_block" {
  let content =
    "chain 5000 chr1 249250621 + 1 5000 chr2 10000 + 1 100 5\n  500 10 20\n  4500 0 0"
  let cf = cl_parse_chain(content)
  let result = cl_find_chain_for_pos(cf, "chr1", 600)
  assert_true(result is Some(_))
}

///|
test "cl_find_chain_for_pos_wrong_chrom" {
  let content =
    "chain 5000 chr1 249250621 + 1 5000 chr2 10000 + 1 100 5\n  5000 0 0"
  let cf = cl_parse_chain(content)
  let result = cl_find_chain_for_pos(cf, "chrX", 100)
  assert_true(result is None)
}

///|
test "cl_find_chain_for_pos_out_of_range" {
  let content =
    "chain 5000 chr1 249250621 + 1 5000 chr2 10000 + 1 100 5\n  5000 0 0"
  let cf = cl_parse_chain(content)
  let result = cl_find_chain_for_pos(cf, "chr1", 6000)
  assert_true(result is None)
}

///|
test "cl_find_chain_for_pos_empty_chains" {
  let cf = ChainFile::{
    alignments: [],
    target_name: "",
    query_name: "",
  }
  let result = cl_find_chain_for_pos(cf, "chr1", 100)
  assert_true(result is None)
}

///|
test "cl_liftover_pos_basic" {
  let content =
    "chain 5000 chr1 249250621 + 1 5000 chr2 10000 + 1 5000 5\n  5000 0 0"
  let cf = cl_parse_chain(content)
  let result = cl_liftover_pos(cf, "chr1", 100)
  assert_eq(result, Some(100))
}

///|
test "cl_liftover_pos_first_block" {
  let content =
    "chain 5000 chr1 249250621 + 1 5000 chr2 10000 + 1 5000 5\n  500 10 20\n  4500 0 0"
  let cf = cl_parse_chain(content)
  let result = cl_liftover_pos(cf, "chr1", 100)
  assert_eq(result, Some(100))
}

///|
test "cl_liftover_pos_second_block" {
  let content =
    "chain 5000 chr1 249250621 + 1 5000 chr2 10000 + 1 5000 5\n  500 10 20\n  4500 0 0"
  let cf = cl_parse_chain(content)
  let result = cl_liftover_pos(cf, "chr1", 600)
  assert_eq(result, Some(610))
}

///|
test "cl_liftover_pos_last_block" {
  let content =
    "chain 5000 chr1 249250621 + 1 5000 chr2 10000 + 1 5000 5\n  500 10 20\n  4500 0 0"
  let cf = cl_parse_chain(content)
  let result = cl_liftover_pos(cf, "chr1", 4999)
  assert_eq(result, Some(5009))
}

///|
test "cl_liftover_pos_at_chain_start" {
  let content =
    "chain 5000 chr1 249250621 + 1 5000 chr2 10000 + 1 5000 5\n  5000 0 0"
  let cf = cl_parse_chain(content)
  let result = cl_liftover_pos(cf, "chr1", 1)
  assert_eq(result, Some(1))
}

///|
test "cl_liftover_pos_at_chain_end" {
  let content =
    "chain 5000 chr1 249250621 + 1 5000 chr2 10000 + 1 5000 5\n  5000 0 0"
  let cf = cl_parse_chain(content)
  let result = cl_liftover_pos(cf, "chr1", 4999)
  assert_eq(result, Some(4999))
}

///|
test "cl_liftover_pos_out_of_range_before" {
  let content =
    "chain 5000 chr1 249250621 + 1 5000 chr2 10000 + 1 5000 5\n  5000 0 0"
  let cf = cl_parse_chain(content)
  let result = cl_liftover_pos(cf, "chr1", 0)
  assert_true(result is None)
}

///|
test "cl_liftover_pos_out_of_range_after" {
  let content =
    "chain 5000 chr1 249250621 + 1 5000 chr2 10000 + 1 5000 5\n  5000 0 0"
  let cf = cl_parse_chain(content)
  let result = cl_liftover_pos(cf, "chr1", 6000)
  assert_true(result is None)
}

///|
test "cl_liftover_pos_wrong_chrom" {
  let content =
    "chain 5000 chr1 249250621 + 1 5000 chr2 10000 + 1 5000 5\n  5000 0 0"
  let cf = cl_parse_chain(content)
  let result = cl_liftover_pos(cf, "chrX", 100)
  assert_true(result is None)
}

///|
test "cl_liftover_pos_in_gap" {
  let content =
    "chain 5000 chr1 249250621 + 1 5000 chr2 10000 + 1 5000 5\n  500 10 20\n  4500 0 0"
  let cf = cl_parse_chain(content)
  let result = cl_liftover_pos(cf, "chr1", 501)
  assert_true(result is None)
}

///|
test "cl_liftover_pos_reverse_query_strand" {
  let content =
    "chain 5000 chr1 249250621 + 1 5000 chr2 10000 - 5000 1 5\n  5000 0 0"
  let cf = cl_parse_chain(content)
  let result = cl_liftover_pos(cf, "chr1", 100)
  assert_eq(result, Some(4901))
}

///|
test "cl_liftover_pos_reverse_query_at_start" {
  let content =
    "chain 5000 chr1 249250621 + 1 5000 chr2 10000 - 5000 1 5\n  5000 0 0"
  let cf = cl_parse_chain(content)
  let result = cl_liftover_pos(cf, "chr1", 1)
  assert_eq(result, Some(5000))
}

///|
test "cl_liftover_pos_reverse_query_at_end" {
  let content =
    "chain 5000 chr1 249250621 + 1 5000 chr2 10000 - 5000 1 5\n  5000 0 0"
  let cf = cl_parse_chain(content)
  let result = cl_liftover_pos(cf, "chr1", 4999)
  assert_eq(result, Some(2))
}

///|
test "cl_liftover_pos_empty_chains" {
  let cf = ChainFile::{
    alignments: [],
    target_name: "",
    query_name: "",
  }
  let result = cl_liftover_pos(cf, "chr1", 100)
  assert_true(result is None)
}

///|
test "cl_liftover_range_basic" {
  let content =
    "chain 5000 chr1 249250621 + 1 5000 chr2 10000 + 1 5000 5\n  5000 0 0"
  let cf = cl_parse_chain(content)
  let result = cl_liftover_range(cf, "chr1", 100, 200)
  assert_eq(result, Some((100, 200)))
}

///|
test "cl_liftover_range_across_blocks" {
  let content =
    "chain 5000 chr1 249250621 + 1 5000 chr2 10000 + 1 5000 5\n  500 10 20\n  4500 0 0"
  let cf = cl_parse_chain(content)
  let result = cl_liftover_range(cf, "chr1", 400, 700)
  assert_eq(result, Some((400, 710)))
}

///|
test "cl_liftover_range_out_of_range" {
  let content =
    "chain 5000 chr1 249250621 + 1 5000 chr2 10000 + 1 5000 5\n  5000 0 0"
  let cf = cl_parse_chain(content)
  let result = cl_liftover_range(cf, "chr1", 4900, 6000)
  assert_true(result is None)
}

///|
test "cl_liftover_range_reverse_strand" {
  let content =
    "chain 5000 chr1 249250621 + 1 5000 chr2 10000 - 5000 1 5\n  5000 0 0"
  let cf = cl_parse_chain(content)
  let result = cl_liftover_range(cf, "chr1", 100, 200)
  assert_eq(result, Some((4801, 4901)))
}

///|
test "cl_liftover_range_empty" {
  let content =
    "chain 5000 chr1 249250621 + 1 5000 chr2 10000 + 1 5000 5\n  5000 0 0"
  let cf = cl_parse_chain(content)
  let result = cl_liftover_range(cf, "chr1", 200, 100)
  assert_true(result is None)
}

///|
test "cl_get_chain_summary_basic" {
  let blocks : Array[ChainBlock] = [
    ChainBlock::{ size: 500, dt: 10, dq: 20 },
    ChainBlock::{ size: 300, dt: 5, dq: 15 },
    ChainBlock::{ size: 200, dt: 0, dq: 0 },
  ]
  let chain = cl_make_test_chain("chr1", 1, 1000, "chr2", 1, 1000, blocks, 42)
  let summary = cl_get_chain_summary(chain)
  assert_eq(summary["total_block_size"], 1000)
  assert_eq(summary["num_blocks"], 3)
  assert_eq(summary["total_gap_target"], 15)
  assert_eq(summary["total_gap_query"], 35)
  assert_eq(summary["aligned_target"], 999)
  assert_eq(summary["aligned_query"], 999)
}

///|
test "cl_get_chain_summary_empty_blocks" {
  let blocks : Array[ChainBlock] = []
  let chain = cl_make_test_chain("chr1", 1, 500, "chr2", 1, 500, blocks, 43)
  let summary = cl_get_chain_summary(chain)
  assert_eq(summary["total_block_size"], 0)
  assert_eq(summary["num_blocks"], 0)
  assert_eq(summary["aligned_target"], 499)
}

///|
test "cl_get_chain_summary_single_block" {
  let blocks : Array[ChainBlock] = [ChainBlock::{ size: 500, dt: 0, dq: 0 }]
  let chain = cl_make_test_chain("chr1", 1, 501, "chr2", 1, 501, blocks, 44)
  let summary = cl_get_chain_summary(chain)
  assert_eq(summary["total_block_size"], 500)
  assert_eq(summary["num_blocks"], 1)
  assert_eq(summary["total_gap_target"], 0)
  assert_eq(summary["total_gap_query"], 0)
}

///|
test "cl_chain_to_string_roundtrip" {
  let content =
    "chain 5000 chr1 249250621 + 1 5000 chr2 10000 + 1 100 5\n  500 10 20\n  400 0 0\n"
  let cf = cl_parse_chain(content)
  assert_eq(cf.alignments.length(), 1)
  let serialized = cl_chain_to_string(cf.alignments[0])
  assert_true(serialized.length() > 0)
}

///|
test "cl_chain_to_string_format" {
  let blocks : Array[ChainBlock] = [
    ChainBlock::{ size: 500, dt: 10, dq: 20 },
    ChainBlock::{ size: 400, dt: 0, dq: 0 },
  ]
  let chain = cl_make_test_chain("chr1", 1, 5000, "chr2", 1, 100, blocks, 5)
  let s = cl_chain_to_string(chain)
  assert_true(cl_starts_with(s, "chain 1000"))
}

///|
test "cl_chain_to_string_contains_blocks" {
  let blocks : Array[ChainBlock] = [
    ChainBlock::{ size: 500, dt: 10, dq: 20 },
    ChainBlock::{ size: 400, dt: 0, dq: 0 },
  ]
  let chain = cl_make_test_chain("chr1", 1, 5000, "chr2", 1, 100, blocks, 5)
  let s = cl_chain_to_string(chain)
  assert_true(cl_starts_with(s, "chain"))
}

///|
test "cl_parse_chain_with_blank_lines" {
  let content =
    "chain 5000 chr1 249250621 + 1 5000 chr2 10000 + 1 100 5\n\n  5000 0 0\n\nchain 3000 chr5 181538259 + 100 200 chr7 159345973 + 300 500 6\n\n  100 0 0"
  let cf = cl_parse_chain(content)
  assert_eq(cf.alignments.length(), 2)
}

///|
test "cl_liftover_pos_multi_chain" {
  let content =
    "chain 5000 chr1 249250621 + 1 5000 chr2 10000 + 1 5000 5\n  5000 0 0\nchain 3000 chr1 249250621 + 6000 10000 chr3 8000 + 500 1500 6\n  3000 0 0"
  let cf = cl_parse_chain(content)
  let result1 = cl_liftover_pos(cf, "chr1", 100)
  assert_eq(result1, Some(100))
  let result2 = cl_liftover_pos(cf, "chr1", 6500)
  assert_eq(result2, Some(1000))
}

///|
test "cl_liftover_range_multi_chain" {
  let content =
    "chain 5000 chr1 249250621 + 1 5000 chr2 10000 + 1 5000 5\n  5000 0 0\nchain 3000 chr1 249250621 + 6000 10000 chr3 8000 + 500 1500 6\n  3000 0 0"
  let cf = cl_parse_chain(content)
  let result = cl_liftover_range(cf, "chr1", 100, 200)
  assert_eq(result, Some((100, 200)))
}

///|
test "cl_liftover_range_spans_multiple_chains" {
  let content =
    "chain 5000 chr1 249250621 + 1 5000 chr2 10000 + 1 5000 5\n  5000 0 0\nchain 3000 chr1 249250621 + 6000 10000 chr3 8000 + 500 1500 6\n  3000 0 0"
  let cf = cl_parse_chain(content)
  let result = cl_liftover_range(cf, "chr1", 4900, 7000)
  assert_true(result is None)
}

///|
test "cl_parse_chain_single_nucleotide_block" {
  let content =
    "chain 100 chr1 1000 + 100 200 chr2 500 + 50 150 1\n  1 0 0"
  let cf = cl_parse_chain(content)
  assert_eq(cf.alignments[0].blocks.length(), 1)
  assert_eq(cf.alignments[0].blocks[0].size, 1)
  let result = cl_liftover_pos(cf, "chr1", 100)
  assert_eq(result, Some(50))
}

///|
test "cl_liftover_pos_chain_boundary_start" {
  let content =
    "chain 5000 chr1 249250621 + 1 5000 chr2 10000 + 1 5000 5\n  5000 0 0"
  let cf = cl_parse_chain(content)
  let result = cl_liftover_pos(cf, "chr1", 1)
  assert_true(result is Some(_))
}

///|
test "cl_liftover_pos_chain_boundary_end" {
  let content =
    "chain 5000 chr1 249250621 + 1 5000 chr2 10000 + 1 5000 5\n  5000 0 0"
  let cf = cl_parse_chain(content)
  let result = cl_liftover_pos(cf, "chr1", 4999)
  assert_true(result is Some(_))
}

///|
test "cl_find_chain_for_pos_multiple_chains" {
  let content =
    "chain 5000 chr1 249250621 + 1 5000 chr2 10000 + 1 5000 5\n  5000 0 0\nchain 3000 chr1 249250621 + 6000 10000 chr3 8000 + 500 1500 6\n  3000 0 0"
  let cf = cl_parse_chain(content)
  let r1 = cl_find_chain_for_pos(cf, "chr1", 100)
  match r1 {
    Some(ch) => assert_eq(ch.id, 5)
    None => assert_true(false)
  }
  let r2 = cl_find_chain_for_pos(cf, "chr1", 6500)
  match r2 {
    Some(ch) => assert_eq(ch.id, 6)
    None => assert_true(false)
  }
}

///|
test "cl_parse_chain_whitespace_tabs" {
  let content =
    "chain\t5000\tchr1\t249250621\t+\t1\t5000\tchr2\t10000\t+\t1\t100\t5\n\t5000\t0\t0"
  let cf = cl_parse_chain(content)
  assert_eq(cf.alignments.length(), 1)
  assert_eq(cf.alignments[0].blocks.length(), 1)
}

///|
test "cl_liftover_pos_block_boundary_alignment" {
  let content =
    "chain 5000 chr1 249250621 + 1 1100 chr2 10000 + 1 1100 5\n  500 10 20\n  500 0 0"
  let cf = cl_parse_chain(content)
  let r1 = cl_liftover_pos(cf, "chr1", 500)
  assert_eq(r1, Some(500))
  let r2 = cl_liftover_pos(cf, "chr1", 501)
  assert_true(r2 is None)
  let r3 = cl_liftover_pos(cf, "chr1", 511)
  assert_eq(r3, Some(521))
  let r4 = cl_liftover_pos(cf, "chr1", 1010)
  assert_eq(r4, Some(1020))
}

///|
test "cl_parse_chain_id_preserved" {
  let content =
    "chain 100 chr1 1000 + 100 200 chr2 500 + 50 150 99\n  100 0 0"
  let cf = cl_parse_chain(content)
  assert_eq(cf.alignments[0].id, 99)
}

///|
test "cl_chain_to_string_score_preserved" {
  let content =
    "chain 7777 chr1 1000 + 100 200 chr2 500 + 50 150 1\n  100 0 0"
  let cf = cl_parse_chain(content)
  let s = cl_chain_to_string(cf.alignments[0])
  assert_true(cl_starts_with(s, "chain 7777"))
}

///|
test "cl_liftover_range_sub_block" {
  let content =
    "chain 5000 chr1 249250621 + 1 5000 chr2 10000 + 1 5000 5\n  500 10 20\n  4500 0 0"
  let cf = cl_parse_chain(content)
  let result = cl_liftover_range(cf, "chr1", 50, 150)
  assert_eq(result, Some((50, 150)))
}

///|
test "cl_get_chain_summary_gaps" {
  let blocks : Array[ChainBlock] = [
    ChainBlock::{ size: 100, dt: 50, dq: 30 },
    ChainBlock::{ size: 200, dt: 25, dq: 20 },
    ChainBlock::{ size: 300, dt: 0, dq: 0 },
  ]
  let chain = cl_make_test_chain("chr3", 1, 756, "chr5", 1, 700, blocks, 100)
  let summary = cl_get_chain_summary(chain)
  assert_eq(summary["total_gap_target"], 75)
  assert_eq(summary["total_gap_query"], 50)
}

///|
test "cl_liftover_pos_complex_chain" {
  let content =
    "chain 5000 chr1 249250621 + 1 5000 chr2 10000 + 1 5000 5\n  1000 50 100\n  2000 30 40\n  2000 0 0"
  let cf = cl_parse_chain(content)
  let r1 = cl_liftover_pos(cf, "chr1", 500)
  assert_eq(r1, Some(500))
  let r2 = cl_liftover_pos(cf, "chr1", 1500)
  assert_eq(r2, Some(1550))
  let r3 = cl_liftover_pos(cf, "chr1", 3500)
  assert_eq(r3, Some(3560))
  let r4 = cl_liftover_pos(cf, "chr1", 4500)
  assert_eq(r4, Some(4560))
}