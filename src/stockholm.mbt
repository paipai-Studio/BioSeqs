///|
/// Stockholm format parser and writer, mirroring Biopython Bio.Stockholm.
///
/// The Stockholm format is used by Pfam/Rfam for multiple sequence alignments
/// with secondary structure annotation.
///
/// Format example:
///   # STOCKHOLM 1.0
///   #=GF ID   entry_name
///   #=GF AC   accession
///   #=GF DE   description
///   #=GC SS_cons    ........((((...)))).....
///   seq1/1-100    ACGTACGTACGTACGT
///   seq2/1-95     ACGTACGTACGTACG-
///   #=GR seq1    ........((((...)))).....
///   //

///|
/// A single sequence in a Stockholm alignment block.
pub struct StockholmSequence {
  name : String
  accession : String
  description : String
  aligned_seq : String
  start : Int
  end : Int
}

///|
/// Create a new StockholmSequence.
pub fn StockholmSequence::new(
  name : String,
  accession : String,
  description : String,
  aligned_seq : String,
  start : Int,
  end : Int,
) -> StockholmSequence {
  StockholmSequence::{ name, accession, description, aligned_seq, start, end }
}

///|
/// A block within a Stockholm alignment.
/// Blocks are separated by empty lines and contain sequences
/// and annotations for a specific region.
pub struct StockholmBlock {
  sequences : Array[StockholmSequence]
  secondary_structure : String
  consensus : String
  annotation_lines : Array[String]
  gc_annotations : Map[String, String]
  gr_annotations : Map[String, String]
}

///|
/// Create a new StockholmBlock.
pub fn StockholmBlock::new(
  sequences : Array[StockholmSequence],
  secondary_structure : String,
  consensus : String,
  annotation_lines : Array[String],
  gc_annotations : Map[String, String],
  gr_annotations : Map[String, String],
) -> StockholmBlock {
  StockholmBlock::{
    sequences,
    secondary_structure,
    consensus,
    annotation_lines,
    gc_annotations,
    gr_annotations,
  }
}

///|
/// Full Stockholm alignment with blocks and annotations.
pub struct StockholmAlignment {
  version : String
  blocks : Array[StockholmBlock]
  gf_annotations : Map[String, String]
  gs_annotations : Map[String, String]
  markup : Array[String]
}

///|
/// Create a new StockholmAlignment.
pub fn StockholmAlignment::new(
  version : String,
  blocks : Array[StockholmBlock],
  gf_annotations : Map[String, String],
  gs_annotations : Map[String, String],
  markup : Array[String],
) -> StockholmAlignment {
  StockholmAlignment::{ version, blocks, gf_annotations, gs_annotations, markup }
}

///|
/// Parse Stockholm format string into a StockholmAlignment.
pub fn stockholm_parse(content : String) -> StockholmAlignment {
  let lines = content.split("\n").to_array()
  let n = lines.length()

  let version = "1.0"
  let gf_annotations : Map[String, String] = Map([], capacity=16)
  let gs_annotations : Map[String, String] = Map([], capacity=16)
  let markup : Array[String] = Array::new()

  let blocks : Array[StockholmBlock] = Array::new()

  let mut i = 0
  let mut current_seqs : Array[StockholmSequence] = Array::new()
  let mut current_gc : Map[String, String] = Map([], capacity=8)
  let mut current_gr : Map[String, String] = Map([], capacity=8)
  let mut current_annotations : Array[String] = Array::new()
  let mut has_block_content = false
  let mut done = false

  while i < n && !done {
    let line = lines[i].to_string()
    let trimmed = trim(line)

    if trimmed.length() == 0 {
      if has_block_content {
        let block = stockholm_build_block(
          current_seqs, current_gc, current_gr, current_annotations,
        )
        blocks.push(block)
        current_seqs = Array::new()
        current_gc = Map([], capacity=8)
        current_gr = Map([], capacity=8)
        current_annotations = Array::new()
        has_block_content = false
      }
      i = i + 1
      continue
    }

    if starts_with(trimmed, 0, "//") {
      if has_block_content {
        let block = stockholm_build_block(
          current_seqs, current_gc, current_gr, current_annotations,
        )
        blocks.push(block)
        current_seqs = Array::new()
        current_gc = Map([], capacity=8)
        current_gr = Map([], capacity=8)
        current_annotations = Array::new()
        has_block_content = false
      }
      done = true
      i = i + 1
      continue
    }

    if starts_with(trimmed, 0, "# STOCKHOLM") {
      i = i + 1
      continue
    }

    if starts_with(trimmed, 0, "#=GF") {
      let rest = trim(trimmed[4:].to_string())
      let parts = split_by_char(rest, ' '.to_int())
      let filtered = stockholm_filter_parts(parts)
      if filtered.length() >= 2 {
        let field = filtered[0]
        let mut value = ""
        let mut pi = 1
        while pi < filtered.length() {
          if pi > 1 { value = value + " " }
          value = value + filtered[pi]
          pi = pi + 1
        }
        match gf_annotations.get(field) {
          Some(existing) => gf_annotations[field] = existing + " " + value
          None => gf_annotations[field] = value
        }
      } else if filtered.length() == 1 {
        match gf_annotations.get(filtered[0]) {
          Some(existing) => gf_annotations[filtered[0]] = existing
          None => gf_annotations[filtered[0]] = ""
        }
      }
      has_block_content = true
      i = i + 1
      continue
    }

    if starts_with(trimmed, 0, "#=GC") {
      let rest = trim(trimmed[4:].to_string())
      let gc_parts = split_by_char(rest, ' '.to_int())
      let gc_filtered = stockholm_filter_parts(gc_parts)
      if gc_filtered.length() >= 2 {
        let field = gc_filtered[0]
        let mut value = ""
        let mut pi = 1
        while pi < gc_filtered.length() {
          if pi > 1 { value = value + " " }
          value = value + gc_filtered[pi]
          pi = pi + 1
        }
        match current_gc.get(field) {
          Some(existing) => current_gc[field] = existing + value
          None => current_gc[field] = value
        }
        current_annotations.push(trimmed)
        has_block_content = true
      }
      i = i + 1
      continue
    }

    if starts_with(trimmed, 0, "#=GR") {
      let rest = trim(trimmed[4:].to_string())
      let gr_parts = split_by_char(rest, ' '.to_int())
      let gr_filtered = stockholm_filter_parts(gr_parts)
      if gr_filtered.length() >= 2 {
        let name = gr_filtered[0]
        let mut value = ""
        let mut pi = 1
        while pi < gr_filtered.length() {
          if pi > 1 { value = value + " " }
          value = value + gr_filtered[pi]
          pi = pi + 1
        }
        match current_gr.get(name) {
          Some(existing) => current_gr[name] = existing + value
          None => current_gr[name] = value
        }
        current_annotations.push(trimmed)
        has_block_content = true
      }
      i = i + 1
      continue
    }

    if starts_with(trimmed, 0, "#=GS") {
      let rest = trim(trimmed[4:].to_string())
      let gs_parts = split_by_char(rest, ' '.to_int())
      let gs_filtered = stockholm_filter_parts(gs_parts)
      if gs_filtered.length() >= 2 {
        let name = gs_filtered[0]
        let mut value = ""
        let mut pi = 1
        while pi < gs_filtered.length() {
          if pi > 1 { value = value + " " }
          value = value + gs_filtered[pi]
          pi = pi + 1
        }
        match gs_annotations.get(name) {
          Some(existing) => gs_annotations[name] = existing + " " + value
          None => gs_annotations[name] = value
        }
        has_block_content = true
      }
      i = i + 1
      continue
    }

    if starts_with(trimmed, 0, "#=GA") {
      current_annotations.push(trimmed)
      has_block_content = true
      i = i + 1
      continue
    }

    if starts_with(trimmed, 0, "#") {
      markup.push(trimmed)
      has_block_content = true
      i = i + 1
      continue
    }

    let seq = stockholm_parse_sequence_line(trimmed)
    match seq {
      Some(s) => {
        current_seqs.push(s)
        has_block_content = true
      }
      None => {
        markup.push(trimmed)
      }
    }
    i = i + 1
  }

  if has_block_content {
    let block = stockholm_build_block(
      current_seqs, current_gc, current_gr, current_annotations,
    )
    blocks.push(block)
  }

  StockholmAlignment::{
    version,
    blocks,
    gf_annotations,
    gs_annotations,
    markup,
  }
}

///|
/// Write StockholmAlignment to Stockholm format string.
pub fn stockholm_write(ali : StockholmAlignment) -> String {
  let buf = StringBuilder::new()
  buf.write_string("# STOCKHOLM ")
  buf.write_string(ali.version)
  buf.write_string("\n")

  let gf_keys = ali.gf_annotations.keys().collect()
  let mut gi = 0
  while gi < gf_keys.length() {
    let key = gf_keys[gi]
    let value = ali.gf_annotations[key]
    if value.length() > 0 {
      buf.write_string("#=GF ")
      buf.write_string(key)
      buf.write_string("   ")
      buf.write_string(value)
      buf.write_string("\n")
    } else {
      buf.write_string("#=GF ")
      buf.write_string(key)
      buf.write_string("\n")
    }
    gi = gi + 1
  }

  let gs_keys = ali.gs_annotations.keys().collect()
  let mut gs = 0
  while gs < gs_keys.length() {
    let key = gs_keys[gs]
    let value = ali.gs_annotations[key]
    buf.write_string("#=GS ")
    buf.write_string(key)
    buf.write_string(" ")
    buf.write_string(value)
    buf.write_string("\n")
    gs = gs + 1
  }

  let mut bi = 0
  while bi < ali.blocks.length() {
    let block = ali.blocks[bi]

    let gc_keys = block.gc_annotations.keys().collect()
    let mut gci = 0
    while gci < gc_keys.length() {
      let key = gc_keys[gci]
      let value = block.gc_annotations[key]
      buf.write_string("#=GC ")
      buf.write_string(key)
      buf.write_string("    ")
      buf.write_string(value)
      buf.write_string("\n")
      gci = gci + 1
    }

    let seq_names : Array[String] = Array::new()
    let mut si = 0
    while si < block.sequences.length() {
      let seq = block.sequences[si]
      seq_names.push(seq.name)
      si = si + 1
    }

    let max_name_len = stockholm_max_name_length(block.sequences)

    let mut sj = 0
    while sj < block.sequences.length() {
      let seq = block.sequences[sj]
      buf.write_string(seq.name)
      buf.write_string("/")
      buf.write_string(seq.start.to_string())
      buf.write_string("-")
      buf.write_string(seq.end.to_string())
      let padding = max_name_len - seq.name.length() + 4
      let mut p = 0
      while p < padding {
        buf.write_string(" ")
        p = p + 1
      }
      buf.write_string(seq.aligned_seq)
      buf.write_string("\n")
      sj = sj + 1
    }

    let gr_keys = block.gr_annotations.keys().collect()
    let mut gri = 0
    while gri < gr_keys.length() {
      let key = gr_keys[gri]
      let value = block.gr_annotations[key]
      buf.write_string("#=GR ")
      buf.write_string(key)
      let gr_padding = max_name_len - key.length() + 4
      let mut gp = 0
      while gp < gr_padding {
        buf.write_string(" ")
        gp = gp + 1
      }
      buf.write_string(value)
      buf.write_string("\n")
      gri = gri + 1
    }

    if bi < ali.blocks.length() - 1 {
      buf.write_string("\n")
    }
    bi = bi + 1
  }

  buf.write_string("//\n")
  buf.to_string()
}

///|
/// Calculate percent identity for a Stockholm block.
/// Returns the average pairwise percent identity, excluding gap-only positions.
pub fn stockholm_percent_identity(block : StockholmBlock) -> Double {
  let seqs = block.sequences
  let n = seqs.length()
  if n < 2 { return 0.0 }

  let len = seqs[0].aligned_seq.length()
  if len == 0 { return 0.0 }

  let mut total_identity = 0.0
  let mut total_comparisons = 0.0

  let mut i = 0
  while i < n {
    let mut j = i + 1
    while j < n {
      let mut matches = 0
      let mut comparable = 0
      let mut k = 0
      while k < len {
        let c1 = char_at(seqs[i].aligned_seq, k)
        let c2 = char_at(seqs[j].aligned_seq, k)
        let is_gap1 = c1 == "-" || c1 == "."
        let is_gap2 = c2 == "-" || c2 == "."
        if !is_gap1 && !is_gap2 {
          comparable = comparable + 1
          if c1 == c2 {
            matches = matches + 1
          }
        }
        k = k + 1
      }
      if comparable > 0 {
        total_identity = total_identity + matches.to_double() / comparable.to_double()
        total_comparisons = total_comparisons + 1.0
      }
      j = j + 1
    }
    i = i + 1
  }

  if total_comparisons > 0.0 {
    total_identity / total_comparisons
  } else {
    0.0
  }
}

///|
/// Calculate per-column conservation for a Stockholm block.
/// Returns an array of Double values representing conservation at each column.
/// Uses frequency of the most common residue (excluding gaps).
pub fn stockholm_conservation(block : StockholmBlock) -> Array[Double] {
  let result : Array[Double] = Array::new()
  let seqs = block.sequences
  let n = seqs.length()
  if n == 0 { return result }

  let len = seqs[0].aligned_seq.length()
  if len == 0 { return result }

  let mut col = 0
  while col < len {
    let counts : Map[String, Int] = Map([], capacity=20)
    let mut non_gap = 0

    let mut si = 0
    while si < n {
      let c = char_at(seqs[si].aligned_seq, col)
      let is_gap = c == "-" || c == "."
      if !is_gap {
        non_gap = non_gap + 1
        match counts.get(c) {
          Some(existing) => counts[c] = existing + 1
          None => counts[c] = 1
        }
      }
      si = si + 1
    }

    if non_gap == 0 {
      result.push(0.0)
    } else {
      let keys = counts.keys().collect()
      let mut max_count = 0
      let mut ki = 0
      while ki < keys.length() {
        if counts[keys[ki]] > max_count {
          max_count = counts[keys[ki]]
        }
        ki = ki + 1
      }
      result.push(max_count.to_double() / non_gap.to_double())
    }

    col = col + 1
  }

  result
}

///|
/// Convert Stockholm alignment to FASTA format string.
pub fn stockholm_to_fasta(ali : StockholmAlignment) -> String {
  let buf = StringBuilder::new()

  let seen : Map[String, Bool] = Map([], capacity=32)

  let mut bi = 0
  while bi < ali.blocks.length() {
    let block = ali.blocks[bi]
    let mut si = 0
    while si < block.sequences.length() {
      let seq = block.sequences[si]
      let display_name = stockholm_display_name(seq)

      match seen.get(seq.name) {
        Some(_) => {
          buf.write_string(seq.aligned_seq)
          buf.write_string("\n")
        }
        None => {
          buf.write_string(">")
          buf.write_string(display_name)
          buf.write_string("\n")
          buf.write_string(seq.aligned_seq)
          buf.write_string("\n")
          seen[seq.name] = true
        }
      }
      si = si + 1
    }
    bi = bi + 1
  }

  buf.to_string()
}

///|
/// Merge all blocks in a Stockholm alignment into a single block.
/// Sequences with the same name are concatenated across blocks.
pub fn stockholm_merge_blocks(ali : StockholmAlignment) -> StockholmAlignment {
  if ali.blocks.length() <= 1 {
    return ali
  }

  let merged_seqs : Array[StockholmSequence] = Array::new()
  let seq_map : Map[String, StringBuilder] = Map([], capacity=32)
  let seq_info : Map[String, StockholmSequence] = Map([], capacity=32)

  let mut bi = 0
  while bi < ali.blocks.length() {
    let block = ali.blocks[bi]
    let mut si = 0
    while si < block.sequences.length() {
      let seq = block.sequences[si]
      match seq_map.get(seq.name) {
        Some(sb) => {
          sb.write_string(seq.aligned_seq)
        }
        None => {
          let sb = StringBuilder::new()
          sb.write_string(seq.aligned_seq)
          seq_map[seq.name] = sb
          seq_info[seq.name] = seq
        }
      }
      si = si + 1
    }
    bi = bi + 1
  }

  let names = seq_map.keys().collect()
  let mut ni = 0
  while ni < names.length() {
    let name = names[ni]
    let sb = seq_map[name]
    let info = seq_info[name]
    merged_seqs.push(StockholmSequence::{
      name: info.name,
      accession: info.accession,
      description: info.description,
      aligned_seq: sb.to_string(),
      start: info.start,
      end: info.end,
    })
    ni = ni + 1
  }

  let merged_gc : Map[String, String] = Map([], capacity=8)
  let merged_gr : Map[String, String] = Map([], capacity=8)
  let merged_annotations : Array[String] = Array::new()

  let gc_keys : Array[String] = Array::new()
  let gr_keys : Array[String] = Array::new()

  let mut ci = 0
  while ci < ali.blocks.length() {
    let block = ali.blocks[ci]

    let block_gc_keys = block.gc_annotations.keys().collect()
    let mut gci = 0
    while gci < block_gc_keys.length() {
      let key = block_gc_keys[gci]
      let value = block.gc_annotations[key]
      match merged_gc.get(key) {
        Some(existing) => merged_gc[key] = existing + value
        None => {
          merged_gc[key] = value
          gc_keys.push(key)
        }
      }
      gci = gci + 1
    }

    let block_gr_keys = block.gr_annotations.keys().collect()
    let mut gri = 0
    while gri < block_gr_keys.length() {
      let key = block_gr_keys[gri]
      let value = block.gr_annotations[key]
      match merged_gr.get(key) {
        Some(existing) => merged_gr[key] = existing + value
        None => {
          merged_gr[key] = value
          gr_keys.push(key)
        }
      }
      gri = gri + 1
    }

    let mut ai = 0
    while ai < block.annotation_lines.length() {
      merged_annotations.push(block.annotation_lines[ai])
      ai = ai + 1
    }

    ci = ci + 1
  }

  let merged_block = StockholmBlock::{
    sequences: merged_seqs,
    secondary_structure: match merged_gc.get("SS_cons") {
      Some(ss) => ss
      None => ""
    },
    consensus: "",
    annotation_lines: merged_annotations,
    gc_annotations: merged_gc,
    gr_annotations: merged_gr,
  }

  StockholmAlignment::{
    version: ali.version,
    blocks: [merged_block],
    gf_annotations: ali.gf_annotations,
    gs_annotations: ali.gs_annotations,
    markup: ali.markup,
  }
}

// --- Internal helper functions ---

fn stockholm_filter_parts(parts : Array[String]) -> Array[String] {
  let result : Array[String] = Array::new()
  let mut i = 0
  while i < parts.length() {
    if parts[i].length() > 0 {
      result.push(parts[i])
    }
    i = i + 1
  }
  result
}

fn stockholm_parse_sequence_line(line : String) -> StockholmSequence? {
  if line.length() == 0 { return None }

  let sp = ' '.to_int().to_uint16()
  let mut space_pos = -1
  let mut i = 0
  while i < line.length() {
    if line.unsafe_get(i) == sp {
      space_pos = i
      break
    }
    i = i + 1
  }

  if space_pos <= 0 { return None }

  let name_part = line[0:space_pos].to_owned()
  let seq_part_start = space_pos + 1
  let seq_part = if seq_part_start < line.length() {
    trim(line[seq_part_start:].to_owned())
  } else {
    ""
  }

  if seq_part.length() == 0 { return None }

  let slash_pos = stockholm_find_slash(name_part)
  let (name, start, end) = if slash_pos >= 0 {
    let n = name_part[0:slash_pos].to_owned()
    let range = name_part[slash_pos + 1:].to_owned()
    let dash_pos = stockholm_find_dash(range)
    if dash_pos >= 0 {
      (n, parse_int(range[0:dash_pos].to_owned()), parse_int(range[dash_pos + 1:].to_owned()))
    } else {
      (n, 1, parse_int(range))
    }
  } else {
    (name_part, 1, seq_part.length())
  }

  Some(StockholmSequence::{
    name,
    accession: "",
    description: "",
    aligned_seq: seq_part,
    start,
    end,
  })
}

fn stockholm_find_slash(s : String) -> Int {
  let mut i = 0
  while i < s.length() {
    if s.unsafe_get(i).to_int() == '/'.to_int() {
      return i
    }
    i = i + 1
  }
  -1
}

fn stockholm_find_dash(s : String) -> Int {
  let mut i = 0
  while i < s.length() {
    if s.unsafe_get(i).to_int() == '-'.to_int() {
      return i
    }
    i = i + 1
  }
  -1
}

fn stockholm_build_block(
  seqs : Array[StockholmSequence],
  gc : Map[String, String],
  gr : Map[String, String],
  annotations : Array[String],
) -> StockholmBlock {
  let ss = match gc.get("SS_cons") {
    Some(s) => s
    None => ""
  }

  StockholmBlock::{
    sequences: seqs,
    secondary_structure: ss,
    consensus: "",
    annotation_lines: annotations,
    gc_annotations: gc,
    gr_annotations: gr,
  }
}

fn stockholm_max_name_length(seqs : Array[StockholmSequence]) -> Int {
  let mut max_len = 0
  let mut i = 0
  while i < seqs.length() {
    let name_len = seqs[i].name.length()
    if name_len > max_len {
      max_len = name_len
    }
    i = i + 1
  }
  max_len
}

fn stockholm_display_name(seq : StockholmSequence) -> String {
  if seq.accession.length() > 0 {
    seq.accession
  } else if seq.description.length() > 0 {
    seq.name + " " + seq.description
  } else {
    seq.name
  }
}

///|
/// Create a sample Stockholm alignment content for testing.
pub fn sample_stockholm_content() -> String {
  "# STOCKHOLM 1.0\n" +
  "#=GF ID   TestFamily\n" +
  "#=GF AC   PF00001.20\n" +
  "#=GF DE   A test protein family\n" +
  "#=GC SS_cons    ........((((...)))).....\n" +
  "seq1/1-50    ACGTACGTACGTACGTACGTACGTACGTACGTACGT\n" +
  "seq2/1-50    ACGTACGTACGTACGTACGTACGTACGTACGTACGT\n" +
  "seq3/1-50    ACGTACGTACGTACGTACGTACGTACGTACGTACGT\n" +
  "#=GR seq1    ........((((...)))).....\n" +
  "//\n"
}

///|
/// Create a sample Stockholm alignment for testing.
pub fn sample_stockholm_alignment() -> StockholmAlignment {
  stockholm_parse(sample_stockholm_content())
}