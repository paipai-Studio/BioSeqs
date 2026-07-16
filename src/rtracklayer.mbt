///|
/// rtracklayer - Track format import/export module inspired by Bioconductor rtracklayer.
/// Supports BED, WIG, BEDGraph, and GFF/GTF formats.

///|
/// BED record.
pub struct BedRecord {
  chrom : String
  chrom_start : Int
  chrom_end : Int
  name : String
  score : Int?
  strand : Strand
  thick_start : Int?
  thick_end : Int?
  item_rgb : String
  block_count : Int?
  block_sizes : Array[Int]
  block_starts : Array[Int]
}

///|
/// Parse BED format (3-12 columns).
pub fn parse_bed(content : String) -> Array[BedRecord] {
  let lines = split_by_char(content, '\n'.to_int())
  let records : Array[BedRecord] = Array::new()
  
  let mut i = 0
  while i < lines.length() {
    let line = lines[i]
    if line.length() == 0 || line.unsafe_get(0).unsafe_to_char() == '#' {
      i = i + 1
      continue
    }
    
    let fields = split_by_char(line, '\t'.to_int())
    if fields.length() < 3 {
      i = i + 1
      continue
    }
    
    let name = if fields.length() > 3 { fields[3] } else { "" }
    let score = if fields.length() > 4 { Some(parse_int(fields[4])) } else { None }
    let strand = if fields.length() > 5 {
      match fields[5] {
        "+" => Strand::Plus
        "-" => Strand::Minus
        _ => Strand::Star
      }
    } else { Strand::Star }
    let thick_start = if fields.length() > 6 { Some(parse_int(fields[6])) } else { None }
    let thick_end = if fields.length() > 7 { Some(parse_int(fields[7])) } else { None }
    let item_rgb = if fields.length() > 8 { fields[8] } else { "" }
    let block_count = if fields.length() > 9 { Some(parse_int(fields[9])) } else { None }
    
    let block_sizes = if fields.length() > 10 {
      let sizes = split_by_char(fields[10], ','.to_int())
      let sizes_arr : Array[Int] = Array::new()
      let mut j = 0
      while j < sizes.length() {
        if sizes[j].length() > 0 {
          sizes_arr.push(parse_int(sizes[j]))
        }
        j = j + 1
      }
      sizes_arr
    } else { [] }
    
    let block_starts = if fields.length() > 11 {
      let starts = split_by_char(fields[11], ','.to_int())
      let starts_arr : Array[Int] = Array::new()
      let mut j = 0
      while j < starts.length() {
        if starts[j].length() > 0 {
          starts_arr.push(parse_int(starts[j]))
        }
        j = j + 1
      }
      starts_arr
    } else { [] }
    
    records.push(BedRecord::{
      chrom: fields[0],
      chrom_start: parse_int(fields[1]),
      chrom_end: parse_int(fields[2]),
      name,
      score,
      strand,
      thick_start,
      thick_end,
      item_rgb,
      block_count,
      block_sizes,
      block_starts,
    })
    
    i = i + 1
  }
  
  records
}

///|
/// Convert BedRecord array to GRanges.
pub fn bed_to_granges(records : Array[BedRecord]) -> GRanges {
  let n = records.length()
  let seqnames : Array[String] = Array::make(n, "")
  let ranges : Array[(Int, Int)] = Array::make(n, (0, 0))
  let strands : Array[Strand] = Array::make(n, Strand::Star)
  
  let mut i = 0
  while i < n {
    let r = records[i]
    seqnames[i] = r.chrom
    ranges[i] = (r.chrom_start, r.chrom_end)
    strands[i] = r.strand
    i = i + 1
  }
  
  granges(seqnames, ranges, strands)
}

///|
/// Write BedRecord array to BED format string.
pub fn write_bed(records : Array[BedRecord]) -> String {
  let mut result = ""
  
  let mut i = 0
  while i < records.length() {
    let r = records[i]
    let mut line = r.chrom + "\t" + r.chrom_start.to_string() + "\t" + r.chrom_end.to_string()
    
    if r.name.length() > 0 {
      line = line + "\t" + r.name
    }
    match r.score {
      Some(s) => line = line + "\t" + s.to_string()
      None => ()
    }
    match r.strand {
      Strand::Plus => line = line + "\t+"
      Strand::Minus => line = line + "\t-"
      Strand::Star => line = line + "\t."
    }
    match r.thick_start {
      Some(s) => line = line + "\t" + s.to_string()
      None => ()
    }
    match r.thick_end {
      Some(e) => line = line + "\t" + e.to_string()
      None => ()
    }
    if r.item_rgb.length() > 0 {
      line = line + "\t" + r.item_rgb
    }
    match r.block_count {
      Some(bc) => {
        line = line + "\t" + bc.to_string()
        let mut sizes = ""
        let mut j = 0
        while j < r.block_sizes.length() {
          if j > 0 { sizes = sizes + "," }
          sizes = sizes + r.block_sizes[j].to_string()
          j = j + 1
        }
        line = line + "\t" + sizes
        let mut starts = ""
        j = 0
        while j < r.block_starts.length() {
          if j > 0 { starts = starts + "," }
          starts = starts + r.block_starts[j].to_string()
          j = j + 1
        }
        line = line + "\t" + starts
      }
      None => ()
    }
    
    if i > 0 { result = result + "\n" }
    result = result + line
    i = i + 1
  }
  
  result
}

///|
/// WIG record types.
pub enum WigType {
  VariableStep(String, Int?)
  FixedStep(String, Int, Int)
}

///|
/// WIG record.
pub struct WigRecord {
  wig_type : WigType
  positions : Array[Int]
  values : Array[Double]
}

///|
/// Parse WIG format.
pub fn parse_wig(content : String) -> Array[WigRecord] {
  let lines = split_by_char(content, '\n'.to_int())
  let records : Array[WigRecord] = Array::new()
  
  let mut current_chrom = ""
  let mut current_start = 0
  let mut current_step = 1
  let mut current_span : Int? = None
  let mut is_fixed_step = false
  let mut positions : Array[Int] = Array::new()
  let mut values : Array[Double] = Array::new()
  
  let mut i = 0
  while i < lines.length() {
    let line = lines[i]
    if line.length() == 0 || line.unsafe_get(0).unsafe_to_char() == '#' {
      i = i + 1
      continue
    }
    
    let trimmed = trim(line)
    if starts_with(trimmed, 0, "variableStep") || starts_with(trimmed, 0, "VariableStep") {
      if positions.length() > 0 {
        let wig_type = if is_fixed_step {
          WigType::FixedStep(current_chrom, current_start, current_step)
        } else {
          WigType::VariableStep(current_chrom, current_span)
        }
        records.push(WigRecord::{ wig_type, positions: positions.clone(), values: values.clone() })
        positions = Array::new()
        values = Array::new()
      }
      current_chrom = extract_param(trimmed, "chrom=")
      let span_str = extract_param(trimmed, "span=")
      current_span = if span_str.length() > 0 { Some(parse_int(span_str)) } else { None }
      is_fixed_step = false
    } else if starts_with(trimmed, 0, "fixedStep") || starts_with(trimmed, 0, "FixedStep") {
      if positions.length() > 0 {
        let wig_type = if is_fixed_step {
          WigType::FixedStep(current_chrom, current_start, current_step)
        } else {
          WigType::VariableStep(current_chrom, current_span)
        }
        records.push(WigRecord::{ wig_type, positions: positions.clone(), values: values.clone() })
        positions = Array::new()
        values = Array::new()
      }
      current_chrom = extract_param(trimmed, "chrom=")
      let start_str = extract_param(trimmed, "start=")
      let step_str = extract_param(trimmed, "step=")
      current_start = if start_str.length() > 0 { parse_int(start_str) } else { 1 }
      current_step = if step_str.length() > 0 { parse_int(step_str) } else { 1 }
      is_fixed_step = true
    } else if starts_with(trimmed, 0, "track") || starts_with(trimmed, 0, "browser") {
      i = i + 1
      continue
    } else {
      let fields = split_by_char(trimmed, '\t'.to_int())
      if fields.length() >= 2 {
        positions.push(parse_int(fields[0]))
        match parse_double(fields[1]) {
          Some(v) => values.push(v)
          None => values.push(0.0)
        }
      } else if fields.length() == 1 {
        positions.push(current_start)
        match parse_double(fields[0]) {
          Some(v) => values.push(v)
          None => values.push(0.0)
        }
        current_start = current_start + current_step
      }
    }
    
    i = i + 1
  }
  
  if positions.length() > 0 {
    let wig_type = if is_fixed_step {
      WigType::FixedStep(current_chrom, current_start - current_step, current_step)
    } else {
      WigType::VariableStep(current_chrom, current_span)
    }
    records.push(WigRecord::{ wig_type, positions, values })
  }
  
  records
}

///|
/// Write WIG format.
pub fn write_wig(records : Array[WigRecord]) -> String {
  let mut result = ""
  
  let mut i = 0
  while i < records.length() {
    let r = records[i]
    
    match r.wig_type {
      WigType::VariableStep(chrom, span) => {
        result = result + "variableStep chrom=" + chrom
        match span {
          Some(s) => result = result + " span=" + s.to_string()
          None => ()
        }
        result = result + "\n"
      }
      WigType::FixedStep(chrom, start, step) => {
        result = result + "fixedStep chrom=" + chrom + " start=" + start.to_string() + " step=" + step.to_string() + "\n"
      }
    }
    
    let mut j = 0
    while j < r.positions.length() {
      result = result + r.positions[j].to_string() + "\t" + r.values[j].to_string() + "\n"
      j = j + 1
    }
    
    i = i + 1
  }
  
  result
}

///|
/// BEDGraph record.
pub struct BedGraphRecord {
  chrom : String
  chrom_start : Int
  chrom_end : Int
  value : Double
}

///|
/// Parse BEDGraph format.
pub fn parse_bedgraph(content : String) -> Array[BedGraphRecord] {
  let lines = split_by_char(content, '\n'.to_int())
  let records : Array[BedGraphRecord] = Array::new()
  
  let mut i = 0
  while i < lines.length() {
    let line = lines[i]
    if line.length() == 0 || line.unsafe_get(0).unsafe_to_char() == '#' || starts_with(line, 0, "track") {
      i = i + 1
      continue
    }
    
    let fields = split_by_char(line, '\t'.to_int())
    if fields.length() >= 4 {
      let value = match parse_double(fields[3]) {
        Some(v) => v
        None => 0.0
      }
      records.push(BedGraphRecord::{
        chrom: fields[0],
        chrom_start: parse_int(fields[1]),
        chrom_end: parse_int(fields[2]),
        value,
      })
    }
    
    i = i + 1
  }
  
  records
}

///|
/// Write BEDGraph format.
pub fn write_bedgraph(records : Array[BedGraphRecord]) -> String {
  let mut result = "track type=bedGraph\n"
  
  let mut i = 0
  while i < records.length() {
    let r = records[i]
    let line = r.chrom + "\t" + r.chrom_start.to_string() + "\t" + r.chrom_end.to_string() + "\t" + r.value.to_string()
    if i > 0 { result = result + "\n" }
    result = result + line
    i = i + 1
  }
  
  result
}

///|
/// GFF/GTF record.
pub struct GffRecord {
  seqid : String
  source : String
  feature_type : String
  start : Int
  end : Int
  score : Double?
  strand : Strand
  phase : Int?
  attributes : Map[String, String]
}

///|
/// Parse GFF/GTF format (GFF3 and GTF).
pub fn parse_gff(content : String) -> Array[GffRecord] {
  let lines = split_by_char(content, '\n'.to_int())
  let records : Array[GffRecord] = Array::new()
  
  let mut i = 0
  while i < lines.length() {
    let line = lines[i]
    if line.length() == 0 || line.unsafe_get(0).unsafe_to_char() == '#' {
      i = i + 1
      continue
    }
    
    let fields = split_by_char(line, '\t'.to_int())
    if fields.length() < 8 {
      i = i + 1
      continue
    }
    
    let strand = match fields[6] {
      "+" => Strand::Plus
      "-" => Strand::Minus
      "." => Strand::Star
      _ => Strand::Star
    }
    
    let score = if fields[5] != "." {
      match parse_double(fields[5]) {
        Some(s) => Some(s)
        None => None
      }
    } else { None }
    
    let phase = if fields[7] != "." { Some(parse_int(fields[7])) } else { None }
    
    let attributes = if fields.length() > 8 {
      parse_gff_attributes(fields[8])
    } else { Map([], capacity=8) }
    
    records.push(GffRecord::{
      seqid: fields[0],
      source: fields[1],
      feature_type: fields[2],
      start: parse_int(fields[3]),
      end: parse_int(fields[4]),
      score,
      strand,
      phase,
      attributes,
    })
    
    i = i + 1
  }
  
  records
}

///|
/// Parse GFF attributes string.
fn parse_gff_attributes(attr_str : String) -> Map[String, String] {
  let attrs = Map([], capacity=8)
  let parts = split_by_char(attr_str, ';'.to_int())
  
  let mut i = 0
  while i < parts.length() {
    let part = trim(parts[i])
    if part.length() == 0 {
      i = i + 1
      continue
    }
    
    let eq_pos = find_char(part, '=')
    if eq_pos >= 0 {
      let key = substring(part, 0, eq_pos)
      let value_len = part.length() - eq_pos - 1
      let mut value = substring(part, eq_pos + 1, value_len)
      if value.length() >= 2 && value.unsafe_get(0).unsafe_to_char() == '"' && value.unsafe_get(value.length() - 1).unsafe_to_char() == '"' {
        value = substring(value, 1, value.length() - 2)
      }
      attrs[key] = value
    }
    
    i = i + 1
  }
  
  attrs
}

///|
/// Write GFF format.
pub fn write_gff(records : Array[GffRecord]) -> String {
  let mut result = "##gff-version 3\n"
  
  let mut i = 0
  while i < records.length() {
    let r = records[i]
    let strand_str = match r.strand {
      Strand::Plus => "+"
      Strand::Minus => "-"
      Strand::Star => "."
    }
    let score_str = match r.score {
      Some(s) => s.to_string()
      None => "."
    }
    let phase_str = match r.phase {
      Some(p) => p.to_string()
      None => "."
    }
    
    let mut attrs = ""
    for key, value in r.attributes {
      if attrs.length() > 0 { attrs = attrs + "; " }
      attrs = attrs + key + "=\"" + value + "\""
    }
    
    let line = r.seqid + "\t" + r.source + "\t" + r.feature_type + "\t" + r.start.to_string() + "\t" + r.end.to_string() + "\t" + score_str + "\t" + strand_str + "\t" + phase_str + "\t" + attrs
    if i > 0 { result = result + "\n" }
    result = result + line
    i = i + 1
  }
  
  result
}

///|
/// Convert GffRecord array to GRanges.
pub fn gff_to_granges(records : Array[GffRecord]) -> GRanges {
  let n = records.length()
  let seqnames : Array[String] = Array::make(n, "")
  let ranges : Array[(Int, Int)] = Array::make(n, (0, 0))
  let strands : Array[Strand] = Array::make(n, Strand::Star)
  
  let mut i = 0
  while i < n {
    let r = records[i]
    seqnames[i] = r.seqid
    ranges[i] = (r.start, r.end)
    strands[i] = r.strand
    i = i + 1
  }
  
  granges(seqnames, ranges, strands)
}

///|
/// Extract parameter from WIG header line.
fn extract_param(line : String, prefix : String) -> String {
  let parts = split_by_char(line, ' '.to_int())
  
  let mut i = 0
  while i < parts.length() {
    let part = parts[i]
    if starts_with(part, 0, prefix) {
      let value_start = prefix.length()
      return substring(part, value_start, part.length() - value_start)
    }
    i = i + 1
  }
  ""
}

///|
/// Trim whitespace from string.
fn trim(s : String) -> String {
  let mut start = 0
  while start < s.length() && is_whitespace(s.unsafe_get(start).unsafe_to_char()) {
    start = start + 1
  }
  
  let mut end = s.length()
  while end > start && is_whitespace(s.unsafe_get(end - 1).unsafe_to_char()) {
    end = end - 1
  }
  
  substring(s, start, end - start)
}

///|
fn is_whitespace(c : Char) -> Bool {
  c == ' ' || c == '\t' || c == '\r'
}

///|
fn find_char(s : String, c : Char) -> Int {
  let mut i = 0
  while i < s.length() {
    if s.unsafe_get(i).unsafe_to_char() == c {
      return i
    }
    i = i + 1
  }
  -1
}

///|
fn find_char_from(s : String, c : Char, start : Int) -> Int {
  let mut i = start
  while i < s.length() {
    if s.unsafe_get(i).unsafe_to_char() == c {
      return i
    }
    i = i + 1
  }
  -1
}