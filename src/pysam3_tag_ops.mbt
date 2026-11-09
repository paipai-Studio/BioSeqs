// Pysam third-round: Read and tag operations (items 25-28)
// B array tags, CIGAR =/X operators, RG→sample mapping, strand-specific counts

///|
pub enum Pysam3TagValue {
  Pysam3TagChar(Char)
  Pysam3TagInt(Int32)
  Pysam3TagFloat(Double)
  Pysam3TagString(String)
  Pysam3TagHex(String)
  Pysam3TagArray(Pysam3TagArrayType, Array[Int32])
}

///|
pub enum Pysam3TagArrayType {
  Pysam3TagArrayInt8
  Pysam3TagArrayUInt8
  Pysam3TagArrayInt16
  Pysam3TagArrayUInt16
  Pysam3TagArrayInt32
  Pysam3TagArrayUInt32
  Pysam3TagArrayFloat
}

///|
pub fn pysam3_tag_array_type_char(t : Pysam3TagArrayType) -> Char {
  match t {
    Pysam3TagArrayInt8 => 'c'
    Pysam3TagArrayUInt8 => 'C'
    Pysam3TagArrayInt16 => 's'
    Pysam3TagArrayUInt16 => 'S'
    Pysam3TagArrayInt32 => 'i'
    Pysam3TagArrayUInt32 => 'I'
    Pysam3TagArrayFloat => 'f'
  }
}

///|
pub fn pysam3_parse_tag_array(value : String) -> Pysam3TagValue {
  // Format: B:type,v1,v2,v3,...
  if value.length() < 3 || value[0].to_int() != 66 {
    return Pysam3TagString(value)
  }
  let type_char = value[1]
  let rest = value[3:value.length()].to_owned()
  let parts = rest.split(",").to_array()
  let arr : Array[Int32] = []
  for i in 0..<parts.length() {
    let s = parts[i].to_owned().trim().to_owned()
    let v = parse_int(s)
    arr.push(v)
  }
  let arr_type = match type_char {
    'c' => Pysam3TagArrayInt8
    'C' => Pysam3TagArrayUInt8
    's' => Pysam3TagArrayInt16
    'S' => Pysam3TagArrayUInt16
    'i' => Pysam3TagArrayInt32
    'I' => Pysam3TagArrayUInt32
    'f' => Pysam3TagArrayFloat
    _ => Pysam3TagArrayInt32
  }
  Pysam3TagArray(arr_type, arr)
}

///|
pub fn pysam3_format_tag_array(
  arr_type : Pysam3TagArrayType,
  values : Array[Int32],
) -> String {
  let sb = StringBuilder()
  let _ = sb.write("B,")
  let _ = sb.write_char(pysam3_tag_array_type_char(arr_type))
  for i in 0..<values.length() {
    let _ = sb.write(",")
    let _ = sb.write(values[i].to_string())
  }
  sb.to_string()
}

///|
pub fn pysam3_parse_tag_typed(
  _tag : String,
  type_char : Char,
  value : String,
) -> Pysam3TagValue {
  match type_char {
    'A' =>
      if value.length() > 0 {
        Pysam3TagChar(Int32::unsafe_to_char(value[0].to_int()))
      } else {
        Pysam3TagChar('?')
      }
    'i' => Pysam3TagInt(parse_int(value))
    'f' => Pysam3TagFloat(parse_double(value).unwrap_or(0.0))
    'Z' => Pysam3TagString(value)
    'H' => Pysam3TagHex(value)
    'B' => pysam3_parse_tag_array(value)
    _ => Pysam3TagString(value)
  }
}

// --- CIGAR =/X operators ---

///|
pub fn pysam3_cigar_has_seq_match(cigar : Array[(Int32, Int32)]) -> Bool {
  for i in 0..<cigar.length() {
    let (op, _) = cigar[i]
    if op == 7 || op == 8 {
      return true
    }
  }
  false
}

///|
pub fn pysam3_cigar_expand_eq_x(
  cigar : Array[(Int32, Int32)],
  ref_seq : String,
  read_seq : String,
  ref_start : Int32,
) -> Array[(Int32, Int32)] {
  // Convert M operations to = (match) or X (mismatch) based on sequences
  let result : Array[(Int32, Int32)] = []
  let mut ref_idx = ref_start
  let mut seq_idx = 0
  for i in 0..<cigar.length() {
    let (op, length) = cigar[i]
    if op == 0 {
      // M: expand into = and X
      let mut run_op = -1
      let mut run_len = 0
      for k in 0..<length {
        let r = if ref_idx + k < ref_seq.length() {
          ref_seq[ref_idx + k].to_int()
        } else {
          78
        }
        let s = if seq_idx + k < read_seq.length() {
          read_seq[seq_idx + k].to_int()
        } else {
          78
        }
        let cur_op = if r == s { 7 } else { 8 }
        if cur_op == run_op {
          run_len = run_len + 1
        } else {
          if run_len > 0 {
            result.push((run_op, run_len))
          }
          run_op = cur_op
          run_len = 1
        }
      }
      if run_len > 0 {
        result.push((run_op, run_len))
      }
      ref_idx = ref_idx + length
      seq_idx = seq_idx + length
    } else {
      result.push((op, length))
      if op == 0 || op == 7 || op == 8 {
        ref_idx = ref_idx + length
        seq_idx = seq_idx + length
      } else if op == 1 {
        seq_idx = seq_idx + length
      } else if op == 2 || op == 3 {
        ref_idx = ref_idx + length
      } else if op == 4 {
        seq_idx = seq_idx + length
      }
    }
  }
  result
}

///|
pub fn pysam3_cigar_collapse_eq_x(
  cigar : Array[(Int32, Int32)],
) -> Array[(Int32, Int32)] {
  // Collapse = and X back to M
  let result : Array[(Int32, Int32)] = []
  for i in 0..<cigar.length() {
    let (op, length) = cigar[i]
    if op == 7 || op == 8 {
      // Check if previous result is M
      if result.length() > 0 {
        let (prev_op, prev_len) = result[result.length() - 1]
        if prev_op == 0 {
          result[result.length() - 1] = (0, prev_len + length)
        } else {
          result.push((0, length))
        }
      } else {
        result.push((0, length))
      }
    } else {
      result.push((op, length))
    }
  }
  result
}

///|
pub fn pysam3_cigar_count_matches_mismatches(
  cigar : Array[(Int32, Int32)],
) -> (Int32, Int32) {
  // (matches, mismatches)
  let mut matches = 0
  let mut mismatches = 0
  for i in 0..<cigar.length() {
    let (op, length) = cigar[i]
    if op == 7 {
      matches = matches + length
    } else if op == 8 {
      mismatches = mismatches + length
    } else if op == 0 {
      // M: count as matches (unknown)
      matches = matches + length
    }
  }
  (matches, mismatches)
}

// --- RG→Sample mapping ---

///|
pub struct Pysam3ReadGroup {
  id_ : String
  sample : String
  library : String
  platform : String
  platform_unit : String
}

///|
pub fn pysam3_parse_rg_line(line : String) -> Pysam3ReadGroup {
  // Parse @RG\tID:xxx\tSM:yyy\tLB:zzz\tPL:illumina\tPU:unit1
  let mut id_ = ""
  let mut sample = ""
  let mut library = ""
  let mut platform = ""
  let mut platform_unit = ""
  let parts = line.split("\t").to_array()
  for i in 0..<parts.length() {
    let field = parts[i].to_owned()
    if i == 0 {
      continue
    }
    let mut colon_idx : Int32 = -1
    for ci in 0..<field.length() {
      if field[ci].to_int() == 58 {
        colon_idx = ci
        break
      }
    }
    if colon_idx >= 0 {
      let key = field[0:colon_idx].to_owned()
      let val = field[colon_idx + 1:field.length()].to_owned()
      match key {
        "ID" => id_ = val
        "SM" => sample = val
        "LB" => library = val
        "PL" => platform = val
        "PU" => platform_unit = val
        _ => ()
      }
    }
  }
  { id_, sample, library, platform, platform_unit, }
}

///|
pub fn pysam3_build_rg_sample_map(
  header_lines : Array[String],
) -> Array[(String, String)] {
  // (RG_ID, sample_name)
  let result : Array[(String, String)] = []
  for i in 0..<header_lines.length() {
    let line = header_lines[i]
    if line.has_prefix("@RG") {
      let rg = pysam3_parse_rg_line(line)
      if rg.id_.length() > 0 {
        result.push((rg.id_, rg.sample))
      }
    }
  }
  result
}

///|
pub fn pysam3_rg_to_sample(
  rg_map : Array[(String, String)],
  rg_id : String,
) -> String {
  for i in 0..<rg_map.length() {
    let (id_, sample) = rg_map[i]
    if id_ == rg_id {
      return sample
    }
  }
  ""
}

// --- Strand-specific counts ---

///|
pub struct Pysam3StrandCounts {
  forward : Int32
  reverse : Int32
}

///|
pub fn pysam3_strand_counts_empty() -> Pysam3StrandCounts {
  { forward: 0, reverse: 0, }
}

///|
pub fn pysam3_strand_counts_add(
  counts : Pysam3StrandCounts,
  is_reverse : Bool,
) -> Pysam3StrandCounts {
  if is_reverse {
    { ..counts, reverse: counts.reverse + 1, }
  } else {
    { ..counts, forward: counts.forward + 1, }
  }
}

///|
pub fn pysam3_strand_bias_ratio(counts : Pysam3StrandCounts) -> Double {
  let total = counts.forward + counts.reverse
  if total == 0 {
    return 0.5
  }
  counts.forward.to_double() / total.to_double()
}

///|
pub fn pysam3_strand_counts_per_position(
  alignments : Array[(Int32, Int32, Bool)],
  start : Int32,
  end_ : Int32,
) -> Array[(Int32, Pysam3StrandCounts)] {
  // (ref_start, read_length, is_reverse)
  let result : Array[(Int32, Pysam3StrandCounts)] = []
  for pos in start..<end_ {
    let mut counts = pysam3_strand_counts_empty()
    for i in 0..<alignments.length() {
      let (ref_start, read_len, is_reverse) = alignments[i]
      if pos >= ref_start && pos < ref_start + read_len {
        counts = pysam3_strand_counts_add(counts, is_reverse)
      }
    }
    result.push((pos, counts))
  }
  result
}
