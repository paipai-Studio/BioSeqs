/// pysam_bam_tag_typed.mbt — Typed tag manipulation for BAM records
/// Inspired by pysam AlignedSegment: get_tag/set_tag with types

///|
/// Tag value types in BAM/SAM.
pub enum PysamTagType {
  Char(String) // A: single character
  Int(Int32) // i: integer
  Float(Double) // f: float
  String(String) // Z: string
  ByteArray(Array[Int32]) // H: hex byte array
  IntArray(Array[Int32]) // B: integer array
  FloatArray(Array[Double]) // B: float array
}

///|
/// Parse a SAM tag field into (tag_name, tag_type, value).
pub fn pysam_parse_tag(tag_field : String) -> (String, String, PysamTagType) {
  // Format: XX:Y:value where XX is 2-char tag name, Y is type
  if tag_field.length() < 5 {
    return ("", "", PysamTagType::String(""))
  }
  let tag_name = tag_field[0:2].to_owned()
  let tag_type = tag_field[3:4].to_owned()
  let value = tag_field[5:].to_owned()
  let typed_value = match tag_type {
    "A" => PysamTagType::Char(value)
    "i" => PysamTagType::Int(parse_int(value))
    "f" => {
      let d = parse_double(value)
      match d {
        Some(v) => PysamTagType::Float(v)
        None => PysamTagType::Float(0.0)
      }
    }
    "Z" => PysamTagType::String(value)
    "H" => PysamTagType::String(value) // hex string
    _ => PysamTagType::String(value)
  }
  (tag_name, tag_type, typed_value)
}

///|
/// Format a typed tag for SAM output.
pub fn pysam_format_tag(tag_name : String, value : PysamTagType) -> String {
  match value {
    PysamTagType::Char(c) => tag_name + ":A:" + c
    PysamTagType::Int(i) => tag_name + ":i:" + i.to_string()
    PysamTagType::Float(f) => tag_name + ":f:" + f.to_string()
    PysamTagType::String(s) => tag_name + ":Z:" + s
    PysamTagType::ByteArray(_) => tag_name + ":H:" // simplified
    PysamTagType::IntArray(arr) => {
      let sb = StringBuilder()
      let _ = sb.write(tag_name + ":B:i")
      for i = 0; i < arr.length(); i = i + 1 {
        let _ = sb.write("," + arr[i].to_string())
      }
      sb.to_string()
    }
    PysamTagType::FloatArray(arr) => {
      let sb = StringBuilder()
      let _ = sb.write(tag_name + ":B:f")
      for i = 0; i < arr.length(); i = i + 1 {
        let _ = sb.write("," + arr[i].to_string())
      }
      sb.to_string()
    }
  }
}

///|
/// Get a tag from SAM fields by name.
pub fn pysam_get_tag(
  fields : Array[StringView],
  tag_name : String,
) -> PysamTagType? {
  for i = 11; i < fields.length(); i = i + 1 {
    let tag = fields[i].to_owned()
    if tag.length() >= 2 && tag[0:2].to_owned() == tag_name {
      let (_, _, value) = pysam_parse_tag(tag)
      return Some(value)
    }
  }
  None
}

///|
/// Set a tag in SAM line.
pub fn pysam_set_tag(
  sam_line : String,
  tag_name : String,
  value : PysamTagType,
) -> String {
  let fields = sam_line.split("\t").to_array()
  let sb = StringBuilder()
  let mut replaced = false
  for i = 0; i < fields.length(); i = i + 1 {
    if i > 0 {
      let _ = sb.write_char('\t')
    }
    let field = fields[i].to_owned()
    if i >= 11 && field.length() >= 2 && field[0:2].to_owned() == tag_name {
      let _ = sb.write(pysam_format_tag(tag_name, value))
      replaced = true
    } else {
      let _ = sb.write(field)
    }
  }
  if !replaced {
    let _ = sb.write_char('\t')
    let _ = sb.write(pysam_format_tag(tag_name, value))
  }
  sb.to_string()
}

///|
/// Remove a tag from SAM line.
pub fn pysam_remove_tag(sam_line : String, tag_name : String) -> String {
  let fields = sam_line.split("\t").to_array()
  let sb = StringBuilder()
  let mut first = true
  for i = 0; i < fields.length(); i = i + 1 {
    let field = fields[i].to_owned()
    if i >= 11 && field.length() >= 2 && field[0:2].to_owned() == tag_name {
      continue // skip this tag
    }
    if !first {
      let _ = sb.write_char('\t')
    }
    let _ = sb.write(field)
    first = false
  }
  sb.to_string()
}

///|
/// List all tags in a SAM line.
pub fn pysam_list_tags(sam_line : String) -> Array[(String, String)] {
  let fields = sam_line.split("\t").to_array()
  let result : Array[(String, String)] = Array(capacity=0)
  for i = 11; i < fields.length(); i = i + 1 {
    let tag = fields[i].to_owned()
    if tag.length() >= 5 {
      result.push((tag[0:2].to_owned(), tag[3:4].to_owned()))
    }
  }
  result
}

///|
/// Suppress unused constructor warnings.
pub fn pysam_bam_unused_suppress_4() -> Unit {
  ignore(PysamTagType::ByteArray([]))
  ignore(PysamTagType::IntArray([]))
  ignore(PysamTagType::FloatArray([]))
}
