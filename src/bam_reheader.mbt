/// Replace a BAM file's header without rewriting all records.
///
/// pysam/samtools equivalent: `samtools reheader`

///|
/// Replace the header of a BamFile with a new BamHeader.
/// Does NOT validate reference compatibility — use bam_reheader_safe for that.
pub fn bam_reheader(bam : BamFile, new_header : BamHeader) -> BamFile {
  { header: new_header, records: bam.records, }
}

///|
/// Replace header with validation that reference names and lengths match.
pub fn bam_reheader_safe(
  bam : BamFile,
  new_header : BamHeader,
) -> BamFile raise BamError {
  let old_names = bam.get_reference_names()
  let new_names : Array[String] = Array(capacity=0)
  for ref_ in new_header.references {
    new_names.push(ref_.name)
  }
  if old_names.length() != new_names.length() {
    raise BamError("Reference count mismatch")
  }
  for i = 0; i < old_names.length(); i = i + 1 {
    if old_names[i] != new_names[i] {
      raise BamError("Reference name mismatch")
    }
    match bam.get_reference_length(old_names[i]) {
      Some(old_len) =>
        if new_header.references[i].length != old_len {
          raise BamError("Reference length mismatch")
        }
      None => ()
    }
  }
  { header: new_header, records: bam.records, }
}

///|
/// Replace just the header text (HD/PG/CO/RG lines) while keeping references.
pub fn bam_reheader_text(bam : BamFile, new_text : String) -> BamFile {
  let new_header = { text: new_text, references: bam.header.references, }
  { header: new_header, records: bam.records, }
}

///|
/// Parse a SAM header text into a BamHeader, extracting @SQ lines.
pub fn bam_parse_sam_header(header_text : String) -> BamHeader {
  let refs : Array[BamReference] = Array(capacity=0)
  let lines = header_text.split("\n").to_array()
  for line in lines {
    if line.length() >= 3 && line[0:3].to_owned() == "@SQ" {
      let mut name = ""
      let mut length : Int32 = 0
      let fields = line.split("\t").to_array()
      for field in fields {
        if field.length() >= 3 &&
          field[0:2].to_owned() == "SN" &&
          field[2:3].to_owned() == ":" {
          name = field[3:].to_owned()
        } else if field.length() >= 3 &&
          field[0:2].to_owned() == "LN" &&
          field[2:3].to_owned() == ":" {
          let len_str = field[3:].to_owned()
          length = parse_int(len_str)
        }
      }
      if name != "" {
        refs.push({ name, length, })
      }
    }
  }
  { text: header_text, references: FixedArray::from_array(refs), }
}
