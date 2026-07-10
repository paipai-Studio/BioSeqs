// FASTQ format parser and writer, mirroring BioPython's
// Bio.SeqIO.QualityIO.
//
// FASTQ format (Sanger, Phred+33):
//   @sequence_id optional_description
//   ACGTACGT...
//   +optional_repeat_of_header
//   !"#$%&'...
//
// Line 1: @ + title (first word = id, whole title = description)
// Line 2: raw sequence
// Line 3: + (optionally followed by the title again)
// Line 4: quality string (same length as sequence)
//
// Quality scores: Q = ord(char) - 33 (Phred+33 / Sanger encoding)

///|
/// Parse FASTQ format text (Sanger/Phred+33) and return an array of
/// SeqRecord objects with phred_quality letter annotations.
///
/// Mirrors BioPython's SeqIO.parse(handle, "fastq").
pub fn parse_fastq(content : String) -> Array[SeqRecord] raise SeqIOError {
  let records = Array::new()
  let lines = content.split("\n").to_array()
  let n = lines.length()
  let mut i = 0
  let at : UInt16 = '@'.to_int().to_uint16()
  let plus : UInt16 = '+'.to_int().to_uint16()
  let cr : UInt16 = '\r'.to_int().to_uint16()

  while i + 3 < n {
    let line0 = lines[i]
    let line1 = lines[i + 1]
    let line2 = lines[i + 2]
    let line3 = lines[i + 3]

    let len0 = line0.length()
    let mut trimmed0 = len0
    if len0 > 0 && line0.unsafe_get(len0 - 1) == cr {
      trimmed0 = len0 - 1
    }
    if trimmed0 == 0 && i + 4 >= n {
      break
    }
    if trimmed0 == 0 || line0.unsafe_get(0) != at {
      let line0_str = line0.to_owned()
      raise SeqIOError(
        "Expected '@' at start of FASTQ header, got: " + line0_str,
      )
    }

    let len2 = line2.length()
    let mut trimmed2 = len2
    if len2 > 0 && line2.unsafe_get(len2 - 1) == cr {
      trimmed2 = len2 - 1
    }
    if trimmed2 == 0 || line2.unsafe_get(0) != plus {
      let line2_str = line2.to_owned()
      raise SeqIOError(
        "Expected '+' at start of FASTQ separator, got: " + line2_str,
      )
    }

    let len1 = line1.length()
    let mut trimmed1 = len1
    if len1 > 0 && line1.unsafe_get(len1 - 1) == cr {
      trimmed1 = len1 - 1
    }

    let len3 = line3.length()
    let mut trimmed3 = len3
    if len3 > 0 && line3.unsafe_get(len3 - 1) == cr {
      trimmed3 = len3 - 1
    }

    if trimmed1 != trimmed3 {
      raise SeqIOError(
        "Sequence length " +
        trimmed1.to_string() +
        " and quality length " +
        trimmed3.to_string() +
        " don't match in FASTQ record",
      )
    }

    let title = if trimmed0 > 1 { line0[1:trimmed0].to_owned() } else { "" }
    let (id, description) = match title.find(" ") {
      Some(pos) => (title[0:pos].to_owned(), title)
      None => (title, title)
    }

    let seq_str = if trimmed1 > 0 { line1[0:trimmed1].to_owned() } else { "" }

    let quals : FixedArray[Int] = FixedArray::make(trimmed3, 0)
    for j = 0; j < trimmed3; j = j + 1 {
      quals[j] = line3.unsafe_get(j).to_int() - 33
    }
    let quals_arr = Array::new()
    for j = 0; j < trimmed3; j = j + 1 {
      quals_arr.push(quals[j])
    }

    let letter_ann : Map[String, Array[Int]] = Map([], capacity=1)
    letter_ann["phred_quality"] = quals_arr
    records.push(SeqRecord::{
      seq: Seq::new(seq_str),
      id,
      name: id,
      description,
      letter_annotations: letter_ann,
    })
    i = i + 4
  }
  records
}

///|
/// Write an array of SeqRecord objects to FASTQ format text (Sanger/Phred+33).
///
/// Mirrors BioPython's SeqIO.write(records, handle, "fastq").
/// - Records must have phred_quality letter annotations.
/// - Quality scores are encoded as chr(Q + 33).
pub fn write_fastq(records : Array[SeqRecord]) -> String raise SeqIOError {
  let buf = StringBuilder::new()
  for rec in records {
    let quals = match rec.letter_annotations.get("phred_quality") {
      Some(q) => q
      None =>
        raise SeqIOError(
          "Missing phred_quality annotation for record " + rec.id,
        )
    }
    let seq = rec.seq.to_string()
    if seq.length() != quals.length() {
      raise SeqIOError(
        "Sequence length " +
        seq.length().to_string() +
        " and quality length " +
        quals.length().to_string() +
        " don't match in record " +
        rec.id,
      )
    }
    // Title line: "@" + description (matching BioPython)
    buf.write_char('@')
    buf.write_string(build_fasta_title(rec.id, rec.description))
    buf.write_char('\n')
    // Sequence line
    buf.write_string(seq)
    buf.write_char('\n')
    // Separator line
    buf.write_string("+\n")
    // Quality line: encode each Q as chr(Q + 33)
    let qbuf : FixedArray[UInt16] = FixedArray::make(quals.length(), 0)
    for j = 0; j < quals.length(); j = j + 1 {
      qbuf[j] = (quals[j] + 33).to_uint16()
    }
    buf.write_string(unsafe_from_uint16_fixedarray(qbuf))
    buf.write_char('\n')
  }
  buf.to_string()
}
