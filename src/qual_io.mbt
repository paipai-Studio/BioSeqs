// QUAL format parser and writer, mirroring BioPython's
// Bio.SeqIO.QualityIO QUAL support.
//
// QUAL format:
//   >sequence_id optional_description
//   10 20 30 40 10 20 30 40 ...
//
// QUAL files contain PHRED quality scores (integers from 0 to ~90) as
// space-separated text, in a FASTA-like format. The corresponding sequences
// are stored in a separate FASTA file.
//
// In BioSeqs, the QUAL record has an undefined sequence (length = number of
// quality scores) with phred_quality letter annotations.

///|
/// Parse QUAL format text and return an array of SeqRecord objects.
///
/// QUAL format is similar to FASTA, but contains space-separated PHRED
/// quality scores instead of sequence characters. The resulting SeqRecord
/// objects have:
///   - seq: undefined sequence of length = number of quality scores
///   - letter_annotations["phred_quality"]: array of quality scores
///
/// Mirrors BioPython's SeqIO.parse(handle, "qual").
pub fn parse_qual(content : String) -> Array[SeqRecord] {
  let records = Array::new()
  let lines = content.split("\n").to_array()
  let mut current_title : String? = None
  let mut current_quals = Array::new()  // Array[Int]
  let gt : UInt16 = '>'.to_int().to_uint16()
  let cr : UInt16 = '\r'.to_int().to_uint16()
  let sp : UInt16 = ' '.to_int().to_uint16()
  let tab : UInt16 = '\t'.to_int().to_uint16()

  for line in lines {
    let n = line.length()
    let mut trimmed_len = n
    if n > 0 && line.unsafe_get(n - 1) == cr {
      trimmed_len = n - 1
    }
    if trimmed_len == 0 {
      continue
    }
    if line.unsafe_get(0) == gt {
      match current_title {
        Some(title) =>
          records.push(qual_build_record(title, current_quals))
        None => ()
      }
      current_title = Some(
        if trimmed_len > 1 {
          line[1:trimmed_len].to_owned()
        } else {
          ""
        },
      )
      current_quals = Array::new()
    } else {
      // Parse space-separated integer tokens from the line
      let mut num_start = -1
      for i = 0; i < trimmed_len; i = i + 1 {
        let c = line.unsafe_get(i)
        if c == sp || c == tab {
          if num_start >= 0 {
            let token = line[num_start:i].to_owned()
            current_quals.push(parse_int(token))
            num_start = -1
          }
        } else if num_start < 0 {
          num_start = i
        }
      }
      if num_start >= 0 {
        let token = line[num_start:trimmed_len].to_owned()
        current_quals.push(parse_int(token))
      }
    }
  }
  match current_title {
    Some(title) =>
      records.push(qual_build_record(title, current_quals))
    None => ()
  }
  records
}

///|
/// Convert a QUAL header and quality scores array into a SeqRecord.
fn qual_build_record(title : String, quals : Array[Int]) -> SeqRecord {
  let (id, description) = match title.find(" ") {
    Some(pos) => (title[0:pos].to_owned(), title)
    None => (title, title)
  }
  let letter_ann : Map[String, Array[Int]] = Map([], capacity=1)
  letter_ann["phred_quality"] = quals
  SeqRecord::{
    seq: Seq::new(make_n_string(quals.length())),
    id,
    name: id,
    description,
    letter_annotations: letter_ann,
  }
}

///|
/// Write an array of SeqRecord objects to QUAL format text.
///
/// QUAL format is FASTA-like with space-separated PHRED quality scores:
///   >id description
///   10 20 30 40 ...
///
/// Scores are wrapped at 60 characters per line (matching Biopython).
/// Records must have phred_quality letter annotations.
pub fn write_qual(records : Array[SeqRecord]) -> String raise SeqIOError {
  let buf = StringBuilder::new()
  for rec in records {
    let quals = match rec.letter_annotations.get("phred_quality") {
      Some(q) => q
      None =>
        raise SeqIOError(
          "Missing phred_quality annotation for record " + rec.id,
        )
    }
    // Title line: ">" + description (matching BioPython)
    buf.write_char('>')
    buf.write_string(build_fasta_title(rec.id, rec.description))
    buf.write_char('\n')
    // Quality scores: space-separated, wrapped at ~60 chars per line
    // Biopython wraps at 60 characters (not scores), so roughly 30 scores per line
    let line_len = 60
    let mut line_buf = StringBuilder::new()
    let mut line_pos = 0
    for j = 0; j < quals.length(); j = j + 1 {
      let token = quals[j].to_string()
      let token_len = token.length()
      if line_pos > 0 && line_pos + 1 + token_len > line_len {
        buf.write_string(line_buf.to_string())
        buf.write_char('\n')
        line_buf = StringBuilder::new()
        line_pos = 0
      }
      if line_pos > 0 {
        line_buf.write_char(' ')
        line_pos = line_pos + 1
      }
      line_buf.write_string(token)
      line_pos = line_pos + token_len
    }
    if line_pos > 0 {
      buf.write_string(line_buf.to_string())
      buf.write_char('\n')
    }
  }
  buf.to_string()
}

///|
/// Create a string of 'N' characters of the given length.
/// Used for QUAL records where the sequence is undefined.
fn make_n_string(length : Int) -> String {
  if length <= 0 {
    ""
  } else {
    let buf = FixedArray::make(length, 'N'.to_int().to_uint16())
    unsafe_from_uint16_fixedarray(buf)
  }
}

///|
/// Pair a FASTA file (containing sequences) with a QUAL file (containing
/// quality scores) to produce SeqRecord objects with both sequence and
/// quality information.
///
/// Mirrors BioPython's Bio.SeqIO.QualityIO.PairedFastaQualIterator.
///
/// The FASTA and QUAL records are matched by record id (first word of the
/// title line). If a record is present in one file but not the other, it
/// is skipped. Returns records in the order they appear in the FASTA file.
pub fn paired_fasta_qual_iterator(
  fasta_content : String,
  qual_content : String,
) -> Array[SeqRecord] {
  let fasta_records = parse_fasta(fasta_content)
  let qual_records = parse_qual(qual_content)
  // Build a dictionary of QUAL records by id
  let qual_by_id : Map[String, SeqRecord] = Map([], capacity=qual_records.length())
  for qrec in qual_records {
    qual_by_id[qrec.id] = qrec
  }
  // Match and combine
  let result = Array::new()
  for frec in fasta_records {
    match qual_by_id.get(frec.id) {
      Some(qrec) =>
        match qrec.letter_annotations.get("phred_quality") {
          Some(quals) =>
            if frec.seq.length() == quals.length() {
              let ann = frec.letter_annotations
              ann["phred_quality"] = quals
              result.push(SeqRecord::{
                seq: frec.seq,
                id: frec.id,
                name: frec.name,
                description: frec.description,
                letter_annotations: ann,
              })
            }
          None => result.push(frec)
        }
      None => result.push(frec)
    }
  }
  result
}