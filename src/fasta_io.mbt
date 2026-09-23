// FASTA format parser and writer, mirroring BioPython's
// Bio.SeqIO.FastaIO.
//
// FASTA format:
//   >header_line
//   sequence_line_1
//   sequence_line_2
//   ...
//   >next_header
//   ...
//
// The header line (after '>') is the title. The first whitespace-delimited
// token becomes the record id; the entire title becomes the description.

///|
/// Default line width for wrapping FASTA sequence output (matches BioPython).
let fasta_line_width : Int = 60

///|
/// Parse FASTA format text and return an array of SeqRecord objects.
///
/// Optimized: single pass over the content, minimal allocations.
/// Uses StringBuilder for sequences instead of repeated String operations.
pub fn parse_fasta(content : String) -> Array[SeqRecord] {
  let records = Array(capacity=0)
  let lines = content.split("\n").to_array()
  let mut current_title : String? = None
  let mut current_seq = StringBuilder()
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
          records.push(fasta_title_to_record(title, current_seq.to_string()))
        None => ()
      }
      current_title = Some(
        if trimmed_len > 1 {
          line[1:trimmed_len].to_owned()
        } else {
          ""
        },
      )
      current_seq = StringBuilder()
    } else {
      for i = 0; i < trimmed_len; i = i + 1 {
        let c = line.unsafe_get(i)
        if c != sp && c != tab {
          current_seq.write_char(c.unsafe_to_char())
        }
      }
    }
  }
  match current_title {
    Some(title) =>
      records.push(fasta_title_to_record(title, current_seq.to_string()))
    None => ()
  }
  records
}

///|
/// Build a SeqRecord from a FASTA title and sequence string.
/// The title is split at the first whitespace: first token = id, whole
/// title = description (matching BioPython).
fn fasta_title_to_record(title : String, seq_str : String) -> SeqRecord {
  let (id, description) = match title.find(" ") {
    Some(pos) => (title[0:pos].to_owned(), title)
    None => (title, title)
  }
  SeqRecord::{
    seq: Seq::new(seq_str),
    id,
    name: id,
    description,
    letter_annotations: Map([], capacity=0),
    annotations: Map([], capacity=0),
    features: [],
  }
}

///|
/// Write an array of SeqRecord objects to FASTA format text.
///
/// Mirrors BioPython's FastaIO writer:
/// - The title line is ">" + description (if description starts with id)
///   or ">" + id + " " + description (otherwise).
/// - Sequences are wrapped at 60 characters per line.
pub fn write_fasta(records : Array[SeqRecord]) -> String {
  let buf = StringBuilder()
  for rec in records {
    // Determine the title line, matching BioPython's FastaIO logic:
    //   - If description starts with id: title = description
    //   - If description is non-empty but doesn't start with id:
    //     title = id + " " + description
    //   - If description is empty: title = id
    let title = build_fasta_title(rec.id, rec.description)
    buf.write_char('>')
    buf.write_string(title)
    buf.write_char('\n')
    // Write sequence wrapped at fasta_line_width
    let seq = rec.seq.to_string()
    let n = seq.length()
    let mut i = 0
    while i < n {
      let end = if i + fasta_line_width < n { i + fasta_line_width } else { n }
      buf.write_string(seq[i:end].to_owned())
      buf.write_char('\n')
      i = end
    }
  }
  buf.to_string()
}

///|
/// Build the FASTA title line from id and description.
/// If description starts with id, use description as-is.
/// Otherwise, prepend id + " ".
fn build_fasta_title(id : String, description : String) -> String {
  if description.length() == 0 {
    id
  } else if description.has_prefix(id) {
    description
  } else {
    id + " " + description
  }
}

///|
/// SimpleFastaParser: low-level FASTA parser returning (title, sequence)
/// tuples, mirroring Bio.SeqIO.FastaIO.SimpleFastaParser.
///
/// Empty line handling (mirroring Biopython):
///   - Lines with only whitespace (spaces, tabs) are treated as blank
///   - Blank lines between records are skipped silently
///   - CRLF line endings are supported
///
/// Example:
///   ```moonbit nocheck
///   let records = simple_fasta_parser(">seq1\nATCG\n>seq2\nGCTA\n")
///   for (title, seq) in records {
///     println("Title: \(title), Seq: \(seq)")
///   }
///   ```
pub fn simple_fasta_parser(content : String) -> Array[(String, String)] {
  let records = Array(capacity=0)
  let lines = content.split("\n").to_array()
  let mut current_title = ""
  let mut current_seq = StringBuilder()
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

    // Check if line is empty or contains only whitespace
    let mut is_blank = true
    for i = 0; i < trimmed_len; i = i + 1 {
      let c = line.unsafe_get(i)
      if c != sp && c != tab {
        is_blank = false
        break
      }
    }
    if is_blank {
      continue
    }

    if line.unsafe_get(0) == gt {
      // Save previous record if one is open
      if current_title.length() > 0 {
        records.push((current_title, current_seq.to_string()))
        current_title = ""
        current_seq = StringBuilder()
      }
      current_title = if trimmed_len > 1 {
        line[1:trimmed_len].to_owned()
      } else {
        ""
      }
    } else {
      for i = 0; i < trimmed_len; i = i + 1 {
        let c = line.unsafe_get(i)
        if c != sp && c != tab {
          current_seq.write_char(c.unsafe_to_char())
        }
      }
    }
  }
  if current_title.length() > 0 {
    records.push((current_title, current_seq.to_string()))
  }
  records
}

///|
/// FastqGeneralIterator: low-level FASTQ parser returning
/// (title, sequence, qualities) tuples, mirroring
/// Bio.SeqIO.QualityIO.FastqGeneralIterator.
///
/// Qualities are Phred scores decoded from ASCII-33. Records whose four
/// lines are incomplete at the end of input are skipped. CRLF line
/// endings are supported.
///
/// Example:
///   ```moonbit nocheck
///   let records = fastq_general_parser("@seq1\nACGT\n+\n!''*\n")
///   for (title, seq, quals) in records {
///     println("Title: \(title), Seq: \(seq), Quals: \(quals)")
///   }
///   ```
pub fn fastq_general_parser(
  content : String,
) -> Array[(String, String, Array[Int])] {
  let records = Array(capacity=0)
  let lines = content.split("\n").to_array()
  let n_lines = lines.length()
  let at : UInt16 = '@'.to_int().to_uint16()
  let cr : UInt16 = '\r'.to_int().to_uint16()

  let mut index = 0
  while index + 3 < n_lines {
    let header = lines[index]
    let header_len = header.length()
    let mut h_end = header_len
    if h_end > 0 && header.unsafe_get(h_end - 1) == cr {
      h_end = h_end - 1
    }
    if h_end < 1 || header.unsafe_get(0) != at {
      index = index + 1
      continue
    }
    let title = if h_end > 1 { header[1:h_end].to_owned() } else { "" }
    let seq_line = lines[index + 1]
    let seq_len = seq_line.length()
    let mut s_end = seq_len
    if s_end > 0 && seq_line.unsafe_get(s_end - 1) == cr {
      s_end = s_end - 1
    }
    let seq = seq_line[0:s_end].to_owned()
    let quals_line = lines[index + 3]
    let quals_len = quals_line.length()
    let mut q_end = quals_len
    if q_end > 0 && quals_line.unsafe_get(q_end - 1) == cr {
      q_end = q_end - 1
    }
    let quals : Array[Int] = Array(capacity=0)
    for i = 0; i < q_end; i = i + 1 {
      quals.push(quals_line.unsafe_get(i).to_int() - 33)
    }
    records.push((title, seq, quals))
    index = index + 4
  }
  records
}

///|
/// FastaTwoLineParser: parse two-line ("compact") FASTA where each record
/// is a single header line ">id sequence" with no continuation lines.
///
/// Lines that do not start with '>' are ignored. If a header line contains
/// no space, the whole text after '>' is the title and the sequence is
/// empty.
///
/// Example:
///   ```moonbit nocheck
///   let records = fasta_two_line_parser(">seq1 ACGT\n>seq2 GCTA\n")
///   ```
pub fn fasta_two_line_parser(content : String) -> Array[(String, String)] {
  let records = Array(capacity=0)
  let lines = content.split("\n").to_array()
  let gt : UInt16 = '>'.to_int().to_uint16()
  let cr : UInt16 = '\r'.to_int().to_uint16()
  let sp : UInt16 = ' '.to_int().to_uint16()

  for line in lines {
    let n = line.length()
    let mut trimmed_len = n
    if n > 0 && line.unsafe_get(n - 1) == cr {
      trimmed_len = n - 1
    }
    if trimmed_len < 1 || line.unsafe_get(0) != gt {
      continue
    }
    // Split at the first space after '>'
    let mut space_pos = -1
    for i = 1; i < trimmed_len; i = i + 1 {
      if line.unsafe_get(i) == sp {
        space_pos = i
        break
      }
    }
    if space_pos == -1 {
      records.push((line[1:trimmed_len].to_owned(), ""))
    } else {
      let title = line[1:space_pos].to_owned()
      let seq = line[space_pos + 1:trimmed_len].to_owned()
      records.push((title, seq))
    }
  }
  records
}

///|
/// FastqTwoLineParser: parse two-line ("compact") FASTQ where each record
/// is "@title sequence" followed by a quality line (Phred+33 encoded).
///
/// Lines not starting with '@' are ignored.
///
/// Example:
///   ```moonbit nocheck
///   let records = fastq_two_line_parser("@seq1 ACGT\n!''*\n")
///   ```
pub fn fastq_two_line_parser(
  content : String,
) -> Array[(String, String, Array[Int])] {
  let records = Array(capacity=0)
  let lines = content.split("\n").to_array()
  let at : UInt16 = '@'.to_int().to_uint16()
  let cr : UInt16 = '\r'.to_int().to_uint16()
  let sp : UInt16 = ' '.to_int().to_uint16()

  let mut i = 0
  while i < lines.length() {
    let line0 = lines[i]
    let n0 = line0.length()
    let mut trimmed0 = n0
    if n0 > 0 && line0.unsafe_get(n0 - 1) == cr {
      trimmed0 = n0 - 1
    }
    if trimmed0 < 1 || line0.unsafe_get(0) != at {
      i = i + 1
      continue
    }
    // Split at the first space after '@'
    let mut space_pos = -1
    for j = 1; j < trimmed0; j = j + 1 {
      if line0.unsafe_get(j) == sp {
        space_pos = j
        break
      }
    }
    let title = if space_pos > 0 {
      line0[1:space_pos].to_owned()
    } else {
      line0[1:trimmed0].to_owned()
    }
    let seq = if space_pos > 0 {
      line0[space_pos + 1:trimmed0].to_owned()
    } else {
      ""
    }
    let quals : Array[Int] = Array(capacity=0)
    if i + 1 < lines.length() {
      let line1 = lines[i + 1]
      let n1 = line1.length()
      let mut trimmed1 = n1
      if n1 > 0 && line1.unsafe_get(n1 - 1) == cr {
        trimmed1 = n1 - 1
      }
      for j = 0; j < trimmed1; j = j + 1 {
        quals.push(line1.unsafe_get(j).to_int() - 33)
      }
    }
    records.push((title, seq, quals))
    i = i + 2
  }
  records
}

///|
/// FastaIterator: parse FASTA into SeqRecord objects, mirroring
/// Bio.SeqIO.FastaIO.FastaIterator (the "fasta" format reader).
pub fn fasta_iterator(content : String) -> Array[SeqRecord] {
  parse_fasta(content)
}

///|
/// FastaTwoLineIterator: parse two-line ("compact") FASTA into
/// SeqRecord objects, mirroring Bio.SeqIO.FastaIO.FastaTwoLineIterator.
pub fn fasta_two_line_iterator(content : String) -> Array[SeqRecord] {
  let records = Array(capacity=0)
  for pair in fasta_two_line_parser(content) {
    records.push(fasta_title_to_record(pair.0, pair.1))
  }
  records
}

///|
/// Write a single SeqRecord in FASTA format with the given line width.
///
/// The title line follows BioPython's FastaIO logic:
///   - description starts with id (or is empty): use description/id as-is
///   - otherwise: id + " " + description
pub fn write_fasta_record(
  record : SeqRecord,
  width? : Int = fasta_line_width,
) -> String {
  let buf = StringBuilder()
  buf.write_char('>')
  buf.write_string(build_fasta_title(record.id, record.description))
  buf.write_char('\n')
  let seq = record.seq.to_string()
  let n = seq.length()
  let mut i = 0
  while i < n {
    let end = if i + width < n { i + width } else { n }
    buf.write_string(seq[i:end].to_owned())
    buf.write_char('\n')
    i = end
  }
  buf.to_string()
}

///|
/// Write an array of SeqRecord objects to FASTA text with the given
/// line width (see `write_fasta` for the default-width variant).
pub fn write_fasta_records(
  records : Array[SeqRecord],
  width? : Int = fasta_line_width,
) -> String {
  let buf = StringBuilder()
  for rec in records {
    buf.write_string(write_fasta_record(rec, width~))
  }
  buf.to_string()
}
