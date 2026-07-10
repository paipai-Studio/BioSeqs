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
  let records = Array::new()
  let lines = content.split("\n").to_array()
  let mut current_title : String? = None
  let mut current_seq = StringBuilder::new()
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
      current_seq = StringBuilder::new()
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
  let buf = StringBuilder::new()
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
