// Tab (tab-separated) format parser and writer, mirroring BioPython's
// Bio.SeqIO.TabIO.
//
// Tab format is the simplest plain-text sequence format:
//
//   ID123456\tCATCNAGATGACACTACGACTACGACTCAGACTAC
//   ID123457\tACACTACGACTACGACTCAGACTACAAN
//
// Each line has exactly one tab, separating:
//   - field 1: identifier (becomes record.id and record.name)
//   - field 2: sequence (whitespace ignored within sequence)
//
// Blank lines are ignored. Description is not stored (set to record.id).

///|
/// Parse tab-delimited sequence format and return an array of SeqRecord
/// objects.
///
/// Mirrors BioPython's SeqIO.parse(handle, "tab").
/// Raises SeqIOError if a line has anything other than exactly one tab.
pub fn parse_tab(content : String) -> Array[SeqRecord] raise SeqIOError {
  let records = Array::new()
  let lines = content.split("\n").to_array()
  let cr : UInt16 = '\r'.to_int().to_uint16()
  let tab : UInt16 = '\t'.to_int().to_uint16()
  let sp : UInt16 = ' '.to_int().to_uint16()

  for line in lines {
    let n = line.length()
    let mut trimmed = n
    if n > 0 && line.unsafe_get(n - 1) == cr {
      trimmed = n - 1
    }
    // Skip blank lines
    let mut is_blank = true
    for j = 0; j < trimmed; j = j + 1 {
      if line.unsafe_get(j) != sp {
        is_blank = false
        break
      }
    }
    if is_blank {
      continue
    }

    // Count tabs
    let mut tab_count = 0
    let mut tab_pos = -1
    for j = 0; j < trimmed; j = j + 1 {
      if line.unsafe_get(j) == tab {
        tab_count = tab_count + 1
        tab_pos = j
      }
    }

    if tab_count != 1 {
      raise SeqIOError(
        "Each line should have exactly one tab separating title and sequence, this line has " +
        tab_count.to_string() +
        " tabs: " +
        line[0:trimmed].to_owned(),
      )
    }

    let id = line[0:tab_pos].to_owned()
    let seq_part = line[tab_pos + 1:trimmed].to_owned()
    // Strip whitespace from the sequence
    let seq_buf = StringBuilder::new()
    let sp_part = seq_part.length()
    for j = 0; j < sp_part; j = j + 1 {
      let c = seq_part.unsafe_get(j)
      if c != sp && c != tab {
        seq_buf.write_char(c.unsafe_to_char())
      }
    }

    records.push(SeqRecord::new(
      Seq::new(seq_buf.to_string()),
      id=id,
      name=id,
      description=id,
    ))
  }
  records
}

///|
/// Write an array of SeqRecord objects to tab-delimited format.
///
/// Mirrors BioPython's TabIO writer:
///   - Each line: record.id + "\t" + sequence
///   - Description/name and other annotations are NOT preserved
pub fn write_tab(records : Array[SeqRecord]) -> String {
  let buf = StringBuilder::new()
  for rec in records {
    buf.write_string(rec.id)
    buf.write_char('\t')
    buf.write_string(rec.seq.to_string())
    buf.write_char('\n')
  }
  buf.to_string()
}
