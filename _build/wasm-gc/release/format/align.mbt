// Bio.Align.MultipleSeqAlignment implementation for MoonBit.
//
// MultipleSeqAlignment represents a collection of sequences (usually shown as
// rows) which are all the same length (usually with gap characters for
// insertions or padding). The data can be regarded as a matrix of letters,
// with well defined columns.
//
// Mirrors BioPython's Bio.Align.MultipleSeqAlignment class.

///|
/// Error raised by AlignIO operations.
pub suberror AlignIOError {
  AlignIOError(String)
}

///|
/// Represents a multiple sequence alignment (MSA).
///
/// All sequences in the alignment must have the same length. Can be accessed
/// like an array of SeqRecords.
pub struct MultipleSeqAlignment {
  records : Array[SeqRecord]
  annotations : Map[String, String]
  column_annotations : Map[String, Array[String]]
} derive(Eq)

///|
/// Create a new MultipleSeqAlignment from an array of SeqRecords.
/// All sequences must have the same length.
pub fn MultipleSeqAlignment::new(
  records : Array[SeqRecord],
  annotations? : Map[String, String] = Map([], capacity=0),
  column_annotations? : Map[String, Array[String]] = Map([], capacity=0),
) -> MultipleSeqAlignment raise AlignIOError {
  if records.length() > 0 {
    let first_len = records[0].length()
    for rec in records {
      if rec.length() != first_len {
        raise AlignIOError(
          "All sequences must have the same length in an alignment",
        )
      }
    }
  }
  MultipleSeqAlignment::{ records, annotations, column_annotations }
}

///|
/// Return the number of sequences (rows) in the alignment.
pub fn MultipleSeqAlignment::num_records(self : MultipleSeqAlignment) -> Int {
  self.records.length()
}

///|
/// Return the number of columns in the alignment.
pub fn MultipleSeqAlignment::get_alignment_length(
  self : MultipleSeqAlignment,
) -> Int {
  if self.records.length() == 0 {
    0
  } else {
    self.records[0].length()
  }
}

///|
/// Access the SeqRecord at index `i`.
pub fn MultipleSeqAlignment::get(
  self : MultipleSeqAlignment,
  i : Int,
) -> SeqRecord? {
  if i >= 0 && i < self.records.length() {
    Some(self.records[i])
  } else {
    None
  }
}

///|
/// Get the column at position `col` as a String.
pub fn MultipleSeqAlignment::column(
  self : MultipleSeqAlignment,
  col : Int,
) -> String? {
  if self.records.length() == 0 {
    return None
  }
  if col < 0 || col >= self.get_alignment_length() {
    return None
  }
  let buf = StringBuilder::new(size_hint=self.records.length())
  for rec in self.records {
    match rec.seq.get_char(col) {
      Some(c) => buf.write_char(c)
      None => buf.write_char('?')
    }
  }
  Some(buf.to_string())
}

///|
/// Slice the alignment: align[start_row:end_row, start_col:end_col]
/// Returns a new MultipleSeqAlignment with the specified rows and columns.
pub fn MultipleSeqAlignment::slice(
  self : MultipleSeqAlignment,
  row_start? : Int = 0,
  row_end? : Int = -1,
  col_start? : Int = 0,
  col_end? : Int = -1,
) -> MultipleSeqAlignment {
  let n_rows = self.records.length()
  let n_cols = self.get_alignment_length()

  let rs = if row_start < 0 {
    0
  } else if row_start >= n_rows {
    n_rows
  } else {
    row_start
  }
  let re = if row_end < 0 {
    n_rows
  } else if row_end > n_rows {
    n_rows
  } else {
    row_end
  }

  let cs = if col_start < 0 {
    0
  } else if col_start >= n_cols {
    n_cols
  } else {
    col_start
  }
  let ce = if col_end < 0 {
    n_cols
  } else if col_end > n_cols {
    n_cols
  } else {
    col_end
  }

  let new_records = Array::new()
  for i = rs; i < re; i = i + 1 {
    let rec = self.records[i]
    let sliced_seq = rec.seq.slice(start=cs, end=ce)
    let new_rec = SeqRecord::{ ..rec, seq: sliced_seq }
    new_records.push(new_rec)
  }

  MultipleSeqAlignment::{
    records: new_records,
    annotations: Map([], capacity=0),
    column_annotations: Map([], capacity=0),
  }
}

///|
/// Add a sequence to the alignment.
pub fn MultipleSeqAlignment::add_sequence(
  self : MultipleSeqAlignment,
  seq : SeqRecord,
) -> MultipleSeqAlignment raise AlignIOError {
  let align_len = self.get_alignment_length()
  if align_len > 0 && seq.length() != align_len {
    raise AlignIOError("Sequence length does not match alignment length")
  }
  let new_records = Array::new()
  for r in self.records {
    new_records.push(r)
  }
  new_records.push(seq)
  MultipleSeqAlignment::{ ..self, records: new_records }
}

///|
/// Convert the alignment to a string representation (similar to BioPython's __str__).
pub fn MultipleSeqAlignment::repr(self : MultipleSeqAlignment) -> String {
  let buf = StringBuilder::new()
  let n_records = self.records.length()
  let n_cols = self.get_alignment_length()
  buf.write_string("Alignment with ")
  buf.write_string(n_records.to_string())
  buf.write_string(" rows and ")
  buf.write_string(n_cols.to_string())
  buf.write_string(" columns")
  if n_records > 0 && n_cols > 0 {
    buf.write_char('\n')
    for rec in self.records {
      let seq_str = rec.seq.to_string()
      let display_seq = if seq_str.length() > 50 {
        seq_str[0:50].to_owned() + "..."
      } else {
        seq_str
      }
      buf.write_string(display_seq)
      buf.write_char(' ')
      buf.write_string(rec.id)
      buf.write_char('\n')
    }
  }
  buf.to_string()
}

///|
/// Return all SeqRecords in the alignment.
pub fn MultipleSeqAlignment::iter(
  self : MultipleSeqAlignment,
) -> Array[SeqRecord] {
  self.records
}
