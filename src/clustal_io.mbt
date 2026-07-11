// ClustalW format parser and writer, mirroring BioPython's
// Bio.AlignIO.ClustalIO.
//
// ClustalW format:
//   CLUSTAL X (1.81) multiple sequence alignment
//   
//   AF191659.1      TATACATTAAAGAAGGGGGATGCGGATAAATGGAAAGGCGAAAG
//   AF191658.1      TATACATTAAAGAAGGGGGATGCGGATAAATGGAAAGGCGAAAG
//   AF191661.1      TATACATTAAAGAAGGGGGATGCGGATAAATGGAAAGGCGAAAG
//   AF191660.1      TATACATAAAAGAAGGGGGATGCGGATAAATGGAAAGGCGAAAG
//   
//   AF191659.1      AGA
//   AF191658.1      AGA
//   AF191661.1      AGA
//   AF191660.1      AGA

///|
/// Parse ClustalW format text and return a MultipleSeqAlignment.
pub fn parse_clustal(
  content : String,
) -> MultipleSeqAlignment raise AlignIOError {
  let lines = content.split("\n").to_array()
  let n = lines.length()
  let mut i = 0
  let records : Map[String, StringBuilder] = Map([], capacity=16)
  let gt : UInt16 = '>'.to_int().to_uint16()
  let cr : UInt16 = '\r'.to_int().to_uint16()
  let sp : UInt16 = ' '.to_int().to_uint16()

  while i < n {
    let line = lines[i]
    let len = line.length()
    let mut trimmed = len
    if len > 0 && line.unsafe_get(len - 1) == cr {
      trimmed = len - 1
    }
    if trimmed == 0 {
      i = i + 1
      continue
    }
    let first = line.unsafe_get(0)
    if first == gt {
      break
    }
    if first == '*'.to_int().to_uint16() ||
      first == ':'.to_int().to_uint16() ||
      first == '.'.to_int().to_uint16() ||
      first == ' '.to_int().to_uint16() {
      i = i + 1
      continue
    }
    if trimmed >= 7 &&
      line.unsafe_get(0) == 'C'.to_int().to_uint16() &&
      line.unsafe_get(1) == 'L'.to_int().to_uint16() &&
      line.unsafe_get(2) == 'U'.to_int().to_uint16() &&
      line.unsafe_get(3) == 'S'.to_int().to_uint16() &&
      line.unsafe_get(4) == 'T'.to_int().to_uint16() &&
      line.unsafe_get(5) == 'A'.to_int().to_uint16() &&
      line.unsafe_get(6) == 'L'.to_int().to_uint16() {
      i = i + 1
      continue
    }
    let mut j = 0
    while j < trimmed && line.unsafe_get(j) != sp {
      j = j + 1
    }
    if j == 0 {
      i = i + 1
      continue
    }
    let name = line[0:j].to_owned()
    let mut k = j
    while k < trimmed && line.unsafe_get(k) == sp {
      k = k + 1
    }
    let seq_part = line[k:trimmed].to_owned()
    match records.get(name) {
      Some(sb) => sb.write_string(seq_part)
      None => {
        let sb = StringBuilder::new()
        sb.write_string(seq_part)
        records[name] = sb
      }
    }
    i = i + 1
  }

  let seq_records = Array::new()
  for name, sb in records {
    let seq = Seq::new(sb.to_string())
    seq_records.push(SeqRecord::{
      seq,
      id: name,
      name,
      description: name,
      letter_annotations: Map([], capacity=0),
    })
  }
  MultipleSeqAlignment::new(seq_records)
}

///|
/// Write a MultipleSeqAlignment to ClustalW format text.
pub fn write_clustal(align : MultipleSeqAlignment) -> String {
  let buf = StringBuilder::new()
  buf.write_string("CLUSTAL X (1.81) multiple sequence alignment\n\n\n")

  let records = align.iter()
  let n_cols = align.get_alignment_length()
  let line_width = 60
  let id_padding = 36

  let mut pos = 0
  while pos < n_cols {
    let end = if pos + line_width < n_cols { pos + line_width } else { n_cols }
    for rec in records {
      buf.write_string(rec.id)
      let padding = id_padding - rec.id.length()
      if padding > 0 {
        buf.write_string(" ".repeat(padding))
      } else {
        buf.write_char(' ')
      }
      buf.write_string(rec.seq.to_string()[pos:end].to_owned())
      buf.write_char('\n')
    }
    buf.write_char('\n')
    pos = end
  }
  buf.write_char('\n')
  buf.to_string()
}
