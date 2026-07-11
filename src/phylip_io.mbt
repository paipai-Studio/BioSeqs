// PHYLIP format parser and writer, mirroring BioPython's
// Bio.AlignIO.PhylipIO.
//
// PHYLIP interlaced format:
//   4 384
//   CYS1_DICDI -----MKVILLFVLAVFTVFVSS---------------RGIPPE
//   ALEU_HORVU MAHARVLLLALAVLATAAVAVASSSSFADSNPIRPVTDRAASTL
//   CATH_HUMAN ------MWATLPLLCAGAWLLGV--------PVCGAAELSVNSL
//   CYS1_DICDI ...I--
//   ALEU_HORVU ...VAA
//   CATH_HUMAN ...PLV
//
// PHYLIP sequential format:
//   4 384
//   CYS1_DICDI -----MKVILLFVLAVFTVFVSS---------------RGIPPE...I--
//   ALEU_HORVU MAHARVLLLALAVLATAAVAVASSSSFADSNPIRPVTDRAASTL...VAA
//   CATH_HUMAN ------MWATLPLLCAGAWLLGV--------PVCGAAELSVNSL...PLV

///|
/// Parse an integer from a string.

///|
/// Parse PHYLIP format text (interlaced) and return a MultipleSeqAlignment.
pub fn parse_phylip(
  content : String,
) -> MultipleSeqAlignment raise AlignIOError {
  parse_phylip_generic(content, false)
}

///|
/// Parse PHYLIP sequential format text and return a MultipleSeqAlignment.
pub fn parse_phylip_sequential(
  content : String,
) -> MultipleSeqAlignment raise AlignIOError {
  parse_phylip_generic(content, true)
}

///|
/// Generic PHYLIP parser that handles both interlaced and sequential formats.
fn parse_phylip_generic(
  content : String,
  sequential : Bool,
) -> MultipleSeqAlignment raise AlignIOError {
  let lines = content.split("\n").to_array()
  let n = lines.length()
  let mut i = 0
  let cr : UInt16 = '\r'.to_int().to_uint16()
  let sp : UInt16 = ' '.to_int().to_uint16()

  while i < n {
    let line = lines[i]
    let len = line.length()
    let mut trimmed = len
    if len > 0 && line.unsafe_get(len - 1) == cr {
      trimmed = len - 1
    }
    if trimmed > 0 {
      break
    }
    i = i + 1
  }

  if i >= n {
    raise AlignIOError("Empty PHYLIP file")
  }

  let header = lines[i]
  let hlen = header.length()
  let mut htrimmed = hlen
  if hlen > 0 && header.unsafe_get(hlen - 1) == cr {
    htrimmed = hlen - 1
  }

  let mut space_pos = 0
  while space_pos < htrimmed && header.unsafe_get(space_pos) != sp {
    space_pos = space_pos + 1
  }
  let num_seqs_str = header[0:space_pos].to_owned()
  let num_seqs = parse_int(num_seqs_str)

  let mut col_pos = space_pos
  while col_pos < htrimmed && header.unsafe_get(col_pos) == sp {
    col_pos = col_pos + 1
  }
  let num_cols_str = header[col_pos:htrimmed].to_owned()
  let _num_cols = parse_int(num_cols_str)

  i = i + 1
  let records : Map[String, StringBuilder] = Map([], capacity=num_seqs)
  let current_names = Array::new()

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
    let mut j = 0
    while j < trimmed && line.unsafe_get(j) != sp {
      j = j + 1
    }
    if j == 0 {
      i = i + 1
      continue
    }
    let name = if j > 10 { line[0:10].to_owned() } else { line[0:j].to_owned() }
    let mut k = j
    while k < trimmed && line.unsafe_get(k) == sp {
      k = k + 1
    }
    let seq_part = line[k:trimmed].replace(old=" ", new="").to_owned()

    match records.get(name) {
      Some(sb) => sb.write_string(seq_part)
      None => {
        let sb = StringBuilder::new()
        sb.write_string(seq_part)
        records[name] = sb
        current_names.push(name)
      }
    }
    if sequential && current_names.length() >= num_seqs {
      break
    }
    i = i + 1
  }

  let seq_records = Array::new()
  for name in current_names {
    match records.get(name) {
      Some(sb) => {
        let seq = Seq::new(sb.to_string())
        seq_records.push(SeqRecord::{
          seq,
          id: name,
          name,
          description: name,
          letter_annotations: Map([], capacity=0),
        })
      }
      None => ()
    }
  }
  MultipleSeqAlignment::new(seq_records)
}

///|
/// Write a MultipleSeqAlignment to PHYLIP interlaced format text.
pub fn write_phylip(align : MultipleSeqAlignment) -> String {
  let buf = StringBuilder::new()
  let records = align.iter()
  let n_seqs = records.length()
  let n_cols = align.get_alignment_length()

  buf.write_char(' ')
  buf.write_string(n_seqs.to_string())
  buf.write_char(' ')
  buf.write_string(n_cols.to_string())
  buf.write_char('\n')

  let line_width = 60
  let name_width = 10
  let group_size = 10

  let mut pos = 0
  while pos < n_cols {
    let end = if pos + line_width < n_cols { pos + line_width } else { n_cols }
    for rec in records {
      let name = if rec.id.length() >= name_width {
        rec.id[0:name_width].to_owned()
      } else {
        rec.id + " ".repeat(name_width - rec.id.length())
      }
      buf.write_string(name)
      buf.write_char(' ')
      let seq_str = rec.seq.to_string()
      let mut i = pos
      while i < end {
        let gend = if i + group_size < end { i + group_size } else { end }
        buf.write_string(seq_str[i:gend].to_owned())
        if gend < end {
          buf.write_char(' ')
        }
        i = gend
      }
      buf.write_char(' ')
      buf.write_char('\n')
    }
    if pos + line_width < n_cols {
      buf.write_char('\n')
    }
    pos = end
  }
  buf.to_string()
}

///|
/// Write a MultipleSeqAlignment to PHYLIP sequential format text.
pub fn write_phylip_sequential(align : MultipleSeqAlignment) -> String {
  let buf = StringBuilder::new()
  let records = align.iter()
  let n_seqs = records.length()
  let n_cols = align.get_alignment_length()

  buf.write_string(n_seqs.to_string())
  buf.write_char(' ')
  buf.write_string(n_cols.to_string())
  buf.write_char('\n')

  let name_width = 10
  for rec in records {
    let name = if rec.id.length() >= name_width {
      rec.id[0:name_width].to_owned()
    } else {
      rec.id + " ".repeat(name_width - rec.id.length())
    }
    buf.write_string(name)
    buf.write_char(' ')
    buf.write_string(rec.seq.to_string())
    buf.write_char('\n')
  }
  buf.write_char('\n')
  buf.to_string()
}
