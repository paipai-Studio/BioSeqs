// GenBank format parser, mirroring BioPython's Bio.SeqIO.GenBankIO.
//
// GenBank is a complex multi-section format:
//
//   LOCUS       NM_001     5028 bp    DNA    linear   PRI 21-JUN-1999
//   DEFINITION  Description text here.
//   ACCESSION   NM_001
//   VERSION     NM_001.1  GI:12345
//   SOURCE      Organism name
//   FEATURES    ...
//   ORIGIN
//           1 atgcatgcat gcatgcatgc atgcatgcat gcatgcatgc atgcatgcat gcatgcatgc
//          61 ...
//   //
//
// Key fields extracted:
//   - id:          VERSION (first token, e.g. "NM_001.1")
//   - name:        LOCUS name (first token after LOCUS)
//   - description: DEFINITION text (joined, trailing '.' stripped)
//   - seq:         ORIGIN section (numbers and spaces stripped)

///|
/// Parse GenBank format text and return an array of SeqRecord objects.
///
/// Mirrors BioPython's SeqIO.parse(handle, "genbank").
/// Extracts: id (VERSION), name (LOCUS), description (DEFINITION),
/// and sequence (ORIGIN section). Multiple records separated by '//'.

///|
/// Precomputed uppercase conversion table for ASCII letters.
/// 0 means "pass through unchanged".
let uppercase_table : FixedArray[UInt16] = {
  let t : FixedArray[UInt16] = FixedArray::make(128, 0)
  for c = 'a'.to_int(); c <= 'z'.to_int(); c = c + 1 {
    t[c] = (c - 32).to_uint16()
  }
  t
}

///|
pub fn parse_genbank(content : String) -> Array[SeqRecord] {
  let records = Array::new()
  let lines = content.split("\n").to_array()
  let n = lines.length()
  let mut i = 0
  let cr : UInt16 = '\r'.to_int().to_uint16()
  let sp : UInt16 = ' '.to_int().to_uint16()
  let dot : UInt16 = '.'.to_int().to_uint16()

  while i < n {
    let line = lines[i]
    let len = line.length()
    let mut trimmed = len
    if len > 0 && line.unsafe_get(len - 1) == cr {
      trimmed = len - 1
    }
    if trimmed >= 5 &&
      line.unsafe_get(0) == 'L'.to_int().to_uint16() &&
      line.unsafe_get(1) == 'O'.to_int().to_uint16() &&
      line.unsafe_get(2) == 'C'.to_int().to_uint16() &&
      line.unsafe_get(3) == 'U'.to_int().to_uint16() &&
      line.unsafe_get(4) == 'S'.to_int().to_uint16() {
      let mut name = "<unknown>"
      let mut id = "<unknown id>"
      let mut description = "<unknown description>"
      let mut seq_buf : FixedArray[UInt16] = FixedArray::make(0, 0)
      let mut seq_len = 0
      let mut in_origin = false

      let locus_rest = line[5:trimmed]
      let mut j = 0
      let rest_len = locus_rest.length()
      while j < rest_len && locus_rest.unsafe_get(j) == sp {
        j = j + 1
      }
      let mut name_end = j
      while name_end < rest_len && locus_rest.unsafe_get(name_end) != sp {
        name_end = name_end + 1
      }
      name = locus_rest[j:name_end].to_owned()

      i = i + 1
      while i < n {
        let l = lines[i]
        let ll = l.length()
        let mut lt = ll
        if ll > 0 && l.unsafe_get(ll - 1) == cr {
          lt = ll - 1
        }
        if lt >= 2 &&
          l.unsafe_get(0) == '/'.to_int().to_uint16() &&
          l.unsafe_get(1) == '/'.to_int().to_uint16() {
          break
        }
        if lt >= 10 &&
          l.unsafe_get(0) == 'D'.to_int().to_uint16() &&
          l.unsafe_get(1) == 'E'.to_int().to_uint16() &&
          l.unsafe_get(2) == 'F'.to_int().to_uint16() &&
          l.unsafe_get(3) == 'I'.to_int().to_uint16() &&
          l.unsafe_get(4) == 'N'.to_int().to_uint16() &&
          l.unsafe_get(5) == 'I'.to_int().to_uint16() &&
          l.unsafe_get(6) == 'T'.to_int().to_uint16() &&
          l.unsafe_get(7) == 'I'.to_int().to_uint16() &&
          l.unsafe_get(8) == 'O'.to_int().to_uint16() &&
          l.unsafe_get(9) == 'N'.to_int().to_uint16() {
          let def_part = l[10:lt]
          let mut dj = 0
          let dpl = def_part.length()
          while dj < dpl && def_part.unsafe_get(dj) == sp {
            dj = dj + 1
          }
          let def_buf = StringBuilder::new()
          def_buf.write_string(def_part[dj:dpl].to_owned())
          i = i + 1
          while i < n {
            let cont = lines[i]
            let cl = cont.length()
            let mut ct = cl
            if cl > 0 && cont.unsafe_get(cl - 1) == cr {
              ct = cl - 1
            }
            if ct == 0 || cont.unsafe_get(0) != sp {
              break
            }
            def_buf.write_char(' ')
            let mut cj = 0
            while cj < ct && cont.unsafe_get(cj) == sp {
              cj = cj + 1
            }
            def_buf.write_string(cont[cj:ct].to_owned())
            i = i + 1
          }
          description = def_buf.to_string()
          if description.length() > 0 &&
            description.unsafe_get(description.length() - 1) == dot {
            description = description[0:description.length() - 1].to_owned()
          }
          continue
        }
        if lt >= 7 &&
          l.unsafe_get(0) == 'V'.to_int().to_uint16() &&
          l.unsafe_get(1) == 'E'.to_int().to_uint16() &&
          l.unsafe_get(2) == 'R'.to_int().to_uint16() &&
          l.unsafe_get(3) == 'S'.to_int().to_uint16() &&
          l.unsafe_get(4) == 'I'.to_int().to_uint16() &&
          l.unsafe_get(5) == 'O'.to_int().to_uint16() &&
          l.unsafe_get(6) == 'N'.to_int().to_uint16() {
          let ver_part = l[7:lt]
          let mut vj = 0
          let vpl = ver_part.length()
          while vj < vpl && ver_part.unsafe_get(vj) == sp {
            vj = vj + 1
          }
          let mut vid_end = vj
          while vid_end < vpl && ver_part.unsafe_get(vid_end) != sp {
            vid_end = vid_end + 1
          }
          id = ver_part[vj:vid_end].to_owned()
        }
        if lt >= 6 &&
          l.unsafe_get(0) == 'O'.to_int().to_uint16() &&
          l.unsafe_get(1) == 'R'.to_int().to_uint16() &&
          l.unsafe_get(2) == 'I'.to_int().to_uint16() &&
          l.unsafe_get(3) == 'G'.to_int().to_uint16() &&
          l.unsafe_get(4) == 'I'.to_int().to_uint16() &&
          l.unsafe_get(5) == 'N'.to_int().to_uint16() {
          in_origin = true
          i = i + 1
          continue
        }
        if in_origin {
          let est_cap = seq_len + lt
          if seq_buf.length() < est_cap {
            let new_buf : FixedArray[UInt16] = FixedArray::make(est_cap * 2, 0)
            for k = 0; k < seq_len; k = k + 1 {
              new_buf[k] = seq_buf[k]
            }
            seq_buf = new_buf
          }
          for j = 0; j < lt; j = j + 1 {
            let c = l.unsafe_get(j)
            if c >= 'a'.to_int().to_uint16() && c <= 'z'.to_int().to_uint16() {
              seq_buf[seq_len] = uppercase_table[c.to_int()]
              seq_len = seq_len + 1
            } else if c >= 'A'.to_int().to_uint16() &&
              c <= 'Z'.to_int().to_uint16() {
              seq_buf[seq_len] = c
              seq_len = seq_len + 1
            }
          }
        }
        i = i + 1
      }
      let final_seq_buf : FixedArray[UInt16] = FixedArray::make(seq_len, 0)
      for k = 0; k < seq_len; k = k + 1 {
        final_seq_buf[k] = seq_buf[k]
      }
      records.push(SeqRecord::{
        seq: Seq::new(unsafe_from_uint16_fixedarray(final_seq_buf)),
        id,
        name,
        description,
        letter_annotations: Map([], capacity=0),
      })
    }
    i = i + 1
  }
  records
}

// --- GenBank Writer ---

///|
/// Write an array of SeqRecord objects to GenBank format text.
///
/// Mirrors BioPython's GenBankIO writer:
///   - LOCUS line: name + length + "bp" + molecule type + topology + division + date
///   - DEFINITION, ACCESSION, VERSION lines
///   - ORIGIN section with sequence wrapped at 60 chars, grouped in 10s, position labels
///   - End with "//"
pub fn write_genbank(records : Array[SeqRecord]) -> String {
  let buf = StringBuilder::new()
  for rec in records {
    // --- LOCUS line ---
    // LOCUS       NM_001     5028 bp    DNA    linear   PRI 21-JUN-1999
    // Fields (fixed positions): 0-5 LOCUS, 6-11 spaces, 12-27 locus name, 28-39 length+bp,
    // 40-47 mol_type, 48-54 topology, 55-57 division, 58- date
    let seq_len = rec.seq.length()
    let locus_name = if rec.name.length() > 16 { rec.name[0:16].to_owned() } else { rec.name }
    buf.write_string("LOCUS       ")
    buf.write_string(locus_name)
    // Pad to position 28 (12 + 16 = 28)
    let pad1 = 28 - (12 + locus_name.length())
    for _j = 0; pad1 > 0 && _j < pad1; _j = _j + 1 {
      buf.write_char(' ')
    }
    // Length + bp (right-aligned in 11 chars, positions 28-39)
    let len_str = seq_len.to_string() + " bp"
    let pad_len = 11 - len_str.length()
    for _j = 0; pad_len > 0 && _j < pad_len; _j = _j + 1 {
      buf.write_char(' ')
    }
    buf.write_string(len_str)
    buf.write_string(" DNA ")
    // Topology
    buf.write_string("linear  ")
    // Division + date (use "UNK" division and a placeholder date)
    buf.write_string("UNK 01-JAN-2024\n")

    // --- DEFINITION ---
    let mut def_text = rec.description
    if def_text.length() == 0 {
      def_text = rec.id
    }
    // Wrap DEFINITION at ~58 chars after the keyword
    buf.write_string("DEFINITION  ")
    let wrapped_def = build_wrapped(def_text, 12, 58)
    buf.write_string(wrapped_def)
    buf.write_char('.')
    buf.write_char('\n')

    // --- ACCESSION ---
    buf.write_string("ACCESSION   ")
    let acc = if rec.id.length() == 0 { rec.name } else { rec.id }
    buf.write_string(acc)
    buf.write_char('\n')

    // --- VERSION ---
    buf.write_string("VERSION     ")
    buf.write_string(acc)
    buf.write_char('\n')

    // --- ORIGIN ---
    buf.write_string("ORIGIN      \n")
    let seq_section = build_genbank_sequence(rec.seq.to_string())
    buf.write_string(seq_section)

    // End record
    buf.write_string("//\n")
  }
  buf.to_string()
}

///|
/// Build text wrapped into lines, returning the wrapped string.
/// `indent` = chars already written (so we know where we are)
/// `width`  = max chars per line
fn build_wrapped(text : String, indent : Int, width : Int) -> String {
  let buf = StringBuilder::new()
  let total = text.length()
  let mut pos = 0
  let mut cur_indent = indent
  while pos < total {
    let remaining = total - pos
    let avail = width - cur_indent
    if remaining <= avail {
      buf.write_string(text[pos:total].to_owned())
      pos = total
    } else {
      // Find a space to break at, within avail
      let mut break_at = -1
      let check_end = if pos + avail < total { pos + avail } else { total }
      for k = pos; k < check_end; k = k + 1 {
        if text.unsafe_get(k) == ' '.to_int().to_uint16() {
          break_at = k
        }
      }
      if break_at == -1 {
        // No space - hard break
        break_at = pos + avail
      }
      buf.write_string(text[pos:break_at].to_owned())
      buf.write_char('\n')
      // Continuation indent
      buf.write_string("            ")
      cur_indent = 12
      // Skip space if we broke at one
      if break_at < total && text.unsafe_get(break_at) == ' '.to_int().to_uint16() {
        pos = break_at + 1
      } else {
        pos = break_at
      }
    }
  }
  buf.to_string()
}

///|
/// Build the ORIGIN sequence section string.
/// Format: "      1 atgcatgcat gcatgcatgc atgcatgcat gcatgcatgc atgcatgcat gcatgcatgc"
/// Position label: right-aligned in 9 chars
/// Sequence: 6 groups of 10 bases per line
fn build_genbank_sequence(seq : String) -> String {
  let buf = StringBuilder::new()
  let total = seq.length()
  let mut i = 0
  let mut line_start = 0
  while i < total {
    // Position label (1-indexed)
    let pos_label = (line_start + 1).to_string()
    let pad_pos = 9 - pos_label.length()
    for _j = 0; pad_pos > 0 && _j < pad_pos; _j = _j + 1 {
      buf.write_char(' ')
    }
    buf.write_string(pos_label)
    buf.write_char(' ')

    // Write 6 groups of 10 bases
    let mut group = 0
    while group < 6 && i < total {
      let group_end = if i + 10 < total { i + 10 } else { total }
      buf.write_string(seq[i:group_end].to_owned())
      i = group_end
      buf.write_char(' ')
      group = group + 1
    }
    buf.write_char('\n')
    line_start = i
  }
  buf.to_string()
}
