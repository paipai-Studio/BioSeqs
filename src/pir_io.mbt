// PIR/NBRF format parser, mirroring BioPython's Bio.SeqIO.PirIO.
//
// PIR (Protein Information Resource) / NBRF format:
//
//   >P1;ID123456
//   entry description text here
//   MVHYRACEVICDK*
//   >P1;ID123457
//   another entry description
//   MVLSQDEVVCFSNQLPRQGLDAA*
//
// Structure:
//   - First line: ">" + entry_code + ";" + name
//     Entry codes: P1 = protein (full), P2 = protein (fragment),
//                  F1 = fragment, DL = DNA linear, DC = DNA circular,
//                  RL = RNA linear, RC = RNA circular, N1 = other nucleic,
//                  N3 = tRNA, N4 = other RNA
//   - Second line: free-text description
//   - Remaining lines: sequence, ending with "*"
//
// Key fields extracted:
//   - id:   name (after semicolon on first line)
//   - name: name
//   - description: second line text
//   - seq:  sequence characters (spaces ignored, trailing "*" stripped)

///|
/// Parse PIR/NBRF format text and return an array of SeqRecord objects.
///
/// Mirrors BioPython's SeqIO.parse(handle, "pir").
pub fn parse_pir(content : String) -> Array[SeqRecord] {
  let records = Array::new()
  let lines = content.split("\n").to_array()
  let n = lines.length()
  let mut i = 0
  let cr : UInt16 = '\r'.to_int().to_uint16()
  let gt : UInt16 = '>'.to_int().to_uint16()
  let semi : UInt16 = ';'.to_int().to_uint16()
  let star : UInt16 = '*'.to_int().to_uint16()
  let sp : UInt16 = ' '.to_int().to_uint16()

  while i < n {
    let line = lines[i]
    let len = line.length()
    let mut trimmed = len
    if len > 0 && line.unsafe_get(len - 1) == cr {
      trimmed = len - 1
    }
    if trimmed > 1 && line.unsafe_get(0) == gt {
      // Parse header ">P1;NAME"
      let header = line[1:trimmed].to_owned()
      let hl = header.length()
      let mut semi_pos : Int? = None
      for j = 0; j < hl; j = j + 1 {
        if header.unsafe_get(j) == semi {
          semi_pos = Some(j)
          break
        }
      }
      let name = match semi_pos {
        Some(p) => if p + 1 < hl { header[p + 1:hl].to_owned() } else { "" }
        None => header
      }

      i = i + 1
      // Description line (second line of entry)
      let mut description = ""
      if i < n {
        let dl = lines[i]
        let dll = dl.length()
        let mut dlt = dll
        if dll > 0 && dl.unsafe_get(dll - 1) == cr {
          dlt = dll - 1
        }
        description = dl[0:dlt].to_owned()
        i = i + 1
      }

      // Sequence lines until "*" (included) or next ">"
      let seq_buf = StringBuilder::new()
      while i < n {
        let sl = lines[i]
        let sll = sl.length()
        let mut slt = sll
        if sll > 0 && sl.unsafe_get(sll - 1) == cr {
          slt = sll - 1
        }
        // Check for next header
        if slt > 1 && sl.unsafe_get(0) == gt {
          break
        }
        // Process sequence chars
        for j = 0; j < slt; j = j + 1 {
          let c = sl.unsafe_get(j)
          if c == star {
            // End of sequence
          } else if c != sp {
            seq_buf.write_char(c.unsafe_to_char())
          }
        }
        // Check if this line ended with * (end of entry)
        let mut ends_star = false
        for j = 0; j < slt; j = j + 1 {
          if sl.unsafe_get(j) == star {
            ends_star = true
            break
          }
        }
        i = i + 1
        if ends_star {
          break
        }
      }

      records.push(SeqRecord::new(
        Seq::new(seq_buf.to_string()),
        id=name,
        name=name,
        description=description,
      ))
      continue
    }

    i = i + 1
  }
  records
}
