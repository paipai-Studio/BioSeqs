// EMBL format parser, mirroring BioPython's Bio.SeqIO.EmblIO.
//
// EMBL is the European nucleotide archive format, similar to GenBank but
// with different two-letter field codes:
//
//   ID   HSBGLOD; SV 1; linear; genomic DNA; STD; HUM; 3876 BP.
//   XX
//   AC   M16496;
//   XX
//   DT   21-JUN-1999 (Rel. 38, Created)
//   DT   21-JUN-1999 (Rel. 38, Last updated, Version 1)
//   XX
//   DE   Homo sapiens beta-globin gene region, 5' regulatory and coding
//   DE   sequences, clone BG17.
//   XX
//   KW   .
//   XX
//   OS   Homo sapiens (human)
//   OC   Eukaryota; Metazoa; Chordata; Craniata; Vertebrata; Euteleostomi;
//   OC   Mammalia; Eutheria; Euarchontoglires; Primates; Haplorrhini;
//   OC   Simiiformes; Catarrhini; Hominoidea; Hominidae; Homo.
//   XX
//   SQ   Sequence 3876 BP; 924 A; 1024 C; 980 G; 948 T; 0 other;
//        atgcatgcat gcatgcatgc atgcatgcat gcatgcatgc atgcatgcat gcatgcatgc
//       61 ...
//   //
//
// Key fields extracted:
//   - id:          AC (first accession number) or ID locus name
//   - name:        ID locus name (first token before ';' on ID line)
//   - description: DE lines (joined)
//   - seq:         SQ section (numbers and spaces stripped)

///|
/// Parse EMBL format text and return an array of SeqRecord objects.
///
/// Mirrors BioPython's SeqIO.parse(handle, "embl").
/// Extracts: id (AC or ID), name (ID locus), description (DE lines),
/// and sequence (SQ section). Multiple records separated by '//'.
// Uppercase a string (copy of the one defined in genbank_io.mbt; kept
// self-contained so the file doesn't depend on implicit ordering.)
fn uppercase_string(s : String) -> String {
  let buf = StringBuilder::new()
  let n = s.length()
  let a : Int = 'a'.to_int()
  let z : Int = 'z'.to_int()
  for i = 0; i < n; i = i + 1 {
    let c = s.unsafe_get(i)
    let ci = c.to_int()
    if ci >= a && ci <= z {
      buf.write_char((ci - 32).to_uint16().unsafe_to_char())
    } else {
      buf.write_char(c.unsafe_to_char())
    }
  }
  buf.to_string()
}

pub fn parse_embl(content : String) -> Array[SeqRecord] {
  let records = Array::new()
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
    if trimmed >= 2 &&
      line.unsafe_get(0) == 'I'.to_int().to_uint16() &&
      line.unsafe_get(1) == 'D'.to_int().to_uint16() &&
      (trimmed == 2 || line.unsafe_get(2) == sp) {
      let mut name = "<unknown>"
      let mut id = "<unknown id>"
      let mut description = "<unknown description>"
      let seq_buf = StringBuilder::new()
      let mut in_sq = false
      let mut have_ac = false

      // Parse ID line: "ID   name; SV ver; ..."
      if trimmed > 2 {
        let id_rest = line[2:trimmed]
        let mut j = 0
        let rl = id_rest.length()
        while j < rl && id_rest.unsafe_get(j) == sp {
          j = j + 1
        }
        let mut name_end = j
        let semi : UInt16 = ';'.to_int().to_uint16()
        while name_end < rl && id_rest.unsafe_get(name_end) != semi && id_rest.unsafe_get(name_end) != sp {
          name_end = name_end + 1
        }
        name = id_rest[j:name_end].to_owned()
        id = name
      }

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

        // Skip XX separator lines
        if lt >= 2 &&
          l.unsafe_get(0) == 'X'.to_int().to_uint16() &&
          l.unsafe_get(1) == 'X'.to_int().to_uint16() &&
          (lt == 2 || l.unsafe_get(2) == sp) {
          i = i + 1
          continue
        }

        // AC line: "AC   ACC1; ACC2;"
        if !have_ac && lt >= 2 &&
          l.unsafe_get(0) == 'A'.to_int().to_uint16() &&
          l.unsafe_get(1) == 'C'.to_int().to_uint16() &&
          (lt == 2 || l.unsafe_get(2) == sp) {
          if lt > 2 {
            let ac_part = l[2:lt]
            let mut aj = 0
            let apl = ac_part.length()
            while aj < apl && ac_part.unsafe_get(aj) == sp {
              aj = aj + 1
            }
            let mut ac_end = aj
            let semi : UInt16 = ';'.to_int().to_uint16()
            while ac_end < apl && ac_part.unsafe_get(ac_end) != semi && ac_part.unsafe_get(ac_end) != sp {
              ac_end = ac_end + 1
            }
            if ac_end > aj {
              id = ac_part[aj:ac_end].to_owned()
              have_ac = true
            }
          }
          i = i + 1
          continue
        }

        // DE line: "DE   description text"
        if lt >= 2 &&
          l.unsafe_get(0) == 'D'.to_int().to_uint16() &&
          l.unsafe_get(1) == 'E'.to_int().to_uint16() &&
          (lt == 2 || l.unsafe_get(2) == sp) {
          let def_buf = StringBuilder::new()
          if lt > 2 {
            let def_part = l[2:lt]
            let mut dj = 0
            let dpl = def_part.length()
            while dj < dpl && def_part.unsafe_get(dj) == sp {
              dj = dj + 1
            }
            def_buf.write_string(def_part[dj:dpl].to_owned())
          }
          i = i + 1
          while i < n {
            let cont = lines[i]
            let cl = cont.length()
            let mut ct = cl
            if cl > 0 && cont.unsafe_get(cl - 1) == cr {
              ct = cl - 1
            }
            if ct < 2 ||
              !(cont.unsafe_get(0) == 'D'.to_int().to_uint16() &&
                cont.unsafe_get(1) == 'E'.to_int().to_uint16() &&
                (ct == 2 || cont.unsafe_get(2) == sp)) {
              break
            }
            def_buf.write_char(' ')
            if ct > 2 {
              let mut cj = 0
              while cj < ct && cont.unsafe_get(cj) == sp {
                cj = cj + 1
              }
              def_buf.write_string(cont[cj:ct].to_owned())
            }
            i = i + 1
          }
          description = def_buf.to_string()
          continue
        }

        // SQ section start: "SQ   Sequence N BP; ..."
        if lt >= 2 &&
          l.unsafe_get(0) == 'S'.to_int().to_uint16() &&
          l.unsafe_get(1) == 'Q'.to_int().to_uint16() &&
          (lt == 2 || l.unsafe_get(2) == sp) {
          in_sq = true
          i = i + 1
          continue
        }

        // Parse SQ sequence lines: start with spaces, then digits optional, then bases
        if in_sq {
          let l_start = if lt > 5 { 5 } else { lt }
          // Check if this is a sequence line (not all whitespace, contains bases)
          let mut has_base = false
          let mut k = 0
          while k < lt {
            let c = l.unsafe_get(k)
            let ci = c.to_int()
            let is_base =
              (ci >= 'A'.to_int() && ci <= 'Z'.to_int()) ||
              (ci >= 'a'.to_int() && ci <= 'z'.to_int()) ||
              c == '-'.to_int().to_uint16()
            if is_base {
              has_base = true
              break
            }
            k = k + 1
          }
          if has_base {
            for j = 0; j < lt; j = j + 1 {
              let c = l.unsafe_get(j)
              let ci = c.to_int()
              let is_base =
                (ci >= 'A'.to_int() && ci <= 'Z'.to_int()) ||
                (ci >= 'a'.to_int() && ci <= 'z'.to_int()) ||
                c == '-'.to_int().to_uint16()
              if is_base {
                seq_buf.write_char(c.unsafe_to_char())
              }
            }
          }
        }

        i = i + 1
      }

      // Strip trailing period from description (mirror GenBank parser behavior)
      let dot : UInt16 = '.'.to_int().to_uint16()
      if description.length() > 0 &&
        description.unsafe_get(description.length() - 1) == dot {
        description = description[0:description.length() - 1].to_owned()
      }
      // Uppercase the sequence for consistency with GenBank parser
      let seq_str = seq_buf.to_string()
      let seq_uppercased = uppercase_string(seq_str)
      records.push(SeqRecord::new(
        Seq::new(seq_uppercased),
        id=id,
        name=name,
        description=description,
      ))
      continue
    }

    i = i + 1
  }
  records
}
