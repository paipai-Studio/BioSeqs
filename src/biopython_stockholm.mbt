/// biopython_stockholm.mbt — Stockholm alignment format I/O
/// Inspired by Biopython Bio.AlignIO.StockholmIO

///|
/// Stockholm alignment record.
pub struct BiopythonStockholmAlignment {
  seqs : Array[(String, String)] // (name, sequence) pairs
  gf : Array[(String, String)] // GF lines (annotation)
  gs : Array[(String, String, String)] // (seqname, feature, value)
  gr : Array[(String, String, String)] // (seqname, feature, value)
}

///|
/// Parse a Stockholm format alignment.
pub fn biopython_parse_stockholm(text : String) -> BiopythonStockholmAlignment? {
  let lines = text.split("\n").to_array()
  let mut found_header = false
  let mut found_end = false
  let seqs : Array[(String, String)] = Array(capacity=0)
  let gf : Array[(String, String)] = Array(capacity=0)
  let gs : Array[(String, String, String)] = Array(capacity=0)
  let gr : Array[(String, String, String)] = Array(capacity=0)
  // Accumulate sequences (may be split across multiple lines)
  let seq_names : Array[String] = Array(capacity=0)
  let seq_data : Array[StringBuilder] = Array(capacity=0)
  for i = 0; i < lines.length(); i = i + 1 {
    let line = lines[i].trim().to_owned()
    if line.length() == 0 {
      continue
    }
    if line.has_prefix("# STOCKHOLM") {
      found_header = true
      continue
    }
    if line == "//" {
      found_end = true
      continue
    }
    if !found_header {
      continue
    }
    if line.has_prefix("#=GF") {
      // #=GF <feature> <value>
      let rest = line[4:].to_owned().trim().to_owned()
      let sp = biopython_sth_split2(rest)
      if sp.0.length() > 0 && sp.1.length() > 0 {
        gf.push((sp.0, sp.1))
      }
      continue
    }
    if line.has_prefix("#=GS") {
      let rest = line[4:].to_owned().trim().to_owned()
      let sp1 = biopython_sth_split2(rest)
      let sp2 = biopython_sth_split2(sp1.1)
      if sp1.0.length() > 0 && sp2.0.length() > 0 {
        gs.push((sp1.0, sp2.0, sp2.1))
      }
      continue
    }
    if line.has_prefix("#=GR") {
      let rest = line[4:].to_owned().trim().to_owned()
      let sp1 = biopython_sth_split2(rest)
      let sp2 = biopython_sth_split2(sp1.1)
      if sp1.0.length() > 0 && sp2.0.length() > 0 {
        gr.push((sp1.0, sp2.0, sp2.1))
      }
      continue
    }
    if line.has_prefix("#") {
      continue // other comment
    }
    // Sequence line: <name> <sequence>
    let sp = biopython_sth_split2(line)
    if sp.0.length() > 0 && sp.1.length() > 0 {
      let name = sp.0
      let seq_part = sp.1
      let mut found = false
      for j = 0; j < seq_names.length(); j = j + 1 {
        if seq_names[j] == name {
          let _ = seq_data[j].write(seq_part)
          found = true
          break
        }
      }
      if !found {
        seq_names.push(name)
        let sb = StringBuilder()
        let _ = sb.write(seq_part)
        seq_data.push(sb)
      }
    }
  }
  if !found_header || !found_end {
    return None
  }
  for j = 0; j < seq_names.length(); j = j + 1 {
    seqs.push((seq_names[j], seq_data[j].to_string()))
  }
  Some({ seqs, gf, gs, gr, })
}

///|
/// Split a string into first token and the rest.
fn biopython_sth_split2(s : String) -> (String, String) {
  let trimmed = s.trim().to_owned()
  let n = trimmed.length()
  let mut space_pos = -1
  for i = 0; i < n; i = i + 1 {
    if trimmed[i].to_int() == 32 { // space
      space_pos = i
      break
    }
  }
  if space_pos < 0 {
    (trimmed, "")
  } else {
    (trimmed[0:space_pos].to_owned(), trimmed[space_pos:].trim().to_owned())
  }
}

///|
/// Write a Stockholm format alignment.
pub fn biopython_write_stockholm(aln : BiopythonStockholmAlignment) -> String {
  let sb = StringBuilder()
  let _ = sb.write("# STOCKHOLM 1.0\n")
  let _ = sb.write("\n")
  // GF annotations
  for i = 0; i < aln.gf.length(); i = i + 1 {
    let (feature, value) = aln.gf[i]
    let _ = sb.write("#=GF " + feature + " " + value + "\n")
  }
  if aln.gf.length() > 0 {
    let _ = sb.write("\n")
  }
  // GS annotations
  for i = 0; i < aln.gs.length(); i = i + 1 {
    let (name, feature, value) = aln.gs[i]
    let _ = sb.write("#=GS " + name + " " + feature + " " + value + "\n")
  }
  if aln.gs.length() > 0 {
    let _ = sb.write("\n")
  }
  // Find max name length for alignment
  let mut max_name = 0
  for i = 0; i < aln.seqs.length(); i = i + 1 {
    let (name, _) = aln.seqs[i]
    if name.length() > max_name {
      max_name = name.length()
    }
  }
  // Sequences
  for i = 0; i < aln.seqs.length(); i = i + 1 {
    let (name, seq) = aln.seqs[i]
    let padding = biopython_sth_pad(name, max_name)
    let _ = sb.write(padding + " " + seq + "\n")
  }
  let _ = sb.write("\n")
  // GR annotations
  for i = 0; i < aln.gr.length(); i = i + 1 {
    let (name, feature, value) = aln.gr[i]
    let padding = biopython_sth_pad("#=GR " + name, max_name + 4)
    let _ = sb.write(padding + " " + feature + " " + value + "\n")
  }
  if aln.gr.length() > 0 {
    let _ = sb.write("\n")
  }
  let _ = sb.write("//\n")
  sb.to_string()
}

///|
/// Pad a string to given width with spaces.
fn biopython_sth_pad(s : String, width : Int32) -> String {
  let result = StringBuilder()
  let _ = result.write(s)
  for _i = s.length(); _i < width; _i = _i + 1 {
    let _ = result.write_char(' ')
  }
  result.to_string()
}
