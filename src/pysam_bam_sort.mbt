/// pysam_bam_sort.mbt — Sort SAM/BAM records
/// Inspired by pysam/samtools sort

///|
/// Sort order type.
pub enum PysamSortOrder {
  Coordinate // sort by (rname, pos)
  QueryName // sort by qname
  Tag(String) // sort by tag value
}

///|
/// Sort SAM text by coordinate or query name.
pub fn pysam_sort_sam(
  sam_text : String,
  order~ : PysamSortOrder = PysamSortOrder::Coordinate,
) -> String {
  let lines = sam_text.split("\n").to_array()
  let headers : Array[String] = Array(capacity=0)
  let records : Array[String] = Array(capacity=0)
  for i = 0; i < lines.length(); i = i + 1 {
    let line = lines[i].to_owned()
    if line.length() == 0 {
      continue
    }
    if line[0].to_int() == 64 {
      headers.push(line)
    } else {
      records.push(line)
    }
  }
  // Sort records
  match order {
    PysamSortOrder::Coordinate =>
      records.sort_by(fn(a, b) {
        let fa = a.split("\t").to_array()
        let fb = b.split("\t").to_array()
        if fa.length() < 3 || fb.length() < 3 {
          return 0
        }
        let ra = fa[2].to_owned()
        let rb = fb[2].to_owned()
        if ra < rb {
          -1
        } else if ra > rb {
          1
        } else {
          let pa = parse_int(fa[3].to_owned())
          let pb = parse_int(fb[3].to_owned())
          if pa < pb {
            -1
          } else if pa > pb {
            1
          } else {
            0
          }
        }
      })
    PysamSortOrder::QueryName =>
      records.sort_by(fn(a, b) {
        let fa = a.split("\t").to_array()
        let fb = b.split("\t").to_array()
        if fa.length() < 1 || fb.length() < 1 {
          return 0
        }
        let na = fa[0].to_owned()
        let nb = fb[0].to_owned()
        if na < nb {
          -1
        } else if na > nb {
          1
        } else {
          0
        }
      })
    PysamSortOrder::Tag(_) =>
      // Sort by tag value (simplified: no-op)
      ()
  }
  // Build output
  let sb = StringBuilder()
  for i = 0; i < headers.length(); i = i + 1 {
    let _ = sb.write(headers[i] + "\n")
  }
  for i = 0; i < records.length(); i = i + 1 {
    let _ = sb.write(records[i] + "\n")
  }
  sb.to_string()
}

///|
/// Check if SAM is sorted by coordinate.
pub fn pysam_is_sorted(sam_text : String) -> Bool {
  let lines = sam_text.split("\n").to_array()
  let mut prev_rname = ""
  let mut prev_pos = -1
  for i = 0; i < lines.length(); i = i + 1 {
    let line = lines[i].to_owned()
    if line.length() == 0 || line[0].to_int() == 64 {
      continue
    }
    let fields = line.split("\t").to_array()
    if fields.length() < 4 {
      continue
    }
    let rname = fields[2].to_owned()
    let pos = parse_int(fields[3].to_owned())
    if prev_pos >= 0 {
      if rname < prev_rname {
        return false
      }
      if rname == prev_rname && pos < prev_pos {
        return false
      }
    }
    prev_rname = rname
    prev_pos = pos
  }
  true
}

///|
/// Suppress unused constructor warnings.
pub fn pysam_bam_unused_suppress_3() -> Unit {
  ignore(PysamSortOrder::QueryName)
  ignore(PysamSortOrder::Tag(""))
}
