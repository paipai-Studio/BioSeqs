/// pysam_bam_flag_utils.mbt — Flag conversion, CIGAR simplification, SAM validation
/// Items 18-20

///|
/// Convert BAM flag to human-readable string.
pub fn pysam_flag_to_str(flag : Int32) -> String {
  let parts : Array[String] = Array(capacity=0)
  if (flag & 1) != 0 {
    parts.push("PAIRED")
  }
  if (flag & 2) != 0 {
    parts.push("PROPER_PAIR")
  }
  if (flag & 4) != 0 {
    parts.push("UNMAPPED")
  }
  if (flag & 8) != 0 {
    parts.push("MATE_UNMAPPED")
  }
  if (flag & 16) != 0 {
    parts.push("REVERSE")
  }
  if (flag & 32) != 0 {
    parts.push("MATE_REVERSE")
  }
  if (flag & 64) != 0 {
    parts.push("READ1")
  }
  if (flag & 128) != 0 {
    parts.push("READ2")
  }
  if (flag & 256) != 0 {
    parts.push("SECONDARY")
  }
  if (flag & 512) != 0 {
    parts.push("QCFAIL")
  }
  if (flag & 1024) != 0 {
    parts.push("DUP")
  }
  if (flag & 2048) != 0 {
    parts.push("SUPPLEMENTARY")
  }
  parts.join(",")
}

///|
/// Convert human-readable flag string to integer.
pub fn pysam_str_to_flag(s : String) -> Int32 {
  let mut flag = 0
  let parts = s.split(",").to_array()
  for i = 0; i < parts.length(); i = i + 1 {
    let part = parts[i].to_owned().to_upper()
    match part {
      "PAIRED" => flag = flag | 1
      "PROPER_PAIR" => flag = flag | 2
      "UNMAPPED" => flag = flag | 4
      "MATE_UNMAPPED" => flag = flag | 8
      "REVERSE" => flag = flag | 16
      "MATE_REVERSE" => flag = flag | 32
      "READ1" => flag = flag | 64
      "READ2" => flag = flag | 128
      "SECONDARY" => flag = flag | 256
      "QCFAIL" => flag = flag | 512
      "DUP" => flag = flag | 1024
      "SUPPLEMENTARY" => flag = flag | 2048
      _ => ()
    }
  }
  flag
}

///|
/// Simplify CIGAR: merge adjacent same ops, remove zero-length ops.
pub fn pysam_simplify_cigar(cigar : String) -> String {
  // Parse into elements
  let ops : Array[(Int32, Int32)] = Array(capacity=0) // (length, op_char)
  let mut num = 0
  for i = 0; i < cigar.length(); i = i + 1 {
    let c = cigar[i].to_int()
    if c >= 48 && c <= 57 {
      num = num * 10 + (c - 48)
    } else {
      if num > 0 {
        ops.push((num, c))
      }
      num = 0
    }
  }
  // Merge adjacent same ops
  let merged : Array[(Int32, Int32)] = Array(capacity=0)
  for i = 0; i < ops.length(); i = i + 1 {
    let (len, op) = ops[i]
    if merged.length() > 0 {
      let (prev_len, prev_op) = merged[merged.length() - 1]
      if prev_op == op {
        merged[merged.length() - 1] = (prev_len + len, op)
      } else {
        merged.push((len, op))
      }
    } else {
      merged.push((len, op))
    }
  }
  // Build output
  let sb = StringBuilder()
  for i = 0; i < merged.length(); i = i + 1 {
    let (len, op) = merged[i]
    let _ = sb.write(len.to_string())
    let _ = sb.write_char(op.unsafe_to_char())
  }
  sb.to_string()
}

///|
/// Validate a SAM line.
/// Returns list of error messages (empty = valid).
pub fn pysam_validate_sam_line(line : String) -> Array[String] {
  let errors : Array[String] = Array(capacity=0)
  let fields = line.split("\t").to_array()
  // Must have at least 11 mandatory fields
  if fields.length() < 11 {
    errors.push("Insufficient fields: " + fields.length().to_string() + " < 11")
    return errors
  }
  // QNAME: non-empty
  if fields[0].length() == 0 {
    errors.push("Empty QNAME")
  }
  // FLAG: numeric, 0-65535
  let flag = parse_int(fields[1].to_owned())
  if flag < 0 || flag > 65535 {
    errors.push("Invalid FLAG: " + flag.to_string())
  }
  // RNAME: * or non-empty
  let rname = fields[2].to_owned()
  if rname.length() == 0 {
    errors.push("Empty RNAME")
  }
  // POS: 0 or positive
  let pos = parse_int(fields[3].to_owned())
  if pos < 0 {
    errors.push("Invalid POS: " + pos.to_string())
  }
  // MAPQ: 0-255
  let mapq = parse_int(fields[4].to_owned())
  if mapq < 0 || mapq > 255 {
    errors.push("Invalid MAPQ: " + mapq.to_string())
  }
  // CIGAR: * or valid pattern
  let cigar = fields[5].to_owned()
  if cigar != "*" {
    if !pysam_valid_cigar(cigar) {
      errors.push("Invalid CIGAR: " + cigar)
    }
  }
  // RNEXT: * or = or same as RNAME
  let rnext = fields[6].to_owned()
  if rnext.length() == 0 {
    errors.push("Empty RNEXT")
  }
  // PNEXT: 0 or positive
  let pnext = parse_int(fields[7].to_owned())
  if pnext < 0 {
    errors.push("Invalid PNEXT: " + pnext.to_string())
  }
  // TLEN: any integer
  // SEQ: * or valid sequence
  let seq = fields[9].to_owned()
  if seq != "*" {
    for i = 0; i < seq.length(); i = i + 1 {
      let c = seq[i].to_int()
      let valid = c == 65 ||
        c == 67 ||
        c == 71 ||
        c == 84 ||
        c == 78 ||
        c == 97 ||
        c == 99 ||
        c == 103 ||
        c == 116 ||
        c == 110 ||
        c == 46 // . for padding
      if !valid {
        errors.push("Invalid base in SEQ: " + c.unsafe_to_char().to_string())
        break
      }
    }
  }
  // QUAL: * or same length as SEQ
  let qual = fields[10].to_owned()
  if qual != "*" && seq != "*" {
    if qual.length() != seq.length() {
      errors.push(
        "QUAL length " +
        qual.length().to_string() +
        " != SEQ length " +
        seq.length().to_string(),
      )
    }
  }
  // Flag consistency checks
  if (flag & 4) != 0 && rname != "*" {
    errors.push("UNMAPPED flag set but RNAME != *")
  }
  if (flag & 4) != 0 && pos != 0 {
    errors.push("UNMAPPED flag set but POS != 0")
  }
  if (flag & 4) == 0 && rname == "*" {
    errors.push("Mapped read but RNAME == *")
  }
  errors
}

///|
/// Check if CIGAR string is valid.
fn pysam_valid_cigar(cigar : String) -> Bool {
  if cigar.length() == 0 {
    return false
  }
  let mut has_digit = false
  for i = 0; i < cigar.length(); i = i + 1 {
    let c = cigar[i].to_int()
    if c >= 48 && c <= 57 {
      has_digit = true
    } else {
      // Must be a valid CIGAR op
      let valid_op = c == 77 ||
        c == 73 ||
        c == 68 ||
        c == 78 ||
        c == 83 ||
        c == 72 ||
        c == 80 ||
        c == 61 ||
        c == 88
      if !valid_op {
        return false
      }
      if !has_digit {
        return false // op without preceding number
      }
      has_digit = false
    }
  }
  !has_digit // should end with an op, not a digit
}

///|
/// Validate entire SAM text.
pub fn pysam_validate_sam(sam_text : String) -> Array[(Int32, Array[String])] {
  let lines = sam_text.split("\n").to_array()
  let errors : Array[(Int32, Array[String])] = Array(capacity=0)
  for i = 0; i < lines.length(); i = i + 1 {
    let line = lines[i].to_owned()
    if line.length() == 0 || line[0].to_int() == 64 {
      continue
    }
    let line_errors = pysam_validate_sam_line(line)
    if line_errors.length() > 0 {
      errors.push((i, line_errors))
    }
  }
  errors
}
