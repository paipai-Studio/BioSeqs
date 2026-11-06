/// pysam_bam_fastx.mbt — FASTA/FASTQ file iteration and writing
/// Items 11-13: FastxFile迭代器, FASTQ写入器, FASTA spliced fetch

///|
/// A FASTA/FASTQ record.
pub struct PysamFastxRecord {
  name : String
  seq : String
  qual : String // empty for FASTA
  is_fastq : Bool
}

///|
/// Iterate FASTA/FASTQ text, auto-detecting format.
/// Returns records one by one.
pub fn pysam_fastx_iter(text : String) -> Array[PysamFastxRecord] {
  let lines = text.split("\n").to_array()
  let records : Array[PysamFastxRecord] = Array(capacity=0)
  let i = 0
  let n = lines.length()
  // Auto-detect: first non-empty line starts with @ (fastq) or > (fasta)
  let mut idx = i
  while idx < n {
    let line = lines[idx].to_owned()
    if line.length() == 0 {
      idx = idx + 1
      continue
    }
    if line[0].to_int() == 64 { // '@' = FASTQ
      // Name line
      let name = line[1:].to_owned()
      idx = idx + 1
      if idx >= n {
        break
      }
      let seq = lines[idx].to_owned()
      idx = idx + 1
      if idx >= n {
        break
      }
      // Plus line
      idx = idx + 1
      if idx >= n {
        break
      }
      let qual = lines[idx].to_owned()
      idx = idx + 1
      records.push({ name, seq, qual, is_fastq: true, })
    } else if line[0].to_int() == 62 { // '>' = FASTA
      let name = line[1:].to_owned()
      idx = idx + 1
      let seq_sb = StringBuilder()
      while idx < n {
        let seq_line = lines[idx].to_owned()
        if seq_line.length() == 0 {
          idx = idx + 1
          continue
        }
        if seq_line[0].to_int() == 62 || seq_line[0].to_int() == 64 {
          break
        }
        let _ = seq_sb.write(seq_line)
        idx = idx + 1
      }
      records.push({ name, seq: seq_sb.to_string(), qual: "", is_fastq: false, })
    } else {
      idx = idx + 1
    }
  }
  records
}

///|
/// FASTQ writer: format a record to FASTQ text.
pub fn pysam_fastq_format(
  name : String,
  seq : String,
  qual : String,
  phred~ : Int32 = 33,
) -> String {
  ignore(phred)
  "@" + name + "\n" + seq + "\n+\n" + qual + "\n"
}

///|
/// Write multiple records to FASTQ text.
pub fn pysam_fastq_write(records : Array[PysamFastxRecord]) -> String {
  let sb = StringBuilder()
  for i = 0; i < records.length(); i = i + 1 {
    let r = records[i]
    if r.is_fastq {
      let _ = sb.write("@" + r.name + "\n")
      let _ = sb.write(r.seq + "\n")
      let _ = sb.write("+\n")
      let _ = sb.write(r.qual + "\n")
    } else {
      let _ = sb.write(">" + r.name + "\n")
      let _ = sb.write(r.seq + "\n")
    }
  }
  sb.to_string()
}

///|
/// Write FASTA with wrapped sequence (default 60 chars per line).
pub fn pysam_fasta_write(
  name : String,
  seq : String,
  wrap~ : Int32 = 60,
) -> String {
  let sb = StringBuilder()
  let _ = sb.write(">" + name + "\n")
  let mut i = 0
  while i < seq.length() {
    let end = if i + wrap < seq.length() { i + wrap } else { seq.length() }
    let _ = sb.write(seq[i:end].to_owned() + "\n")
    i = end
  }
  sb.to_string()
}

///|
/// FASTA spliced fetch: extract spliced sequence using exon coordinates.
/// exons: array of (start, end) in reference coordinates (0-based, exclusive end)
pub fn pysam_fasta_spliced_fetch(
  ref_seq : String,
  exons : Array[(Int32, Int32)],
  strand~ : Bool = true,
) -> String {
  let sb = StringBuilder()
  for i = 0; i < exons.length(); i = i + 1 {
    let (start, end) = exons[i]
    let s = if start < 0 { 0 } else { start }
    let e = if end > ref_seq.length() { ref_seq.length() } else { end }
    if s < e {
      let _ = sb.write(ref_seq[s:e].to_owned())
    }
  }
  let result = sb.to_string()
  if strand {
    result
  } else {
    // Reverse complement
    pysam_spliced_revcomp(result)
  }
}

///|
/// Reverse complement helper.
fn pysam_spliced_revcomp(seq : String) -> String {
  let n = seq.length()
  let result : Array[Int32] = Array::make(n, 0)
  for i = 0; i < n; i = i + 1 {
    let comp = match9(seq[i].to_int())
    result[n - 1 - i] = comp
  }
  let sb = StringBuilder()
  for i = 0; i < n; i = i + 1 {
    let _ = sb.write_char(result[i].unsafe_to_char())
  }
  sb.to_string()
}

///|
/// Match for complement base.
fn match9(c : Int32) -> Int32 {
  match c {
    65 => 84 // A->T
    84 => 65 // T->A
    71 => 67 // G->C
    67 => 71 // C->G
    97 => 116 // a->t
    116 => 97 // t->a
    103 => 99 // g->c
    99 => 103 // c->g
    x => x
  }
}

///|
/// Convert SAM records to FASTQ.
pub fn pysam_sam_to_fastq(sam_text : String) -> (String, String) {
  let lines = sam_text.split("\n").to_array()
  let r1_sb = StringBuilder()
  let r2_sb = StringBuilder()
  for i = 0; i < lines.length(); i = i + 1 {
    let line = lines[i].to_owned()
    if line.length() == 0 || line[0].to_int() == 64 {
      continue
    }
    let fields = line.split("\t").to_array()
    if fields.length() < 11 {
      continue
    }
    let qname = fields[0].to_owned()
    let flag = parse_int(fields[1].to_owned())
    let seq = fields[9].to_owned()
    let qual = fields[10].to_owned()
    let is_read1 = (flag & 64) != 0
    let is_read2 = (flag & 128) != 0
    if is_read1 {
      let _ = r1_sb.write(pysam_fastq_format(qname, seq, qual))
    } else if is_read2 {
      let _ = r2_sb.write(pysam_fastq_format(qname, seq, qual))
    } else {
      let _ = r1_sb.write(pysam_fastq_format(qname, seq, qual))
    }
  }
  (r1_sb.to_string(), r2_sb.to_string())
}
