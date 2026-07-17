///|
/// BLAST result parsing, inspired by Biopython Bio.Blast and Bio.SearchIO.
///
/// Supports parsing BLAST tabular format (-outfmt 6 and -outfmt 7) and basic
/// BLAST XML format.

///|
/// High-scoring Segment Pair (HSP) - represents a local alignment between
/// query and subject sequences.
pub struct BlastHsp {
  query_start : Int
  query_end : Int
  subject_start : Int
  subject_end : Int
  query_strand : Strand
  subject_strand : Strand
  identity : Double
  positive : Int
  gaps : Int
  align_length : Int
  bitscore : Double
  evalue : Double
  query_seq : String
  subject_seq : String
  match_seq : String
}

///|
/// Create a BlastHsp.
pub fn BlastHsp::new(
  query_start : Int,
  query_end : Int,
  subject_start : Int,
  subject_end : Int,
  query_strand : Strand,
  subject_strand : Strand,
  identity : Double,
  positive : Int,
  gaps : Int,
  align_length : Int,
  bitscore : Double,
  evalue : Double,
  query_seq : String,
  subject_seq : String,
  match_seq : String,
) -> BlastHsp {
  BlastHsp::{
    query_start,
    query_end,
    subject_start,
    subject_end,
    query_strand,
    subject_strand,
    identity,
    positive,
    gaps,
    align_length,
    bitscore,
    evalue,
    query_seq,
    subject_seq,
    match_seq,
  }
}

///|
/// Get query start position (1-based).
pub fn BlastHsp::query_start(self : BlastHsp) -> Int {
  self.query_start
}

///|
/// Get query end position (1-based).
pub fn BlastHsp::query_end(self : BlastHsp) -> Int {
  self.query_end
}

///|
/// Get subject start position (1-based).
pub fn BlastHsp::subject_start(self : BlastHsp) -> Int {
  self.subject_start
}

///|
/// Get subject end position (1-based).
pub fn BlastHsp::subject_end(self : BlastHsp) -> Int {
  self.subject_end
}

///|
/// Get query strand.
pub fn BlastHsp::query_strand(self : BlastHsp) -> Strand {
  self.query_strand
}

///|
/// Get subject strand.
pub fn BlastHsp::subject_strand(self : BlastHsp) -> Strand {
  self.subject_strand
}

///|
/// Get percent identity.
pub fn BlastHsp::identity(self : BlastHsp) -> Double {
  self.identity
}

///|
/// Get number of positive (similar) residues.
pub fn BlastHsp::positive(self : BlastHsp) -> Int {
  self.positive
}

///|
/// Get number of gaps.
pub fn BlastHsp::gaps(self : BlastHsp) -> Int {
  self.gaps
}

///|
/// Get alignment length.
pub fn BlastHsp::align_length(self : BlastHsp) -> Int {
  self.align_length
}

///|
/// Get bit score.
pub fn BlastHsp::bitscore(self : BlastHsp) -> Double {
  self.bitscore
}

///|
/// Get E-value.
pub fn BlastHsp::evalue(self : BlastHsp) -> Double {
  self.evalue
}

///|
/// Get query sequence in alignment.
pub fn BlastHsp::query_seq(self : BlastHsp) -> String {
  self.query_seq
}

///|
/// Get subject sequence in alignment.
pub fn BlastHsp::subject_seq(self : BlastHsp) -> String {
  self.subject_seq
}

///|
/// Get match sequence (| for match, + for positive, space otherwise).
pub fn BlastHsp::match_seq(self : BlastHsp) -> String {
  self.match_seq
}

///|
/// A BLAST hit representing a subject sequence with one or more HSPs.
pub struct BlastHit {
  subject_id : String
  subject_def : String
  subject_length : Int
  hsps : Array[BlastHsp]
}

///|
/// Create a BlastHit.
pub fn BlastHit::new(
  subject_id : String,
  subject_def : String,
  subject_length : Int,
  hsps : Array[BlastHsp],
) -> BlastHit {
  BlastHit::{ subject_id, subject_def, subject_length, hsps }
}

///|
/// Get subject ID.
pub fn BlastHit::subject_id(self : BlastHit) -> String {
  self.subject_id
}

///|
/// Get subject definition/description.
pub fn BlastHit::subject_def(self : BlastHit) -> String {
  self.subject_def
}

///|
/// Get subject length.
pub fn BlastHit::subject_length(self : BlastHit) -> Int {
  self.subject_length
}

///|
/// Get HSPs for this hit.
pub fn BlastHit::hsps(self : BlastHit) -> Array[BlastHsp] {
  self.hsps.copy()
}

///|
/// Get number of HSPs.
pub fn BlastHit::num_hsps(self : BlastHit) -> Int {
  self.hsps.length()
}

///|
/// Get best (highest bit score) HSP.
pub fn BlastHit::best_hsp(self : BlastHit) -> BlastHsp? {
  if self.hsps.length() == 0 {
    return None
  }
  let mut best = self.hsps[0]
  let mut i = 1
  while i < self.hsps.length() {
    if self.hsps[i].bitscore() > best.bitscore() {
      best = self.hsps[i]
    }
    i = i + 1
  }
  Some(best)
}

///|
/// A complete BLAST search result containing one query and its hits.
pub struct BlastRecord {
  query_id : String
  query_def : String
  query_length : Int
  hits : Array[BlastHit]
  program : String
  version : String
  database : String
}

///|
/// Create a BlastRecord.
pub fn BlastRecord::new(
  query_id : String,
  query_def : String,
  query_length : Int,
  hits : Array[BlastHit],
  program : String,
  version : String,
  database : String,
) -> BlastRecord {
  BlastRecord::{
    query_id,
    query_def,
    query_length,
    hits,
    program,
    version,
    database,
  }
}

///|
/// Get query ID.
pub fn BlastRecord::query_id(self : BlastRecord) -> String {
  self.query_id
}

///|
/// Get query definition/description.
pub fn BlastRecord::query_def(self : BlastRecord) -> String {
  self.query_def
}

///|
/// Get query length.
pub fn BlastRecord::query_length(self : BlastRecord) -> Int {
  self.query_length
}

///|
/// Get hits for this query.
pub fn BlastRecord::hits(self : BlastRecord) -> Array[BlastHit] {
  self.hits.copy()
}

///|
/// Get number of hits.
pub fn BlastRecord::num_hits(self : BlastRecord) -> Int {
  self.hits.length()
}

///|
/// Get BLAST program name.
pub fn BlastRecord::program(self : BlastRecord) -> String {
  self.program
}

///|
/// Get BLAST version.
pub fn BlastRecord::version(self : BlastRecord) -> String {
  self.version
}

///|
/// Get database name.
pub fn BlastRecord::database(self : BlastRecord) -> String {
  self.database
}

///|
/// Get best hit (highest bit score among all HSPs).
pub fn BlastRecord::best_hit(self : BlastRecord) -> BlastHit? {
  if self.hits.length() == 0 {
    return None
  }
  let mut best_hit = self.hits[0]
  let mut best_score = match best_hit.best_hsp() {
    Some(hsp) => hsp.bitscore()
    None => 0.0
  }
  let mut i = 1
  while i < self.hits.length() {
    match self.hits[i].best_hsp() {
      Some(hsp) =>
        if hsp.bitscore() > best_score {
          best_score = hsp.bitscore()
          best_hit = self.hits[i]
        }
      None => ()
    }
    i = i + 1
  }
  Some(best_hit)
}

///|
/// Filter hits by E-value threshold.
pub fn BlastRecord::filter_by_evalue(
  self : BlastRecord,
  threshold : Double,
) -> BlastRecord {
  let filtered_hits : Array[BlastHit] = []
  let mut i = 0
  while i < self.hits.length() {
    let hit = self.hits[i]
    let filtered_hsps : Array[BlastHsp] = []
    let mut j = 0
    while j < hit.hsps().length() {
      if hit.hsps()[j].evalue() <= threshold {
        filtered_hsps.push(hit.hsps()[j])
      }
      j = j + 1
    }
    if filtered_hsps.length() > 0 {
      filtered_hits.push(
        BlastHit::new(
          hit.subject_id(),
          hit.subject_def(),
          hit.subject_length(),
          filtered_hsps,
        ),
      )
    }
    i = i + 1
  }
  BlastRecord::new(
    self.query_id(),
    self.query_def(),
    self.query_length(),
    filtered_hits,
    self.program(),
    self.version(),
    self.database(),
  )
}

///|
/// Filter hits by percent identity threshold.
pub fn BlastRecord::filter_by_identity(
  self : BlastRecord,
  threshold : Double,
) -> BlastRecord {
  let filtered_hits : Array[BlastHit] = []
  let mut i = 0
  while i < self.hits.length() {
    let hit = self.hits[i]
    let filtered_hsps : Array[BlastHsp] = []
    let mut j = 0
    while j < hit.hsps().length() {
      if hit.hsps()[j].identity() >= threshold {
        filtered_hsps.push(hit.hsps()[j])
      }
      j = j + 1
    }
    if filtered_hsps.length() > 0 {
      filtered_hits.push(
        BlastHit::new(
          hit.subject_id(),
          hit.subject_def(),
          hit.subject_length(),
          filtered_hsps,
        ),
      )
    }
    i = i + 1
  }
  BlastRecord::new(
    self.query_id(),
    self.query_def(),
    self.query_length(),
    filtered_hits,
    self.program(),
    self.version(),
    self.database(),
  )
}

///|
/// Parse BLAST tabular format (-outfmt 6 or -outfmt 7).
///
/// Standard fields for -outfmt 6:
/// qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore
///
/// -outfmt 7 is the same as -outfmt 6 but with comment lines starting with '#'.
pub fn parse_blast_tabular(content : String) -> BlastRecord? {
  let lines = split_by_char(content, '\n'.to_int())
  if lines.length() == 0 {
    return None
  }

  let mut query_id = "unknown"
  let mut query_def = ""
  let query_length = 0
  let hits : Array[BlastHit] = Array::new()

  let mut i = 0
  while i < lines.length() {
    let line = lines[i].trim().to_string()
    if line.length() == 0 {
      i = i + 1
      continue
    }

    if starts_with(line, 0, "# Query: ") {
      query_id = line[9:].to_string()
      query_def = query_id
    } else if !starts_with(line, 0, "#") {
      let fields = split_by_char(line, 9)
      if fields.length() >= 12 {
        if query_id == "unknown" {
          query_id = fields[0]
          query_def = fields[0]
        }
        let sseqid = fields[1]

        let pident = match parse_double(fields[2]) {
          Some(v) => v
          None => 0.0
        }
        let align_length = match parse_double(fields[3]) {
          Some(v) => v.to_int()
          None => 0
        }
        let gaps = match parse_double(fields[5]) {
          Some(v) => v.to_int()
          None => 0
        }
        let qstart = match parse_double(fields[6]) {
          Some(v) => v.to_int()
          None => 0
        }
        let qend = match parse_double(fields[7]) {
          Some(v) => v.to_int()
          None => 0
        }
        let sstart = match parse_double(fields[8]) {
          Some(v) => v.to_int()
          None => 0
        }
        let send = match parse_double(fields[9]) {
          Some(v) => v.to_int()
          None => 0
        }
        let evalue = match parse_double(fields[10]) {
          Some(v) => v
          None => 0.0
        }
        let bitscore = match parse_double(fields[11]) {
          Some(v) => v
          None => 0.0
        }

        let query_strand = if qstart <= qend {
          Strand::Plus
        } else {
          Strand::Minus
        }
        let subject_strand = if sstart <= send {
          Strand::Plus
        } else {
          Strand::Minus
        }

        let hsp = BlastHsp::new(
          if qstart <= qend {
            qstart
          } else {
            qend
          },
          if qstart <= qend {
            qend
          } else {
            qstart
          },
          if sstart <= send {
            sstart
          } else {
            send
          },
          if sstart <= send {
            send
          } else {
            sstart
          },
          query_strand,
          subject_strand,
          pident,
          (pident / 100.0 * align_length.to_double()).to_int(),
          gaps,
          align_length,
          bitscore,
          evalue,
          "",
          "",
          "",
        )

        let mut found = false
        let mut j = 0
        while j < hits.length() {
          if hits[j].subject_id() == sseqid {
            let current_hsps = hits[j].hsps()
            let new_hsps = current_hsps.copy()
            new_hsps.push(hsp)
            hits[j] = BlastHit::new(
              hits[j].subject_id(),
              hits[j].subject_def(),
              hits[j].subject_length(),
              new_hsps,
            )
            found = true
            break
          }
          j = j + 1
        }

        if !found {
          let hsps : Array[BlastHsp] = Array::new()
          hsps.push(hsp)
          hits.push(BlastHit::new(sseqid, sseqid, 0, hsps))
        }
      }
    }

    i = i + 1
  }

  if hits.length() == 0 {
    return None
  }

  Some(
    BlastRecord::new(
      query_id, query_def, query_length, hits, "blast", "unknown", "unknown",
    ),
  )
}

///|
/// Parse BLAST XML format (basic parser for standard NCBI BLAST XML output).
pub fn parse_blast_xml(content : String) -> BlastRecord? {
  let lines = split_by_char(content, '\n'.to_int())
  let mut i = 0
  let mut query_id = "unknown"
  let mut query_def = ""
  let mut query_length = 0
  let hits : Array[BlastHit] = []
  let mut in_hit = false
  let mut in_hsp = false
  let mut current_hit_id = ""
  let mut current_hit_def = ""
  let mut current_hit_length = 0
  let mut current_hsps : Array[BlastHsp] = []
  let mut hsp_fields = Map([], capacity=20)

  while i < lines.length() {
    let line = lines[i].trim().to_owned()
    if line.length() == 0 {
      i = i + 1
      continue
    }

    if starts_with(line, 0, "<Query_id>") {
      query_id = extract_between(line, "<Query_id>", "</Query_id>")
    } else if starts_with(line, 0, "<Query_def>") {
      query_def = extract_between(line, "<Query_def>", "</Query_def>")
    } else if starts_with(line, 0, "<Query_len>") {
      let len_str = extract_between(line, "<Query_len>", "</Query_len>")
      let len_opt = parse_double(len_str)
      if len_opt.is_some() {
        query_length = len_opt.unwrap().to_int()
      }
    } else if starts_with(line, 0, "<Hit>") {
      in_hit = true
      current_hsps = []
    } else if starts_with(line, 0, "</Hit>") {
      if in_hit && current_hit_id.length() > 0 {
        hits.push(
          BlastHit::new(
            current_hit_id, current_hit_def, current_hit_length, current_hsps,
          ),
        )
      }
      in_hit = false
      current_hit_id = ""
    } else if in_hit && starts_with(line, 0, "<Hit_id>") {
      current_hit_id = extract_between(line, "<Hit_id>", "</Hit_id>")
    } else if in_hit && starts_with(line, 0, "<Hit_def>") {
      current_hit_def = extract_between(line, "<Hit_def>", "</Hit_def>")
    } else if in_hit && starts_with(line, 0, "<Hit_len>") {
      let len_str = extract_between(line, "<Hit_len>", "</Hit_len>")
      let len_opt = parse_double(len_str)
      if len_opt.is_some() {
        current_hit_length = len_opt.unwrap().to_int()
      }
    } else if in_hit && starts_with(line, 0, "<Hsp>") {
      in_hsp = true
      hsp_fields = Map([], capacity=20)
    } else if in_hit && in_hsp && starts_with(line, 0, "</Hsp>") {
      let hsp = build_hsp_from_fields(hsp_fields)
      if hsp.is_some() {
        current_hsps.push(hsp.unwrap())
      }
      in_hsp = false
    } else if in_hit && in_hsp {
      let key_val = parse_xml_tag(line)
      if key_val.is_some() {
        let (key, val) = key_val.unwrap()
        hsp_fields[key] = val
      }
    }

    i = i + 1
  }

  if query_id == "unknown" && hits.length() == 0 {
    return None
  }

  Some(
    BlastRecord::new(
      query_id, query_def, query_length, hits, "blast", "unknown", "unknown",
    ),
  )
}

///|
/// Extract text between two delimiters.
fn extract_between(s : String, start : String, end : String) -> String {
  let start_idx = s.find(start)
  if start_idx.is_none() {
    return ""
  }
  let start_pos = start_idx.unwrap() + start.length()
  let end_part = s[start_pos:].to_string()
  let end_idx = end_part.find(end)
  if end_idx.is_none() {
    return ""
  }
  end_part[0:end_idx.unwrap()].to_string()
}

///|
/// Parse an XML tag into (key, value).
fn parse_xml_tag(line : String) -> (String, String)? {
  let start_idx = line.find("<")
  let end_open_idx = line.find(">")
  let close_start_idx = line.find("</")
  let close_end_idx = line.find(">")

  if start_idx.is_none() ||
    end_open_idx.is_none() ||
    close_start_idx.is_none() ||
    close_end_idx.is_none() {
    return None
  }

  let tag_name = line[start_idx.unwrap() + 1:end_open_idx.unwrap()].to_string()
  let value = line[end_open_idx.unwrap() + 1:close_start_idx.unwrap()].to_string()
  Some((tag_name, value))
}

///|
/// Build an HSP from parsed XML fields.
fn build_hsp_from_fields(fields : Map[String, String]) -> BlastHsp? {
  let get_int = fn(key : String) -> Int {
    if fields.contains(key) {
      let opt = parse_double(fields[key])
      if opt.is_some() {
        opt.unwrap().to_int()
      } else {
        0
      }
    } else {
      0
    }
  }

  let get_double = fn(key : String) -> Double {
    if fields.contains(key) {
      let opt = parse_double(fields[key])
      if opt.is_some() {
        opt.unwrap()
      } else {
        0.0
      }
    } else {
      0.0
    }
  }

  let qstart = get_int("Hsp_query-from")
  let qend = get_int("Hsp_query-to")
  let sstart = get_int("Hsp_hit-from")
  let send = get_int("Hsp_hit-to")

  let query_strand = if qstart <= qend { Strand::Plus } else { Strand::Minus }
  let subject_strand = if sstart <= send { Strand::Plus } else { Strand::Minus }

  Some(
    BlastHsp::new(
      if qstart <= qend {
        qstart
      } else {
        qend
      },
      if qstart <= qend {
        qend
      } else {
        qstart
      },
      if sstart <= send {
        sstart
      } else {
        send
      },
      if sstart <= send {
        send
      } else {
        sstart
      },
      query_strand,
      subject_strand,
      get_double("Hsp_identity") / get_double("Hsp_align-len") * 100.0,
      get_int("Hsp_positive"),
      get_int("Hsp_gaps"),
      get_int("Hsp_align-len"),
      get_double("Hsp_bit-score"),
      get_double("Hsp_evalue"),
      fields.get("Hsp_qseq").unwrap_or(""),
      fields.get("Hsp_hseq").unwrap_or(""),
      fields.get("Hsp_midline").unwrap_or(""),
    ),
  )
}

///|
/// Format a BlastRecord as a human-readable string.
pub fn BlastRecord::to_string(self : BlastRecord) -> String {
  let sb = StringBuilder::new()
  sb.write_string("Query: " + self.query_id())
  if self.query_def().length() > 0 {
    sb.write_string(" (" + self.query_def() + ")")
  }
  sb.write_string("\nLength: " + self.query_length().to_string())
  sb.write_string("\nProgram: " + self.program())
  sb.write_string("\nVersion: " + self.version())
  sb.write_string("\nDatabase: " + self.database())
  sb.write_string("\n\nHits (" + self.num_hits().to_string() + "):\n")

  let mut i = 0
  while i < self.hits.length() {
    let hit = self.hits[i]
    sb.write_string("  " + (i + 1).to_string() + ". " + hit.subject_id())
    if hit.subject_def().length() > 0 {
      sb.write_string(" (" + hit.subject_def() + ")")
    }
    sb.write_string("\n     Length: " + hit.subject_length().to_string())
    sb.write_string(", HSPs: " + hit.num_hsps().to_string())

    let mut j = 0
    while j < hit.hsps().length() {
      let hsp = hit.hsps()[j]
      sb.write_string("\n     HSP " + (j + 1).to_string() + ":")
      sb.write_string(" bitscore=" + hsp.bitscore().to_string())
      sb.write_string(", evalue=" + hsp.evalue().to_string())
      sb.write_string(", identity=" + hsp.identity().to_string() + "%")
      sb.write_string(
        "\n         Query: " +
        hsp.query_start().to_string() +
        "-" +
        hsp.query_end().to_string(),
      )
      sb.write_string(" (" + strand_to_string(hsp.query_strand()) + ")")
      sb.write_string(
        "\n         Subject: " +
        hsp.subject_start().to_string() +
        "-" +
        hsp.subject_end().to_string(),
      )
      sb.write_string(" (" + strand_to_string(hsp.subject_strand()) + ")")
      j = j + 1
    }
    sb.write_string("\n")
    i = i + 1
  }

  sb.to_string()
}

///|
/// Convert strand to string.
fn strand_to_string(s : Strand) -> String {
  match s {
    Strand::Plus => "+"
    Strand::Minus => "-"
    Strand::Star => "?"
  }
}
