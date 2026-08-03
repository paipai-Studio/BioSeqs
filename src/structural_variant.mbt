// Bioconductor StructuralVariantAnnotation - Structural variant annotation.
//
// Inspired by the Bioconductor StructuralVariantAnnotation package
// (Cmero et al.), which provides infrastructure for parsing, classifying,
// and annotating structural variants (SVs) from VCF files. This pure-MoonBit
// module implements:
// - VCF breakend (BND) ALT field parsing
// - SV type inference (DEL, DUP, INV, INS, TRA, BND) from INFO fields
// - Partner breakend pairing (matching BND mates)
// - Gene/exon overlap annotation using genomic intervals
// - SV filtering and summarization
//
// Reference:
//   https://bioconductor.org/packages/StructuralVariantAnnotation/

// ============================================================================
// SV type enum
// ============================================================================

///|
/// The type of a structural variant.
pub enum SvType {
  DEL // deletion
  DUP // duplication
  INV // inversion
  INS // insertion
  TRA // translocation (inter-chromosomal)
  BND // breakend (unresolved junction)
}

///|
/// Construct a DEL type tag.
pub fn sv_type_del() -> SvType {
  SvType::DEL
}

///|
/// Construct a DUP type tag.
pub fn sv_type_dup() -> SvType {
  SvType::DUP
}

///|
/// Construct an INV type tag.
pub fn sv_type_inv() -> SvType {
  SvType::INV
}

///|
/// Construct an INS type tag.
pub fn sv_type_ins() -> SvType {
  SvType::INS
}

///|
/// Construct a TRA type tag.
pub fn sv_type_tra() -> SvType {
  SvType::TRA
}

///|
/// Construct a BND type tag.
pub fn sv_type_bnd() -> SvType {
  SvType::BND
}

///|
/// Parse a type string into an SvType. Returns BND for unknown types.
pub fn sv_parse_type(type_str : String) -> SvType {
  let up = type_str.to_upper()
  match up {
    "DEL" => SvType::DEL
    "DUP" => SvType::DUP
    "INV" => SvType::INV
    "INS" => SvType::INS
    "TRA" => SvType::TRA
    "BND" => SvType::BND
    _ => SvType::BND
  }
}

///|
/// String label for an SvType.
pub fn sv_type_label(t : SvType) -> String {
  match t {
    DEL => "DEL"
    DUP => "DUP"
    INV => "INV"
    INS => "INS"
    TRA => "TRA"
    BND => "BND"
  }
}

// ============================================================================
// Breakend representation
// ============================================================================

///|
/// A breakend: a single end of a structural variant junction. The breakend
/// has a position (chrom, pos) and an optional mate location
/// (mate_chrom, mate_pos) for paired BND records.
pub struct SvBreakend {
  chrom : String
  pos : Int
  strand : String // "+" or "-"
  mate_chrom : String // empty if no mate
  mate_pos : Int // -1 if no mate
  mate_strand : String // "+" or "-" or "" if no mate
  mate_id : String // ID of the partner breakend (empty if unknown)
}

///|
/// Construct an SvBreakend with mate information.
pub fn SvBreakend::new(
  chrom : String,
  pos : Int,
  strand : String,
  mate_chrom : String,
  mate_pos : Int,
  mate_strand : String
) -> SvBreakend {
  SvBreakend::{ chrom, pos, strand, mate_chrom, mate_pos, mate_strand, mate_id: "" }
}

///|
pub fn SvBreakend::chrom(self : SvBreakend) -> String {
  self.chrom
}

///|
pub fn SvBreakend::pos(self : SvBreakend) -> Int {
  self.pos
}

///|
pub fn SvBreakend::strand(self : SvBreakend) -> String {
  self.strand
}

///|
pub fn SvBreakend::mate_chrom(self : SvBreakend) -> String {
  self.mate_chrom
}

///|
pub fn SvBreakend::mate_pos(self : SvBreakend) -> Int {
  self.mate_pos
}

///|
pub fn SvBreakend::mate_strand(self : SvBreakend) -> String {
  self.mate_strand
}

///|
pub fn SvBreakend::mate_id(self : SvBreakend) -> String {
  self.mate_id
}

///|
/// Whether this breakend has a mate (i.e., is part of a paired BND).
pub fn SvBreakend::has_mate(self : SvBreakend) -> Bool {
  self.mate_chrom.length() > 0
}

///|
/// Whether this breakend is inter-chromosomal (mate on a different
/// chromosome).
pub fn SvBreakend::is_inter_chromosomal(self : SvBreakend) -> Bool {
  self.has_mate() && self.mate_chrom != self.chrom
}

// ============================================================================
// SV record
// ============================================================================

///|
/// A structural variant record parsed from a VCF line. Contains the SV
/// type, position, end coordinate, SV length, and (for BND records) the
/// breakend information.
pub struct SvRecord {
  id : String
  chrom : String
  pos : Int
  sv_type : SvType
  end : Int // END coordinate from INFO
  svlen : Int // SVLEN from INFO (negative for DEL, positive for DUP)
  breakend : SvBreakend? // present for BND records
  quality : Double // QUAL field (NaN if missing)
  filter : String // FILTER field
}

///|
/// Construct an SvRecord.
pub fn SvRecord::new(
  id : String,
  chrom : String,
  pos : Int,
  sv_type : SvType,
  end_pos : Int,
  svlen : Int,
  breakend : SvBreakend?,
  quality : Double,
  filter : String
) -> SvRecord {
  SvRecord::{ id, chrom, pos, sv_type, end: end_pos, svlen, breakend, quality, filter }
}

///|
pub fn SvRecord::id(self : SvRecord) -> String {
  self.id
}

///|
pub fn SvRecord::chrom(self : SvRecord) -> String {
  self.chrom
}

///|
pub fn SvRecord::pos(self : SvRecord) -> Int {
  self.pos
}

///|
pub fn SvRecord::sv_type(self : SvRecord) -> SvType {
  self.sv_type
}

///|
pub fn SvRecord::end(self : SvRecord) -> Int {
  self.end
}

///|
pub fn SvRecord::svlen(self : SvRecord) -> Int {
  self.svlen
}

///|
pub fn SvRecord::breakend(self : SvRecord) -> SvBreakend? {
  self.breakend
}

///|
pub fn SvRecord::quality(self : SvRecord) -> Double {
  self.quality
}

///|
pub fn SvRecord::filter(self : SvRecord) -> String {
  self.filter
}

///|
/// Whether this is a BND (breakend) record.
pub fn SvRecord::is_bnd(self : SvRecord) -> Bool {
  self.sv_type is BND
}

///|
/// Whether this SV is inter-chromosomal (BND with mate on a different
/// chromosome, or TRA type).
pub fn SvRecord::is_inter_chromosomal(self : SvRecord) -> Bool {
  match self.breakend {
    Some(b) => b.is_inter_chromosomal()
    None => self.sv_type is TRA
  }
}

///|
/// Size of the SV (absolute value of svlen). Returns 0 for BND/TRA
/// records without a defined length.
pub fn SvRecord::size(self : SvRecord) -> Int {
  self.svlen.abs()
}

// ============================================================================
// BND ALT field parsing
// ============================================================================

///|
/// Trim leading/trailing whitespace from a string.
fn sv_trim(s : String) -> String {
  let n = s.length()
  let mut start = 0
  while start < n {
    let ch = s[start : start + 1]
    if ch == " " || ch == "\t" || ch == "\r" || ch == "\n" {
      start = start + 1
    } else {
      break
    }
  }
  let mut end = n
  while end > start {
    let ch = s[end - 1 : end]
    if ch == " " || ch == "\t" || ch == "\r" || ch == "\n" {
      end = end - 1
    } else {
      break
    }
  }
  s[start:end].to_owned()
}

///|
/// Reverse find: return the index of the last occurrence of `ch` in `s`.
fn sv_rfind(s : String, ch : String) -> Int? {
  let n = s.length()
  let mut i = n - 1
  while i >= 0 {
    if s[i : i + 1] == ch {
      return Some(i)
    }
    i = i - 1
  }
  None
}

///|
/// Parse a BND ALT field into a breakend.
///
/// The BND ALT notation (from the VCF 4.2 specification) encodes a breakend
/// junction:
///   `t[p[`  — extends right to join with position p (forward strand)
///   `t]p]`  — extends right to join with position p (reverse strand)
///   `]p]t`  — extends left to join with position p (reverse strand)
///   `[p[t`  — extends left to join with position p (forward strand)
/// where `t` is the joining nucleotide and `p` is `chrom:pos`.
///
/// Returns Some(SvBreakend) if parsing succeeds, None otherwise.
pub fn sv_parse_bnd_alt(
  alt : String,
  chrom : String,
  pos : Int
) -> SvBreakend? {
  let trimmed = sv_trim(alt)
  let n = trimmed.length()
  if n < 4 {
    return None
  }
  // Determine orientation by examining brackets.
  // Case 1: t[p[ — starts with a base, then [, then chrom:pos, then [
  // Case 2: t]p] — starts with a base, then ], then chrom:pos, then ]
  // Case 3: ]p]t — starts with ], then chrom:pos, then ], then base
  // Case 4: [p[t — starts with [, then chrom:pos, then [, then base
  let first = trimmed[0 : 1]
  if first == "[" || first == "]" {
    // Case 3 or 4: ]p]t or [p[t
    let open_bracket = first
    // Find the matching closing bracket (the second occurrence).
    // open_bracket is at index 0; search for the next occurrence starting at 1.
    let mut close_idx = -1
    for k in 1..<n {
      if trimmed[k : k + 1] == open_bracket {
        close_idx = k
        break
      }
    }
    if close_idx == -1 {
      return None
    }
    // close_idx is the index of the second occurrence of the bracket.
    // The mate location is between index 1 and close_idx.
    if close_idx <= 1 {
      return None
    }
    let mate_loc = trimmed[1 : close_idx]
    let (mate_chrom, mate_pos) = sv_parse_location(mate_loc.to_owned())
    if mate_chrom.length() == 0 {
      return None
    }
    let strand = if open_bracket == "[" { "+" } else { "-" }
    let mate_strand = if open_bracket == "[" { "+" } else { "-" }
    Some(SvBreakend::new(chrom, pos, strand, mate_chrom, mate_pos, mate_strand))
  } else {
    // Case 1 or 2: t[p[ or t]p]
    // The joining base is the first character.
    let bracket_char = trimmed[1 : 2]
    if bracket_char != "[" && bracket_char != "]" {
      return None
    }
    // Find the matching closing bracket (the second occurrence).
    // bracket_char is at index 1; search for the next occurrence starting at 2.
    let mut close_idx = -1
    for k in 2..<n {
      if trimmed[k : k + 1] == bracket_char {
        close_idx = k
        break
      }
    }
    if close_idx == -1 {
      return None
    }
    // The mate location is between index 2 and close_idx.
    if close_idx <= 2 {
      return None
    }
    let mate_loc = trimmed[2 : close_idx]
    let (mate_chrom, mate_pos) = sv_parse_location(mate_loc.to_owned())
    if mate_chrom.length() == 0 {
      return None
    }
    let strand = if bracket_char == "[" { "+" } else { "-" }
    let mate_strand = if bracket_char == "[" { "+" } else { "-" }
    Some(SvBreakend::new(chrom, pos, strand, mate_chrom, mate_pos, mate_strand))
  }
}

///|
/// Parse a "chrom:pos" location string.
fn sv_parse_location(loc : String) -> (String, Int) {
  // Find the last colon (chrom may contain colons in some notations,
  // but typically not).
  let colon_pos = sv_rfind(loc, ":")
  match colon_pos {
    Some(idx) => {
      let chrom = loc[0 : idx].to_owned()
      let pos_str = loc[idx + 1 :]
      let pos = match parse_double(pos_str.to_owned()) {
        Some(v) => v.to_int()
        None => 0
      }
      (chrom, pos)
    }
    None => ("", 0)
  }
}

// ============================================================================
// INFO field parsing
// ============================================================================

///|
/// Parse a VCF INFO field string into a Map.
/// INFO format: "KEY1=VALUE1;KEY2=VALUE2;KEY3;..."
/// Keys without "=" are treated as flags with value "true".
fn sv_parse_info(info : String) -> Map[String, String] {
  let result : Map[String, String] = Map([], capacity=16)
  let fields = info.split(";").to_array()
  for field in fields {
    let eq_pos = field.find("=")
    match eq_pos {
      Some(idx) => {
        let key = field[0 : idx].to_owned()
        let value = field[idx + 1 :].to_owned()
        result[key] = value
      }
      None => result[field.to_owned()] = "true"
    }
  }
  result
}

// ============================================================================
// VCF line parsing
// ============================================================================

///|
/// Parse a single VCF data line into an SvRecord.
///
/// The VCF line format (tab-separated):
///   CHROM POS ID REF ALT QUAL FILTER INFO
/// For SV records, the ALT field is either:
///   - An SV type string: <DEL>, <DUP>, <INV>, <INS>, <TRA>
///   - A BND notation: G]chr:pos], C[chr:pos[, etc.
/// The INFO field must contain SVTYPE, and may contain END, SVLEN.
pub fn sv_parse_vcf_line(line : String) -> SvRecord? {
  let trimmed = sv_trim(line)
  if trimmed.length() == 0 {
    return None
  }
  // Skip header lines.
  if trimmed[0 : 1] == "#" {
    return None
  }
  let fields = trimmed.split("\t").to_array()
  if fields.length() < 8 {
    return None
  }
  let chrom = fields[0].to_owned()
  let pos = match parse_double(fields[1].to_owned()) {
    Some(v) => v.to_int()
    None => return None
  }
  let id = fields[2].to_owned()
  let alt = fields[4].to_owned()
  let qual_str = fields[5].to_owned()
  let filter = fields[6].to_owned()
  let info = fields[7].to_owned()
  let quality = match parse_double(qual_str) {
    Some(v) => v
    None => 0.0
  }
  // Parse INFO.
  let info_map = sv_parse_info(info)
  // Determine SV type.
  let svtype_str = info_map.get_or_default("SVTYPE", "BND")
  let sv_type = sv_parse_type(svtype_str)
  // Parse END.
  let end_pos = match info_map.get("END") {
    Some(v) =>
      match parse_double(v) {
        Some(d) => d.to_int()
        None => pos
      }
    None => pos
  }
  // Parse SVLEN.
  let svlen = match info_map.get("SVLEN") {
    Some(v) =>
      match parse_double(v) {
        Some(d) => d.to_int()
        None => 0
      }
    None => 0
  }
  // Parse BND breakend if applicable.
  let breakend = if sv_type is BND || alt.contains("[") || alt.contains("]") {
    sv_parse_bnd_alt(alt, chrom, pos)
  } else {
    None
  }
  Some(
    SvRecord::new(id, chrom, pos, sv_type, end_pos, svlen, breakend, quality, filter),
  )
}

///|
/// Parse multiple VCF lines into an array of SvRecords. Header lines
/// (starting with #) are skipped.
pub fn sv_parse_vcf(content : String) -> Array[SvRecord] {
  let records = Array::new()
  let lines = content.split("\n").to_array()
  for line in lines {
    match sv_parse_vcf_line(line.to_owned()) {
      Some(rec) => records.push(rec)
      None => continue
    }
  }
  records
}

// ============================================================================
// Partner pairing
// ============================================================================

///|
/// Find partner breakends among a set of BND records. Two BND records are
/// partners if record A's breakend mate location matches record B's
/// location (and vice versa). Returns an array of (id_a, id_b) pairs.
pub fn sv_find_partners(records : Array[SvRecord]) -> Array[(String, String)] {
  let pairs : Array[(String, String)] = Array::new()
  let n = records.length()
  let matched : Array[Bool] = Array::make(n, false)
  for i in 0..<n {
    if matched[i] {
      continue
    }
    let ri = records[i]
    match ri.breakend {
      Some(bi) => {
        // Look for a record j whose chrom/pos matches bi's mate location.
        for j in (i + 1)..<n {
          if matched[j] {
            continue
          }
          let rj = records[j]
          match rj.breakend {
            Some(bj) => {
              // Check if i's mate matches j's location and vice versa.
              if bi.mate_chrom == rj.chrom &&
                bi.mate_pos == rj.pos &&
                bj.mate_chrom == ri.chrom &&
                bj.mate_pos == ri.pos {
                pairs.push((ri.id, rj.id))
                matched[i] = true
                matched[j] = true
                break
              }
            }
            None => continue
          }
        }
      }
      None => continue
    }
  }
  pairs
}

///|
/// Check whether two SV records are partners (mutually-referencing BND
/// breakends).
pub fn sv_are_partners(a : SvRecord, b : SvRecord) -> Bool {
  match (a.breakend, b.breakend) {
    (Some(ba), Some(bb)) =>
      ba.mate_chrom == b.chrom &&
      ba.mate_pos == b.pos &&
      bb.mate_chrom == a.chrom &&
      bb.mate_pos == a.pos
    _ => false
  }
}

// ============================================================================
// Gene region annotation
// ============================================================================

///|
/// A gene region for overlap annotation.
pub struct SvGeneRegion {
  gene : String
  chrom : String
  start : Int
  end : Int
  strand : String
}

///|
/// Construct an SvGeneRegion.
pub fn SvGeneRegion::new(
  gene : String,
  chrom : String,
  start : Int,
  end : Int,
  strand : String
) -> SvGeneRegion {
  SvGeneRegion::{ gene, chrom, start, end, strand }
}

///|
pub fn SvGeneRegion::gene(self : SvGeneRegion) -> String {
  self.gene
}

///|
pub fn SvGeneRegion::chrom(self : SvGeneRegion) -> String {
  self.chrom
}

///|
pub fn SvGeneRegion::start(self : SvGeneRegion) -> Int {
  self.start
}

///|
pub fn SvGeneRegion::end(self : SvGeneRegion) -> Int {
  self.end
}

///|
pub fn SvGeneRegion::strand(self : SvGeneRegion) -> String {
  self.strand
}

///|
/// Check if a position interval [qstart, qend] on `qchrom` overlaps this
/// gene region.
fn SvGeneRegion::overlaps(
  self : SvGeneRegion,
  qchrom : String,
  qstart : Int,
  qend : Int
) -> Bool {
  self.chrom == qchrom && self.start <= qend && self.end >= qstart
}

// ============================================================================
// SV annotation
// ============================================================================

///|
/// An annotated SV: the original record plus overlapping genes.
pub struct SvAnnotation {
  record : SvRecord
  genes : Array[String] // names of overlapping genes
}

///|
/// Construct an SvAnnotation.
pub fn SvAnnotation::new(
  record : SvRecord,
  genes : Array[String]
) -> SvAnnotation {
  SvAnnotation::{ record, genes }
}

///|
pub fn SvAnnotation::record(self : SvAnnotation) -> SvRecord {
  self.record
}

///|
pub fn SvAnnotation::genes(self : SvAnnotation) -> Array[String] {
  self.genes
}

///|
/// Whether the SV overlaps any gene.
pub fn SvAnnotation::has_gene_overlap(self : SvAnnotation) -> Bool {
  self.genes.length() > 0
}

///|
/// Annotate SV records with overlapping genes.
///
/// For each SV, the interval [pos, end] (or [pos, mate_pos] for BND) is
/// checked against all gene regions. Genes that overlap the SV interval
/// are listed in the annotation.
pub fn sv_annotate(
  records : Array[SvRecord],
  genes : Array[SvGeneRegion]
) -> Array[SvAnnotation] {
  let annotations : Array[SvAnnotation] = Array::new()
  for rec in records {
    let sv_start = rec.pos
    let sv_end = if rec.end > rec.pos {
      rec.end
    } else {
      match rec.breakend {
        Some(b) => if b.mate_pos > sv_start { b.mate_pos } else { sv_start }
        None => sv_start
      }
    }
    let sv_chrom = rec.chrom
    let overlapping : Array[String] = Array::new()
    for gene in genes {
      if gene.overlaps(sv_chrom, sv_start, sv_end) {
        overlapping.push(gene.gene)
      }
    }
    // Also check mate chromosome for BND records.
    match rec.breakend {
      Some(b) => {
        if b.mate_chrom.length() > 0 && b.mate_chrom != sv_chrom {
          for gene in genes {
            if gene.chrom == b.mate_chrom &&
              gene.start <= b.mate_pos &&
              gene.end >= b.mate_pos &&
              !overlapping.contains(gene.gene) {
              overlapping.push(gene.gene)
            }
          }
        }
      }
      None => ()
    }
    annotations.push(SvAnnotation::new(rec, overlapping))
  }
  annotations
}

// ============================================================================
// Filtering
// ============================================================================

///|
/// Filter SV records by minimum quality score. Returns records with
/// quality >= min_qual.
pub fn sv_filter_by_quality(
  records : Array[SvRecord],
  min_qual : Double
) -> Array[SvRecord] {
  let result : Array[SvRecord] = Array::new()
  for rec in records {
    if rec.quality >= min_qual {
      result.push(rec)
    }
  }
  result
}

///|
/// Filter SV records by minimum size. Returns records with
/// size >= min_size. BND/TRA records (size 0) are always excluded.
pub fn sv_filter_by_size(
  records : Array[SvRecord],
  min_size : Int
) -> Array[SvRecord] {
  let result : Array[SvRecord] = Array::new()
  for rec in records {
    if rec.size() >= min_size {
      result.push(rec)
    }
  }
  result
}

///|
/// Filter SV records by type. Returns records matching the given type.
pub fn sv_filter_by_type(
  records : Array[SvRecord],
  t : SvType
) -> Array[SvRecord] {
  let result : Array[SvRecord] = Array::new()
  for rec in records {
    if sv_type_label(rec.sv_type) == sv_type_label(t) {
      result.push(rec)
    }
  }
  result
}

// ============================================================================
// Summary
// ============================================================================

///|
/// Count SV records by type. Returns a map from type label to count.
pub fn sv_count_by_type(
  records : Array[SvRecord]
) -> Map[String, Int] {
  let counts : Map[String, Int] = Map([], capacity=6)
  for rec in records {
    let label = sv_type_label(rec.sv_type)
    counts[label] = counts.get_or_default(label, 0) + 1
  }
  counts
}

///|
/// Text summary of a single SV record.
pub fn sv_record_summary(rec : SvRecord) -> String {
  let sb = StringBuilder::new()
  sb.write_string("ID: ")
  sb.write_string(rec.id)
  sb.write_char('\n')
  sb.write_string("  Type: ")
  sb.write_string(sv_type_label(rec.sv_type))
  sb.write_char('\n')
  sb.write_string("  Location: ")
  sb.write_string(rec.chrom)
  sb.write_char(':')
  sb.write_string(rec.pos.to_string())
  if rec.end != rec.pos {
    sb.write_char('-')
    sb.write_string(rec.end.to_string())
  }
  sb.write_char('\n')
  if rec.svlen != 0 {
    sb.write_string("  SVLEN: ")
    sb.write_string(rec.svlen.to_string())
    sb.write_char('\n')
  }
  sb.write_string("  Quality: ")
  sb.write_string(rec.quality.to_string())
  sb.write_char('\n')
  sb.write_string("  Filter: ")
  sb.write_string(rec.filter)
  sb.write_char('\n')
  if rec.is_inter_chromosomal() {
    sb.write_string("  Inter-chromosomal: true\n")
  }
  match rec.breakend {
    Some(b) => {
      sb.write_string("  Breakend mate: ")
      sb.write_string(b.mate_chrom)
      sb.write_char(':')
      sb.write_string(b.mate_pos.to_string())
      sb.write_string(" (")
      sb.write_string(b.mate_strand)
      sb.write_string(")\n")
    }
    None => ()
  }
  sb.to_string()
}

///|
/// Text summary of an SV annotation.
pub fn sv_annotation_summary(ann : SvAnnotation) -> String {
  let sb = StringBuilder::new()
  sb.write_string(sv_record_summary(ann.record))
  if ann.genes.length() > 0 {
    sb.write_string("  Overlapping genes: ")
    for i in 0..<ann.genes.length() {
      if i > 0 {
        sb.write_string(", ")
      }
      sb.write_string(ann.genes[i])
    }
    sb.write_char('\n')
  } else {
    sb.write_string("  Overlapping genes: (none)\n")
  }
  sb.to_string()
}

// ============================================================================
// Sample data
// ============================================================================

///|
/// Sample VCF text with structural variants: a DEL, DUP, INV, INS, and
/// two paired BND records (translocation).
pub fn sv_sample_vcf() -> String {
  let sb = StringBuilder::new()
  let lines = [
    "##fileformat=VCFv4.2",
    "##fileDate=20260101",
    "##ALT=<ID=DEL,Description=\"Deletion\">",
    "##ALT=<ID=DUP,Description=\"Duplication\">",
    "##ALT=<ID=INV,Description=\"Inversion\">",
    "##ALT=<ID=INS,Description=\"Insertion\">",
    "##ALT=<ID=BND,Description=\"Breakend\">",
    "##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type of structural variant\">",
    "##INFO=<ID=END,Number=1,Type=Integer,Description=\"End position\">",
    "##INFO=<ID=SVLEN,Number=1,Type=Integer,Description=\"Length of SV\">",
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
    "chr1\t1000\tsv_del_1\tA\t<DEL>\t500.0\tPASS\tSVTYPE=DEL;END=2000;SVLEN=-1000",
    "chr1\t3000\tsv_dup_1\tA\t<DUP>\t450.0\tPASS\tSVTYPE=DUP;END=5000;SVLEN=2000",
    "chr2\t1500\tsv_inv_1\tA\t<INV>\t300.0\tPASS\tSVTYPE=INV;END=2500;SVLEN=1000",
    "chr2\t4000\tsv_ins_1\tA\t<INS>\t200.0\tLOWQUAL\tSVTYPE=INS;SVLEN=500",
    "chr1\t6000\tsv_bnd_1\tA\tA]chr2:7000]\t600.0\tPASS\tSVTYPE=BND",
    "chr2\t7000\tsv_bnd_2\tA\tA]chr1:6000]\t600.0\tPASS\tSVTYPE=BND",
  ]
  for line in lines {
    sb.write_string(line)
    sb.write_char('\n')
  }
  sb.to_string()
}

///|
/// Sample gene regions for annotation testing.
pub fn sv_sample_genes() -> Array[SvGeneRegion] {
  [
    SvGeneRegion::new("BRCA1", "chr1", 800, 2500, "+"),
    SvGeneRegion::new("BRCA2", "chr1", 2800, 5200, "-"),
    SvGeneRegion::new("TP53", "chr2", 1200, 2800, "+"),
    SvGeneRegion::new("EGFR", "chr2", 5500, 8000, "+"),
    SvGeneRegion::new("KRAS", "chr1", 5800, 6500, "-"),
  ]
}

///|
/// Sample SV records parsed from the sample VCF.
pub fn sv_sample_records() -> Array[SvRecord] {
  sv_parse_vcf(sv_sample_vcf())
}
