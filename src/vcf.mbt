///|
/// VCF format parsing module, matching pysam VCF specification.

///|
/// VCF meta header record (## lines).
pub struct VcfMetaHeader {
  key : String
  value : String
  tags : Map[String, String]
}

///|
/// VCF fixed header record (# line).
pub struct VcfFixedHeader {
  columns : Array[String]
}

///|
/// VCF record.
pub struct VcfRecord {
  chrom : String
  pos : Int32
  mut id : Array[String]
  ref_seq : String
  alt : Array[String]
  mut qual : Double?
  mut filter : Array[String]
  info : Map[String, String]
  format : Array[String]?
  samples : Map[String, Array[String]]
}

///|
/// VCF file representation.
pub struct VcfFile {
  meta_headers : Array[VcfMetaHeader]
  fixed_header : VcfFixedHeader?
  records : Array[VcfRecord]
}

///|
/// Parse INFO field into HashMap.
fn parse_info(info_str : String) -> Map[String, String] {
  let result = Map([], capacity=0)
  if info_str == "." || info_str.length() == 0 {
    return result
  }

  let entries = parse_semicolon_list(info_str)
  let mut i = 0
  while i < entries.length() {
    let entry = entries[i]
    let eq_pos = entry.find("=")
    if eq_pos is Some(_) {
      let pos = eq_pos.unwrap()
      let key = entry[0:pos].to_owned()
      let value = if pos + 1 < entry.length() {
        entry[pos + 1:entry.length()].to_owned()
      } else {
        ""
      }
      result.set(key, value)
    } else {
      result.set(entry, "")
    }
    i = i + 1
  }

  result
}

///|
/// Parse FORMAT field into array.
fn parse_format(format_str : String) -> Array[String]? {
  if format_str == "." || format_str.length() == 0 {
    return None
  }
  Some(parse_comma_list(format_str))
}

///|
/// Parse VCF meta header line (##...).
fn parse_vcf_meta_header(line : String) -> VcfMetaHeader? {
  if line.length() < 3 ||
    line.unsafe_get(0).to_int() != '#'.to_int() ||
    line.unsafe_get(1).to_int() != '#'.to_int() {
    return None
  }

  let content = line[2:line.length()].to_owned()
  let eq_pos = content.find("=")
  if eq_pos is None {
    return None
  }

  let pos = eq_pos.unwrap()
  let key = content[0:pos].to_owned()
  let value_part = content[pos + 1:content.length()].to_owned()

  let tags = Map([], capacity=0)
  let mut value = value_part

  if value_part.length() > 0 &&
    value_part.unsafe_get(0).to_int() == '<'.to_int() {
    let n = value_part.length()
    if n > 1 && value_part.unsafe_get(n - 1).to_int() == '>'.to_int() {
      let inner = value_part[1:n - 1].to_owned()
      let fields = parse_comma_list(inner)
      let mut i = 0
      while i < fields.length() {
        let field = fields[i]
        let feq_pos = field.find("=")
        if feq_pos is Some(_) {
          let fpos = feq_pos.unwrap()
          let fkey = field[0:fpos].to_owned()
          let fvalue = if fpos + 1 < field.length() {
            field[fpos + 1:field.length()].to_owned()
          } else {
            ""
          }
          tags.set(fkey, fvalue)
        }
        i = i + 1
      }
      value = inner
    }
  }

  Some(VcfMetaHeader::{ key, value, tags, })
}

///|
/// Parse VCF fixed header line (#...).
fn parse_vcf_fixed_header(line : String) -> VcfFixedHeader? {
  if line.length() < 2 ||
    line.unsafe_get(0).to_int() != '#'.to_int() ||
    line.unsafe_get(1).to_int() == '#'.to_int() {
    return None
  }

  let content = line[1:line.length()].to_owned()
  let columns = split_by_tab(content)
  Some(VcfFixedHeader::{ columns, })
}

///|
/// Parse VCF record line.
fn parse_vcf_record(line : String, sample_names : Array[String]) -> VcfRecord? {
  let parts = split_by_tab(line)
  if parts.length() < 8 {
    return None
  }

  let chrom = parts[0]
  let pos = parse_int(parts[1])
  let id = parse_comma_list(parts[2])
  let ref_seq = parts[3]
  let alt = parse_comma_list(parts[4])
  let qual = parse_double(parts[5])
  let filter = parse_comma_list(parts[6])
  let info = parse_info(parts[7])

  let format = if parts.length() > 8 { parse_format(parts[8]) } else { None }

  let samples = Map([], capacity=0)
  if parts.length() > 9 && format is Some(_) {
    let _format_cols = format.unwrap()
    let mut i = 9
    while i < parts.length() {
      let sample_name = if i - 9 < sample_names.length() {
        sample_names[i - 9]
      } else {
        (i - 9).to_string()
      }
      let sample_data = parse_comma_list(parts[i])
      samples.set(sample_name, sample_data)
      i = i + 1
    }
  }

  Some(VcfRecord::{
    chrom,
    pos,
    id,
    ref_seq,
    alt,
    qual,
    filter,
    info,
    format,
    samples,
  })
}

///|
/// Parse VCF content from string.
pub fn parse_vcf(content : String) -> VcfFile {
  let meta_headers = Array(capacity=0)
  let mut fixed_header : VcfFixedHeader? = None
  let records = Array(capacity=0)

  let lines = content.split("\n").to_array()
  let mut i = 0
  let n = lines.length()
  while i < n {
    let line = lines[i].to_owned()
    let line_len = line.length()
    if line_len == 0 {
      i = i + 1
      continue
    }

    if line.unsafe_get(0).to_int() == '#'.to_int() {
      if line.length() > 1 && line.unsafe_get(1).to_int() == '#'.to_int() {
        let meta_opt = parse_vcf_meta_header(line)
        if meta_opt is Some(_) {
          meta_headers.push(meta_opt.unwrap())
        }
      } else {
        let fixed_opt = parse_vcf_fixed_header(line)
        if fixed_opt is Some(_) {
          fixed_header = fixed_opt
        }
      }
    } else {
      let sample_names = match fixed_header {
        Some(h) => {
          let names = Array(capacity=0)
          let mut j = 9
          while j < h.columns.length() {
            names.push(h.columns[j])
            j = j + 1
          }
          names
        }
        None => []
      }
      let record_opt = parse_vcf_record(line, sample_names)
      if record_opt is Some(_) {
        records.push(record_opt.unwrap())
      }
    }
    i = i + 1
  }

  VcfFile::{ meta_headers, fixed_header, records, }
}

///|
/// Get number of records in VCF file.
pub fn VcfFile::len(self : VcfFile) -> Int32 {
  self.records.length()
}

///|
/// Get number of meta headers.
pub fn VcfFile::num_meta_headers(self : VcfFile) -> Int32 {
  self.meta_headers.length()
}

///|
/// Get sample names from fixed header.
pub fn VcfFile::get_sample_names(self : VcfFile) -> Array[String] {
  match self.fixed_header {
    Some(h) => {
      let result = Array(capacity=0)
      let mut i = 9
      while i < h.columns.length() {
        result.push(h.columns[i])
        i = i + 1
      }
      result
    }
    None => []
  }
}

///|
/// Get VCF version from meta headers.
pub fn VcfFile::get_version(self : VcfFile) -> String? {
  let mut i = 0
  while i < self.meta_headers.length() {
    let header = self.meta_headers[i]
    if header.key == "fileformat" {
      return Some(header.value)
    }
    i = i + 1
  }
  None
}

///|
/// Check if variant is a SNP (single nucleotide polymorphism).
pub fn VcfRecord::is_snp(self : VcfRecord) -> Bool {
  self.ref_seq.length() == 1 &&
  self.alt.length() == 1 &&
  self.alt[0].length() == 1
}

///|
/// Check if variant is an insertion.
pub fn VcfRecord::is_insertion(self : VcfRecord) -> Bool {
  if self.alt.length() == 0 {
    return false
  }
  self.ref_seq.length() < self.alt[0].length()
}

///|
/// Check if variant is a deletion.
pub fn VcfRecord::is_deletion(self : VcfRecord) -> Bool {
  if self.alt.length() == 0 {
    return false
  }
  self.ref_seq.length() > self.alt[0].length()
}

///|
/// Get genotype for a sample.
pub fn VcfRecord::get_genotype(
  self : VcfRecord,
  sample_name : String,
) -> Array[String]? {
  match self.samples.get(sample_name) {
    Some(data) => Some(data.copy())
    None => None
  }
}

///|
fn join_str_array(arr : Array[String], sep : String) -> String {
  if arr.length() == 0 {
    return "."
  }
  let buf = StringBuilder()
  let mut first = true
  for v in arr {
    if !first {
      buf.write_string(sep)
    }
    buf.write_string(v)
    first = false
  }
  buf.to_string()
}

///|
fn vcf_format_qual(q : Double?) -> String {
  match q {
    None => "."
    Some(v) => v.to_string()
  }
}

///|
/// Serialize INFO back to `k=v;k2=v2` (flag fields without `=` become bare `k`).
fn vcf_format_info(info : Map[String, String]) -> String {
  if info.length() == 0 {
    return "."
  }
  let buf = StringBuilder()
  let mut first = true
  for k in info.keys() {
    if !first {
      buf.write_string(";")
    }
    buf.write_string(k)
    match info.get(k) {
      None => ()
      Some(v) =>
        if v != "" {
          buf.write_string("=")
          buf.write_string(v)
        }
    }
    first = false
  }
  buf.to_string()
}

///|
fn vcf_format_record(rec : VcfRecord, sample_names : Array[String]) -> String {
  let parts : Array[String] = Array(capacity=0)
  parts.push(rec.chrom)
  parts.push(rec.pos.to_string())
  parts.push(join_str_array(rec.id, ","))
  parts.push(rec.ref_seq)
  parts.push(join_str_array(rec.alt, ","))
  parts.push(vcf_format_qual(rec.qual))
  parts.push(join_str_array(rec.filter, ","))
  parts.push(vcf_format_info(rec.info))
  match rec.format {
    None => ()
    Some(fmt) => {
      parts.push(join_str_array(fmt, ","))
      for name in sample_names {
        match rec.samples.get(name) {
          None => parts.push(".")
          Some(data) => parts.push(join_str_array(data, ","))
        }
      }
    }
  }
  let buf = StringBuilder()
  let mut first = true
  for p in parts {
    if !first {
      buf.write_string("\t")
    }
    buf.write_string(p)
    first = false
  }
  buf.to_string()
}

///|
/// Serialize a `VcfFile` back to VCF text (`VariantFile.write` equivalent),
/// including `##` meta headers, the `#CHROM ...` fixed header, and records.
/// Round-trips with `parse_vcf`.
pub fn VcfFile::to_string(self : VcfFile) -> String {
  let lines : Array[String] = Array(capacity=0)
  for mh in self.meta_headers {
    let line = if mh.tags.length() == 0 {
      "##" + mh.key + "=" + mh.value
    } else {
      "##" + mh.key + "=<" + mh.value + ">"
    }
    lines.push(line)
  }
  match self.fixed_header {
    None => ()
    Some(fh) => {
      let mut hline = "#"
      let mut first = true
      for col in fh.columns {
        if !first {
          hline = hline + "\t"
        }
        hline = hline + col
        first = false
      }
      lines.push(hline)
      let sample_names : Array[String] = Array(capacity=0)
      for j = 9; j < fh.columns.length(); j = j + 1 {
        sample_names.push(fh.columns[j])
      }
      for rec in self.records {
        lines.push(vcf_format_record(rec, sample_names))
      }
    }
  }
  let buf = StringBuilder()
  let mut first = true
  for l in lines {
    if !first {
      buf.write_string("\n")
    }
    buf.write_string(l)
    first = false
  }
  buf.to_string()
}

// ============================================================================
// P1-3: Variant 完整 API — VcfFile / VcfRecord / header 构建
// ============================================================================

///|
/// Allele index array `[0, 1, ... n_alt]`. Index 0 is REF, remaining indices
/// are the ALT alleles. Matches Pysam `VariantRecord.allele_indices`.
pub fn VcfRecord::allele_indices(self : VcfRecord) -> Array[Int32] {
  let n = self.alt.length() + 1
  let out = Array::make(n, 0)
  let mut i = 0
  while i < n {
    out[i] = i
    i = i + 1
  }
  out
}

///|
/// Length of the variant: `|ref_len - max(alt_len)|` for indels; otherwise
/// the ref length (SNPs). Matches Pysam `VariantRecord.rlen` /
/// `VariantRecord.rlen_sv` simplified semantics.
pub fn VcfRecord::rlen(self : VcfRecord) -> Int32 {
  let rl = self.ref_seq.length()
  if self.alt.length() == 0 {
    return rl
  }
  let mut max_alt = self.alt[0].length()
  let mut i = 1
  while i < self.alt.length() {
    let len = self.alt[i].length()
    if len > max_alt {
      max_alt = len
    }
    i = i + 1
  }
  if self.is_insertion() || self.is_deletion() {
    let diff = rl - max_alt
    if diff < 0 { -diff } else { diff }
  } else {
    rl
  }
}

///|
/// The variant's end position (1-based inclusive), computed as
/// `pos + rlen - 1`. Useful for range queries / BED-style overlap tests.
pub fn VcfRecord::pos_end(self : VcfRecord) -> Int32 {
  self.pos + self.rlen() - 1
}

///|
/// True when at least one filter is non-PASS (i.e. record is filtered).
pub fn VcfRecord::is_filtered(self : VcfRecord) -> Bool {
  for f in self.filter {
    if f != "PASS" && f != "." && f != "" {
      return true
    }
  }
  false
}

///|
/// Get an INFO field value by key, or None if absent.
pub fn VcfRecord::info_get(self : VcfRecord, key : String) -> String? {
  self.info.get(key)
}

///|
/// Set an INFO field value (flag fields use `""` as value). Returns a new
/// record (immutable convention, consistent with `BamRecord::set_tag`).
pub fn VcfRecord::info_set(self : VcfRecord, key : String, value : String) -> VcfRecord {
  let info = self.info.copy()
  info.set(key, value)
  let new_format : Array[String]? = match self.format {
    Some(f) => Some(f.copy())
    None => None
  }
  VcfRecord::{
    chrom: self.chrom,
    pos: self.pos,
    id: self.id.copy(),
    ref_seq: self.ref_seq,
    alt: self.alt.copy(),
    qual: self.qual,
    filter: self.filter.copy(),
    info,
    format: new_format,
    samples: self.samples.copy(),
  }
}

///|
/// Get a FORMAT column value for `sample_name`. Returns None when the
/// sample / format column is missing.
pub fn VcfRecord::sample_format(
  self : VcfRecord,
  sample_name : String,
  format_col : String,
) -> String? {
  let fmt = match self.format {
    Some(f) => f
    None => return None
  }
  let mut col_idx = -1
  let mut i = 0
  while i < fmt.length() {
    if fmt[i] == format_col {
      col_idx = i
      break
    }
    i = i + 1
  }
  if col_idx < 0 {
    return None
  }
  let data = match self.samples.get(sample_name) {
    Some(d) => d
    None => return None
  }
  if col_idx >= data.length() {
    None
  } else {
    Some(data[col_idx])
  }
}

///|
/// Contig names collected from `##contig=<...>` meta-headers. Empty when none
/// are declared. Matches Pysam `VariantHeader.contigs`.
pub fn VcfFile::get_contigs(self : VcfFile) -> Array[String] {
  let out : Array[String] = Array(capacity=0)
  for mh in self.meta_headers {
    if mh.key == "contig" {
      match mh.tags.get("ID") {
        Some(name) => out.push(name)
        None => ()
      }
    }
  }
  out
}

///|
/// Look up the declared length of `contig_name` from meta-headers; None if
/// absent or length tag missing. Matches Pysam
/// `VariantFile.header.contigs[].length`.
pub fn VcfFile::get_contig_length(
  self : VcfFile,
  contig_name : String,
) -> Int32? {
  let none : Int32? = None
  for mh in self.meta_headers {
    if mh.key == "contig" {
      match mh.tags.get("ID") {
        Some(id) =>
          if id == contig_name {
            match mh.tags.get("length") {
              Some(l) => return Some(parse_int(l))
              None => return none
            }
          }
        None => ()
      }
    }
  }
  none
}

///|
/// Contig ID → index (0-based) lookup, consistent with Pysam
/// `VariantFile.header.get_tid(name)`. Returns -1 when absent.
pub fn VcfFile::get_tid(self : VcfFile, contig_name : String) -> Int32 {
  let contigs = self.get_contigs()
  let mut i = 0
  while i < contigs.length() {
    if contigs[i] == contig_name {
      return i
    }
    i = i + 1
  }
  -1
}

///|
/// Fetch records overlapping `[start, stop)` on `contig`. Start/stop are
/// **1-based inclusive/exclusive** (VCF convention). Matches Pysam
/// `VariantFile.fetch(contig, start, stop)`.
pub fn VcfFile::fetch_region(
  self : VcfFile,
  contig : String,
  start : Int32,
  stop : Int32,
) -> Array[VcfRecord] {
  let out : Array[VcfRecord] = Array(capacity=0)
  for rec in self.records {
    if rec.chrom != contig {
      continue
    }
    let rec_end = rec.pos_end()
    if rec.pos < stop && rec_end >= start {
      out.push(rec)
    }
  }
  out
}

///|
/// Fetch records by samtools-style region string (`chr1:100-200` etc). Uses
/// the same parser as `parse_region` (BAM 模块).
pub fn VcfFile::fetch_region_str(self : VcfFile, region : String) -> Array[VcfRecord] {
  let (contig, start0, stop0) = parse_region(region)
  // VCF uses 1-based start, so +1
  let start = match start0 {
    Some(s) => s + 1
    None => 1
  }
  let stop = match stop0 {
    Some(e) => e
    None => 0 // no upper bound → use a very large value
  }
  if stop == 0 {
    let out : Array[VcfRecord] = Array(capacity=0)
    for rec in self.records {
      if rec.chrom == contig && rec.pos >= start {
        out.push(rec)
      }
    }
    out
  } else {
    self.fetch_region(contig, start, stop)
  }
}

///|
/// Build a minimal VcfFile from an array of contig names and lengths (for
/// writing). Matches Pysam `VariantFile.header =
/// VariantHeader().add_contig(...); VariantFile.new_record()`.
pub fn vcf_from_contigs(
  contig_names : Array[String],
  contig_lengths : Array[Int32],
  sample_names : Array[String]?,
) -> VcfFile {
  let meta = Array(capacity=0)
  meta.push(VcfMetaHeader::{
    key: "fileformat",
    value: "VCFv4.2",
    tags: Map([], capacity=0),
  })
  let n = contig_names.length()
  let mut i = 0
  while i < n {
    let tags = Map([], capacity=0)
    tags.set("ID", contig_names[i])
    tags.set("length", contig_lengths[i].to_string())
    meta.push(VcfMetaHeader::{
      key: "contig",
      value: "ID=\{contig_names[i]},length=\{contig_lengths[i]}",
      tags,
    })
    i = i + 1
  }
  let columns : Array[String] = Array(capacity=0)
  columns.push("#CHROM")
  columns.push("POS")
  columns.push("ID")
  columns.push("REF")
  columns.push("ALT")
  columns.push("QUAL")
  columns.push("FILTER")
  columns.push("INFO")
  match sample_names {
    None => ()
    Some(names) => {
      columns.push("FORMAT")
      for name in names {
        columns.push(name)
      }
    }
  }
  VcfFile::{
    meta_headers: meta,
    fixed_header: Some(VcfFixedHeader::{ columns, }),
    records: [],
  }
}

///|
/// Create a blank VcfRecord template, optionally pre-populating contig + pos.
/// The caller is expected to fill REF / ALT / FORMAT / samples.
pub fn vcf_new_record(
  vcf : VcfFile,
  chrom~ : String = "",
  pos~ : Int32 = 0,
) -> VcfRecord {
  // Build empty samples dict from fixed_header's sample names
  let samples = Map([], capacity=0)
  let fmt : Array[String]? = match vcf.fixed_header {
    None => None
    Some(fh) => {
      if fh.columns.length() > 8 && fh.columns[8] == "FORMAT" {
        let names = Array(capacity=0)
        let mut i = 9
        while i < fh.columns.length() {
          names.push(fh.columns[i])
          samples.set(fh.columns[i], [])
          i = i + 1
        }
        Some(names)
      } else {
        None
      }
    }
  }
  VcfRecord::{
    chrom,
    pos,
    id: [],
    ref_seq: "",
    alt: [],
    qual: None,
    filter: [],
    info: Map([], capacity=0),
    format: fmt,
    samples,
  }
}

///|
/// Append a contig declaration to `vcf.meta_headers`. Convenience wrapper for
/// the writer workflow.
pub fn VcfFile::add_contig(
  self : VcfFile,
  name : String,
  length : Int32,
) -> VcfFile {
  let meta = self.meta_headers.copy()
  let tags = Map([], capacity=0)
  tags.set("ID", name)
  tags.set("length", length.to_string())
  meta.push(VcfMetaHeader::{
    key: "contig",
    value: "ID=\{name},length=\{length}",
    tags,
  })
  VcfFile::{
    meta_headers: meta,
    fixed_header: self.fixed_header,
    records: self.records,
  }
}

///|
/// Append an arbitrary meta header (##KEY=VALUE) to the VCF.
pub fn VcfFile::add_meta(
  self : VcfFile,
  key : String,
  value : String,
) -> VcfFile {
  let meta = self.meta_headers.copy()
  meta.push(VcfMetaHeader::{
    key,
    value,
    tags: Map([], capacity=0),
  })
  VcfFile::{
    meta_headers: meta,
    fixed_header: self.fixed_header,
    records: self.records,
  }
}
