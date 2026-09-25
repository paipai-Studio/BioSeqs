///|
/// bcftools-style variant operations: filter, merge, isec, normalize, concat,
/// view (query/select). These operate on `VcfFile` from `vcf.mbt`.

// ============================================================
// bcftools filter — select/exclude records by expression
// ============================================================

///|
/// Filter type: PASS or a named filter.
pub fn VcfRecord::has_filter(self : VcfRecord, filter_name : String) -> Bool {
  for f in self.filter {
    if f == filter_name {
      return true
    }
  }
  false
}

///|
/// Check if a record passed all filters (FILTER is "." or "PASS").
pub fn VcfRecord::is_pass(self : VcfRecord) -> Bool {
  if self.filter.length() == 0 {
    return true
  }
  for f in self.filter {
    if f != "." && f != "PASS" {
      return false
    }
  }
  true
}

///|
/// bcftools filter: keep only records where FILTER matches (PASS or specific).
/// If `filter_name` is "PASS", keeps records with no filtering failure.
pub fn vcf_filter_pass(file : VcfFile) -> VcfFile {
  let result : Array[VcfRecord] = []
  for rec in file.records {
    if rec.is_pass() {
      result.push(rec)
    }
  }
  VcfFile::{
    meta_headers: file.meta_headers,
    fixed_header: file.fixed_header,
    records: result,
  }
}

///|
/// bcftools filter -i: keep records where QUAL >= threshold.
pub fn vcf_filter_min_qual(file : VcfFile, min_qual : Double) -> VcfFile {
  let result : Array[VcfRecord] = []
  for rec in file.records {
    match rec.qual {
      Some(q) => {
        if q >= min_qual {
          result.push(rec)
        }
      }
      None => ()
    }
  }
  VcfFile::{
    meta_headers: file.meta_headers,
    fixed_header: file.fixed_header,
    records: result,
  }
}

///|
/// bcftools filter -i: keep records where a specific FILTER tag is present.
pub fn vcf_filter_tag(file : VcfFile, tag : String) -> VcfFile {
  let result : Array[VcfRecord] = []
  for rec in file.records {
    if rec.has_filter(tag) {
      result.push(rec)
    }
  }
  VcfFile::{
    meta_headers: file.meta_headers,
    fixed_header: file.fixed_header,
    records: result,
  }
}

///|
/// bcftools filter -e: exclude records with a specific FILTER tag.
pub fn vcf_exclude_tag(file : VcfFile, tag : String) -> VcfFile {
  let result : Array[VcfRecord] = []
  for rec in file.records {
    if !(rec.has_filter(tag)) {
      result.push(rec)
    }
  }
  VcfFile::{
    meta_headers: file.meta_headers,
    fixed_header: file.fixed_header,
    records: result,
  }
}

// ============================================================
// bcftools view — select by variant type
// ============================================================

///|
/// bcftools view -v snps: keep only SNP records.
pub fn vcf_view_snps(file : VcfFile) -> VcfFile {
  let result : Array[VcfRecord] = []
  for rec in file.records {
    if rec.is_snp() {
      result.push(rec)
    }
  }
  VcfFile::{
    meta_headers: file.meta_headers,
    fixed_header: file.fixed_header,
    records: result,
  }
}

///|
/// bcftools view -v indels: keep only insertion/deletion records.
pub fn vcf_view_indels(file : VcfFile) -> VcfFile {
  let result : Array[VcfRecord] = []
  for rec in file.records {
    if rec.is_insertion() || rec.is_deletion() {
      result.push(rec)
    }
  }
  VcfFile::{
    meta_headers: file.meta_headers,
    fixed_header: file.fixed_header,
    records: result,
  }
}

// ============================================================
// bcftools merge — combine records from multiple VCF files
// ============================================================

///|
/// bcftools merge: combine multiple VCF files into one.
/// Records are concatenated in order. Headers from the first file are kept;
/// sample columns from subsequent files are appended.
pub fn vcf_merge(files : Array[VcfFile]) -> VcfFile {
  if files.length() == 0 {
    return VcfFile::{
      meta_headers: [],
      fixed_header: None,
      records: [],
    }
  }
  let first = files[0]
  let all_records : Array[VcfRecord] = []
  for rec in first.records {
    all_records.push(rec)
  }
  for i = 1; i < files.length(); i = i + 1 {
    for rec in files[i].records {
      all_records.push(rec)
    }
  }
  VcfFile::{
    meta_headers: first.meta_headers,
    fixed_header: first.fixed_header,
    records: all_records,
  }
}

// ============================================================
// bcftools concat — concatenate VCF files (same samples)
// ============================================================

///|
/// bcftools concat: concatenate VCF files with the same sample set.
/// Records are sorted by (chrom, pos) after merging.
pub fn vcf_concat(files : Array[VcfFile]) -> VcfFile {
  let merged = vcf_merge(files)
  vcf_sort(merged)
}

// ============================================================
// bcftools isec — intersection/complement of VCF files
// ============================================================

///|
/// bcftools isec -n=2: intersection of two VCF files (records present in both,
/// matched by chrom:pos:ref:alt).
pub fn vcf_isec(file_a : VcfFile, file_b : VcfFile) -> VcfFile {
  let set_b = vcf_record_keyset(file_b)
  let result : Array[VcfRecord] = []
  for rec in file_a.records {
    let key = vcf_record_key(rec)
    if set_b.contains(key) {
      result.push(rec)
    }
  }
  VcfFile::{
    meta_headers: file_a.meta_headers,
    fixed_header: file_a.fixed_header,
    records: result,
  }
}

///|
/// bcftools isec -n=2 -C: complement (records in A not in B).
pub fn vcf_isec_complement(file_a : VcfFile, file_b : VcfFile) -> VcfFile {
  let set_b = vcf_record_keyset(file_b)
  let result : Array[VcfRecord] = []
  for rec in file_a.records {
    let key = vcf_record_key(rec)
    if !(set_b.contains(key)) {
      result.push(rec)
    }
  }
  VcfFile::{
    meta_headers: file_a.meta_headers,
    fixed_header: file_a.fixed_header,
    records: result,
  }
}

///|
/// bcftools isec -n=2: union (records in either A or B, deduplicated by key).
pub fn vcf_union(file_a : VcfFile, file_b : VcfFile) -> VcfFile {
  let seen : Map[String, Bool] = Map([], capacity=0)
  let result : Array[VcfRecord] = []
  for rec in file_a.records {
    let key = vcf_record_key(rec)
    if !(seen.contains(key)) {
      seen.set(key, true)
      result.push(rec)
    }
  }
  for rec in file_b.records {
    let key = vcf_record_key(rec)
    if !(seen.contains(key)) {
      seen.set(key, true)
      result.push(rec)
    }
  }
  VcfFile::{
    meta_headers: file_a.meta_headers,
    fixed_header: file_a.fixed_header,
    records: result,
  }
}

// ============================================================
// bcftools norm — left-align and normalize variants
// ============================================================

///|
/// bcftools norm -m -: split multi-allelic records into separate records.
/// Each ALT allele becomes its own record (same chrom:pos:ref).
pub fn vcf_norm_split(file : VcfFile) -> VcfFile {
  let result : Array[VcfRecord] = []
  for rec in file.records {
    if rec.alt.length() <= 1 {
      result.push(rec)
    } else {
      for i = 0; i < rec.alt.length(); i = i + 1 {
        result.push(
          VcfRecord::{
            chrom: rec.chrom,
            pos: rec.pos,
            id: rec.id,
            ref_seq: rec.ref_seq,
            alt: [rec.alt[i]],
            qual: rec.qual,
            filter: rec.filter,
            info: rec.info,
            format: rec.format,
            samples: rec.samples,
          },
        )
      }
    }
  }
  VcfFile::{
    meta_headers: file.meta_headers,
    fixed_header: file.fixed_header,
    records: result,
  }
}

///|
/// bcftools norm -m +: merge records with the same chrom:pos:ref into one
/// multi-allelic record (ALT alleles concatenated).
pub fn vcf_norm_merge(file : VcfFile) -> VcfFile {
  let groups : Map[String, Array[VcfRecord]] = Map([], capacity=0)
  let order : Array[String] = []
  for rec in file.records {
    let key = rec.chrom + ":" + rec.pos.to_string() + ":" + rec.ref_seq
    match groups.get(key) {
      None => {
        groups.set(key, [rec])
        order.push(key)
      }
      Some(arr) => arr.push(rec)
    }
  }
  let result : Array[VcfRecord] = []
  for key in order {
    let group = groups[key]
    let first = group[0]
    let alt : Array[String] = []
    for rec in group {
      for a in rec.alt {
        // dedup
        let mut found = false
        for existing in alt {
          if existing == a {
            found = true
            break
          }
        }
        if !(found) {
          alt.push(a)
        }
      }
    }
    result.push(
      VcfRecord::{
        chrom: first.chrom,
        pos: first.pos,
        id: first.id,
        ref_seq: first.ref_seq,
        alt,
        qual: first.qual,
        filter: first.filter,
        info: first.info,
        format: first.format,
        samples: first.samples,
      },
    )
  }
  VcfFile::{
    meta_headers: file.meta_headers,
    fixed_header: file.fixed_header,
    records: result,
  }
}

// ============================================================
// bcftools sort — sort records by chrom then position
// ============================================================

///|
/// Sort VCF records by (chrom, pos) in ascending order.
pub fn vcf_sort(file : VcfFile) -> VcfFile {
  let records = file.records.copy()
  records.sort_by(fn(a, b) {
    if a.chrom != b.chrom {
      if a.chrom > b.chrom {
        return 1
      }
      return -1
    }
    a.pos - b.pos
  })
  VcfFile::{
    meta_headers: file.meta_headers,
    fixed_header: file.fixed_header,
    records,
  }
}

// ============================================================
// Helper: record key for isec/dedup
// ============================================================

///|
/// Generate a unique key for a VCF record: chrom:pos:ref:alt1,alt2,...
fn vcf_record_key(rec : VcfRecord) -> String {
  let alt_str = join_str_array(rec.alt, ",")
  rec.chrom + ":" + rec.pos.to_string() + ":" + rec.ref_seq + ":" + alt_str
}

///|
/// Build a set of record keys for fast lookup.
fn vcf_record_keyset(file : VcfFile) -> Map[String, Bool] {
  let m : Map[String, Bool] = Map([], capacity=0)
  for rec in file.records {
    m.set(vcf_record_key(rec), true)
  }
  m
}

// ============================================================
// bcftools query — extract fields
// ============================================================

///|
/// bcftools query -f: extract CHROM,POS,REF,ALT from records.
/// Returns one line per record, tab-separated.
pub fn vcf_query_basic(file : VcfFile) -> String {
  let sb = StringBuilder()
  for rec in file.records {
    sb.write_string(rec.chrom)
    sb.write_char('\t')
    sb.write_string(rec.pos.to_string())
    sb.write_char('\t')
    sb.write_string(rec.ref_seq)
    sb.write_char('\t')
    sb.write_string(join_str_array(rec.alt, ","))
    sb.write_char('\n')
  }
  sb.to_string()
}

///|
/// bcftools query -f: extract INFO field value for a specific key.
/// Returns one line per record that has the key (records missing it are skipped).
pub fn vcf_query_info(file : VcfFile, key : String) -> String {
  let sb = StringBuilder()
  for rec in file.records {
    match rec.info.get(key) {
      Some(v) => {
        sb.write_string(v)
        sb.write_char('\n')
      }
      None => ()
    }
  }
  sb.to_string()
}

// ============================================================
// bcftools stats — basic statistics
// ============================================================

///|
/// Basic VCF statistics (similar to `bcftools stats`).
pub struct VcfStats {
  total_records : Int
  snps : Int
  insertions : Int
  deletions : Int
  mnps : Int
  multi_allelic : Int
  pass_records : Int
  failed_records : Int
  qual_mean : Double
}

///|
/// Compute basic statistics for a VCF file.
pub fn vcf_stats(file : VcfFile) -> VcfStats {
  let n = file.records.length()
  let mut snps = 0
  let mut ins = 0
  let mut dels = 0
  let mut mnps = 0
  let mut multi = 0
  let mut pass_n = 0
  let mut fail_n = 0
  let mut qual_sum = 0.0
  let mut qual_count = 0
  for rec in file.records {
    if rec.is_snp() {
      snps = snps + 1
    } else if rec.is_insertion() {
      ins = ins + 1
    } else if rec.is_deletion() {
      dels = dels + 1
    } else {
      // Multi-nucleotide variant (MNP) or complex
      if rec.ref_seq.length() > 1 || (rec.alt.length() > 0 && rec.alt[0].length() > 1) {
        mnps = mnps + 1
      }
    }
    if rec.alt.length() > 1 {
      multi = multi + 1
    }
    if rec.is_pass() {
      pass_n = pass_n + 1
    } else {
      fail_n = fail_n + 1
    }
    match rec.qual {
      Some(q) => {
        qual_sum = qual_sum + q
        qual_count = qual_count + 1
      }
      None => ()
    }
  }
  VcfStats::{
    total_records: n,
    snps,
    insertions: ins,
    deletions: dels,
    mnps,
    multi_allelic: multi,
    pass_records: pass_n,
    failed_records: fail_n,
    qual_mean: if qual_count > 0 { qual_sum / qual_count.to_double() } else { 0.0 },
  }
}