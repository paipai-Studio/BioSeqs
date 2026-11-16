// ============================================================================
// FastaFile — in-memory FASTA with fast random access, aligned with Pysam's
// `pysam.FastaFile`. Supports optional `.fai` index metadata for contig
// ordering / total length (produced by samtools faidx).
//
// Inherits all parsing from `fasta_io::parse_fasta`; adds contig-lookup,
// region fetch, UCSC/samtools region string parsing.
// ============================================================================

///|
/// One contig record inside a FastaFile.
#derive(Eq)
pub struct FastaContig {
  /// Contig name (first whitespace-delimited token of the FASTA header).
  name : String
  /// Full header title (everything after '>' in the header line).
  description : String
  /// Sequence string (uppercase, whitespace removed).
  seq : String
  /// Length of the contig (cached for convenience).
  length : Int32
}

///|
pub extend FastaContig with Eq::{equal, not_equal}

///|
/// In-memory FASTA container with O(n) contig lookup (suitable for genomes
/// with a reasonable number of chromosomes). Matches Pysam's `FastaFile`.
#derive(Eq)
pub struct FastaFile {
  contigs : Array[FastaContig]
}

///|
pub extend FastaFile with Eq::{equal, not_equal}

///|
/// Build a FastaFile from FASTA text. Reuses the project's existing FASTA
/// parser to construct SeqRecords, then re-indexes by contig name.
pub fn FastaFile::from_text(text : String) -> FastaFile {
  let records = parse_fasta(text)
  let contigs : Array[FastaContig] = Array(capacity=records.length())
  for r in records {
    let seq_str = to_uppercase(r.seq.to_string())
    contigs.push(FastaContig::{
      name: r.id,
      description: r.description,
      seq: seq_str,
      length: seq_str.length(),
    })
  }
  FastaFile::{ contigs, }
}

///|
/// Look up a contig by name; returns None when absent. Linear scan (acceptable
/// for typical genome contig counts). Matches Pysam
/// `FastaFile.references`.
pub fn FastaFile::get_contig(self : FastaFile, name : String) -> FastaContig? {
  for c in self.contigs {
    if c.name == name {
      return Some(c)
    }
  }
  None
}

///|
/// Contig names in declaration order (matches Pysam
/// `FastaFile.references`).
pub fn FastaFile::references(self : FastaFile) -> Array[String] {
  let out : Array[String] = Array(capacity=self.contigs.length())
  for c in self.contigs {
    out.push(c.name)
  }
  out
}

///|
/// Number of contigs.
pub fn FastaFile::n_references(self : FastaFile) -> Int32 {
  self.contigs.length()
}

///|
/// Total genome size (sum of contig lengths).
pub fn FastaFile::total_length(self : FastaFile) -> Int64 {
  let mut total : Int64 = 0
  for c in self.contigs {
    total = total + c.length.to_int64()
  }
  total
}

///|
/// Length of `contig_name`; returns None when not present. Matches Pysam
/// `FastaFile.lengths[contig_name]`.
pub fn FastaFile::length_of(self : FastaFile, contig_name : String) -> Int32? {
  for c in self.contigs {
    if c.name == contig_name {
      return Some(c.length)
    }
  }
  None
}

///|
/// Fetch a sub-sequence from `contig_name` in `[start, stop)` (0-based
/// half-open). Returns empty string on out-of-bounds / missing contig.
/// Matches Pysam `FastaFile.fetch(contig, start, stop)`.
pub fn FastaFile::fetch(
  self : FastaFile,
  contig_name : String,
  start : Int32,
  stop : Int32,
) -> String {
  let contig = match self.get_contig(contig_name) {
    None => return ""
    Some(c) => c
  }
  // Clip to valid range
  let s = if start < 0 { 0 } else { start }
  let e = if stop > contig.length { contig.length } else { stop }
  if e <= s || s < 0 {
    return ""
  }
  contig.seq[s:e].to_owned()
}

///|
/// Fetch by samtools-style region string (`chr1:100-200`, `chr1:100+50`,
/// `chr1`, `chr1:100..200`). Region uses **1-based** start and **1-based
/// inclusive** end; we shift internally. Matches Pysam
/// `FastaFile.fetch("chr1:100-200")`.
pub fn FastaFile::fetch_region(self : FastaFile, region : String) -> String {
  let (contig, start0, stop0) = parse_region(region)
  // parse_region returns 0-based half-open; we need 0-based half-open for
  // FastaFile.fetch, so pass through directly.
  match (start0, stop0) {
    (Some(s), Some(e)) => self.fetch(contig, s, e)
    (Some(s), None) =>
      match self.length_of(contig) {
        Some(len) => self.fetch(contig, s, len)
        None => ""
      }
    (None, None) =>
      match self.get_contig(contig) {
        Some(c) => c.seq
        None => ""
      }
    (None, Some(e)) => self.fetch(contig, 0, e)
  }
}

///|
/// Content of a whole contig by name; returns empty string when absent.
pub fn FastaFile::get_seq(self : FastaFile, contig_name : String) -> String {
  match self.get_contig(contig_name) {
    Some(c) => c.seq
    None => ""
  }
}

// ============================================================================
// FASTA indexing helpers (samtools faidx `.fai` format)
// ============================================================================

///|
/// Single row of a `.fai` index file.
pub struct FaiEntry {
  name : String
  length : Int32
  offset : Int64
  bases_per_line : Int32
  bytes_per_line : Int32
}

///|
/// Parse a `.fai` index file string (one entry per line, tab-separated).
/// Returns an empty array when input is empty or malformed. Named
/// `parse_fai_index` to avoid collision with the `faidx.mbt` parser which
/// raises a typed error on malformed input.
pub fn parse_fai_index(content : String) -> Array[FaiEntry] {
  let out : Array[FaiEntry] = Array(capacity=0)
  if content.length() == 0 {
    return out
  }
  let lines = content.split("\n").to_array()
  for line in lines {
    let tokens = line.split("\t").to_array()
    if tokens.length() < 5 {
      continue
    }
    let name = tokens[0].to_owned()
    let length = parse_int(tokens[1].to_owned())
    // offset may not be present in synthetic fasta; tolerate default 0
    let offset = if tokens[2].length() > 0 {
      parse_int(tokens[2].to_owned()).to_int64()
    } else {
      0L
    }
    let bpl = parse_int(tokens[3].to_owned())
    let btl = parse_int(tokens[4].to_owned())
    out.push(FaiEntry::{
      name,
      length,
      offset,
      bases_per_line: bpl,
      bytes_per_line: btl,
    })
  }
  out
}

///|
/// Build a `.fai`-style index string from a FastaFile. Useful for samtools
/// faidx-like workflows that need the canonical index representation.
pub fn FastaFile::build_fai(self : FastaFile) -> String {
  let sb = StringBuilder()
  for c in self.contigs {
    sb.write_string(c.name)
    sb.write_char('\t')
    sb.write_string(c.length.to_string())
    sb.write_char('\t')
    // synthetic offset placeholder (0s); real offsets depend on file layout
    sb.write_string("0")
    sb.write_char('\t')
    sb.write_string(c.length.to_string()) // bases_per_line = full contig
    sb.write_char('\t')
    sb.write_string((c.length + 1).to_string()) // bytes_per_line incl. newline
    sb.write_char('\n')
  }
  sb.to_string()
}
