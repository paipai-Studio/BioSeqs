// FormatRegistry: central registry of sequence/alignment format names,
// aliases, and capability flags for the unified SeqIO entry points.
//
// This is the single source of truth required by spec 4.4.3 ("the format
// alias table must be defined in one place"): every entry point normalizes
// its format argument through the alias tables below before dispatching.
//
// Normalization is a two-level constant-time match (spec 4.1.1):
//   1. the caller's string (lowercased) is matched to the canonical name;
//   2. the entry point matches the canonical name to an implementation.

///|
/// Normalize a sequence format identifier to its canonical name.
///
/// The lookup is case-insensitive (the input is lowercased first) and
/// covers every alias registered for the unified SeqIO entry points
/// (seqio_parse / seqio_write / seqio_index). Unregistered identifiers
/// raise `Unknown format: <original>` with the original input echoed back,
/// without scanning the content.
pub fn seqio_normalize_format(format : String) -> String raise SeqIOError {
  match format.to_lower() {
    "fasta" => "fasta"
    "fastq" | "fastq-sanger" => "fastq-sanger"
    "fastq-illumina" => "fastq-illumina"
    "fastq-solexa" => "fastq-solexa"
    "qual" => "qual"
    "genbank" | "gb" => "genbank"
    "embl" => "embl"
    "pir" | "nbrf" => "pir"
    "tab" => "tab"
    "swiss" => "swiss"
    "uniprot-xml" => "uniprot-xml"
    "seqxml" => "seqxml"
    "imgt" => "imgt"
    "ig" => "ig"
    "insdc" => "insdc"
    "phd" => "phd"
    "ace" => "ace"
    "gfa" => "gfa"
    "pdb-seqres" => "pdb-seqres"
    "pdb-atom" => "pdb-atom"
    "2bit" => "2bit"
    "sff" => "sff"
    "sff-trim" => "sff-trim"
    "snapgene" => "snapgene"
    "gck" => "gck"
    "xdna" => "xdna"
    "nib" => "nib"
    "fasta-2line" => "fasta-2line"
    "cif-seqres" => "cif-seqres"
    "cif-atom" => "cif-atom"
    _ => raise SeqIOError("Unknown format: " + format)
  }
}

///|
/// Normalize an alignment format identifier to its canonical name.
///
/// Mirrors seqio_normalize_format for the AlignIO entry points
/// (alignio_parse / alignio_write). Unregistered identifiers raise
/// `Unknown format: <original>` with the original input echoed back.
pub fn alignio_normalize_format(format : String) -> String raise AlignIOError {
  match format.to_lower() {
    "clustal" => "clustal"
    "phylip" => "phylip"
    "phylip-sequential" => "phylip-sequential"
    "fasta" => "fasta"
    "stockholm" | "sto" => "stockholm"
    "nexus" => "nexus"
    "maf" => "maf"
    "msf" => "msf"
    "emboss" => "emboss"
    "a2m" => "a2m"
    "sam" => "sam"
    "psl" => "psl"
    "mauve" => "mauve"
    _ => raise AlignIOError("Unknown format: " + format)
  }
}

///|
/// All sequence formats the unified parse entry point can handle
/// (canonical names, aliases resolved).
pub fn seqio_supported_formats() -> Array[String] {
  [
    "fasta", "fastq-sanger", "fastq-illumina", "fastq-solexa", "qual", "genbank",
    "embl", "pir", "tab", "swiss", "uniprot-xml", "seqxml", "imgt", "ig", "insdc",
    "phd", "ace", "gfa", "pdb-seqres", "pdb-atom", "2bit", "sff", "sff-trim", "snapgene",
    "gck", "xdna", "nib", "fasta-2line", "cif-seqres", "cif-atom",
  ]
}

///|
/// All sequence formats the unified write entry point can produce.
pub fn seqio_writable_formats() -> Array[String] {
  [
    "fasta", "fastq-sanger", "fastq-illumina", "fastq-solexa", "qual", "genbank",
    "tab", "embl", "pir", "fasta-2line", "seqxml", "swiss", "ig", "imgt", "phd",
  ]
}

///|
/// All sequence formats the unified index entry point supports.
pub fn seqio_indexable_formats() -> Array[String] {
  [
    "fasta", "fastq-sanger", "fastq-illumina", "fastq-solexa", "genbank", "embl",
    "tab", "qual", "swiss", "seqxml",
  ]
}

///|
/// Formats whose payloads are binary and must be carried as base64 text
/// by the in-memory unified entry point (spec 4.3.1).
pub fn seqio_binary_formats() -> Array[String] {
  ["2bit", "sff", "sff-trim", "snapgene", "gck", "xdna", "nib"]
}
