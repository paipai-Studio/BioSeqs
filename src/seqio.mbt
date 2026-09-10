// SeqIO: unified sequence file I/O, mirroring BioPython's Bio.SeqIO.
//
// Provides:
//   - seqio_parse(content, format)     -> Array[SeqRecord]
//   - seqio_read(content, format)     -> SeqRecord  (exactly one record)
//   - seqio_write(records, format)     -> String
//   - seqio_convert(content, in_fmt, out_fmt) -> String
//   - seqio_to_dict(records)            -> Map[String, SeqRecord]
//
// Supported formats (parse):
//   "fasta", "fastq"/"fastq-sanger", "fastq-illumina", "fastq-solexa",
//   "qual", "genbank"/"gb", "embl", "pir"/"nbrf", "tab",
//   "swiss", "uniprot-xml", "seqxml", "imgt", "ig", "insdc", "phd", "ace",
//   "gfa", "pdb-seqres", "pdb-atom", "fasta-2line", "cif-seqres",
//   "cif-atom", and the base64-carried binary formats "2bit", "sff",
//   "sff-trim", "snapgene", "gck", "xdna", "nib"
// Supported formats (write):
//   "fasta", "fastq"/"fastq-sanger", "fastq-illumina", "fastq-solexa",
//   "qual", "genbank"/"gb", "tab", "embl", "pir", "fasta-2line",
//   "seqxml", "swiss", "ig", "imgt", "phd"

///|
/// Parse sequence file content in the given format and return all records.
///
/// Mirrors BioPython's `Bio.SeqIO.parse(handle, format)`.
pub fn seqio_parse(
  content : String,
  format : String,
) -> Array[SeqRecord] raise SeqIOError {
  match format {
    "fasta" => parse_fasta(content)
    "fastq" | "fastq-sanger" => parse_fastq(content)
    "fastq-illumina" => parse_fastq_illumina(content)
    "fastq-solexa" => parse_fastq_solexa(content)
    "qual" => parse_qual(content)
    "genbank" | "gb" => parse_genbank(content)
    "embl" => parse_embl(content)
    "pir" | "nbrf" => parse_pir(content)
    "tab" => parse_tab(content)
    // Text/container formats unified through the seqio_unify adapters.
    "swiss" => parse_swiss_records(content)
    "uniprot-xml" => parse_uniprot_xml_records(content)
    "seqxml" => parse_seqxml_records(content)
    "imgt" => parse_imgt_records(content)
    "ig" => parse_ig_records(content)
    "insdc" => parse_insdc_records(content)
    "phd" => parse_phd_records(content)
    "ace" => parse_ace_records(content)
    "gfa" => parse_gfa_records(content)
    "pdb-seqres" => parse_pdb_seqres_records(content)
    "pdb-atom" => parse_pdb_atom_records(content)
    // New formats producing SeqRecord directly.
    "fasta-2line" => parse_fasta_2line(content)
    "cif-seqres" => parse_cif_seqres(content)
    "cif-atom" => parse_cif_atom(content)
    // Binary formats carried as base64 through the seqio_binary_unify
    // adapters (strict base64 validation, spec 4.3.1).
    "2bit" => parse_twobit_records(content)
    "sff" => parse_sff_records(content)
    "sff-trim" => parse_sff_trim_records(content)
    "snapgene" => parse_snapgene_records(content)
    "gck" => parse_gck_records(content)
    "xdna" => parse_xdna_records(content)
    "nib" => parse_nib_records(content)
    _ => raise SeqIOError("Unknown format: " + format)
  }
}

///|
/// Read exactly one record from the content. Raises SeqIOError if the
/// content contains zero or more than one record.
///
/// Mirrors BioPython's `Bio.SeqIO.read(handle, format)`.
pub fn seqio_read(
  content : String,
  format : String,
) -> SeqRecord raise SeqIOError {
  let records = seqio_parse(content, format)
  if records.length() == 0 {
    raise SeqIOError("No records found in handle")
  }
  if records.length() > 1 {
    raise SeqIOError("More than one record found in handle")
  }
  records[0]
}

///|
/// Write records to a string in the given format.
///
/// Mirrors BioPython's `Bio.SeqIO.write(sequences, handle, format)`.
/// Returns the formatted text.
pub fn seqio_write(
  records : Array[SeqRecord],
  format : String,
) -> String raise SeqIOError {
  match format {
    "fasta" => write_fasta(records)
    "fastq" | "fastq-sanger" => write_fastq(records)
    "fastq-illumina" => write_fastq_illumina(records)
    "fastq-solexa" => write_fastq_solexa(records)
    "qual" => write_qual(records)
    "genbank" | "gb" => write_genbank(records)
    "tab" => write_tab(records)
    // Formats unified through the seqio_write_ext adapters (spec 5.3).
    "embl" => write_embl_records(records)
    "pir" => write_pir_records(records)
    "fasta-2line" => write_fasta_2line(records)
    "seqxml" => write_seqxml_records(records)
    "swiss" => write_swiss_records(records)
    "ig" => write_ig_records(records)
    "imgt" => write_imgt_records(records)
    "phd" => write_phd_records(records)
    _ => {
      // Known but non-writable formats (e.g. 2bit, uniprot-xml) list the
      // writable set; truly unknown formats keep the legacy message.
      match (try? seqio_normalize_format(format)) {
        Ok(_) => {
          let supported = seqio_writable_formats().join(", ")
          raise SeqIOError(
            "format does not support writing:\n" +
            format +
            "。Supported: " +
            supported,
          )
        }
        Err(SeqIOError(msg)) => raise SeqIOError(msg)
      }
    }
  }
}

///|
/// Convert sequence file content from one format to another.
///
/// Mirrors BioPython's `Bio.SeqIO.convert(in_handle, in_format, out_handle, out_format)`.
/// Returns the converted text.
/// Raises SeqIOError if either format is unsupported or parsing fails.
pub fn seqio_convert(
  content : String,
  in_format : String,
  out_format : String,
) -> String raise SeqIOError {
  let records = seqio_parse(content, in_format)
  seqio_write(records, out_format)
}

///|
/// Convert an array of records to a dictionary keyed by record id.
///
/// Mirrors BioPython's `Bio.SeqIO.to_dict(sequences)`.
/// If duplicate ids exist, later records overwrite earlier ones.
pub fn seqio_to_dict(records : Array[SeqRecord]) -> Map[String, SeqRecord] {
  let d : Map[String, SeqRecord] = Map([], capacity=16)
  for rec in records {
    d[rec.id] = rec
  }
  d
}

// ============================================================================
// Error handling helpers with line number support
// ============================================================================

///|
/// Create a SeqIOError with line number information.
/// 
/// Mirrors Biopython's approach of including line numbers in error messages.
/// This makes debugging much easier for large files.
pub fn seqio_error_with_line(
  message : String,
  line_number : Int,
) -> SeqIOError {
  SeqIOError(message + " (at line " + line_number.to_string() + ")")
}
