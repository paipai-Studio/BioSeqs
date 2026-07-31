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
//   "fasta", "fastq"/"fastq-sanger", "genbank"/"gb",
//   "embl", "pir"/"nbrf", "tab"
// Supported formats (write):
//   "fasta", "fastq"/"fastq-sanger", "genbank"/"gb", "tab"

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
    "genbank" | "gb" => parse_genbank(content)
    "embl" => parse_embl(content)
    "pir" | "nbrf" => parse_pir(content)
    "tab" => parse_tab(content)
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
    "genbank" | "gb" => write_genbank(records)
    "tab" => write_tab(records)
    _ => raise SeqIOError("Unknown format: " + format)
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
