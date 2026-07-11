// SeqIO: unified sequence file I/O, mirroring BioPython's Bio.SeqIO.
//
// Provides:
//   - seqio_parse(content, format)   -> Array[SeqRecord]
//   - seqio_read(content, format)   -> SeqRecord  (exactly one record)
//   - seqio_write(records, format)   -> String
//   - seqio_to_dict(records)          -> Map[String, SeqRecord]
//
// Supported formats: "fasta", "fastq" (Sanger/Phred+33).

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
    _ => raise SeqIOError("Unknown format: " + format)
  }
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
