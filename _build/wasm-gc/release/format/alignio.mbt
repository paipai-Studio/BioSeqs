// AlignIO module, mirroring BioPython's Bio.AlignIO.
//
// Provides functions for reading and writing multiple sequence alignments
// in various formats: clustal, phylip, phylip-sequential, fasta, etc.

///|
/// Parse alignment format text and return an array of MultipleSeqAlignment objects.
///
/// Mirrors BioPython's AlignIO.parse(handle, format).
pub fn alignio_parse(
  content : String,
  format : String,
) -> Array[MultipleSeqAlignment] raise AlignIOError {
  match format {
    "clustal" => [parse_clustal(content)]
    "phylip" => [parse_phylip(content)]
    "phylip-sequential" => [parse_phylip_sequential(content)]
    "fasta" => [parse_fasta_as_alignment(content)]
    _ => raise AlignIOError("Unknown format: " + format)
  }
}

///|
/// Read a single alignment from format text.
///
/// Mirrors BioPython's AlignIO.read(handle, format).
/// Raises AlignIOError if there is not exactly one alignment.
pub fn alignio_read(
  content : String,
  format : String,
) -> MultipleSeqAlignment raise AlignIOError {
  let alignments = alignio_parse(content, format)
  if alignments.length() == 0 {
    raise AlignIOError("No alignment found in handle")
  }
  if alignments.length() > 1 {
    raise AlignIOError("More than one alignment found in handle")
  }
  alignments[0]
}

///|
/// Write alignments to format text.
///
/// Mirrors BioPython's AlignIO.write(alignments, handle, format).
pub fn alignio_write(
  alignments : Array[MultipleSeqAlignment],
  format : String,
) -> String raise AlignIOError {
  match format {
    "clustal" => write_clustal_multiple(alignments)
    "phylip" => write_phylip_multiple(alignments)
    "phylip-sequential" => write_phylip_sequential_multiple(alignments)
    "fasta" => write_fasta_as_alignment(alignments)
    _ => raise AlignIOError("Unknown format: " + format)
  }
}

///|
/// Parse a FASTA alignment (same as parse_fasta but returns MultipleSeqAlignment).
fn parse_fasta_as_alignment(
  content : String,
) -> MultipleSeqAlignment raise AlignIOError {
  let records = parse_fasta(content)
  MultipleSeqAlignment::new(records)
}

///|
/// Write alignments as FASTA format.
fn write_fasta_as_alignment(alignments : Array[MultipleSeqAlignment]) -> String {
  let all_records = Array::new()
  for align in alignments {
    for rec in align.iter() {
      all_records.push(rec)
    }
  }
  write_fasta(all_records)
}

///|
/// Write multiple alignments as Clustal format (concatenated).
fn write_clustal_multiple(alignments : Array[MultipleSeqAlignment]) -> String {
  let buf = StringBuilder::new()
  for i, align in alignments {
    if i > 0 {
      buf.write_char('\n')
    }
    buf.write_string(write_clustal(align))
  }
  buf.to_string()
}

///|
/// Write multiple alignments as PHYLIP interlaced format (concatenated).
fn write_phylip_multiple(alignments : Array[MultipleSeqAlignment]) -> String {
  let buf = StringBuilder::new()
  for i, align in alignments {
    if i > 0 {
      buf.write_char('\n')
    }
    buf.write_string(write_phylip(align))
  }
  buf.to_string()
}

///|
/// Write multiple alignments as PHYLIP sequential format (concatenated).
fn write_phylip_sequential_multiple(
  alignments : Array[MultipleSeqAlignment],
) -> String {
  let buf = StringBuilder::new()
  for i, align in alignments {
    if i > 0 {
      buf.write_char('\n')
    }
    buf.write_string(write_phylip_sequential(align))
  }
  buf.to_string()
}
