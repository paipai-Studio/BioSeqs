///|
/// Bio.SeqIO.FastaIO - two-line FASTA ("fasta-2line") support.
///
/// Mirrors Biopython's `Bio.SeqIO.FastaIO` fasta-2line variant: each record
/// is written as exactly two lines - the ">" header line and one single
/// line containing the whole sequence (never wrapped).
///
/// Parsing note (spec 5.2.1.1b): parsing reuses the standard FASTA
/// semantics, so a sequence that is unexpectedly split across multiple
/// lines is accepted and concatenated (i.e. fasta-2line input is parsed
/// exactly like multi-line FASTA). The strict two-line property is
/// enforced on writing only.

///|
/// Parse fasta-2line content. Reuses the standard FASTA parser: the
/// identifier is the first whitespace-delimited field of the header line
/// and the description is the remainder. A sequence unexpectedly split
/// over several lines is handled like standard multi-line FASTA (the
/// lines are concatenated).
pub fn parse_fasta_2line(content : String) -> Array[SeqRecord] {
  parse_fasta(content)
}

///|
/// Write records as two-line FASTA: one ">" header line plus one single
/// unwrapped sequence line per record.
pub fn write_fasta_2line(records : Array[SeqRecord]) -> String {
  let buf = StringBuilder()
  for rec in records {
    buf.write_char('>')
    buf.write_string(build_fasta_title(rec.id, rec.description))
    buf.write_char('\n')
    buf.write_string(rec.seq.to_string())
    buf.write_char('\n')
  }
  buf.to_string()
}
