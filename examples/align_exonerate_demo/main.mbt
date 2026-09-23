///|
fn coordinate_text(values : Array[Int]) -> String {
  let output = StringBuilder()
  output.write_char('[')
  for index in 0..<values.length() {
    if index > 0 {
      output.write_string(", ")
    }
    output.write_string(values[index].to_string())
  }
  output.write_char(']')
  output.to_string()
}

///|
fn main {
  println("=== Biopython Bio.Align.exonerate Demo ===")

  let document = @src.align_exonerate_parse(@src.align_exonerate_example_text()) catch {
    AlignExonerateError(message) =>
      abort("Exonerate parsing failed: " + message)
  }
  let spliced = document.alignments[0]
  let translated = document.alignments[1]

  println("\n1. Parse report metadata and alignments")
  println("  " + document.summary())
  println("  command: " + document.metadata.command_line)
  for alignment in document.alignments {
    println("  " + alignment.summary())
  }

  println("\n2. Preserve spliced vulgar operations")
  let spliced_path = spliced.coordinate_path()
  println("  operations:  " + spliced.operation_codes())
  println("  target path: " + coordinate_text(spliced_path[0]))
  println("  query path:  " + coordinate_text(spliced_path[1]))
  let spliced_counts = spliced.counts()
  println(
    "  introns=" +
    spliced_counts.introns.to_string() +
    ", intron units=" +
    spliced_counts.intron_units.to_string() +
    ", matching operations=" +
    spliced_counts.matching_operations.to_string(),
  )

  println("\n3. Map residues across an intron")
  println(
    "  query residue 6 -> target residue " +
    spliced.query_position_to_target(6).unwrap_or(-1).to_string(),
  )
  println(
    "  target residue 108 -> query residue " +
    spliced.target_position_to_query(108).unwrap_or(-1).to_string() +
    " (-1 denotes an intron)",
  )

  println("\n4. Handle protein-to-DNA reverse translation")
  let translated_path = translated.coordinate_path()
  println("  target path: " + coordinate_text(translated_path[0]))
  println("  query path:  " + coordinate_text(translated_path[1]))
  for residue in 0..<translated.query_end {
    println(
      "  protein " +
      residue.to_string() +
      " -> reverse-strand nucleotide " +
      translated.query_position_to_target(residue).unwrap_or(-1).to_string(),
    )
  }

  println("\n5. Write canonical vulgar and cigar reports")
  let vulgar = @src.align_exonerate_write(document, format="vulgar") catch {
    AlignExonerateError(message) =>
      abort("vulgar serialization failed: " + message)
  }
  let cigar = @src.align_exonerate_write(document, format="cigar") catch {
    AlignExonerateError(message) =>
      abort("cigar serialization failed: " + message)
  }
  let vulgar_round_trip = @src.align_exonerate_parse(vulgar) catch {
    AlignExonerateError(message) =>
      abort("vulgar round trip failed: " + message)
  }
  let cigar_round_trip = @src.align_exonerate_parse(cigar) catch {
    AlignExonerateError(message) => abort("cigar round trip failed: " + message)
  }
  println(
    "  vulgar preserves document: " +
    (vulgar_round_trip == document).to_string(),
  )
  println(
    "  cigar preserves first coordinate path: " +
    (cigar_round_trip.alignments[0].coordinate_path() ==
    spliced.coordinate_path()).to_string(),
  )
  println("  vulgar bytes: " + vulgar.length().to_string())
  println("  cigar bytes:  " + cigar.length().to_string())

  println("\n=== Demo Complete ===")
}
