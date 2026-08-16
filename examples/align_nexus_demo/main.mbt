///|
fn nexus_integer_array_text(values : Array[Int]) -> String {
  let output = StringBuilder::new()
  output.write_char('[')
  for index = 0; index < values.length(); index = index + 1 {
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
  println("=== Biopython Bio.Align.nexus Demo ===")

  let alignment = @src.align_nexus_parse(@src.align_nexus_example_text()) catch {
    AlignNexusError(message) => abort("NEXUS parsing failed: " + message)
  }

  println("\n1. Parse an interleaved DATA matrix")
  println("  " + alignment.summary())
  println(
    "  datatype: " +
    alignment.metadata.data_type.code() +
    " (" +
    alignment.metadata.data_type.molecule_type() +
    ")",
  )
  let match_character = match alignment.metadata.match_character {
    Some(value) => value
    None => "none"
  }
  println(
    "  missing=" +
    alignment.metadata.missing_character +
    ", gap=" +
    alignment.metadata.gap_character +
    ", matchchar=" +
    match_character,
  )
  for sequence in alignment.sequences {
    println(
      "  " +
      sequence.id +
      ": " +
      sequence.aligned_sequence +
      " -> " +
      sequence.sequence,
    )
  }

  println("\n2. Inspect coordinate paths and residue mapping")
  let path = alignment.coordinate_path()
  for row = 0; row < alignment.num_sequences(); row = row + 1 {
    println(
      "  " +
      alignment.sequences[row].id +
      ": " +
      nexus_integer_array_text(path[row]),
    )
  }
  println(
    "  reference residue 5 -> query one residue " +
    alignment.map_position(0, 1, 5).unwrap_or(-1).to_string(),
  )
  println(
    "  reference residue 7 -> query one residue " +
    alignment.map_position(0, 1, 7).unwrap_or(-1).to_string() +
    " (-1 denotes a gap)",
  )

  println("\n3. Calculate all-pairs statistics")
  let counts = alignment.counts()
  println(
    "  pairs=" +
    counts.pairs.to_string() +
    ", aligned=" +
    counts.aligned.to_string() +
    ", identities=" +
    counts.identities.to_string() +
    ", gaps=" +
    counts.gap_columns.to_string(),
  )
  let consensus = alignment.consensus(minimum_fraction=0.67) catch {
    AlignNexusError(message) => abort("consensus failed: " + message)
  }
  println("  consensus: " + consensus)
  println(
    "  occupancy at source gap column: " + alignment.occupancy()[4].to_string(),
  )

  println("\n4. Write canonical interleaved NEXUS and parse it again")
  let canonical = @src.align_nexus_write(
    alignment,
    interleave=Some(true),
    block_width=4,
  ) catch {
    AlignNexusError(message) => abort("NEXUS serialization failed: " + message)
  }
  let round_trip = @src.align_nexus_parse(canonical) catch {
    AlignNexusError(message) => abort("round-trip parsing failed: " + message)
  }
  println("  output bytes: " + canonical.length().to_string())
  println(
    "  quoted IDs and rows preserved: " +
    (round_trip.sequences == alignment.sequences).to_string(),
  )
  println(
    "  coordinate path preserved: " +
    (round_trip.coordinate_path() == alignment.coordinate_path()).to_string(),
  )

  println("\n5. Remove all-gap source columns during construction")
  let compact = @src.align_nexus_from_aligned(
    ["reference", "query"],
    ["A-C--GT", "ATC--GT"],
    @src.AlignNexusDna,
  ) catch {
    AlignNexusError(message) =>
      abort("alignment construction failed: " + message)
  }
  println(
    "  source columns=" +
    compact.source_alignment_length().to_string() +
    ", alignment columns=" +
    compact.alignment_length().to_string() +
    ", removed=" +
    compact.metadata.removed_all_gap_columns.to_string(),
  )
  println("  compact reference: " + compact.sequences[0].aligned_sequence)
}
